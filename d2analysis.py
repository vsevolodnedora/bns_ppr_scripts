from __future__ import division
from sys import path
path.append('modules/')
from _curses import raw
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib import ticker
import matplotlib.pyplot as plt
from matplotlib import rc
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
import scivis.units as ut # for tmerg
import statsmodels.formula.api as smf
import scipy.optimize as opt
from math import pi, sqrt
import matplotlib as mpl
import pandas as pd
import numpy as np
import itertools
import os.path
import cPickle
import time
import copy
import h5py
import csv
import os

from scidata.utils import locate
import scidata.carpet.hdf5 as h5
import scidata.xgraph as xg
from scidata.carpet.interp import Interpolator

from scipy import interpolate
cmap = plt.get_cmap("viridis")
# from sklearn.linear_model import LinearRegression-
from scipy.optimize import fmin
from matplotlib.ticker import AutoMinorLocator, FixedLocator, NullFormatter, \
    MultipleLocator
from matplotlib.colors import LogNorm, Normalize

from general import *
from lists import *
from filework import *
from units import time_constant, volume_constant, energy_constant

from math import pi, log10
import time

class LOAD_STORE_DATASETS:
    """
    Allows easy access to a scidata datasets of 2d data of a given simulation,
    loading them only if they are needed, and storing them for future accesses.

    Assumes that the simulation data is stired in /output-xxxx/data/ folders, where
    'xxxx' stands for the number.
    To know, which iterations are stored in what output-xxxx it first loads an ascii
    file <'file_for_it'> which should be present in all output-xxxx and have a
    structure with columns: 1:it 2:time (other columns do not matter)

    The 2d datasets, located in output-xxxx, expected to be named like:
    'rho.xy.h5'. The list of possible variables, (like 'rho') is set in
    <'self.list_v_ns'>. The list of possible planes (like 'xy') in set in
    <'self.list_planes'>.

    Logic:
        Every time when the dataset is requested via <'get_dataset(it, plane, v_n)'>,
        the class do:
            1. Checks what output-xxxx contains iteration 'it'
            2. Checks if the dataset for this output, plane and variable name 'v_n'
                has already been loaded and is present in the storage
                <'self.dataset_matrix[]'>
                If so: it will return the dataset from the storage.
                If not: it will load the required dataset and add it to the storage
                for future uses.

    """

    def __init__(self, sim):

        self.sim = sim
        self.nlevels = 7
        self.gen_set = {'nlevels':7,
                        'sim': sim,
                        'file_for_it': 'H.norm2.asc',
                        'iterations':0,
                        'indir': Paths.gw170817 + sim + '/',
                        'outdir': Paths.ppr_sims + sim + '/res_2d/'}

        # self.output_it_map = {}
        self.output_it_map, self.it_time = \
            set_it_output_map(Paths.gw170817+self.sim+'/')

        self.iterations = np.array(self.it_time[:, 0], dtype=int)
        self.times =  np.array(self.it_time[:, 1], dtype=float)

        self.list_v_ns = ['rho', 'Y_e', 'temperature', 's_phi', 'entropy', 'dens_unbnd']
        self.list_planes=['xy', 'xz', 'xy']

        self.set_use_new_output_if_duplicated = False

        self.dataset_matrix = [[[0
                                  for z in range(len(self.list_v_ns))]
                                  for k in range(len(self.list_planes))]
                                  for s in range(len(self.output_it_map.keys()))]

    # def set_it_output_map(self):
    #     """
    #     Loads set of files that have '1:it 2:time ...' structure to get a map
    #     of what output-xxxx contains what iteration (and time)
    #     """
    #     print('-' * 25 + 'LOADING it list ({})'
    #           .format(self.gen_set['file_for_it']) + '-' * 25)
    #     print("\t loading from: {}".format(self.gen_set['indir']))
    #     files = locate(self.gen_set['file_for_it'], root=self.gen_set["indir"], followlinks=True)
    #     # remove folders like 'collated'
    #     selected = []
    #     for file in files:
    #         if file.__contains__('output-'):
    #             selected.append(file)
    #     # for overall count of number of iterations and files
    #     it_time = np.zeros(2)
    #     for file in selected:
    #         o_name = file.split('/')
    #         o_dir = ''
    #         for o_part in o_name:
    #             if o_part.__contains__('output-'):
    #                 o_dir = o_part
    #         if o_dir == '':
    #             raise NameError("Did not find output-xxxx in {}".format(o_name))
    #         it_time_i = np.loadtxt(file, usecols=(0,1))
    #         self.output_it_map[o_dir] = it_time_i
    #         it_time = np.vstack((it_time, it_time_i))
    #     it_time = np.delete(it_time, 0, 0)
    #     print('outputs:{} iterations:{} [{}->{}]'.format(len(selected),
    #                                                      len(it_time[:,0]),
    #                                                      int(it_time[:,0].min()),
    #                                                      int(it_time[:,0].max())))
    #     print('-' * 30 + '------DONE-----' + '-' * 30)
    #
    #     self.output_it_map, it_time = set_it_output_map(Paths.gw170817+self.sim+'/')
    #
    #     return it_time

    def check_v_n(self, v_n):
        if v_n not in self.list_v_ns:
            raise NameError("v_n:{} not in the v_n list (in the class)\n{}".format(v_n, self.list_v_ns))

    def check_plane(self, plane):
        if plane not in self.list_planes:
            raise NameError("plane:{} not in the plane_list (in the class)\n{}".format(plane, self.list_planes))

    def i_v_n(self, v_n):
        self.check_v_n(v_n)
        return int(self.list_v_ns.index(v_n))
    #
    def i_plane(self, plane):
        self.check_plane(plane)
        return int(self.list_planes.index(plane))

    def load_dataset(self, o_dir, plane, v_n):
        fname = v_n + '.' + plane + '.h5'
        files = locate(fname, root=self.gen_set['indir'] + o_dir +'/', followlinks=False)
        print("\t Loading: {} plane:{} v_n:{} dataset ({} files)"
              .format(o_dir, plane, v_n, len(files)))
        if len(files) > 1:
            raise ValueError("More than 1 file ({}) found. \nFile:{} location:{}"
                             "\nFiles: {}"
                             .format(len(files), fname, o_dir, files))
        if len(files) == 0:
            raise ValueError("NO fils found. \nlocation:{}"
                             .format(fname, o_dir))
        dset = h5.dataset(files)
        # print("\t loading it:{} plane:{} v_n:{} dset:{}"
        #       .format(o_dir, plane, v_n, dset))
        dset.get_grid().mesh()
        self.dataset_matrix[self.i_output(o_dir)][self.i_plane(plane)][self.i_v_n(v_n)] = dset

    def i_output(self, o_dir):
        if o_dir not in self.output_it_map.keys():
            raise NameError("plane:{} not in the plane_list (in the class)\n{}"
                            .format(o_dir, self.output_it_map.keys()))

        return int(self.output_it_map.keys().index(o_dir))

    def is_dataset_loaded(self, o_dir, plane, v_n):
        if isinstance(self.dataset_matrix[self.i_output(o_dir)][self.i_plane(plane)][self.i_v_n(v_n)], int):
            self.load_dataset(o_dir, plane, v_n)

    def it_to_output_dir(self, it):
        req_output_data_dir = []
        for output_data_dir in self.output_it_map.keys():
            if int(it) in np.array(self.output_it_map[output_data_dir], dtype=int)[:, 0]:
                req_output_data_dir.append(output_data_dir)

        if len(req_output_data_dir) > 1:
            if self.set_use_new_output_if_duplicated:
                print("Warning: it:{} is found in multiple outputs:{}"
                             .format(it, req_output_data_dir))
                return req_output_data_dir[0]

            raise ValueError("it:{} is found in multiple outputs:{}\n"
                             "to overwrite, set 'set_use_new_output_if_duplicated=True' "
                             .format(it, req_output_data_dir))
        elif len(req_output_data_dir) == 0:
            raise ValueError("it:{} not found in a output_it_map:\n{}\n"
                                 .format(it, self.output_it_map.keys()))
        else:
            return req_output_data_dir[0]

    def get_dataset(self, it, plane, v_n):
        o_dir = self.it_to_output_dir(it)
        self.is_dataset_loaded(o_dir, plane, v_n)
        dset = self.dataset_matrix[self.i_output(o_dir)][self.i_plane(plane)][self.i_v_n(v_n)]
        if not it in dset.iterations:
            it__ = int(dset.iterations[find_nearest_index(np.array(dset.iterations), it)])
            raise ValueError("Iteration it:{} (located in {}) \n"
                             "not in the dataset list. Closest:{} Full list:\n{}"
                             .format(it, o_dir, it__, dset.iterations))
        return dset

    def del_dataset(self, it, plane, v_n):
        o_dir = self.it_to_output_dir(it)
        self.dataset_matrix[self.i_output(o_dir)][self.i_plane(plane)][self.i_v_n(v_n)] = 0

    # def get_time(self, it):
    #
    #     time = self.it_time[np.where(self.it_time[:,0] == it), 1]
    #     time = [item for sublist in time for item in sublist]
    #     print(time)
    #     if len(time) == 2:
    #         Printcolor.yellow("for it:{} more than one timestep found {}"
    #                          .format(it, time))
    #         if time[0] == time[1]:
    #             return float(time[0]) * time_constant / 1000
    #         else:
    #             raise ValueError("for it:{} more than one timestep found {}"
    #                              .format(it, time))
    #     if len(time) == 0:
    #         raise ValueError("for it:{} no timesteps found"
    #                          .format(it))
    #     return float(time[0]) * time_constant / 1000

    def load_all(self, plane, v_n):

        print('-' * 25 + 'LOADING ALL DATASETS ({})'
              .format(self.gen_set['file_for_it']) + '-' * 25)
        Printcolor.yellow("Warning: loading all {} datasets "
                          "is a slow process".format(len(self.output_it_map.keys())))
        for o_dir in self.output_it_map.keys():
            self.is_dataset_loaded(o_dir, plane, v_n)

        self.set_all_it_times_from_outputs(plane, v_n)

        print('-' * 30 + '------DONE-----' + '-' * 30)

    def get_all_iterations_times(self, plane, v_n):

        iterations = []
        times = []
        for o_dir in self.output_it_map.keys():
            if isinstance(self.dataset_matrix[self.i_output(o_dir)][self.i_plane(plane)][self.i_v_n(v_n)], int):
                raise ValueError("Not all datasets are loaded. Missing: {}".format(o_dir))
            dset = self.dataset_matrix[self.i_output(o_dir)][self.i_plane(plane)][self.i_v_n(v_n)]
            # iterations.append(dset.iterations)
            for it in dset.iterations:
                iterations.append(it)
                time = dset.get_time(it) * 0.004925794970773136 / 1000
                times.append(time)
                # print("it:{}, time:{}".format(it, time))

        assert len(list(set(iterations))) == len(list(set(times)))

        iterations = np.sort(list(set(iterations)))
        times = np.sort(list(set(times)))

        return iterations, times

    def set_all_it_times_from_outputs(self, plane, v_n):

        self.iterations, self.times = self.get_all_iterations_times(plane, v_n)
        print('\tIterations [{}->{}] and times [{:.3f}->{:.3f}] have been reset.'
              .format(self.iterations[0], self.iterations[-1],
                      self.times[0], self.times[-1]))

        # return list(set([item for sublist in iterations for item in sublist])), \
        #        list(set([item for sublist in iterations for item in sublist]))

    # def get_all_timesteps(self, plane, v_n):
    #
    #     iterations = self.get_all_iterations(plane, v_n)
    #     times = []
    #     for iteration in iterations:
    #         times.append(self.get_time(iteration))
    #     return times


class EXTRACT_STORE_DATA(LOAD_STORE_DATASETS):
    """
    blablabla
    """

    def __init__(self, sim):

        LOAD_STORE_DATASETS.__init__(self, sim)

        # self.gen_set = {'nlevels': 7,
        #                 'file_for_it': 'H.norm2.asc',
        #                 'iterations': 0,
        #                 'indir': Paths.gw170817 + sim + '/',
        #                 'outdir': Paths.ppr_sims + sim + '/2d/'}

        # self.list_v_ns   = ['rho', 'Y_e', 'temperature', 's_phi', 'entropy', 'dens_unbnd']
        # self.list_planes = ['xy', 'xz']
        self.v_n_map = {
            'rho':          "HYDROBASE::rho",
            'Y_e':          "HYDROBASE::Y_e",
            's_phi':        "BNSANALYSIS::s_phi",
            'temperature':  "HYDROBASE::temperature",
            'entropy':      "HYDROBASE::entropy",
            'dens_unb':     "BNSANALYSIS::dens_unbnd"
        }


        # self.output_it_map = {}
        # self.it_time = self.set_it_output_map()

        self.data_matrix = [[[0
                            for z in range(len(self.list_v_ns))]
                            for k in range(len(self.list_planes))]
                            for s in range(len(self.it_time[:,0]))]

        self.grid_matrix = [[[0
                            for z in range(len(self.list_v_ns))]
                            for k in range(len(self.list_planes))]
                            for s in range(len(self.it_time[:,0]))]

    def check_it(self, it):
        if not int(it) in np.array(self.it_time[:,0], dtype=int):
            it_ = int(self.it_time[find_nearest_index(self.it_time[:,0], it), 0])
            raise NameError("it:{} not in the list on iterations: Closest one: {}"
                            .format(it, it_))

        idx = np.where(np.array(self.it_time[:,0], dtype=int) == int(it))

        if len(idx) == 0:
            raise ValueError("For it:{} NO it are found in the it_time[:,0]".format(it))

        if len(idx) > 1:
            raise ValueError("For it:{} multiple it are found in the it_time[:,0]".format(it))

        # print("it:{} idx:{}".format(it, idx))

    def i_it(self, it):
        self.check_it(it)
        idx = list(self.it_time[:,0]).index(int(it))
        return idx
        #
        #
        #
        #
        # self.check_it(it)
        # idx = np.array(np.where(np.array(self.it_time[:,0], dtype=int) == int(it)), dtype=int)
        # idx = list(set([item for sublist in idx for item in sublist]))
        # print("it:{} idx:{}, type:{} len:{}".format(it, idx, type(idx), len(idx)))
        # return int(idx)

    # ---------- GRID

    def extract_grid(self, it, plane, v_n):
        print("\t extracting grid it:{} plane:{} v_n:{}".format(it, plane, v_n))
        dset = self.get_dataset(it, plane, v_n)
        self.grid_matrix[self.i_it(it)][self.i_plane(plane)][self.i_v_n(v_n)] = \
            dset.get_grid(iteration=it)

    def is_grid_extracted(self, it, plane, v_n):

        if isinstance(self.grid_matrix[self.i_it(it)][self.i_plane(plane)][self.i_v_n(v_n)], int):
            self.extract_grid(it, plane, v_n)

    def get_grid(self, it, plane, v_n):

        self.check_plane(plane)
        self.check_v_n(v_n)
        self.is_grid_extracted(it, plane, v_n)

        return self.grid_matrix[self.i_it(it)][self.i_plane(plane)][self.i_v_n(v_n)]

    def del_grid(self, it, plane, v_n):

        self.check_plane(plane)
        self.check_v_n(v_n)

        self.grid_matrix[self.i_it(it)][self.i_plane(plane)][self.i_v_n(v_n)] = 0

    # ---------- DATA

    def extract(self, it, plane, v_n):

        print("\t extracting it:{} plane:{} v_n:{}".format(it, plane, v_n))
        dset = self.get_dataset(it, plane, v_n)
        try:
            data = dset.get_grid_data(self.get_grid(it, plane, v_n),
                                     iteration=it,
                                     variable=self.v_n_map[v_n])
            self.data_matrix[self.i_it(it)][self.i_plane(plane)][self.i_v_n(v_n)] = data
        except KeyError:
            raise KeyError("Wrong Key. Data not found. dset contains:{} attmeped:{} {} {}".format(dset.metadata[0],
                                                                                         self.v_n_map[v_n],
                                                                                                  plane,
                                                                                                  it))

    def is_data_extracted(self, it, plane, v_n):

        if isinstance(self.data_matrix[self.i_it(it)][self.i_plane(plane)][self.i_v_n(v_n)], int):
            self.extract(it, plane, v_n)

    def get_data(self, it, plane, v_n):
        self.check_plane(plane)
        self.check_v_n(v_n)

        self.is_data_extracted(it, plane, v_n)

        return self.data_matrix[self.i_it(it)][self.i_plane(plane)][self.i_v_n(v_n)]

    def del_data(self, it, plane, v_n):

        self.check_plane(plane)
        self.check_v_n(v_n)

        self.data_matrix[self.i_it(it)][self.i_plane(plane)][self.i_v_n(v_n)] = 0

    # ----------- TIME

    # def get_time(self, it):
    #     self.check_it(it)
    #     return float(self.it_time[np.where(self.it_time[:,0] == it), 1][0]) * time_constant / 1000



    # def get_time_(self, it):
    #     return self.get_time__(it)

class EXTRACT_FOR_RL(EXTRACT_STORE_DATA):

    def __init__(self, sim):
        EXTRACT_STORE_DATA.__init__(self, sim)

        self.list_grid_v_ns = ["x", "y", "z", "delta"]


        self.extracted_grid_matrix = [[[[np.zeros(0,)
                            for z in range(len(self.list_grid_v_ns))]
                            for j in range(self.nlevels)]
                            for k in range(len(self.list_planes))]
                            for s in range(len(self.it_time[:,0]))]

        self.extracted_data_matrix = [[[[np.zeros(0,)
                            for z in range(len(self.list_v_ns))]
                            for j in range(self.nlevels)]
                            for k in range(len(self.list_planes))]
                            for s in range(len(self.it_time[:,0]))]

        self.default_v_n = "rho"

    def check_rl(self, rl):
        if rl < 0:
            raise ValueError("Unphysical rl:{} ".format(rl))
        if rl < 0 or rl > self.nlevels:
            raise ValueError("rl is not in limits: {}"
                             .format(rl, self.nlevels))

    def i_rl(self, rl):
        return int(rl)

    def check_grid_v_n(self, v_n):
        if not v_n in self.list_grid_v_ns:
            raise NameError("v_n:{} not in list_grid_v_ns:{}"
                            .format(v_n, self.list_grid_v_ns))

    def i_gr_v_n(self, v_n):
        return int(self.list_grid_v_ns.index(v_n))

    def __extract_grid_data_rl(self, it, plane, rl, grid):

        coords = grid.coordinates()
        mesh = list(np.meshgrid(coords[self.i_rl(rl)], indexing='ij'))
        points = np.column_stack([x.flatten() for x in mesh])
        xyz = []
        delta = grid.levels[self.i_rl(rl)].delta
        # for x in mesh:
        #     xyz.append(x)
        #     print("x: {} ".format(x.shape))
        # print(len(coords))
        # # print(mesh)
        # print(len(mesh))
        # print(len(xyz))
        # print(points.shape)

        self.extracted_grid_matrix[self.i_it(it)][self.i_plane(plane)][self.i_rl(rl)][self.i_gr_v_n("delta")] = delta
        if plane == "xy":
            x, y = grid.mesh()[rl]
            # print(x.shape)
            # print(y.shape)
            self.extracted_grid_matrix[self.i_it(it)][self.i_plane(plane)][self.i_rl(rl)][self.i_gr_v_n("x")] = x
            self.extracted_grid_matrix[self.i_it(it)][self.i_plane(plane)][self.i_rl(rl)][self.i_gr_v_n("y")] = y
        elif plane == "xz":
            x, z = grid.mesh()[rl]
            # print(x.shape)
            # print(z.shape)
            self.extracted_grid_matrix[self.i_it(it)][self.i_plane(plane)][self.i_rl(rl)][self.i_gr_v_n("x")] = x
            self.extracted_grid_matrix[self.i_it(it)][self.i_plane(plane)][self.i_rl(rl)][self.i_gr_v_n("z")] = z
        else:
            raise NameError("Plane: {} is not recognized")
        # exit(0)

    def extract_grid_data_rl(self, it, plane, rl, v_n):

        self.check_grid_v_n(v_n)

        for data_v_n in self.list_v_ns:
            # scrolling through all possible v_ns, looking for the one loaded
            if not isinstance(self.grid_matrix[self.i_it(it)][self.i_plane(plane)][self.i_v_n(data_v_n)], int):
                grid = self.get_grid(it, plane, data_v_n)
                self.__extract_grid_data_rl(it, plane, rl, grid)
                return 0
        print("\tNo pre-loaded data found. Loading default ({}))".format(self.default_v_n))
        grid = self.get_grid(it, plane, self.default_v_n)
        self.__extract_grid_data_rl(it, plane, rl, grid)
        return 0


    def is_grid_extracted_for_rl(self, it, plane, rl, v_n):

        arr = self.extracted_grid_matrix[self.i_it(it)][self.i_plane(plane)][self.i_rl(rl)][self.i_gr_v_n(v_n)]
        if len(arr) == 0:
            self.extract_grid_data_rl(it, plane, rl, v_n)

    def get_grid_v_n_rl(self, it, plane, rl, v_n):

        self.check_grid_v_n(v_n)
        self.check_plane(plane)
        self.check_it(it)
        self.check_rl(rl)

        self.is_grid_extracted_for_rl(it, plane, rl, v_n)

        data = self.extracted_grid_matrix[self.i_it(it)][self.i_plane(plane)][self.i_rl(rl)][self.i_gr_v_n(v_n)]
        return data



    def extract_data_rl(self, it, plane, rl, v_n):

        data = self.get_data(it, plane, v_n)
        data = np.array(data[self.i_rl(rl)])

        self.extracted_data_matrix[self.i_it(it)][self.i_plane(plane)][self.i_rl(rl)][self.i_v_n(v_n)] = data

    def is_data_extracted_for_rl(self, it, plane,rl, v_n):

        arr = self.extracted_data_matrix[self.i_it(it)][self.i_plane(plane)][self.i_rl(rl)][self.i_v_n(v_n)]
        if len(arr) == 0:
            self.extract_data_rl(it, plane, rl, v_n)

    def get_data_rl(self, it, plane, rl, v_n):

        self.check_v_n(v_n)
        self.check_plane(plane)
        self.check_it(it)
        self.check_rl(rl)

        self.is_data_extracted_for_rl(it, plane, rl, v_n)

        data = self.extracted_data_matrix[self.i_it(it)][self.i_plane(plane)][self.i_rl(rl)][self.i_v_n(v_n)]
        return data




class COMPUTE_STORE(EXTRACT_FOR_RL):

    def __init__(self, sim):
        EXTRACT_FOR_RL.__init__(self, sim)

class MASK_STORE(COMPUTE_STORE):

    def __init__(self, fname, symmetry=None):
        COMPUTE_STORE.__init__(self, fname, symmetry)

        rho_const = 6.176269145886162e+17
        self.mask_setup = {'rm_rl': True,  # REMOVE previouse ref. level from the next
                           'rho': [6.e4 / rho_const, 1.e13 / rho_const],  # REMOVE atmo and NS
                           'lapse': [0.15, 1.]} # remove apparent horizon

        self.mask_matrix = [np.ones(0, dtype=bool) for x in range(self.nlevels)]

        self.list_mask_v_n = ["x", "y", "z"]


    def compute_mask(self):

        nlevelist = np.arange(self.nlevels, 0, -1) - 1

        x = []
        y = []
        z = []

        for ii, rl in enumerate(nlevelist):
            x.append(self.get_comp_data(rl, "x")[3:-3, 3:-3, 3:-3])
            y.append(self.get_comp_data(rl, "y")[3:-3, 3:-3, 3:-3])
            z.append(self.get_comp_data(rl, "z")[3:-3, 3:-3, 3:-3])
            mask = np.ones(x[ii].shape, dtype=bool)
            if ii > 0 and self.mask_setup["rm_rl"]:
                x_ = (x[ii][:, :, :] <= x[ii - 1][:, 0, 0].max()) & (
                        x[ii][:, :, :] >= x[ii - 1][:, 0, 0].min())
                y_ = (y[ii][:, :, :] <= y[ii - 1][0, :, 0].max()) & (
                        y[ii][:, :, :] >= y[ii - 1][0, :, 0].min())
                z_ = (z[ii][:, :, :] <= z[ii - 1][0, 0, :].max()) & (
                        z[ii][:, :, :] >= z[ii - 1][0, 0, :].min())
                mask = mask & np.invert((x_ & y_ & z_))

            for v_n in self.mask_setup.keys()[1:]:
                self.check_v_n(v_n)
                if len(self.mask_setup[v_n]) != 2:
                    raise NameError("Error. 2 values are required to set a limit. Give {} for {}"
                                    .format(self.mask_setup[v_n], v_n))
                arr_1 = self.get_comp_data(rl, v_n)[3:-3, 3:-3, 3:-3]
                min_val = float(self.mask_setup[v_n][0])
                max_val = float(self.mask_setup[v_n][1])
                mask_i = (arr_1 > min_val) & (arr_1 < max_val)
                mask = mask & mask_i
                del arr_1
                del mask_i

            self.mask_matrix[rl] = mask

    def is_mask_available(self, rl):
        mask = self.mask_matrix[rl]
        if len(mask) == 0:
            self.compute_mask()

    def get_masked_data(self, rl, v_n):
        self.check_v_n(v_n)
        self.is_available(rl, v_n)
        self.is_mask_available(rl)
        data = np.array(self.get_comp_data(rl, v_n))[3:-3, 3:-3, 3:-3]
        mask = self.mask_matrix[rl]
        return data[mask]

    def __delete__(self, instance):
        instance.dfile.close()
        instance.data_matrix = [[np.zeros(0, )
                                 for x in range(self.nlevels)]
                                 for y in range(len(self.list_all_v_ns))]
        instance.mask_matrix = [np.ones(0, dtype=bool) for x in range(self.nlevels)]

class CYLINDRICAL_GRID:
    """
        Creates a stretched cylindrical grid and allows
        to interpolate any data from carpet grid onto it.
        Stretched means, that the grid consists of 2 parts:
            1) linear distribution in terms of radius (0-15)
            2) logarithmic dist. in terms of radius (15-512)
        Class stores the grid information in its own variables
        that can be accessed directly or through
        `get_new_grid(v_n)` method

        Requirements:
            > dictionary grid_info{} that describes the grid:
            > class `carpet_grid` from scidata
        Usage:
            to access the new grid mesh arrays use:
                get_new_grid(v_n)
            to do the interpolation of arr, use
                get_int_arr(arr)
    """

    def __init__(self, grid_info):

        self.grid_info = grid_info # "n_r", "n_phi", "n_z"
        self.grid_type = "cyl"
        self.list_int_grid_v_ns = ["x_cyl", "y_cyl", "z_cyl",
                          "r_cyl", "phi_cyl",
                          "dr_cyl", "dphi_cyl", "dz_cyl"]

        print('-' * 25 + 'INITIALIZING CYLINDRICAL GRID' + '-' * 25)

        phi_cyl, r_cyl, z_cyl, \
        self.dphi_cyl_3d, self.dr_cyl_3d, self.dz_cyl_3d = self.get_phi_r_z_grid()

        self.r_cyl_3d, self.phi_cyl_3d, self.z_cyl_3d \
            = np.meshgrid(r_cyl, phi_cyl, z_cyl, indexing='ij')
        self.x_cyl_3d = self.r_cyl_3d * np.cos(self.phi_cyl_3d)
        self.y_cyl_3d = self.r_cyl_3d * np.sin(self.phi_cyl_3d)

        # print(self.z_cyl_3d[0, 0, :]); exit(1)

        print("\t GRID: [phi:r:z] = [{}:{}:{}]".format(len(phi_cyl), len(r_cyl), len(z_cyl)))

        print('-' * 30 + '------DONE-----' + '-' * 30)
        print('\n')

    # cylindrical grid
    @staticmethod
    def make_stretched_grid(x0, x1, x2, nlin, nlog):
        assert x1 > 0
        assert x2 > 0
        x_lin_f = np.linspace(x0, x1, nlin)
        x_log_f = 10.0 ** np.linspace(log10(x1), log10(x2), nlog)
        return np.concatenate((x_lin_f, x_log_f))

    def get_phi_r_z_grid(self):

        # extracting grid info
        n_r = self.grid_info["n_r"]
        n_phi = self.grid_info["n_phi"]
        n_z = self.grid_info["n_z"]

        # constracting the grid
        r_cyl_f = self.make_stretched_grid(0., 15., 512., n_r, n_phi)
        z_cyl_f = self.make_stretched_grid(0., 15., 512., n_r, n_phi)
        phi_cyl_f = np.linspace(0, 2 * np.pi, n_phi)

        # edges -> bins (cells)
        r_cyl = 0.5 * (r_cyl_f[1:] + r_cyl_f[:-1])
        z_cyl = 0.5 * (z_cyl_f[1:] + z_cyl_f[:-1])
        phi_cyl = 0.5 * (phi_cyl_f[1:] + phi_cyl_f[:-1])

        # 1D grind -> 3D grid (to mimic the r, z, phi structure)
        dr_cyl = np.diff(r_cyl_f)[:, np.newaxis, np.newaxis]
        dphi_cyl = np.diff(phi_cyl_f)[np.newaxis, :, np.newaxis]
        dz_cyl = np.diff(z_cyl_f)[np.newaxis, np.newaxis, :]

        return phi_cyl, r_cyl, z_cyl, dphi_cyl, dr_cyl, dz_cyl

    def get_xi(self):
        return np.column_stack([self.x_cyl_3d.flatten(),
                                self.y_cyl_3d.flatten(),
                                self.z_cyl_3d.flatten()])

    def get_shape(self):
        return self.x_cyl_3d.shape

    # generic methods to be present in all INTERPOLATION CLASSES
    # def get_int_arr(self, arr_3d):
    #
    #     # if not self.x_cyl_3d.shape == arr_3d.shape:
    #     #     raise ValueError("Passed for interpolation 3d array has wrong shape:\n"
    #     #                      "{} Expected {}".format(arr_3d.shape, self.x_cyl_3d.shape))
    #     xi = np.column_stack([self.x_cyl_3d.flatten(),
    #                           self.y_cyl_3d.flatten(),
    #                           self.z_cyl_3d.flatten()])
    #     F = Interpolator(self.carpet_grid, arr_3d, interp=1)
    #     res_arr_3d = F(xi).reshape(self.x_cyl_3d.shape)
    #     return res_arr_3d

    # def get_int_grid(self, v_n):
    #     if v_n == "x_cyl":
    #         return self.x_cyl_3d
    #     elif v_n == "y_cyl":
    #         return self.y_cyl_3d
    #     elif v_n == "z_cyl":
    #         return self.z_cyl_3d
    #     elif v_n == "r_cyl":
    #         return self.r_cyl_3d
    #     elif v_n == "phi_cyl":
    #         return self.phi_cyl_3d
    #     elif v_n == "dr_cyl":
    #         return self.dr_cyl_3d
    #     elif v_n == "dphi_cyl":
    #         return self.dphi_cyl_3d
    #     elif v_n == "dz_cyl":
    #         return self.dz_cyl_3d
    #     else:
    #         raise NameError("v_n: {} not recogized in grid. Available:{}"
    #                         .format(v_n, self.grid_v_ns))

    def get_int_grid(self, v_n, projection='xyz'):

        if v_n == "x_cyl":
            d3arr = self.x_cyl_3d
        elif v_n == "y_cyl":
            d3arr = self.y_cyl_3d
        elif v_n == "z_cyl":
            d3arr = self.z_cyl_3d
        elif v_n == "r_cyl":
            d3arr = self.r_cyl_3d
        elif v_n == "phi_cyl":
            d3arr = self.phi_cyl_3d
        elif v_n == "dr_cyl":
            d3arr = self.dr_cyl_3d
        elif v_n == "dphi_cyl":
            d3arr = self.dphi_cyl_3d
        elif v_n == "dz_cyl":
            d3arr = self.dz_cyl_3d
        else:
            raise NameError("v_n: {} not recogized in grid. Available:{}"
                            .format(v_n, self.list_int_grid_v_ns))

        if projection == 'xyz':
            return d3arr
        elif projection == 'xy':
            return d3arr[:, :, 0]
        elif projection == 'xz':
            return d3arr[0, :, :]
        elif projection == 'yz':
            return d3arr[:, 0, :]
        else:
            raise NameError("projection: {} not recogized".format(projection))

    def save_grid(self, sim, projection):

        grid_type = self.grid_type

        path = Paths.ppr_sims + sim + "/res_2d/"
        if projection == 'xyz' or projection == None:
            fname = path + self.grid_type  + '_' + 'grid.h5'
        else:
            fname = path + self.grid_type + '_' + 'grid' + '_' + projection + '.h5'

        outfile = h5py.File(fname, "w")

        if not os.path.exists(path):
            os.makedirs(path)

        # print("Saving grid...")
        for v_n in self.list_int_grid_v_ns:
            outfile.create_dataset(v_n, data=self.get_int_grid(v_n, projection))
        outfile.close()


class INTERPOLATE_STORE_2D(EXTRACT_STORE_DATA):
    """
    blablabla
    """

    def __init__(self, sim):



        self.gen_set = {'nlevels':7,
                        'file_for_it': 'H.norm2.asc',
                        'iterations':0,
                        'indir': Paths.gw170817 + sim + '/',
                        'outdir': Paths.ppr_sims + sim + '/res_2d/',
                        }
        self.grid_set = {'type': 'cylindrical', 'n_r': 150, 'n_phi': 150, 'n_z': -100}

        # self.output_it_map = {}

        self.list_v_ns = ['rho', 'Y_e', 'temperature', 's_phi', 'entropy', 'dens_unbnd']
        self.list_planes=['xy', 'xz']


        EXTRACT_STORE_DATA.__init__(self, sim)
        # print(self.i_it(1052672)) WORKS
        # self.data_cl = data_cl


        self.intdata_matrix = [[[0
                                  for z in range(len(self.list_v_ns))]
                                  for k in range(len(self.list_planes))]
                                  for s in range(len(self.it_time[:,0]))]

        # select a grid class

        if self.grid_set['type'] == 'cylindrical':
            self.new_grid_cl = CYLINDRICAL_GRID(self.grid_set)
        else:
            raise NameError("Grid:{} is not recognized")

    # def tmp(self):
    #     print(self.i_it(1052672))
    #     print(self.i_plane('xy'))

    def do_interpolate(self, it, plane, v_n):

        print("\t Interpolating it:{} plane:{} v_n:{} onto {} grid"
              .format(it, plane, v_n, self.grid_set['type']))

        data = self.get_data(it, plane, v_n)
        # print("\t using it:{} plane:{} v_n:{} dset:{}"
        #       .format(it, plane, v_n, data))
        grid = self.get_grid(it, plane, v_n)
        # print(grid)
        v_n_interpolator = Interpolator(grid, data, interp=0)
        # print("data")
        # print(data)
        # print('\n')

        self.coords = grid.coordinates()
        for i in range(len(grid)):
            if i == len(grid)-1:
                mesh = list(np.meshgrid(self.coords[i], indexing='ij'))
                points = np.column_stack([x.flatten() for x in mesh])
                for x in mesh:
                    print("x: {} ".format(x.shape))
                print("x[0]: {}".format(x[0].shape))
                print("x[1]: {}".format(x[1].shape))
                print("data: {} ".format(np.array(data[i]).shape))
                print("data: {}".format(np.array(data[i]).shape))
                print(points.shape)
                print(len(mesh))
                # exit(1)
                data = np.array(data[i])
                from matplotlib import colors
                fig = plt.figure()
                ax = fig.add_subplot(111)
                norm = colors.LogNorm(vmin=data.min(), vmax=data.max())
                ax.pcolormesh(x[0], x[1], data.T, norm=norm, cmap="inferno_r")
                plt.ylim(0, 50)
                plt.title(r"$xz$-plane")
                plt.savefig('{}'.format(Paths.ppr_sims + self.sim + '/res_2d/' + "{}_{}_{}.png"
                                        .format(it, plane, v_n)), bbox_inches='tight', dpi=128)
                print("saved pi_test2.png")
                plt.close()


            # print(data[i])
        exit(1)


        if self.grid_set['type'] == 'cylindrical':

            x_cyl_2d = self.new_grid_cl.get_int_grid('x_cyl', plane)
            y_cyl_2d = self.new_grid_cl.get_int_grid('y_cyl', plane)
            z_cyl_2d = self.new_grid_cl.get_int_grid('z_cyl', plane)

            # print(x_cyl_2d[:, 0])
            # print('\n')
            # print(z_cyl_2d[0, :])
            # print('\n')

            if plane == 'xy':
                xi = np.column_stack([x_cyl_2d.flatten(), y_cyl_2d.flatten()])
            elif plane == 'xz':
                xi = np.column_stack([x_cyl_2d.flatten(), z_cyl_2d.flatten()])
            elif plane == 'yz':
                xi = np.column_stack([y_cyl_2d.flatten(), z_cyl_2d.flatten()])
            else:
                raise NameError("plane:{} not supported in a new grid creation".format(plane))

            # print("\n")
            # for xii in xi:
            #     print(xii)
            # print("\n")
            res_arr_2d = v_n_interpolator(xi).reshape(x_cyl_2d.shape)
            # print(res_arr_2d)
            # exit(1)
        else:
            raise NameError("Grid class for '{}' is not found".format(self.grid_set['type']))


        data = res_arr_2d
        from matplotlib import colors
        fig = plt.figure()
        ax = fig.add_subplot(111)
        norm = colors.LogNorm(vmin=data.min(), vmax=data.max())
        ax.pcolormesh(x_cyl_2d[:, 0], z_cyl_2d[0, :], data[:, :].T, norm=norm, cmap="inferno_r")
        plt.ylim(0,50)
        plt.title(r"$xz$-plane")
        plt.savefig('{}'.format(Paths.ppr_sims+self.sim+'/res_2d/' + "{}_{}_{}.png"
                                .format(it, plane, v_n)), bbox_inches='tight', dpi=128)
        print("saved pi_test2.png")
        plt.close()
        # exit(1)


        self.intdata_matrix[self.i_it(it)][self.i_plane(plane)][self.i_v_n(v_n)] = \
            res_arr_2d

    def is_data_interpolated(self, it, plane, v_n):

        # print(it)
        # print(self.i_plane(plane))
        # print(self.i_it(it))

        if isinstance(self.intdata_matrix[self.i_it(it)][self.i_plane(plane)][self.i_v_n(v_n)], int):
            self.do_interpolate(it, plane, v_n)

    def get_int(self, it, plane, v_n):
        self.check_v_n(v_n)
        self.check_plane(plane)
        self.is_data_interpolated(it, plane, v_n)
        return self.intdata_matrix[self.i_it(it)][self.i_plane(plane)][self.i_v_n(v_n)]

    def get_new_grid(self, plane, v_n):
        if self.grid_set['type'] == 'cylindrical':
            return self.new_grid_cl.get_int_grid(v_n, plane)
        else:
            raise NameError("Grid:{} is not recognized")

    def del_int(self, it, plane, v_n):
        self.check_v_n(v_n)
        self.check_plane(plane)
        self.is_data_interpolated(it, plane, v_n)
        self.intdata_matrix[self.i_it(it)][self.i_plane(plane)][self.i_v_n(v_n)] = 0

    def save_all(self, plane, v_n):

        self.load_all(plane, v_n)
        self.iterations, self.times = self.get_all_iterations_times(plane, v_n)
        # self.iterations = np.array(self.iterations, dtype=int)

        self.set_use_new_output_if_duplicated = True

        path = Paths.ppr_sims + self.sim + "/res_2d/"
        outfname = path + self.new_grid_cl.grid_type + \
                            '_' + v_n + '_' + plane + '.h5'

        if os.path.isfile(outfname):
            print("Rewriting file{}".format(outfname))
            os.remove(outfname)
        else:
            print("Writing file:{}".format(outfname))

        outfile = h5py.File(path + self.new_grid_cl.grid_type + \
                            '_' + v_n + '_' + plane + '.h5', "w")

        if not os.path.exists(path):
            os.makedirs(path)

        outfile.create_dataset("iterations", data=self.iterations)
        outfile.create_dataset("times", data=self.times)


        # print("Saving grid...")
        for it, t in zip(self.iterations, self.times):
            print("Writing: it:{} time:{}".format(it, t))
            outfile.create_dataset("it={}".format(int(it)),
                                   data=self.get_int(it, plane, v_n))
        outfile.close()

    # def get_time(self, it):
    #     return self.get_time_(it)


""" --- --- --- """

class LOAD_INT_DATA_2D:

    def __init__(self, sim):

        self.sim = sim

        self.grid_type = "cyl"

        self.list_v_ns = ['rho', 'Y_e', 'temperature', 's_phi', 'entropy', 'dens_unbnd']
        self.list_planes=['xy', 'xz']

        self.set_use_new_output_if_duplicated = False

        self.output_it_map, self.it_time = \
            set_it_output_map(Paths.gw170817+self.sim+'/')

        self.iterations = list(self.it_time[:, 0])
        self.times = list(self.it_time[:, 1])

        self.data_matrix = [[[np.zeros(0,)
                                  for z in range(len(self.list_v_ns))]
                                  for k in range(len(self.list_planes))]
                                  for s in range(len(self.iterations))]

        self.list_grid_v_ns = []
        self.grid_matrix = [[np.zeros(0, )
                             for z in range(len(self.list_grid_v_ns))]
                            for k in range(len(self.list_planes))]

    def check_v_n(self, v_n):
        if v_n not in self.list_v_ns:
            raise NameError("v_n:{} not in the v_n list (in the class)\n{}".format(v_n, self.list_v_ns))

    def check_plane(self, plane):
        if plane not in self.list_planes:
            raise NameError("plane:{} not in the plane_list (in the class)\n{}".format(plane, self.list_planes))

    def i_v_n(self, v_n):
        self.check_v_n(v_n)
        return int(self.list_v_ns.index(v_n))

    def i_plane(self, plane):
        self.check_plane(plane)
        return int(self.list_planes.index(plane))

    def check_it(self, it):
        if not it in self.iterations:
            raise NameError("it:{} is not in the list of iterations \n{}"
                            .format(it, self.iterations))

    def i_it(self, it):
        return int(self.iterations.index(it))

    def load_iteration_list(self):

        path = Paths.ppr_sims + self.sim + "/res_2d/"
        for v_n in self.list_v_ns:
            for plane in self.list_planes:
                fname = path + self.grid_type + \
                            '_' + v_n + '_' + plane + '.h5'
                if os.path.isfile(fname):
                    dfile = h5py.File(fname, "r")
                    iterations = np.array(dfile["iterations"], dtype=int)
                    times = np.array(dfile["times"], dtype=float)
                    self.iterations = list(iterations)
                    self.times = list(times)
                    self.data_matrix = [[[np.zeros(0, )
                                             for z in range(len(self.list_v_ns))]
                                            for k in range(len(self.list_planes))]
                                           for s in range(len(self.iterations))]
                    print("\titerations ({}) are set".format(len(self.iterations)))
                    return 0
        raise IOError("for sim:{} no data for loading iterations found. "
                      "no data for any of the files in a for of "
                      "fname = path + self.grid_type + '_' + v_n + '_' + plane + '.h5'")

    def load_data(self, plane, v_n):
        path = Paths.ppr_sims + self.sim + "/res_2d/"
        fname = path + self.grid_type + '_' + v_n + '_' + plane + '.h5'

        if not os.path.isfile(fname):
            raise IOError("File: {} not found"
                          .format(fname))

        dfile = h5py.File(fname, "r")
        iterations = np.array(dfile["iterations"], dtype=int)

        for it in iterations:

            if int(it) not in self.iterations:
                raise ValueError("it:{} found in dataset:{} "
                                 "is not in preloaded iteration list"
                                 "\n{}"
                                 .format(it, fname, self.iterations))

            self.data_matrix[self.i_it(it)][self.i_plane(plane)][self.i_v_n(v_n)] = \
                np.array(dfile["it={}".format(int(it))])



    def is_data_loaded(self, it, plane, v_n):

        if len(self.data_matrix[self.i_it(it)][self.i_plane(plane)][self.i_v_n(v_n)]) == 0:
            self.load_data(plane, v_n)

    def is_grid_loaded(self, plane = 'xy'):

        if len(self.list_grid_v_ns) == 0:
            self.load_grid(plane)

    def load_grid(self, plane):

        ffpath = ''

        path = Paths.ppr_sims + self.sim + "/res_2d/"
        fname1 = path + self.grid_type + '_' + 'grid.h5'
        if os.path.isfile(fname1):
            ffpath = fname1
        fname2 = path + self.grid_type + '_' + 'grid' + '_' + plane + '.h5'
        if os.path.isfile(fname2):
            ffpath = fname2
        if ffpath == '':
            raise IOError("Grid not found. Neither {} nor {} is found"
                          .format(fname1,fname2))

        dfile = h5py.File(ffpath, "r")

        self.list_grid_v_ns = []
        for v_n in dfile:
            self.list_grid_v_ns.append(v_n)

        self.grid_matrix = [[np.zeros(0, )
                              for z in range(len(self.list_grid_v_ns))]
                              for k in range(len(self.list_planes))]

        for i_v_n, v_n in enumerate(self.list_grid_v_ns):
            self.grid_matrix[self.i_plane(plane)][i_v_n] = np.array(dfile[v_n])





        if plane == 'xyz' or plane == None:
            fname = path + self.grid_type  + '_' + 'grid.h5'
        else:
            fname = path + self.grid_type + '_' + 'grid' + '_' + plane + '.h5'

    def i_gr_v_n(self, v_n):
        return int(self.list_grid_v_ns.index(v_n))

    def get_int_data(self, it, plane, v_n):
        self.check_v_n(v_n)
        self.check_plane(plane)
        self.is_data_loaded(it, plane, v_n)
        return self.data_matrix[self.i_it(it)][self.i_plane(plane)][self.i_v_n(v_n)]

    def get_int_grid(self, plane, v_n):

        self.is_grid_loaded(plane)
        if not v_n in self.list_grid_v_ns:
            raise NameError("v_n:{} not in the loaded list of grid v_ns:{}"
                            .format(v_n, self.list_grid_v_ns))

        return self.grid_matrix[self.i_plane(plane)][self.i_gr_v_n(v_n)]


class ADD_METHODS_FOR_2DINT_DATA(LOAD_INT_DATA_2D):

    def __init__(self, sim):

        LOAD_INT_DATA_2D.__init__(self, sim)

    def fill_pho0_and_phi2pi(self, phi1d_arr, z2d_arr):
        # adding phi = 360 point *copy of phi = 358(
        phi1d_arr = np.append(phi1d_arr, 2 * np.pi)
        z2d_arr = np.vstack((z2d_arr.T, z2d_arr[:, -1])).T
        # adding phi == 0 point (copy of phi=1)
        phi1d_arr = np.insert(phi1d_arr, 0, 0)
        z2d_arr = np.vstack((z2d_arr[:, 0], z2d_arr.T)).T
        return phi1d_arr, z2d_arr


    def get_modified_2d_data(self, it, plane, v_n_x, v_n_y, v_n, mod):

        x_arr = self.get_int_grid(plane, v_n_y)
        y_arr = self.get_int_grid(plane, v_n_x)
        z_arr = self.get_int_data(it, plane, v_n)


        if mod == 'xy slice':
            return np.array(x_arr[:, 0]), np.array(y_arr[0, :]), np.array(z_arr),

        elif mod == 'fill_phi':
            r2d_arr = np.array(x_arr[:, :])
            phi_arr = np.array(y_arr[0, :])
            z2d_arr = np.array(z_arr[:, :])

            phi_arr, rz2d = self.fill_pho0_and_phi2pi(phi_arr, z2d_arr)

            return np.array(x_arr[:, 0]), phi_arr, rz2d

        elif mod == 'fill_phi *r':
            r2d_arr = np.array(x_arr[:, :])
            phi_arr = np.array(y_arr[0, :])
            z2d_arr = np.array(z_arr[:, :])

            rz2d = r2d_arr * z2d_arr
            phi_arr, rz2d = self.fill_pho0_and_phi2pi(phi_arr, rz2d)

            return np.array(x_arr[:, 0]), phi_arr, rz2d

        elif mod == 'fill_phi *r log':

            r2d_arr = np.array(x_arr[:, :])
            phi_arr = np.array(y_arr[0, :])
            z2d_arr = np.array(z_arr[:, :])

            rz2d = r2d_arr * z2d_arr
            phi_arr, rz2d = self.fill_pho0_and_phi2pi(phi_arr, rz2d)

            return np.array(x_arr[:, 0]), phi_arr, np.log10(rz2d)

        elif mod == 'fill_phi -ave(r)':

            r2d_arr = np.array(x_arr[:, :])
            phi_arr = np.array(y_arr[0, :])
            z2d_arr = np.array(z_arr[:, :])

            for i in range(len(x_arr[:, 0])):
                z2d_arr[i, :] = z2d_arr[i, :] - (np.sum(z2d_arr[i, :]) / len(z2d_arr[i, :]))

            phi_arr, rz2d = self.fill_pho0_and_phi2pi(phi_arr, z2d_arr)

            return np.array(x_arr[:, 0]), phi_arr, rz2d

        else:
            raise NameError("Unknown 'mod' parameter:{} ".format(mod))


class PLOT_MANY:
    """
    This is the class to plot many subplots with prescribed data,
    arranging them and their colorbars from a single dictionary
    for a given plot.

    Class takes a INTERPOLATE_STORE() object with methods:
        get_int(task_dic["it"], task_dic["plane"], task_dic["v_n"])
        get_new_grid(task_dic["plane"], task_dic["v_n_x"])
        get_new_grid(task_dic["plane"], task_dic["v_n_y"])

    For plotting in poalr:
        v_n_x = 'phi_cyl'
        v_n_y = 'r_cyl'

    """

    def __init__(self, data_class):

        self.data = data_class

        self.gen_set = {
            "figdir": self.data.gen_set["outdir"],
            "figname": "tst.png",
            "figsize": (3.5, 3.8),
            "type": "polar",
            "subplots_adjust_h": 0.2,
            "subplots_adjust_w": 0.2
        }

        self.set_plot_dics = []

        plot_dic1 = {
            'position': (1,1), 'title': 'time [ms]', 'cbar': 'right .05 .0',
            'it': 1003520, 'plane': 'xy', 'v_n_x': 'phi_cyl', 'v_n_y': 'r_cyl', 'v_n': 'temperature',
            'vmin': 5., 'vmax': 10., 'rmin': 0., 'rmax': 50., 'mask_below': None, 'mask_above': None,
            'cmap': 'RdBu_r', 'norm': 'log', 'todo': None
        }
        # self.set_plot_dics.append(plot_dic1)

        plot_dic2 = {
            'position': (1,2), 'title': 'time [ms]', 'cbar': None,
            'it': 1024000, 'plane': 'xy', 'v_n_x': 'phi_cyl', 'v_n_y': 'r_cyl', 'v_n': 'temperature',
            'vmin': 5., 'vmax': 10., 'rmin': 0., 'rmax': 50., 'mask_below': None, 'mask_above': None,
            'cmap': 'RdBu_r', 'norm': 'log', 'todo': None
        }
        # self.set_plot_dics.append(plot_dic2)

        plot_dic2 = {
            'position': (1,3), 'title': 'time [ms]', 'cbar': None,
            'it': 1099776, 'plane': 'xy', 'v_n_x': 'phi_cyl', 'v_n_y': 'r_cyl', 'v_n': 'temperature',
            'vmin': 5., 'vmax': 10., 'rmin': 0., 'rmax': 50., 'mask_below': None, 'mask_above': None,
            'cmap': 'RdBu_r', 'norm': 'log', 'todo': None
        }
        # self.set_plot_dics.append(plot_dic2)


    def set_ncols_nrows(self):

        tmp_rows = []
        tmp_cols = []

        for dic in self.set_plot_dics:
            tmp_cols.append(dic['position'][1])
            tmp_rows.append(dic['position'][0])

        max_row = max(tmp_rows)
        max_col = max(tmp_cols)

        for row in range(1, max_row):
            if not row in tmp_rows:
                raise NameError("Please set vertical plot position in a subsequent order: 1,2,3... not 1,3...")

        for col in range(1, max_col):
            if not col in tmp_cols:
                raise NameError("Please set horizontal plot position in a subsequent order: 1,2,3... not 1,3...")

        print("Set {} rows {} columns (total {}) of plots".format(max_row, max_col, len(self.set_plot_dics)))

        return max_row, max_col

    def set_plot_dics_matrix(self):

        plot_dic_matrix = [[0
                             for x in range(self.n_rows)]
                             for y in range(self.n_cols)]

        # get a matrix of dictionaries describing plots (for ease of representation)
        for dic in self.set_plot_dics:
            col, row = int(dic['position'][1]-1), int(dic['position'][0]-1) # -1 as position starts with 1
            # print(col, row)
            for n_row in range(self.n_rows):
                for n_col in range(self.n_cols):
                    if int(col) == int(n_col) and int(row) == int(n_row):
                        plot_dic_matrix[n_col][n_row] = dic
                        # print('adding {} {}'.format(col, row))

            if isinstance(plot_dic_matrix[col][row], int):
                raise ValueError("Dictionary to found for n_row {} n_col {} in "
                                 "creating matrix of dictionaries".format(col, row))

        return plot_dic_matrix

    def set_plot_matrix(self):

        fig = plt.figure(figsize=self.gen_set['figsize'])  # (<->; v)


        if self.gen_set['type'] == 'cartesian':
            # initializing the matrix with dummy axis objects
            sbplot_matrix = [[fig.add_subplot(self.n_rows, self.n_cols, 1)
                                  for x in range(self.n_rows)]
                                  for y in range(self.n_cols)]

            i = 1
            for n_row in range(self.n_rows):
                for n_col in range(self.n_cols):

                    if n_col == 0 and n_row == 0:
                        sbplot_matrix[n_col][n_row] = fig.add_subplot(self.n_rows, self.n_cols, i)
                    elif n_col == 0 and n_row > 0:
                        sbplot_matrix[n_col][n_row] = fig.add_subplot(self.n_rows, self.n_cols, i,
                                                                      sharex=self.sbplot_matrix[n_col][0])
                    elif n_col > 0 and n_row == 0:
                        sbplot_matrix[n_col][n_row] = fig.add_subplot(self.n_rows, self.n_cols, i,
                                                                      sharey=self.sbplot_matrix[0][n_row])
                    else:
                        sbplot_matrix[n_col][n_row] = fig.add_subplot(self.n_rows, self.n_cols, i,
                                                                      sharex=self.sbplot_matrix[n_col][0],
                                                                      sharey=self.sbplot_matrix[0][n_row])

                        # sbplot_matrix[n_col][n_row].axes.get_yaxis().set_visible(False)
                    # sbplot_matrix[n_col][n_row] = fig.add_subplot(n_rows, n_cols, i)
                    i += 1

        elif self.gen_set['type'] == 'polar':
            # initializing the matrix with dummy axis objects
            sbplot_matrix = [[fig.add_subplot(self.n_rows, self.n_cols, 1, projection='polar')
                                  for x in range(self.n_rows)]
                                  for y in range(self.n_cols)]

            i = 1
            for n_row in range(self.n_rows):
                for n_col in range(self.n_cols):

                    if n_col == 0 and n_row == 0:
                        sbplot_matrix[n_col][n_row] = fig.add_subplot(self.n_rows, self.n_cols, i, projection='polar')
                    elif n_col == 0 and n_row > 0:
                        sbplot_matrix[n_col][n_row] = fig.add_subplot(self.n_rows, self.n_cols, i, projection='polar')
                                                                      # sharex=self.sbplot_matrix[n_col][0])
                    elif n_col > 0 and n_row == 0:
                        sbplot_matrix[n_col][n_row] = fig.add_subplot(self.n_rows, self.n_cols, i, projection='polar')
                                                                      # sharey=self.sbplot_matrix[0][n_row])
                    else:
                        sbplot_matrix[n_col][n_row] = fig.add_subplot(self.n_rows, self.n_cols, i, projection='polar')
                                                                      # sharex=self.sbplot_matrix[n_col][0],
                                                                      # sharey=self.sbplot_matrix[0][n_row])

                        # sbplot_matrix[n_col][n_row].axes.get_yaxis().set_visible(False)
                    # sbplot_matrix[n_col][n_row] = fig.add_subplot(n_rows, n_cols, i)
                    i += 1
        else:
            raise NameError("type of the plot is not recognized. Use 'polar' or 'cartesian' ")

        return fig, sbplot_matrix

    @staticmethod
    def plot_one_task(ax, data_cl, task_dic):

        ax.set_rlim(task_dic["rmin"], task_dic["rmax"])

        data = data_cl.get_int(task_dic["it"], task_dic["plane"], task_dic["v_n"])
        phi  = data_cl.get_new_grid(task_dic["plane"], task_dic["v_n_x"])
        r    = data_cl.get_new_grid(task_dic["plane"], task_dic["v_n_y"])

        # set maximum/minimum
        if task_dic["vmin"] == None: task_dic["vmin"] = data.min()
        if task_dic["vmax"] == None: task_dic["vmax"] = data.max()

        # set norm
        if task_dic["norm"] == "norm":
            norm = LogNorm(vmin=task_dic["vmin"], vmax=task_dic["vmax"])
        elif task_dic["norm"] == "log":
            norm = Normalize(vmin=task_dic["vmin"], vmax=task_dic["vmax"])
        else:
            raise NameError("unrecognized norm: {} in task {}".format(task_dic["norm"],
                                                                      task_dic["v_n"]))

        # set masks
        if task_dic["mask_above"] == None and task_dic["mask_below"] == None:
            im2 = ax.pcolormesh(phi, r, data, norm=norm, cmap=task_dic["cmap"])
        elif task_dic["mask_below"] != None and task_dic["mask_above"] == None:
            im2 = ax.pcolormesh(phi, r, np.ma.masked_array(data, data < task_dic["mask_below"]),
                                norm=norm, cmap=task_dic["cmap"])
        elif task_dic["mask_below"] == None and task_dic["mask_above"] != None:
            im2 = ax.pcolormesh(phi, r, np.ma.masked_array(data, data > task_dic["mask_above"]),
                                norm=norm, cmap=task_dic["cmap"])
        else:
            im2 = ax.pcolormesh(phi, r, np.ma.masked_array(data,
                                                   (data > task_dic["mask_above"])
                                                   & (data < task_dic["mask_below"])),
                                                    norm=norm, cmap=task_dic["cmap"])

        im2.set_rasterized(True)
        return im2

    @staticmethod
    def set_plot_title(ax, data_cl, plot_dic):
        if plot_dic["title"] != '' and plot_dic["title"] != None:

            if plot_dic["title"] == 'it':
                title = plot_dic["it"]
            elif plot_dic["title"] == 'time [s]' or \
                plot_dic["title"] == 'time':
                title = "%.3f" % data_cl.get_time(plot_dic["it"]) + " [s]"
            elif plot_dic["title"] == 'time [ms]':
                title = "%.1f" % (data_cl.get_time(plot_dic["it"]) * 1000) + " [ms]"
            else:
                title = plot_dic["title"]
            ax.title.set_text(r'{}'.format(title))

    def plot_images(self):

        # initializing the matrix of images for colorbars (if needed)
        image_matrix = [[0
                        for x in range(self.n_rows)]
                        for y in range(self.n_cols)]

        for n_row in range(self.n_rows):
            for n_col in range(self.n_cols):
                print("Plotting n_row:{} n_col:{}".format(n_row, n_col))
                ax = self.sbplot_matrix[n_col][n_row]
                dic = self.plot_dic_matrix[n_col][n_row]
                if isinstance(dic, int):
                    Printcolor.yellow("Dictionary for row:{} col:{} not set".format(n_row, n_col))
                    self.fig.delaxes(ax) # delets the axis for empty plot
                else:
                    im = self.plot_one_task(ax, self.data, dic)
                    self.set_plot_title(ax, self.data, dic)
                    image_matrix[n_col][n_row] = im

        return image_matrix

    def plot_one_cbar(self, im, dic, n_row, n_col):

        if dic["cbar"] != None and dic["cbar"] != '':

            location = dic["cbar"].split(' ')[0]
            shift_h = float(dic["cbar"].split(' ')[1])
            shift_w = float(dic["cbar"].split(' ')[2])
            cbar_width = 0.02


            if location == 'right':
                ax_to_use = self.sbplot_matrix[-1][n_row]
                pos1 = ax_to_use.get_position()
                pos2 = [pos1.x0 + pos1.width + shift_h,
                        pos1.y0,
                        cbar_width,
                        pos1.height]
            elif location == 'left':
                ax_to_use = self.sbplot_matrix[-1][n_row]
                pos1 = ax_to_use.get_position()
                pos2 = [pos1.x0 - pos1.width - shift_h,
                        pos1.y0,
                        cbar_width,
                        pos1.height]
            elif location == 'bottom':
                ax_to_use = self.sbplot_matrix[n_col][-1]
                pos1 = ax_to_use.get_position()
                pos2 = [pos1.x0,
                        pos1.y0 - pos1.height + shift_w,
                        cbar_width,
                        pos1.height]
            else:
                raise NameError("cbar location {} not recognized. Use 'right' or 'bottom' "
                                .format(location))

            cax1 = self.fig.add_axes(pos2)
            if location == 'right':
                cbar = plt.colorbar(im, cax=cax1, extend='both')#, format='%.1e')
            elif location == 'left':
                cbar = plt.colorbar(im, cax=cax1, extend='both', format='%.1e')
                cax1.yaxis.set_ticks_position('left')
                cax1.yaxis.set_label_position('left')
            else:
                raise NameError("cbar location {} not recognized. Use 'right' or 'bottom' "
                                .format(location))
            cbar.ax.set_title(r"{}".format(str(dic["v_n"]).replace('_', '\_')))


    def plot_colobars(self):

        for n_row in range(self.n_rows):
            for n_col in range(self.n_cols):
                print("Colobar for n_row:{} n_col:{}".format(n_row, n_col))
                # ax  = self.sbplot_matrix[n_col][n_row]
                dic = self.plot_dic_matrix[n_col][n_row]
                im  = self.image_matrix[n_col][n_row]
                if isinstance(dic, int):
                    Printcolor.yellow("Dictionary for row:{} col:{} not set".format(n_row, n_col))
                else:
                    self.plot_one_cbar(im, dic, n_row, n_col)


    def save_plot(self):



        plt.subplots_adjust(hspace=self.gen_set["subplots_adjust_h"])
        plt.subplots_adjust(wspace=self.gen_set["subplots_adjust_w"])
        # plt.tight_layout()
        plt.savefig('{}{}'.format(self.gen_set["figdir"], self.gen_set["figname"]),
                    bbox_inches='tight', dpi=128)
        plt.close()

    def main(self):

        # initializing the n_cols, n_rows
        self.n_rows, self.n_cols = self.set_ncols_nrows()
        # initializing the matrix of dictionaries of the
        self.plot_dic_matrix = self.set_plot_dics_matrix()
        # initializing the axis matrix (for all subplots) and image matrix fo colorbars
        self.fig, self.sbplot_matrix = self.set_plot_matrix()
        # plotting
        self.image_matrix = self.plot_images()
        # adding colobars
        self.plot_colobars()


        # saving the result
        self.save_plot()


# class PLOT_2D:
#
#     def __init__(self, data_cl):
#         self.data_cl = data_cl
#         self.gen_set = {'nlevels': 7,
#                         'file_for_it': 'H.norm2.asc',
#                         'iterations': 0,
#                         'indir': Paths.gw170817 + sim + '/',
#                         'outdir': Paths.ppr_sims + sim + '/2d/',
#                         }
#         self.grid_set = {'type': 'cylindrical', 'n_r': 150, 'n_phi': 150, 'n_z': -100}
#         self.plot_dics = []
#
#         ang_mom = {
#             'it': [254606],
#             'v_n': ['ang_mom'],
#             'vmin': [-5e-4],  # for plotting
#             'vmax': [5e-4],
#             'rmin': [0.],
#             'rmax': [50.],
#             'mask_below': [None],
#             'mask_above': [None],
#             'cmap': ['RdBu_r'],
#             'norm': ['log'],
#             'todo': ['z_int_phi_ave']
#         }
#         self.plot_dics.append(ang_mom)
#
#         ang_mom_flux = {
#             'it': [254606],
#             'v_n': ['ang_mom_flux'],
#             'vmin': [-5e-6],  # for plotting
#             'vmax': [5e-6],
#             'rmin': [0.],
#             'rmax': [50.],
#             'mask_below': [None],
#             'mask_above': [None],
#             'cmap': ['RdBu_r'],
#             'norm': ['log'],
#             'todo': ['z_int_phi_ave']
#         }
#         self.plot_dics.append(ang_mom_flux)
#
#         dens_unb_bern = {
#             'it': [254606],
#             'v_n': ['dens_unb_bern'],
#             'vmin': [1e-7],  # for plotting
#             'vmax': [6e-7],
#             'rmin': [0.],
#             'rmax': [50.],
#             'mask_below': [None],
#             'mask_above': [None],
#             'cmap': ['Greens'],
#             'norm': ['log'],
#             'todo': ['z_int_phi_ave']
#         }
#         self.plot_dics.append(dens_unb_bern)
#
#         temp = {
#             'it': [254606],
#             'v_n': ['temp'],
#             'plane': ['xy'],
#             'vmin': [5.],  # for plotting
#             'vmax': [10.],
#             'rmin': [0.],
#             'rmax': [50.],
#             'mask_below': [None],
#             'mask_above': [None],
#             'cmap': ['Oranges'],
#             'norm': ['log'],
#             'todo': ['slice']
#         }
#         self.plot_dics.append(temp)
#
#         ye = {
#             'it': [254606],
#             'v_n': ['Ye'],
#             'vmin': [0.1],  # for plotting
#             'vmax': [0.4],
#             'rmin': [0.],
#             'rmax': [50.],
#             'mask_below': [None],
#             'mask_above': [None],
#             'cmap': ['Oranges'],
#             'norm': ['norm'],
#             'todo': ['slice']
#         }
#         self.plot_dics.append(ye)
#
#
#     def plot_one_task(self, ax, task_dic, position=0.):
#         ax.set_rlim(task_dic["rmin"], task_dic["rmax"])
#
#         data = self.data_cl.get_int(task_dic["it"], task_dic["plane"], task_dic["v_n"])
#         phi = self.data_cl.get_new_grid("phi_cyl")
#         r = self.data_cl.get_new_grid("r_cyl")
#
#         # set maximum/minimum
#         if task_dic["vmin"] == None: task_dic["vmin"] = data.min()
#         if task_dic["vmax"] == None: task_dic["vmax"] = data.max()
#
#         # set norm
#         if task_dic["norm"] == "norm":
#             norm = LogNorm(vmin=task_dic["vmin"], vmax=task_dic["vmax"])
#         elif task_dic["norm"] == "log":
#             norm = Normalize(vmin=task_dic["vmin"], vmax=task_dic["vmax"])
#         else:
#             raise NameError("unrecognized norm: {} in task {}".format(task_dic["norm"],
#                                                                       task_dic["v_n"]))
#
#         # set masks
#         if task_dic["mask_above"] == None and task_dic["mask_below"] == None:
#             im2 = ax.pcolormesh(phi, r, data, norm=norm, cmap=task_dic["cmap"])
#         elif task_dic["mask_below"] != None and task_dic["mask_above"] == None:
#             im2 = ax.pcolormesh(phi, r, np.ma.masked_array(data, data < task_dic["mask_below"]),
#                                 norm=norm, cmap=task_dic["cmap"])
#         elif task_dic["mask_below"] == None and task_dic["mask_above"] != None:
#             im2 = ax.pcolormesh(phi, r, np.ma.masked_array(data, data > task_dic["mask_above"]),
#                                 norm=norm, cmap=task_dic["cmap"])
#         else:
#             im2 = ax.pcolormesh(data_cl.get_arr("phi_cyl"), data_cl.get_arr("r_cyl"),
#                                 np.ma.masked_array(data,
#                                                    (data > task_dic["mask_above"])
#                                                    & (data < task_dic["mask_below"])),
#                                 norm=norm, cmap=task_dic["cmap"])
#
#         # # if the norm
#         # if task_dic["norm"] == 'norm':
#         #     if task_dic["mask_above"] == None and task_dic["mask_below"] == None:
#         #         im2 = ax.pcolormesh(phi, r, data, norm=norm, cmap=task_dic["cmap"])
#         #     elif task_dic["mask_below"] != None and task_dic["mask_above"] == None:
#         #         im2 = ax.pcolormesh(phi, r, np.ma.masked_array(data, data < task_dic["mask_below"]),
#         #                             norm=norm, cmap=task_dic["cmap"])
#         #     elif task_dic["mask_below"] == None and task_dic["mask_above"] != None:
#         #         im2 = ax.pcolormesh(phi, r, np.ma.masked_array(data, data > task_dic["mask_above"]),
#         #                             norm=norm, cmap=task_dic["cmap"])
#         #     else:
#         #         im2 = ax.pcolormesh(data_cl.get_arr("phi_cyl"), data_cl.get_arr("r_cyl"),
#         #                             np.ma.masked_array(data,
#         #                                                (data > task_dic["mask_above"])
#         #                                                & (data < task_dic["mask_below"])),
#         #                             norm=norm, cmap=task_dic["cmap"])
#         #
#         # elif task_dic["norm"] == 'log':
#         #     if task_dic["mask_above"] == None and task_dic["mask_below"] == None:
#         #         im2 = ax.pcolormesh(phi, r, data, norm=norm, cmap=task_dic["cmap"])
#         #     elif task_dic["mask_below"] != None and task_dic["mask_above"] == None:
#         #         im2 = ax.pcolormesh(phi, r, np.ma.masked_array(data, data < task_dic["mask_below"]),
#         #                             norm=norm, cmap=task_dic["cmap"])
#         #     elif task_dic["mask_below"] == None and task_dic["mask_above"] != None:
#         #         im2 = ax.pcolormesh(phi, r, np.ma.masked_array(data, data > task_dic["mask_above"]),
#         #                             norm=norm, cmap=task_dic["cmap"])
#         #     else:
#         #         im2 = ax.pcolormesh(phi, r,
#         #                             np.ma.masked_array(data,
#         #                                                (data > task_dic["mask_above"])
#         #                                                & (data < task_dic["mask_below"])),
#         #                             norm=norm, cmap=task_dic["cmap"])
#         # else:
#         #     raise NameError("norm:{} is not recognized (task:{})".format(task_dic["norm"],
#         #                                                                  task_dic["v_n"]))
#
#         im2.set_rasterized(True)
#         return im2
#
#     # if task_dic["v_n"] == "ang_mom":
#     #     im = ax.pcolormesh(data_cl.get_slice_z0("phi_cyl"),
#     #                        data_cl.get_slice_z0("r_cyl"),
#     #                        data_cl.get_z_integ__phi_ave(task_dic["v_n"]), # int and ave
#     #                        cmap=task_dic["cmap"],
#     #                        vmin=task_dic["vmin"],
#     #                        vmax=task_dic["vmax"])
#     #     im.set_rasterized(True)
#     #     return im
#     #
#     # elif task_dic["v_n"] == "ang_mom_flux":
#     #     im = ax.pcolormesh(data_cl.get_slice_z0("phi_cyl"),
#     #                        data_cl.get_slice_z0("r_cyl"),
#     #                        data_cl.get_z_integ(task_dic["v_n"]), # int only
#     #                        cmap=task_dic["cmap"],
#     #                        vmin=task_dic["vmin"],
#     #                        vmax=task_dic["vmax"])
#     #     im.set_rasterized(True)
#     #     return im
#     #
#     # elif task_dic["v_n"] in ["dens_unb_bern", "temp", "ye"]:
#     #     data = data_cl.get_z_integ(task_dic["v_n"])
#     #     phi = data_cl.get_slice_z0("phi_cyl")
#     #     r = data_cl.get_slice_z0("r_cyl")
#     #
#     #     if task_dic["vmin"] == None: task_dic["vmin"] = data.min()
#     #     if task_dic["vmax"] == None: task_dic["vmax"] = data.max()
#     #
#     #     if task_dic["norm"] == "norm":
#     #         norm = LogNorm(vmin=task_dic["vmin"], vmax=task_dic["vmax"])
#     #     elif task_dic["norm"] == "log":
#     #         norm = Normalize(vmin=task_dic["vmin"], vmax=task_dic["vmax"])
#     #     else:
#     #         raise NameError("unrecognized norm: {} in task {}".format(task_dic["norm"],
#     #                                                                   task_dic["v_n"]))
#     #     # if the norm
#     #     if task_dic["norm"] == 'norm':
#     #         if task_dic["mask_above"] == None and task_dic["mask_below"] == None:
#     #             im2 = ax.pcolormesh(phi, r, data, norm=norm, cmap=task_dic["cmap"])
#     #         elif task_dic["mask_below"] != None and task_dic["mask_above"] == None:
#     #             im2 = ax.pcolormesh(phi, r, np.ma.masked_array(data, data < task_dic["mask_below"]),
#     #                                 norm=norm, cmap=task_dic["cmap"])
#     #         elif task_dic["mask_below"] == None and task_dic["mask_above"] != None:
#     #             im2 = ax.pcolormesh(phi, r, np.ma.masked_array(data, data > task_dic["mask_above"]),
#     #                                 norm=norm, cmap=task_dic["cmap"])
#     #         else:
#     #             im2 = ax.pcolormesh(data_cl.get_arr("phi_cyl"), data_cl.get_arr("r_cyl"),
#     #                                 np.ma.masked_array(data,
#     #                                                    (data > task_dic["mask_above"])
#     #                                                    & (data < task_dic["mask_below"])),
#     #                                 norm=norm, cmap=task_dic["cmap"])
#     #     elif task_dic["norm"] == 'log':
#     #         if task_dic["mask_above"] == None and task_dic["mask_below"] == None:
#     #             im2 = ax.pcolormesh(phi, r, data, norm=norm, cmap=task_dic["cmap"])
#     #         elif task_dic["mask_below"] != None and task_dic["mask_above"] == None:
#     #             im2 = ax.pcolormesh(phi, r, np.ma.masked_array(data, data < task_dic["mask_below"]),
#     #                                 norm=norm, cmap=task_dic["cmap"])
#     #         elif task_dic["mask_below"] == None and task_dic["mask_above"] != None:
#     #             im2 = ax.pcolormesh(phi, r, np.ma.masked_array(data, data > task_dic["mask_above"]),
#     #                                 norm=norm, cmap=task_dic["cmap"])
#     #         else:
#     #             im2 = ax.pcolormesh(phi, r,
#     #                                 np.ma.masked_array(data, (data > task_dic["mask_above"])
#     #                                                    & (data < task_dic["mask_below"])),
#     #                                 norm=norm, cmap=task_dic["cmap"])
#     #     else:
#     #         raise NameError("norm:{} is not recognized (task:{})".format(task_dic["norm"],
#     #                                                                      task_dic["v_n"]))
#     #
#     #     im2.set_rasterized(True)
#     #     return im2
#     #
#     # else:
#     #     raise NameError("task: {} not recognized".format(task_dic["v_n"]))
#
#
#     def plot_from_dics(self):
#         """
#         1 row per plot
#         :return:
#         """
#
#         # check if all entries in the dictionary have the same length
#         dum_n = len(self.plot_dics[0]['it'])
#         for n_task in range(len(self.plot_dics)):
#             if dum_n != len(self.plot_dics[n_task]['it']):
#                 raise NameError("plot_dic[{}] contains wrong n_of_elements "
#                                 "(should be the same for all)".format(self.plot_dics[n_task]['it']))
#
#         n_rows = len(self.plot_dics)  # |
#         n_cols = len(self.plot_dics[0]['it'])  # <->
#
#         print("Plotting: {} rows {} columns".format(n_rows, n_cols))
#
#         fig = plt.figure(figsize=self.gen_set_dic['figsize'])  # (<->; v)
#         # initializing the matrix with dummy axis objects
#         sbplot_matrix = [[fig.add_subplot(n_rows, n_cols, 1, projection='polar')
#                           for x in range(n_rows)] for y in range(n_cols)]
#
#         # loading data processing classes
#         data_cl = []
#         for it in self.plot_dics[0]['it']:
#             # file type is "12345.cylh5"
#             fname = self.gen_set_dic["path"] + str(int(it)) + self.gen_set_dic["file_extension"]
#             data_cl.append(ANALYZE_CYLINDRICAL(fname))
#
#         i = 1
#         for n_row in range(n_rows):
#             for n_col in range(n_cols):
#
#                 if n_col == 0 and n_row == 0:
#                     sbplot_matrix[n_col][n_row] = fig.add_subplot(n_rows, n_cols, i, projection='polar')
#
#                 elif n_col == 0 and n_row > 0:
#                     sbplot_matrix[n_col][n_row] = fig.add_subplot(n_rows, n_cols, i, projection='polar')  # ,
#                     # sharex=sbplot_matrix[n_col][0])
#
#                 elif n_cols > 0 and n_row == 0:
#                     sbplot_matrix[n_col][n_row] = fig.add_subplot(n_rows, n_cols, i, projection='polar')  # ,
#                     # sharey=sbplot_matrix[0][n_row])
#
#                 else:
#                     sbplot_matrix[n_col][n_row] = fig.add_subplot(n_rows, n_cols, i, projection='polar')  # ,
#                     # sharex=sbplot_matrix[n_col][0],
#                     # sharey=sbplot_matrix[0][n_row])
#
#                     # sbplot_matrix[n_col][n_row].axes.get_yaxis().set_visible(False)
#                 # sbplot_matrix[n_col][n_row] = fig.add_subplot(n_rows, n_cols, i)
#                 i += 1
#
#         # plotting for every point in the grid
#         for n_row in range(n_rows):
#             im = None
#             for n_col in range(n_cols):
#                 ax = sbplot_matrix[n_col][n_row]
#
#                 task_dic = {}
#                 for v_n in self.plot_dics[n_row].keys():
#                     task_dic[v_n] = self.plot_dics[n_row][v_n][n_col]
#                 im = self.plot_one_task(ax, data_cl[n_col], task_dic, n_cols - (n_col + 1))
#
#             # adding axis for colobar and plotting colobars for every row
#             ax_for_row = sbplot_matrix[-1][n_row]
#             pos1 = ax_for_row.get_position()  # get the original position
#             pos2 = [pos1.x0 + pos1.width,
#                     pos1.y0,  # - 0.5 * pos1.height,
#                     0.02, pos1.height]
#             cax1 = fig.add_axes(pos2)
#             cbar = plt.colorbar(im, cax=cax1, extend='both', format='%.1e')
#             cbar.ax.set_title(r"{}".format(str(self.plot_dics[n_row]["v_n"][-1]).replace('_', '\_')))
#
#         plt.subplots_adjust(hspace=0.2)
#         plt.subplots_adjust(wspace=-0.3)
#         # plt.tight_layout()
#         plt.savefig('{}{}'.format(self.gen_set_dic["figir"], self.gen_set_dic["figformat"]),
#                     bbox_inches='tight', dpi=128)
#         plt.close()
#
#     # class Sim_2d_datasets:
#     #     """
#     #     Allows easy access to a scidata datasets of 2d data of a given simulation,
#     #     loading them only if they are needed, and storing them for future accesses.
#     #
#     #     Assumes that the simulation data is stired in /output-xxxx/data/ folders, where
#     #     'xxxx' stands for the number.
#     #     To know, which iterations are stored in what output-xxxx it first loads an ascii
#     #     file <'file_for_it'> which should be present in all output-xxxx and have a
#     #     structure with columns: 1:it 2:time (other columns do not matter)
#     #
#     #     The 2d datasets, located in output-xxxx, expected to be named like:
#     #     'rho.xy.h5'. The list of possible variables, (like 'rho') is set in
#     #     <'self.list_v_ns'>. The list of possible planes (like 'xy') in set in
#     #     <'self.list_planes'>.
#     #
#     #     Logic:
#     #         Every time when the dataset is requested via <'get_dataset(it, plane, v_n)'>,
#     #         the class do:
#     #             1. Checks what output-xxxx contains iteration 'it'
#     #             2. Checks if the dataset for this output, plane and variable name 'v_n'
#     #                 has already been loaded and is present in the storage
#     #                 <'self.dataset_matrix[]'>
#     #                 If so: it will return the dataset from the storage.
#     #                 If not: it will load the required dataset and add it to the storage
#     #                 for future uses.
#     #
#     #     """
#     #
#     #     def __init__(self, sim):
#     #
#     #         self.gen_set = {'nlevels':7,
#     #                         'file_for_it': 'H.norm2.asc',
#     #                         'iterations':0,
#     #                         'indir': Paths.gw170817 + sim + '/',
#     #                         'outdir': Paths.ppr_sims + sim + '/2d/'}
#     #
#     #         self.output_it_map = {}
#     #         self.it_time = self.set_it_output_map()
#     #
#     #         self.list_v_ns = ['rho', 'Y_e', 'temperature', 's_phi', 'entropy', 'dens_unbnd']
#     #         self.list_planes=['xy', 'xz']
#     #
#     #         self.dataset_matrix = [[[0
#     #                                   for z in range(len(self.list_v_ns))]
#     #                                   for k in range(len(self.list_planes))]
#     #                                   for s in range(len(self.output_it_map.keys()))]
#     #
#     #
#     #     def set_it_output_map(self):
#     #         """
#     #         Loads set of files that have '1:it 2:time ...' structure to get a map
#     #         of what output-xxxx contains what iteration (and time)
#     #         """
#     #         print('-' * 25 + 'LOADING it list ({})'
#     #               .format(self.gen_set['file_for_it']) + '-' * 25)
#     #         print(self.gen_set['indir'])
#     #         files = locate(self.gen_set['file_for_it'], root=self.gen_set["indir"], followlinks=True)
#     #         selected = []
#     #         for file in files:
#     #             if file.__contains__('output-'):
#     #                 selected.append(file)
#     #
#     #
#     #         it_time = np.zeros(2)
#     #         for file in selected:
#     #             o_name = file.split('/')
#     #             o_dir = ''
#     #             for o_part in o_name:
#     #                 if o_part.__contains__('output-'):
#     #                     o_dir = o_part
#     #             if o_dir == '':
#     #                 raise NameError("Did not find output-xxxx in {}".format(o_name))
#     #             # print(o_dir)
#     #             # exit(1)
#     #
#     #             it_time_i = np.loadtxt(file, usecols=(0,1))
#     #             # print(it_time_i.shape)
#     #
#     #             # print(file.split(self.gen_set['file_for_it'])[0])
#     #             # exit(1)
#     #
#     #             self.output_it_map[o_dir] = it_time_i
#     #             # print(self.output_it_map[file])
#     #             # for it in it_time_i[:, 0]:
#     #             #     self.output_it_map[it] = file.split(self.gen_set['file_for_it'])[0] # ... ... /sim/output-xxxx/data/
#     #
#     #             it_time = np.vstack(it_time_i)
#     #         # it_time = np.delete(it_time, 0, 0)
#     #         print('outputs:{} iterations:{}'.format(len(selected), len(it_time[:,0])))
#     #         print('-' * 30 + '------DONE-----' + '-' * 30)
#     #         return it_time
#     #
#     #
#     #     def check_v_n(self, v_n):
#     #         if v_n not in self.list_v_ns:
#     #             raise NameError("v_n:{} not in the v_n list (in the class)\n{}".format(v_n, self.list_v_ns))
#     #     def check_plane(self, plane):
#     #         if plane not in self.list_planes:
#     #             raise NameError("plane:{} not in the plane_list (in the class)\n{}".format(plane, self.list_planes))
#     #     # def check_it(self, it):
#     #     #     if it not in self.it_time[:, 0]:
#     #     #         raise NameError("plane:{} not in the plane_list (in the class)\n{}".format(it, self.it_time[:,0]))
#     #     #
#     #     def i_v_n(self, v_n):
#     #         self.check_v_n(v_n)
#     #         return int(self.list_v_ns.index(v_n))
#     #     #
#     #     def i_plane(self, plane):
#     #         self.check_plane(plane)
#     #         return int(self.list_planes.index(plane))
#     #     #
#     #     # # def i_output(self, it):
#     #     # #     for output_data_dir in self.output_it_map.keys():
#     #     # #         if it in np.array(self.output_it_map[output_data_dir])[:, 0]:
#     #     # #             return output_data_dir
#     #     # #
#     #     # #     raise ValueError("it:{} not found in a output_it_map:\n{}".format(it, output_data_dir))
#     #     #
#     #     # def load_h5_dataset(self, v_n='rho', plane='xy'):
#     #     #     """
#     #     #     Usually, v_ns are <rho, Y_e, temperature, s_phi, entropy, dens_unbnd>
#     #     #
#     #     #     :param v_n:
#     #     #     :return:
#     #     #     """
#     #     #
#     #     #     self.check_plane(plane)
#     #     #     self.check_v_n(v_n)
#     #     #
#     #     #     if v_n not in self.list_v_ns:
#     #     #         raise NameError("v_n:{} not in the v_n list (in the class)\n{}".format(v_n, self.list_v_ns))
#     #     #     if plane not in self.list_planes:
#     #     #         raise NameError("plane:{} not in the plane_list (in the class)\n{}".format(plane, self.list_planes))
#     #     #
#     #     #     fname = v_n + '.' + plane + '.h5'
#     #     #
#     #     #     files = locate(fname, root=self.gen_set["indir"], followlinks=True)
#     #     #     dset = h5.dataset(files)
#     #     #
#     #     #     return dset
#     #     #
#     #     #
#     #     # def get_dataset(self, plane, v_n):
#     #     #
#     #     #     if not self.is_loaded(plane, v_n):
#     #     #         self.dataset_matrix[self.i_plane(plane)][self.i_v_n(v_n)] = \
#     #     #             self.load_h5_dataset(plane, v_n)
#     #     #     return self.dataset_matrix[self.i_plane(plane)][self.i_v_n(v_n)]
#     #     #
#     #     #
#     #     # def load_data(self, it, plane, v_n):
#     #     #
#     #     #
#     #     #
#     #     #     files = locate("rho.xy.h5", root=Paths.gw170817 + "LS220_M13641364_M0_SR" + '/output-0029/', followlinks=True)
#     #     #
#     #     #
#     #
#     #
#     #
#     #
#     #
#     #
#     #
#     #
#     #
#     #     def load_dataset(self, o_dir, plane, v_n):
#     #
#     #         fname = v_n + '.' + plane + '.h5'
#     #         # self.gen_set['indir'] + o_dir +'/'
#     #         files = locate("rho.xy.h5", root=self.gen_set['indir'] + o_dir +'/', followlinks=True)
#     #         print("\t Loading: {} plane:{} v_n:{} dataset"
#     #               .format(o_dir, plane, v_n))
#     #         if len(files) > 1:
#     #             raise ValueError("More than 1 file ({}) found. \nFile:{} location:{}"
#     #                              .format(len(files), fname, o_dir))
#     #         if len(files) == 0:
#     #             raise ValueError("NO fils found. \nlocation:{}"
#     #                              .format(fname, o_dir))
#     #         dset = h5.dataset(files)
#     #         self.dataset_matrix[self.i_output(o_dir)][self.i_plane(plane)][self.i_v_n(v_n)] = dset
#     #
#     #
#     #     def i_output(self, o_dir):
#     #
#     #         if o_dir not in self.output_it_map.keys():
#     #             raise NameError("plane:{} not in the plane_list (in the class)\n{}"
#     #                             .format(o_dir, self.output_it_map.keys()))
#     #
#     #         return int(self.output_it_map.keys().index(o_dir))
#     #
#     #     def is_dataset_loaded(self, o_dir, plane, v_n):
#     #
#     #         if isinstance(self.dataset_matrix[self.i_output(o_dir)][self.i_plane(plane)][self.i_v_n(v_n)], int):
#     #             self.load_dataset(o_dir, plane, v_n)
#     #
#     #
#     #
#     #     def it_to_output_dir(self, it):
#     #
#     #         for output_data_dir in self.output_it_map.keys():
#     #             if int(it) in np.array(self.output_it_map[output_data_dir], dtype=int)[:, 0]:
#     #                 return output_data_dir
#     #
#     #         raise ValueError("it:{} not found in a output_it_map:\n{}".format(it, self.output_it_map.keys()))
#     #
#     #     # def is_loaded(self, it, plane, v_n):
#     #     #     data = self.dataset_matrix[self.i_output(it)][self.i_plane(plane)][self.i_v_n(v_n)]
#     #     #     if isinstance(data, int):
#     #     #         self.load_data(it, plane, v_n)
#     #     #     return self.dataset_matrix[self.i_output(it)][self.i_plane(plane)][self.i_v_n(v_n)]
#     #
#     #     def get_dataset(self, it, plane, v_n):
#     #
#     #         o_dir = self.it_to_output_dir(it)
#     #         self.is_dataset_loaded(o_dir, plane, v_n)
#     #
#     #         return self.dataset_matrix[self.i_output(o_dir)][self.i_plane(plane)][self.i_v_n(v_n)]
#     #
#     #         # self.is_loaded(it, plane, v_n)
#     #
#     #         # return
#     #
#     #
#     #     # def get_arr(self, plane, v_n, it, rl):
#     #     #
#     #     #     data = self.data_matrix[self.i_plane(plane)][self.i_v_n(v_n)][self.i_it(it)][rl]
#     #     #
#     #     #     if isinstance(data, int):
#     #     #         self.set_arr(plane, v_n, it, rl)
#     #     #
#     #     #     return data
#     #
#     #
#     # class Sim_2d:
#     #
#     #     def __init__(self, sim):
#     #
#     #         self.gen_set = {'nlevels':7,
#     #                         'file_for_it': 'H.norm2.asc',
#     #                         'iterations':0,
#     #                         'indir': Paths.gw170817 + sim + '/',
#     #                         'outdir': Paths.ppr_sims + sim + '/2d/'}
#     #
#     #
#     #         self.it_output_map = {}
#     #
#     #         self.it_time = self.set_iterations()
#     #         self.gen_set['iterations'] = len(self.it_time[:,0])
#     #
#     #         # self.list_v_ns_planes = ['rho.xy','rho.xz',
#     #         #                         'Y_e.xy', 'Y_e.xz',
#     #         #                         'temperature.xy', 'temperature.xz',
#     #         #                         's_phi.xy', 's_phi.xz',
#     #         #                         'entropy.xy', 'entropy.xz',
#     #         #                         'dens_unbnd.xy', 'dens_unbnd.xz']
#     #
#     #         # contains loaded datasets (so once loaded, they stay here)
#     #         self.dataset_matrix = [[0
#     #                                   for z in range(len(self.list_v_ns))]
#     #                                   for k in range(len(self.list_planes))]
#     #         self.list_v_ns = ['rho', 'Y_e', 'temperature', 's_phi', 'entropy', 'dens_unbnd']
#     #         self.list_planes=['xy', 'xz']
#     #
#     #         # contains actual data, once extracted, -- accessed easily
#     #         self.data_matrix = [[[[0
#     #                                for x in range(self.gen_set["nlevels"])]
#     #                                for y in range(self.gen_set["iterations"])]
#     #                                for z in range(len(self.list_v_ns))]
#     #                                for k in range(len(self.list_planes))]
#     #
#     #         # self.data_matrix = [[[0 for x in range(self.gen_set["nlevels"])]
#     #         #                      for y in range(len(self.gen_set["iterations"]))]
#     #         #                      for z in range(len(self.list_v_ns_planes))]
#     #
#     #
#     #
#     #     def __delete__(self, instance):
#     #         del instance.data_matrix
#     #
#     #     def check_v_n(self, v_n):
#     #         if v_n not in self.list_v_ns:
#     #             raise NameError("v_n:{} not in the v_n list (in the class)\n{}".format(v_n, self.list_v_ns))
#     #     def check_plane(self, plane):
#     #         if plane not in self.list_planes:
#     #             raise NameError("plane:{} not in the plane_list (in the class)\n{}".format(plane, self.list_planes))
#     #     def check_it(self, it):
#     #         if it not in self.it_time[:, 0]:
#     #             raise NameError("plane:{} not in the plane_list (in the class)\n{}".format(it, self.it_time[:,0]))
#     #
#     #     def i_v_n(self, v_n):
#     #         self.check_v_n(v_n)
#     #         return int(self.list_v_ns.index(v_n))
#     #
#     #     def i_plane(self, plane):
#     #         self.check_plane(plane)
#     #         return int(self.list_planes.index(plane))
#     #
#     #     def i_it(self, it):
#     #         self.check_it(it)
#     #         return int(self.it_time[:,0].index(it))
#     #
#     #     def set_iterations(self):
#     #         print('-' * 25 + 'LOADING it list ({})'
#     #               .format(self.gen_set['file_for_it']) + '-' * 25)
#     #         files = locate(self.gen_set['file_for_it'], root=self.gen_set["indir"], followlinks=True)
#     #         it_time = np.zeros(2)
#     #         for file in files:
#     #             it_time_i = np.loadtxt(file,usecols=(1,2))
#     #             it_time = np.vstack(it_time_i)
#     #         it_time = np.delete(it_time, 0, 0)
#     #         print('files:{} iterations:{}'.format(len(files), len(it_time[:,0])))
#     #         print('-' * 30 + '------DONE-----' + '-' * 30)
#     #         return it_time
#     #
#     #     # --------------------------
#     #
#     #     def set_dataset(self, v_n, plane):
#     #         self.check_v_n(v_n)
#     #         self.check_plane(plane)
#     #         if isinstance(self.dataset_matrix[self.i_v_n(v_n)][self.i_plane(plane)], int):
#     #             self.dataset_matrix[self.i_v_n(v_n)][self.i_plane(plane)] = \
#     #                 self.load_h5_dataset(v_n, plane)
#     #
#     #
#     #     def load_h5_dataset(self, v_n='rho', plane='xy'):
#     #         """
#     #         Usually, v_ns are <rho, Y_e, temperature, s_phi, entropy, dens_unbnd>
#     #
#     #         :param v_n:
#     #         :return:
#     #         """
#     #
#     #         self.check_plane(plane)
#     #         self.check_v_n(v_n)
#     #
#     #         if v_n not in self.list_v_ns:
#     #             raise NameError("v_n:{} not in the v_n list (in the class)\n{}".format(v_n, self.list_v_ns))
#     #         if plane not in self.list_planes:
#     #             raise NameError("plane:{} not in the plane_list (in the class)\n{}".format(plane, self.list_planes))
#     #
#     #         fname = v_n + '.' + plane + '.h5'
#     #
#     #         files = locate(fname, root=self.gen_set["indir"], followlinks=True)
#     #         dset = h5.dataset(files)
#     #
#     #
#     #
#     #         return dset
#     #
#     #     def set_arr(self, plane, v_n, it, rl):
#     #
#     #         data = self.data_matrix[self.i_plane(plane)][self.i_v_n(v_n)][self.i_it(it)][rl]
#     #
#     #     # ----------------------
#     #
#     #
#     #     def get_arr(self, plane, v_n, it, rl):
#     #         self.check_v_n(v_n)
#     #         self.check_plane(plane)
#     #         self.check_it(it)
#     #
#     #         data = self.data_matrix[self.i_plane(plane)][self.i_v_n(v_n)][self.i_it(it)][rl]
#     #
#     #         if isinstance(data, int):
#     #             self.set_arr(plane, v_n, it, rl)
#     #
#     #         return data
#     #
#     #
#     # class SIM_2D:
#     #
#     #     def __init__(self, sim):
#     #
#     #         self.gen_set = {'indir': Paths.gw170817 + sim + '/',
#     #                         'outdir': Paths.ppr_sims + sim + '/2d/' }
#     #
#     #         self.list_v_ns = ['rho', 'Y_e', 'temperature', 's_phi', 'entropy', 'dens_unbnd']
#     #         self.list_planes=['xy', 'xz']
#     #
#     #         # data_sets = [np.zeros(0,) for i in range(len(self.list_v_ns))]
#     #         self.data_matrix = [[0 for x in range(len(self.list_v_ns))] for y in range(len(self.data_v_ns))]
#     #
#     #
#     #
#     #     def load_h5_files(self, v_n='rho', plane='xy'):
#     #         """
#     #         Usually, v_ns are <rho, Y_e, temperature, s_phi, entropy, dens_unbnd>
#     #
#     #         :param v_n:
#     #         :return:
#     #         """
#     #         if v_n not in self.list_v_ns:
#     #             raise NameError("v_n:{} not in the v_n list (in the class)\n{}".format(v_n, self.list_v_ns))
#     #         if plane not in self.list_planes:
#     #             raise NameError("plane:{} not in the plane_list (in the class)\n{}".format(plane, self.list_planes))
#     #
#     #         fname = v_n + '.' + plane + '.h5'
#     #         print("\t Loading [{}] files...".format(len(files))),
#     #         files = locate("rho.xy.h5", root=self.gen_set["indir"], followlinks=True)
#     #         self.dset = h5.dataset(files)
#     #         print(" done! (%.2f sec)" % (time.time() - start_t))
#     #         return files
#     #
#     #     def check_v_n(self, v_n):
#     #         if v_n not in self.list_v_ns:
#     #             raise NameError("v_n:{} not in the v_n list (in the class)\n{}".format(v_n, self.list_v_ns))
#     #     def check_plane(self, plane):
#     #         if plane not in self.list_planes:
#     #             raise NameError("plane:{} not in the plane_list (in the class)\n{}".format(plane, self.list_planes))
#     #
#     #     def i_v_n(self, v_n):
#     #         self.check_v_n(v_n)
#     #         return int(self.list_v_ns.index(v_n))
#     #
#     #     def i_plane(self, plane):
#     #         self.check_plane(plane)
#     #         return int(self.list_planes.index(plane))
#     #
#     #     def is_data_loaded(self, v_n, plane):
#     #         if isinstance(self.data_matrix[self.i_v_n(v_n)][self.i_plane(plane)], int):
#     #             self.data_matrix[self.i_v_n(v_n)][self.i_plane(plane)] = self.load_h5_files(v_n, plane)
#     #
#     #     def get_arr(self, v_n, plane):
#     #         self.check_v_n(v_n)
#     #         self.check_plane(plane)

# class COMPUTE_DENSITY_MODES:
#
#     def __init__(self, data_cl):
#
#         self.data = data_cl
#
#
#         # load all data for density (NOTE not a Desity, just rho. as Density = rho * W * vol)
#
#         self.get_set = {
#             'data': '2D', # 2D or 3D.
#             'plane': 'xy',
#             'v_n': 'rho',
#             'modes': [1,2,3,4,5,6],
#             'iterations': 'all',
#             'norm_to_m': 1,
#             'outfname':'rho_modes.h5',
#             'outdir': Paths.ppr_sims + data_cl.sim + '/res_2d/'
#         }
#
#         self.data.load_all(self.get_set['plane'], self.get_set['v_n'])
#
#         self.iterations = self.data.get_all_iterations()
#
#         if isinstance(self.get_set["iterations"], int):
#             self.iterations = self.iterations[::self.get_set["iterations"]]
#
#
#         print('-' * 20 + 'COMPUTING {} modes using {} iterations'
#               .format(len(self.get_set["modes"]), len(self.iterations)) + '-' * 20)
#
#         # initialize datasets
#         self.m1_int_phi = [[np.zeros(0,)
#                              for z in range(len(self.get_set["modes"]))]
#                              for k in range(len(self.iterations))]
#
#         self.m1_int_phi_r = [[np.zeros(0,)
#                              for z in range(len(self.get_set["modes"]))]
#                              for k in range(len(self.iterations))]
#
#         self.out_grid_r = np.zeros(0,)
#
#
#         '''-------------------'''
#         if self.get_set['data'] == '2D':
#             self.compute_2d()
#
#
#     def i_it(self, it):
#         if not int(it) in self.iterations:
#             raise NameError("it:{} not in the list on iterations"
#                             .format(it))
#
#         idx = np.where(self.iterations == int(it))
#         if len(idx) > 1:
#             raise ValueError("For it:{} multiple it are found in the it_time[:,0]".format(it))
#
#         return int(np.where(self.iterations == int(it))[0])
#
#     def i_mode(self, mode):
#         if not int(mode) in self.get_set["modes"]:
#             raise NameError("mode:{} not in the list on iterations"
#                             .format(mode))
#
#         return int(np.where(self.get_set["modes"] == int(mode))[0])
#
#
#
#     def compute_2d(self):
#
#         r_cyl = self.data.get_new_grid(self.get_set['plane'],
#                                        "r_cyl")
#         dr_cyl = self.data.get_new_grid(self.get_set['plane'],
#                                         "dr_cyl")
#
#         for it in self.iterations:
#             density = self.data.get_int(it,
#                                         self.get_set['plane'],
#                                         self.get_set['v_n'])
#             phi_cyl = self.data.get_new_grid(self.get_set['plane'],
#                                          "phi_cyl")
#             for mode in self.get_set["modes"]:
#                 m_int_phi, m_int_phi_r = \
#                     PHYSICS.get_dens_decomp_2d(density,
#                                                phi_cyl, r_cyl, dr_cyl, mode)
#                 if isinstance(self.get_set["norm_to_m"], int):
#                     m_int_phi_norm, m_int_phi_r_norm = \
#                         PHYSICS.get_dens_decomp_2d(density,
#                                                    phi_cyl, r_cyl, dr_cyl,
#                                                    self.get_set["norm_to_m"])
#                     m_int_phi /= m_int_phi_norm
#                     m_int_phi_r /= m_int_phi_r_norm
#                 else:
#                     raise NameError("Wrong dic. setting. Use integer for 'norm_to_m'"
#                                     "Given: {}".format(self.get_set["norm_to_m"]))
#                 self.m1_int_phi[self.i_it(it)][self.i_mode(mode)] = m_int_phi
#                 self.m1_int_phi_r[self.i_it(it)][self.i_mode(mode)] = m_int_phi_r
#
#         self.out_grid_r = r_cyl[:, 0]
#
#
#     def get_mode_int_phi(self, it, mode):
#
#         data = self.m1_int_phi[self.i_it(it)][self.i_mode(mode)]
#         if len(data) > 0:
#             return data
#         else:
#             raise ValueError("mode_int_phi is not computed for it:{} mode:{} ")
#
#         int_phi_reshaped = np.reshape(np.array(self.dataset_matrix[:][self.i_mode(mode)][0]),
#                                       (len(self.iterations), len(self.out_grid_r)))
#
#         return int_phi_reshaped
#
#     def get_mode_int_phi_r(self, mode):
#
#
#     def save_result(self):
#
#         dfile = h5py.File(self.get_set["outdir"] + self.get_set["outfname"], "r")
#
#         dfile.create_dataset("r_cyl", data=self.out_grid_r)
#
#         for i_m, m in enumerate(self.get_set["modes"]):
#
#             int_phi_reshaped = np.reshape(np.array(self.dataset_matrix[:][self.i_mode(m)][0]),
#                                           (len(self.iterations), len(self.out_grid_r)))
#
#             group = dfile.create_group("rl=%d" % m)
#             group["int_phi"] = int_phi_reshaped
#             group["int_phi_r"] = np.array(self.dataset_matrix[:][self.i_mode(m)][1]) # | time   <-> r



# class SIM_2D_old:
#
#     def __init__(self, sim, outdir=''):
#
#         # new grid parameters: polar
#         self.n_phi = 150
#         self.n_r = 140
#
#
#         simdir = Paths.gw170817 + sim + '/'
#
#         if outdir == '':
#             self.outdir = Paths.ppr_sims + sim + '/'
#         else:
#             Printcolor.yellow("Setitng output directory to ./")
#             self.outdir = './'
#
#         files = locate("rho.xy.h5", root=simdir, followlinks=True)
#         start_t = time.time()
#
#         print("\t Loading [{}] files...".format(len(files))),
#         self.dset = h5.dataset(files)
#         print(" done! (%.2f sec)" % (time.time() - start_t))
#
#     @staticmethod
#     def make_stretched_grid(x0, x1, x2, nlin, nlog):
#         assert x1 > 0
#         assert x2 > 0
#         x_lin_f = np.linspace(x0, x1, nlin)
#         x_log_f = 10.0 ** np.linspace(log10(x1), log10(x2), nlog)
#         return np.concatenate((x_lin_f, x_log_f))
#     def get_2d_phi_r_grid(self):
#
#         r_cyl_f = self.make_stretched_grid(0., 15., 512., self.n_r, self.n_phi)
#         phi_cyl_f = np.linspace(0, 2 * np.pi, self.n_phi)
#
#         r_cyl   = 0.5 * (r_cyl_f[1:] + r_cyl_f[:-1])
#         phi_cyl = 0.5 * (phi_cyl_f[1:] + phi_cyl_f[:-1])
#
#         dr_cyl      = np.diff(r_cyl_f)[:, np.newaxis]
#         dphi_cyl    = np.diff(phi_cyl_f)[np.newaxis, :]
#
#         return phi_cyl, r_cyl, dphi_cyl, dr_cyl
#
#     def interp_onto_polar_grid(self, it, variable="HYDROBASE::rho", interp=0):
#
#         from scidata.carpet.interp import Interpolator
#         # assert x_cyl_2d.shape == arr_3d.shape
#
#         phi_cyl, r_cyl, dphi_cyl_2d, dr_cyl_2d = self.get_2d_phi_r_grid()
#
#         r_cyl_2d, phi_cyl_2d = np.meshgrid(r_cyl, phi_cyl, indexing='ij')
#         x_cyl_2d = r_cyl_2d * np.cos(phi_cyl_2d)
#         y_cyl_2d = r_cyl_2d * np.sin(phi_cyl_2d)
#
#         xi = np.column_stack([x_cyl_2d.flatten(), y_cyl_2d.flatten()])
#
#         grid = self.dset.get_grid(iteration=it)
#         rho = self.dset.get_grid_data(grid, iteration=it, variable=variable)
#         irho = Interpolator(grid, rho, interp=interp)
#         res_arr_2d = irho(xi).reshape(x_cyl_2d.shape)
#
#         return res_arr_2d, r_cyl_2d, phi_cyl_2d, dr_cyl_2d, dphi_cyl_2d
#
#     def get_dens_decomp(self, dens_2d, phi_2d, dphi_2d, dr_2d, m=1):
#
#         integ_over_phi = np.sum(dens_2d * np.exp(1j * m * phi_2d) * dphi_2d, axis=1)
#
#         integ_over_phi_r = np.sum(integ_over_phi * dr_2d[:, 0])
#
#         return integ_over_phi, integ_over_phi_r
#
#     def save_densitymode(self):
#
#         m0_integ_over_phi = []
#         m0_integ_over_phi_r = []
#
#         m1_integ_over_phi = []
#         m1_integ_over_phi_r = []
#
#         m2_integ_over_phi = []
#         m2_integ_over_phi_r = []
#
#         r_grid = []
#         times = []
#         for idx, it in enumerate(self.dset.iterations):
#             print("Processing frame {}/{}...".format(idx, len(self.dset.iterations) - 1)),
#             start_t = time.time()
#
#             res_arr_2d, r_cyl_2d, phi_cyl_2d, dr_cyl_2d, dphi_cyl_2d = \
#                 self.interp_onto_polar_grid(it, variable="HYDROBASE::rho", interp=0)
#
#             # --- m0
#
#             m0_integ_over_phi_it, m0_integ_over_phi_r_it = \
#                 self.get_dens_decomp(res_arr_2d, phi_cyl_2d, dphi_cyl_2d, dr_cyl_2d, m=0)
#
#             # append the absolute values (array - for every r) of the complex number
#             m0_integ_over_phi.append(m0_integ_over_phi_it) # np.absolute(
#
#             # appending the absolute value of the final value (summation is done for complex numbers)
#             m0_integ_over_phi_r.append(m0_integ_over_phi_r_it) # np.absolute(
#
#             # --- m1
#
#             m1_integ_over_phi_it, m1_integ_over_phi_r_it = \
#                 self.get_dens_decomp(res_arr_2d, phi_cyl_2d, dphi_cyl_2d, dr_cyl_2d, m=1)
#
#             # append the absolute values (array - for every r) of the complex number
#             m1_integ_over_phi.append(m1_integ_over_phi_it) # np.absolute(
#
#             # appending the absolute value of the final value (summation is done for complex numbers)
#             m1_integ_over_phi_r.append(m1_integ_over_phi_r_it) # np.absolute(
#
#             # --- m2
#
#             m2_integ_over_phi_it, m2_integ_over_phi_r_it = \
#                 self.get_dens_decomp(res_arr_2d, phi_cyl_2d, dphi_cyl_2d, dr_cyl_2d, m=2)
#
#             # append the absolute values (array - for every r) of the complex number
#             m2_integ_over_phi.append(m2_integ_over_phi_it) # np.absolute(
#
#             # appending the absolute value of the final value (summation is done for complex numbers)
#             m2_integ_over_phi_r.append(m2_integ_over_phi_r_it) # np.absolute(
#
#             times.append(self.dset.get_time(it))
#
#             # saving the r_grid
#             if it == self.dset.iterations[-1]:  r_grid = r_cyl_2d[:, 0]
#
#
#
#             del res_arr_2d
#             del r_cyl_2d
#             del phi_cyl_2d
#             del dr_cyl_2d
#             del dphi_cyl_2d
#
#             print("done! (%.2f sec)" % (time.time() - start_t))
#
#         dfile = h5py.File(self.outdir + '/densitymode.h5', "w")
#         # --- m0
#         dfile.create_dataset("m0_integ_over_phi",
#                              data=np.array(m0_integ_over_phi).reshape(len(self.dset.iterations), len(r_grid)))
#         dfile.create_dataset("m0_integ_over_phi_r", data=np.array(m0_integ_over_phi_r))
#         # --- m1
#         dfile.create_dataset("m1_integ_over_phi",
#                              data=np.array(m1_integ_over_phi).reshape(len(self.dset.iterations), len(r_grid)))
#         dfile.create_dataset("m1_integ_over_phi_r", data=np.array(m1_integ_over_phi_r))
#         # --- m2
#         dfile.create_dataset("m2_integ_over_phi",
#                              data=np.array(m2_integ_over_phi).reshape(len(self.dset.iterations), len(r_grid)))
#         dfile.create_dataset("m2_integ_over_phi_r", data=np.array(m2_integ_over_phi_r))
#
#         dfile.create_dataset("r_cyl", data=r_grid)
#         dfile.create_dataset("iterations", data=self.dset.iterations)
#         dfile.create_dataset("time", data=times)
#         dfile.close()

# def plot_dens_modes(sim, outname):
#
#     fname = Paths.ppr_sims + sim + '/densitymode.h5'
#
#     file = h5py.File(fname, 'r')
#     print("Data File: {} [".format(fname)),
#     for i in file:
#         print(i),
#     print(']')
#
#
#     # COMPLEX NUMBERS WARNING
#     # --- m0
#     m0_integ_over_phi = np.array(file["m0_integ_over_phi"])      # 2D
#     m0_integ_over_phi_r = np.array(file["m0_integ_over_phi_r"])  # 1D
#     # --- m1
#     m1_integ_over_phi = np.array(file["m1_integ_over_phi"])      # 2D
#     m1_integ_over_phi_r = np.array(file["m1_integ_over_phi_r"])  # 1D
#     # --- m2
#     m2_integ_over_phi = np.array(file["m2_integ_over_phi"])      # 2D
#     m2_integ_over_phi_r = np.array(file["m2_integ_over_phi_r"])  # 1D
#
#
#     r_cyl = np.array(file["r_cyl"])  # 1D
#     iterations = np.array(file["iterations"])  # 1D
#     time = np.array(file["time"]) * time_constant * 1e-2 # s
#
#     # --- m1
#
#     m1_integ_over_phi_r_phase = np.angle(m1_integ_over_phi_r, deg=True)
#
#     m1_integ_over_phi_r_normed = m1_integ_over_phi_r / m0_integ_over_phi_r
#     m1_integ_over_phi_normed = m1_integ_over_phi / m0_integ_over_phi
#
#     m1_integ_over_phi_r_normed_abs = np.absolute(m1_integ_over_phi_r_normed)
#     m1_integ_over_phi_r_normed_phase = np.angle(m1_integ_over_phi_r_normed, deg=True)
#
#     m1_integ_over_phi_normed_abs = np.absolute(m1_integ_over_phi_normed)
#     m1_integ_over_phi_normed_phase = np.angle(m1_integ_over_phi_normed, deg=True)
#
#     # --- m2
#
#     m2_integ_over_phi_r_normed = m2_integ_over_phi_r / m0_integ_over_phi_r
#     m2_integ_over_phi_r_normed_abs = np.absolute(m2_integ_over_phi_r_normed)
#
#     # fig = plt.figure()
#     # ax = fig.add_subplot(211)
#     # ax.plot(time * 1e2, m1_integ_over_phi_r_normed_phase, '-', color='black', label=r'Phase $C_1 / C_0$')
#     #
#     # ax2 = fig.add_subplot(212)
#     # ax2.plot(time * 1e2, m1_integ_over_phi_r_normed_abs, '-', color='black', label=r'Mag $C_1$')
#     #
#     #
#     # plt.savefig('{}.png'.format('tst_phase'), bbox_inches='tight', dpi=128)
#     # plt.legend()
#     # plt.close()
#
#     # exit(1)
#     # --- plotting
#
#     n_rows = 4
#     n_cols = 1
#
#     fig = plt.figure(figsize=(6.5, 3.5 * 3.6))  # figsize=(4.5, 2.5 * 3.6)  # (<->; v)
#     axs = []
#     for n in range(1, n_rows + 1):
#         if n == 1:
#             axs.append(fig.add_subplot(n_rows, n_cols, n))
#         else:
#             axs.append(fig.add_subplot(n_rows, n_cols, n, sharex=axs[n - 2]))  # sharex=axs[n - 2]))
#
#     i = 0
#
#     divider1 = make_axes_locatable(axs[i])
#     cax = divider1.append_axes("right", size="5%", pad=0.05)
#     cax.xaxis.set_ticklabels([])
#     cax.get_xaxis().set_visible(False)
#     cax.get_yaxis().set_visible(False)
#     # fig.patch.set_visible(False)
#     cax.axis('off')
#     axs[i].plot(time * 1e2, m1_integ_over_phi_r_normed_abs, '-', color='black', label=r'$C_1/C_0$')
#     axs[i].plot(time * 1e2, m2_integ_over_phi_r_normed_abs, '-', color='gray', label=r'$C_2/C_0$')
#     axs[i].set_yscale('log')
#     axs[i].set_ylabel(r'Mag $C_m = \int{\rho\cdot\exp(i m \phi) dr d\phi}$')
#     axs[i].legend()
#
#     i = 1
#
#     divider1 = make_axes_locatable(axs[i])
#     cax = divider1.append_axes("right", size="5%", pad=0.05)
#     cax.xaxis.set_ticklabels([])
#     cax.get_xaxis().set_visible(False)
#     cax.get_yaxis().set_visible(False)
#     # fig.patch.set_visible(False)
#     cax.axis('off')
#     axs[i].plot(time * 1e2, np.unwrap(m1_integ_over_phi_r_normed_phase), '-', color='black', label=r'$C_1/C_0$')
#     # axs[i].plot(time * 1e2, m2_integ_over_phi_r / m0_integ_over_phi_r, '-', color='gray', label=r'$C_2/C_0$')
#     # axs[i].set_yscale('log')
#     axs[i].set_ylabel(r'Phase $C_m = \int{\rho\cdot\exp(i m \phi) dr d\phi}$')
#     axs[i].legend()
#
#     i = 2
#
#
#     from matplotlib import colors
#
#     norm = colors.LogNorm(vmin=1.e-2, vmax=1.)
#     TIME, R = np.meshgrid(time, r_cyl)
#
#
#     im = axs[i].pcolormesh(TIME*1e2, R, m1_integ_over_phi_normed_abs.T, cmap='Blues', norm=norm)#, vmin=0., vmax=0.0002)#, vmin=self.jmin, vmax=self.jmax)
#     im.set_rasterized(True)
#
#     divider1 = make_axes_locatable(axs[i])
#     cax = divider1.append_axes("right", size="5%", pad=0.05)
#     cbar = plt.colorbar(im, cax=cax, norm=norm)
#     cbar.ax.set_title(r"$Mag: C_1(r) / C_0(r)$")
#     cbar.ax.minorticks_off()
#
#     axs[i].set_ylabel(r'$R$')
#     axs[i].set_ylim(0, 50)
#     # axs[i].set_xlabel('time [ms]')
#
#     i = 3
#
#     norm = colors.Normalize(vmin=-180, vmax=180)
#     TIME, R = np.meshgrid(time, r_cyl)
#
#
#     im = axs[i].pcolormesh(TIME*1e2, R, np.unwrap(m1_integ_over_phi_normed_phase.T), cmap='RdBu_r', norm=norm)#, vmin=0., vmax=0.0002)#, vmin=self.jmin, vmax=self.jmax)
#     im.set_rasterized(True)
#
#     divider1 = make_axes_locatable(axs[i])
#     cax = divider1.append_axes("right", size="5%", pad=0.05)
#     cbar = plt.colorbar(im, cax=cax, norm=norm)
#     cbar.ax.set_title(r"$Phase: C_1(r) / C_0(r)$")
#     cbar.ax.minorticks_off()
#
#     axs[i].set_ylabel(r'$R$')
#     axs[i].set_ylim(0, 50)
#     axs[i].set_xlabel('time [ms]')
#
#
#     # plt.savefig('{}.png'.format('tst_densmodes'), bbox_inches='tight', dpi=128)
#     # plt.close()
#     #
#     #
#     #
#     #
#     #
#     #
#     #
#     #
#     #
#     # fig, ax = plt.subplots()
#     # ax.plot(time*1e2, m1_integ_over_phi_r / m0_integ_over_phi_r, '-', color='black', label='m1/m0')
#     # ax.plot(time*1e2, m2_integ_over_phi_r / m0_integ_over_phi_r, '-', color='gray', label='m2/m0')
#     # ax.set_yscale('log')
#     # plt.savefig('{}.png'.format('tst_dens_modes'), bbox_inches='tight', dpi=128)
#     # plt.legend()
#     # plt.close()
#     #
#     #
#     # fig, ax = plt.subplots()
#     # from matplotlib import colors
#     # norm = colors.LogNorm(vmin=1.e-2, vmax=1.)
#     #
#     # TIME, R = np.meshgrid(time, r_cyl)
#     # arr = m1_integ_over_phi / m0_integ_over_phi
#     # im = ax.pcolormesh(TIME*1e2, R, arr.T, cmap='Blues', norm=norm)#, vmin=0., vmax=0.0002)#, vmin=self.jmin, vmax=self.jmax)
#     # im.set_rasterized(True)
#     #
#     # cax1 = fig.add_axes([0.9, 0.1, 0.03, 0.8])
#     # cbar = plt.colorbar(im, cax=cax1)
#     # cbar.ax.set_title(r"$C_1(r) / C_0(r)$")
#     #
#     # ax.set_ylim(0, 50)
#
#     plt.savefig('{}{}.png'.format(Paths.plots, outname), bbox_inches='tight', dpi=128)
#     plt.close()

class COMPUTE_STORE_DESITYMODES:

    def __init__(self, data_cl):

        self.gen_set = {
            'data': '2D', # 2D or 3D.
            'plane': 'xy',
            'v_n': 'rho',
            'v_n_x': 'phi_cyl',
            'v_n_y': 'r_cyl',
            'modes': [1,2,3,4,5,6],
            'iterations': 'all',
            'norm_to_m': 0,
            'save_time': True,
            'outfname':'rho_modes.h5',
            'outdir': Paths.ppr_sims + data_cl.sim + '/res_2d/'
        }

        self.list_v_ns = ["int_phi", "int_phi_r"]

        self.data = data_cl

        self.data.load_all(self.gen_set['plane'], self.gen_set['v_n'])

        self.iterations = self.data.get_all_iterations_times(self.gen_set["plane"],
                                                             self.gen_set["v_n"])
        self.iterations=np.array(self.iterations, dtype=int)
        self.iterations.sort(axis=0)

        if isinstance(self.gen_set["iterations"], int):
            self.iterations = self.iterations[::self.gen_set["iterations"]]

        self.data_matrix = [[[0
                             for k in range(len(self.list_v_ns))]
                             for z in range(len(self.gen_set["modes"]))]
                             for k in range(len(self.iterations))]

        self.gen_set["modes"] = np.array(self.gen_set["modes"], int)

    def check_v_n(self, v_n):
        if v_n not in self.list_v_ns:
            raise NameError("v_n: {} not in the list of v_ns\n{}"
                            .format(v_n, self.list_v_ns))

    def check_mode(self, mode):
        if not int(mode) in self.gen_set["modes"]:
            raise NameError("mode:{} not in the list of modes\n{}"
                            .format(mode, self.gen_set["modes"]))

    def check_it(self, it):
        if not int(it) in self.iterations:
            raise NameError("it:{} not in the list on iterations"
                            .format(it))
        idx = np.where(np.array(self.iterations, dtype=int) == int(it))

        if len(idx) == 0:
            raise ValueError("For it:{} NO it are found in the it_time[:,0]".format(it))
        if len(idx) > 1:
            raise ValueError("For it:{} multiple it are found in the it_time[:,0]".format(it))

    def i_it(self, it):
        self.check_it(it)
        idx = int(np.where(self.iterations == int(it))[0])
        # print(idx)
        return int(idx)

    def i_mode(self, mode):
        self.check_mode(mode)
        return int(np.where(self.gen_set["modes"] == int(mode))[0])

    def i_v_n(self, v_n):
        self.check_v_n(v_n)
        return int(self.list_v_ns.index(v_n))

    def _compute_2d(self, it, mode):

        # getting grid
        r_cyl = self.data.get_new_grid(self.gen_set['plane'],
                                       "r_cyl")
        dr_cyl = self.data.get_new_grid(self.gen_set['plane'],
                                        "dr_cyl")
        phi_cyl = self.data.get_new_grid(self.gen_set['plane'],
                                         "phi_cyl")

        # getting data
        density = self.data.get_int(it,
                                    self.gen_set['plane'],
                                    self.gen_set['v_n'])

        m_int_phi, m_int_phi_r = \
            PHYSICS.get_dens_decomp_2d(density,
                                       phi_cyl, r_cyl, dr_cyl, mode)
        if isinstance(self.gen_set["norm_to_m"], int):
            m_int_phi_norm, m_int_phi_r_norm = \
                PHYSICS.get_dens_decomp_2d(density,
                                           phi_cyl, r_cyl, dr_cyl,
                                           self.gen_set["norm_to_m"])
            m_int_phi /= m_int_phi_norm
            m_int_phi_r /= m_int_phi_r_norm
        else:
            raise NameError("Wrong dic. setting. Use integer for 'norm_to_m'"
                            "Given: {}".format(self.gen_set["norm_to_m"]))


        self.data_matrix[self.i_it(it)][self.i_mode(mode)][self.i_v_n("int_phi")] = \
            m_int_phi
        self.data_matrix[self.i_it(it)][self.i_mode(mode)][self.i_v_n("int_phi_r")] = \
            m_int_phi_r


        # for it in self.iterations:
        #     density = self.data.get_int(it,
        #                                 self.gen_set['plane'],
        #                                 self.gen_set['v_n'])
        #
        #     for mode in self.gen_set["modes"]:
        #         m_int_phi, m_int_phi_r = \
        #             PHYSICS.get_dens_decomp_2d(density,
        #                                        phi_cyl, r_cyl, dr_cyl, mode)
        #         if isinstance(self.gen_set["norm_to_m"], int):
        #             m_int_phi_norm, m_int_phi_r_norm = \
        #                 PHYSICS.get_dens_decomp_2d(density,
        #                                            phi_cyl, r_cyl, dr_cyl,
        #                                            self.gen_set["norm_to_m"])
        #             m_int_phi /= m_int_phi_norm
        #             m_int_phi_r /= m_int_phi_r_norm
        #         else:
        #             raise NameError("Wrong dic. setting. Use integer for 'norm_to_m'"
        #                             "Given: {}".format(self.gen_set["norm_to_m"]))
        #
        #         self.data_matrix[self.i_it(it)][self.i_mode(mode)][self.i_v_n("int_phi")] = \
        #             m_int_phi
        #         self.data_matrix[self.i_it(it)][self.i_mode(mode)][self.i_v_n("int_phi_r")] = \
        #             m_int_phi_r



    def compute_mode(self, it, mode):

        if self.gen_set['data'] == '2D':
            self._compute_2d(it, mode)
        else:
            raise NameError("data format ({}) is not supported yet"
                            .format(self.gen_set['data']))

    def is_computed(self, it, mode, v_n):

        if isinstance(self.data_matrix[self.i_it(it)][self.i_mode(mode)][self.i_v_n(v_n)], int):
            self.compute_mode(it, mode)


    def get_data(self, it, mode, v_n):
        self.check_it(it)
        self.check_mode(mode)
        self.check_v_n(v_n)
        self.is_computed(it, mode, v_n)
        return self.data_matrix[self.i_it(it)][self.i_mode(mode)][self.i_v_n(v_n)]


    def save_all_modes(self):

        print('-' * 20 + 'COMPUTING {} MODES FOR {} ITERATIONS)'
              .format(len(self.gen_set['modes']), len(self.iterations)) + '-' * 20)

        r_cyl = self.data.get_new_grid(self.gen_set['plane'], "r_cyl")
        if self.gen_set['data'] == '2D': r_cyl = r_cyl[:, 0]
        else: raise NameError("data: {} si not supported yet"
                              .format(self.gen_set['data']))
        print("\t saving into {}".format(self.gen_set["outdir"] + self.gen_set["outfname"]))
        dfile = h5py.File(self.gen_set["outdir"] + self.gen_set["outfname"], "w")
        dfile.create_dataset("r_cyl", data=r_cyl)

        for mode in self.gen_set["modes"]:
            int_phi_all = []
            int_phi_r_all = []
            times = []
            iterations = []
            for it in self.iterations:
                try:
                    print("\tComputing it:{}".format(it))
                    int_phi_all.append(self.get_data(it, mode, "int_phi"))
                    int_phi_r_all.append(self.get_data(it, mode, "int_phi_r"))
                    times.append(self.data.get_time(it))
                    iterations.append(it)
                except ValueError:
                    Printcolor.yellow("Warning. it:{} failed with ValueError".format(it))

            int_phi_reshaped = np.reshape(np.array(int_phi_all), (len(iterations), len(r_cyl)))

            group = dfile.create_group("rl=%d" % mode)
            group["int_phi"] = int_phi_reshaped
            group["int_phi_r"] = np.array(int_phi_r_all).flatten()

        dfile.create_dataset("iterations",
                             data=np.array(iterations, dtype=int).flatten())
        if self.gen_set['save_time']:
            dfile.create_dataset("timesteps",
                                 data=np.array(times, dtype=float).flatten())

        print('-' * 30 + '------DONE-----' + '-' * 30)


def plot_density_modes(sim):

    gen_set = {
        # 'data': '2D',  # 2D or 3D.
        # 'plane': 'xy',
        # 'v_n': 'rho',
        # 'v_n_x': 'phi_cyl',
        # 'v_n_y': 'r_cyl',
        'modes': [1, 2, 3, 4, 5, 6],
        # 'iterations': 'all',
        # 'norm_to_m': 1,
        # 'save_time': True,
        'infname': 'rho_modes.h5',
        'outfname': 'rho_modes.png',
        'outdir': Paths.ppr_sims + sim + '/res_2d/'
    }

    dfile = h5py.File(gen_set['outdir'] + gen_set['infname'], "r")

    print('Content: [ ',)
    for var in dfile:
        print(var),
    print(' ]')

    rows = 2
    cols = 1

    fig = plt.figure(figsize=(6.5, 1.5 * 3.6))  # figsize=(4.5, 2.5 * 3.6)  # (<->; v)
    ax_list = []
    for n in range(1, rows + 1):
        if n == 1:
            ax_list.append(fig.add_subplot(rows, cols, n))
        else:
            ax_list.append(fig.add_subplot(rows, cols, n, sharex=ax_list[n - 2]))  # sharex=axs[n - 2]))

    times = np.array(dfile["timesteps"])
    r = np.array(dfile["r_cyl"])

    for m in gen_set['modes']:

        if m == 1: color = 'black'; ls = '-'; lw=1.
        elif m == 2: color = 'gray';  ls = '-.'; lw=1.
        elif m == 3: color = 'blue';  ls = '-.'; lw=0.4
        elif m == 4: color = 'orange';  ls = '-.'; lw=0.4
        elif m == 5: color = 'green';  ls = '-.'; lw=0.4
        elif m == 6: color = 'pink';  ls = '-.'; lw=0.4
        elif m == 7: color = 'purple';  ls = '-.'; lw=0.4
        else: raise ValueError('m is not in color/ls list')

        group = dfile["rl=%d" % m]
        int_phi2d = np.array(group["int_phi"])
        int_phi_r1d = np.array(group["int_phi_r"])  # | time   <-> r


        # phase plot
        ax_list[0].plot(times * 1e3, np.unwrap(np.angle(int_phi_r1d)), ls, lw=lw, color=color, label='m:{}'.format(m))

        ax_list[0].set_ylabel(r'$C_m/C_0$ Phase')
        # ax_list[0].annotate(r'$C_m = \int{\rho W \sqrt{-g}\cdot\exp(i m \phi) dz dr d\phi}$',
        #             xy=(-150 / 180 * np.pi, 50),  # theta, radius
        #             xytext=(0.65, 0.90),  # fraction, fraction
        #             xycoords='axes fraction',
        #             horizontalalignment='center',
        #             verticalalignment='center'
        #             )
        ax_list[0].legend()

        # magnitude plot
        ax_list[1].plot(times * 1e3, np.abs(int_phi_r1d), ls, lw=lw,  color=color, label='m:{}'.format(m))

        ax_list[1].set_yscale('log')
        ax_list[1].set_ylabel(r'$C_m/C_0$ Magnitude')
        ax_list[1].legend()

    plt.savefig(gen_set['outdir'] + gen_set['outfname'], bbox_inches='tight', dpi=128)
    plt.close()



def task_movie1():

    int_ = INTERPOLATE_STORE_2D("LS220_M13641364_M0_SR")
    pl_ = PLOT_MANY(int_)

    pl_.gen_set = {
        "figdir": int_.gen_set["outdir"] + 'movie/',
        "figname": "tst.png",
        "figsize": (3.5, 3.8),
        "type": "polar",
        "subplots_adjust_h": 0.2,
        "subplots_adjust_w": 0.2
    }

    plot_dic1 = {
        'position': (1, 1), 'title': 'time [ms]', 'cbar': 'right .05 .0',
        'it': 1003520, 'plane': 'xy', 'v_n_x': 'phi_cyl', 'v_n_y': 'r_cyl', 'v_n': 'temperature',
        'vmin': 5., 'vmax': 10., 'rmin': 0., 'rmax': 50., 'mask_below': None, 'mask_above': None,
        'cmap': 'RdBu_r', 'norm': 'log', 'todo': None
    }
    pl_.set_plot_dics.append(plot_dic1)

    int_.load_all(plot_dic1['plane'], plot_dic1['v_n'])

    int_.iterations = int_.get_all_iterations_times(plot_dic1["plane"], plot_dic1["v_n"])
    int_.iterations = np.array(int_.iterations, dtype=int)
    int_.iterations.sort(axis=0)

    for it in int_.iterations:
        try:
            pl_.gen_set["figname"] = "{0:07d}.png".format(int(it)) # 7 digit output
            pl_.set_plot_dics[0]["it"] = it
            pl_.main()
        except ValueError:
            Printcolor.yellow("Warning. it:{} failed with ValueError".format(it))

def task_movie2():

    int_ = INTERPOLATE_STORE_2D("LS220_M13641364_M0_SR")
    pl_ = PLOT_MANY(int_)

    pl_.gen_set = {
        "figdir": int_.gen_set["outdir"] + 'movie/',
        "figname": "_tst.png",
        "figsize": (3.5, 3.8),
        "type": "polar",
        "subplots_adjust_h": 0.2,
        "subplots_adjust_w": 0.2
    }

    # RdBu_r
    plot_dic1 = {
        'position': (1, 1), 'title': 'time [ms]', 'cbar': 'left .15 .0',
        'it': 1003520, 'plane': 'xy', 'v_n_x': 'phi_cyl', 'v_n_y': 'r_cyl', 'v_n': 'rho',
        'vmin': 1e-6, 'vmax': 5e-6, 'rmin': 0., 'rmax': 50., 'mask_below': None, 'mask_above': None,
        'cmap': 'inferno', 'norm': 'log', 'todo': None
    }
    pl_.set_plot_dics.append(plot_dic1)

    plot_dic2 = {
        'position': (1, 2), 'title': 'it', 'cbar': 'right .05 .0',
        'it': 1003520, 'plane': 'xy', 'v_n_x': 'phi_cyl', 'v_n_y': 'r_cyl', 'v_n': 'Y_e',
        'vmin': 0.05, 'vmax': 0.45, 'rmin': 0., 'rmax': 50., 'mask_below': None, 'mask_above': None,
        'cmap': 'inferno', 'norm': 'log', 'todo': None
    }
    pl_.set_plot_dics.append(plot_dic2)

    plot_dic2 = {
        'position': (2, 1), 'title': None, 'cbar': 'left .15 .0',
        'it': 1003520, 'plane': 'xy', 'v_n_x': 'phi_cyl', 'v_n_y': 'r_cyl', 'v_n': 'temperature',
        'vmin': 5., 'vmax': 10., 'rmin': 0., 'rmax': 50., 'mask_below': None, 'mask_above': None,
        'cmap': 'inferno', 'norm': 'log', 'todo': None
    }
    pl_.set_plot_dics.append(plot_dic2)

    plot_dic2 = {
        'position': (2, 2), 'title': None, 'cbar': 'right .05 .0',
        'it': 1003520, 'plane': 'xy', 'v_n_x': 'phi_cyl', 'v_n_y': 'r_cyl', 'v_n': 'entropy',
        'vmin': 0, 'vmax': 30, 'rmin': 0., 'rmax': 50., 'mask_below': None, 'mask_above': None,
        'cmap': 'inferno', 'norm': 'log', 'todo': None
    }
    pl_.set_plot_dics.append(plot_dic2)

    int_.load_all(plot_dic1['plane'], plot_dic1['v_n'])

    int_.iterations = int_.get_all_iterations_times(plot_dic1["plane"], plot_dic1["v_n"])
    int_.iterations = np.array(int_.iterations, dtype=int)
    int_.iterations.sort(axis=0)

    for it in int_.iterations:
        try:
            pl_.gen_set["figname"] = "{0:07d}.png".format(int(it)) # 7 digit output
            for dic in pl_.set_plot_dics:
                dic["it"] = it
            pl_.main()
        except ValueError:
            Printcolor.yellow("Warning. it:{} failed with ValueError".format(it))


def task_plot_many(sim):

    int_ = INTERPOLATE_STORE_2D(sim)
    pl_ = PLOT_MANY(int_)

    pl_.gen_set = {
        "figdir": int_.gen_set["outdir"], # ppr_sim + /sim/ + res_2d/
        "figname": "temps4.png",
        "figsize": (8.5, 2.5), # <-> | # for 4 in the row
        # "figsize": (5.5, 5.5),  # <-> |
        "type": "polar",
        "subplots_adjust_h": 1.2,
        "subplots_adjust_w": 0.2
    }
    # RdBu_r
    plot_dic1 = {
        'position': (1, 1), 'title': 'time [ms]', 'cbar': 'left .15 .0',
        'it': 700416, 'plane': 'xy', 'v_n_x': 'phi_cyl', 'v_n_y': 'r_cyl', 'v_n': 'rho',
        'vmin': 1e-6, 'vmax': 5e-6, 'rmin': 0., 'rmax': 50., 'mask_below': None, 'mask_above': None,
        'cmap': 'inferno', 'norm': 'log', 'todo': None
    }
    # pl_.set_plot_dics.append(plot_dic1)

    plot_dic2 = {
        'position': (1, 2), 'title': 'it', 'cbar': 'right .05 .0',
        'it': 700416, 'plane': 'xy', 'v_n_x': 'phi_cyl', 'v_n_y': 'r_cyl', 'v_n': 'Y_e',
        'vmin': 0.05, 'vmax': 0.45, 'rmin': 0., 'rmax': 50., 'mask_below': None, 'mask_above': None,
        'cmap': 'inferno', 'norm': 'log', 'todo': None
    }
    # pl_.set_plot_dics.append(plot_dic2)

    plot_dic2 = {
        'position': (2, 1), 'title': None, 'cbar': 'left .15 .0',
        'it': 700416, 'plane': 'xy', 'v_n_x': 'phi_cyl', 'v_n_y': 'r_cyl', 'v_n': 'temperature',
        'vmin': 5., 'vmax': 10., 'rmin': 0., 'rmax': 50., 'mask_below': None, 'mask_above': None,
        'cmap': 'inferno', 'norm': 'log', 'todo': None
    }
    # pl_.set_plot_dics.append(plot_dic2)

    plot_dic2 = {
        'position': (2, 2), 'title': None, 'cbar': 'right .05 .0',
        'it': 700416, 'plane': 'xy', 'v_n_x': 'phi_cyl', 'v_n_y': 'r_cyl', 'v_n': 'entropy',
        'vmin': 0, 'vmax': 30, 'rmin': 0., 'rmax': 50., 'mask_below': None, 'mask_above': None,
        'cmap': 'inferno', 'norm': 'log', 'todo': None
    }
    # pl_.set_plot_dics.append(plot_dic2)

    '''4 temperatures'''
    plot_dic2 = {
        'position': (1, 1), 'title': 'time [ms]', 'cbar': None,#'left .15 .0',
        'it': 432128, 'plane': 'xy', 'v_n_x': 'phi_cyl', 'v_n_y': 'r_cyl', 'v_n': 'temperature',
        'vmin': 5., 'vmax': 10., 'rmin': 0., 'rmax': 50., 'mask_below': None, 'mask_above': None,
        'cmap': 'inferno', 'norm': 'log', 'todo': None
    }
    pl_.set_plot_dics.append(plot_dic2)
    plot_dic2 = {
        'position': (1, 2), 'title': 'time [ms]', 'cbar': None,#'left .15 .0',
        'it': 649216, 'plane': 'xy', 'v_n_x': 'phi_cyl', 'v_n_y': 'r_cyl', 'v_n': 'temperature',
        'vmin': 5., 'vmax': 10., 'rmin': 0., 'rmax': 50., 'mask_below': None, 'mask_above': None,
        'cmap': 'inferno', 'norm': 'log', 'todo': None
    }
    pl_.set_plot_dics.append(plot_dic2)
    plot_dic2 = {
        'position': (1, 3), 'title': 'time [ms]', 'cbar': None,#'left .15 .0',
        'it': 866304, 'plane': 'xy', 'v_n_x': 'phi_cyl', 'v_n_y': 'r_cyl', 'v_n': 'temperature',
        'vmin': 5., 'vmax': 10., 'rmin': 0., 'rmax': 50., 'mask_below': None, 'mask_above': None,
        'cmap': 'inferno', 'norm': 'log', 'todo': None
    }
    pl_.set_plot_dics.append(plot_dic2)
    plot_dic2 = {
        'position': (1, 4), 'title': 'time [ms]', 'cbar': 'right .05 .0',
        'it': 1083392, 'plane': 'xy', 'v_n_x': 'phi_cyl', 'v_n_y': 'r_cyl', 'v_n': 'temperature',
        'vmin': 5., 'vmax': 10., 'rmin': 0., 'rmax': 50., 'mask_below': None, 'mask_above': None,
        'cmap': 'inferno', 'norm': 'log', 'todo': None
    }
    pl_.set_plot_dics.append(plot_dic2)



    pl_.main()

def get_nearest_iterations_for_times(sim, required_times, file_for_it="dens.norm1.asc"):

    _, it_time = set_it_output_map(Paths.gw170817 + sim + '/', file_for_it)

    it_time[:, 1] *= 0.004925794970773136 * 1e-3  # time is seconds

    print("\n")
    print("    < tmin: {:.3f} tmax: {:.3f} ({}) >".format(it_time[:, 1].min(), it_time[:, 1].max(),
                                                          file_for_it))
    print("    --------------- TASK -------------------")
    print("    t_req  |  t_aval  |  it       ")
    for required_time in required_times:

        if required_time > it_time[:, 1].max():
            raise ValueError("time {:.3f}s beyond the simulation length ({:.3f}s)".format(required_time, it_time[:, 1].max()))
        if required_time < it_time[:, 1].min():
            raise ValueError("time {:.3f}s is too small, minimum is ({:.3f}s)".format(required_time, it_time[:, 1].min()))

        closest_iteration = int(it_time[find_nearest_index(it_time[:, 1], required_time), 0])
        closest_time = it_time[find_nearest_index(it_time[:, 1], required_time), 1]

        print("    {:.3f}  |  {:.3f}   |  {}  ".format(required_time, closest_time, closest_iteration))
    print("    --------------- DONE -------------------")
    print("\n")


if __name__ == '__main__':

    ''' --- DEBUGGING --- '''
    sim = "DD2_M15091235_M0_LK_HR"
    it = 366592

    rl = 3
    v_n = "rho"
    o_rl = EXTRACT_FOR_RL(sim)


    from matplotlib import colors
    fig = plt.figure(figsize=(3., 6.))


    plane = "xy"
    x = o_rl.get_grid_v_n_rl(it, plane, rl, "x")
    y = o_rl.get_grid_v_n_rl(it, plane, rl, "y")
    data_xy = o_rl.get_data_rl(it, plane, rl, v_n)

    ax = fig.add_subplot(212)
    norm = colors.LogNorm(vmin=data_xy.min(), vmax=data_xy.max())
    ax.pcolormesh(x, y, data_xy, norm=norm, cmap="inferno_r")

    ''' --- '''

    plane = "xz"
    x = o_rl.get_grid_v_n_rl(it, plane, rl, "x")
    z = o_rl.get_grid_v_n_rl(it, plane, rl, "z")
    data_xz = o_rl.get_data_rl(it, plane, rl, v_n)

    ax = fig.add_subplot(211)
    norm = colors.LogNorm(vmin=data_xz.min(), vmax=data_xz.max())
    ax.pcolormesh(x, z, data_xz, norm=norm, cmap="inferno_r")


    print("x shape {}".format(x.shape))
    print("y shape {}".format(y.shape))
    print("z shape {}".format(z.shape))
    print("data_xy shape {}".format(data_xy.shape))
    print("data_xz shape {}".format(data_xz.shape))

    # plt.ylim(0, 50)
    # plt.ylim(0, 50)
    plt.title(r"$xz$-plane")
    plt.savefig('{}'.format(Paths.ppr_sims + sim + '/res_2d/' + "{}_{}.png"
                            .format(it, v_n)), bbox_inches='tight', dpi=128)
    print("saved pi_test2.png")
    plt.close()

    exit(0)



    sim = "LS220_M13641364_M0_SR"

    ''' --- ITERATION FINDER --- '''
    # finds interations for given times (s) in the default files, like dens.norm1.asc
    get_nearest_iterations_for_times(sim, [0.020, 0.030, 0.040, 0.050])
    # exit(1)

    ''' --- MOVIE --- '''
    # task_movie1() # for 1 plot of temperature
    # task_movie2() # 4 plots with rho, temp, ye, entropy

    ''' --- PLOT XY PROJECIONS OF MULTIPLE ITERATIONS --- '''
    task_plot_many(sim)

    ''' --- COMPUTE rho MODES --- '''
    # int_ = INTERPOLATE_STORE('DD2_M13641364_M0_SR')
    # dm = COMPUTE_STORE_DESITYMODES(int_); dm.save_all_modes()

    ''' --- PLOT rho MODES --- '''
    # plot_density_modes('DD2_M13641364_M0_SR')

    # ''' TESTING NEW '''
    # data = EXTRACT_STORE_DATA("LS220_M13641364_M0_SR")
    # print(len(data.get_data(1013760, 'xy', 'temperature')))
    # exit(1)

    # int_ = INTERPOLATE_STORE("LS220_M13641364_M0_SR")
    # dm = COMPUTE_STORE_DESITYMODES(int_)
    # dm.save_all_modes()
    # exit(1)
    # int_ = INTERPOLATE_STORE("LS220_M13641364_M0_SR")
    # int_.load_all('xy', 'temperature')
    # print("done")
    # int_.get_all_iterations('xy', 'temperature')
    # exit(1)


    # int_.get_data(1099776, 'xy', 'temperature')
    # pl.get_plot_matrix()
    # exit(1)

    # it=1014272
    # files = locate("rho.xy.h5", root=Paths.gw170817 + "LS220_M13641364_M0_SR" + '/output-0029/', followlinks=True)
    # dset = h5.dataset(files)
    # print(dset.metadata[0])
    # print(dset.iterations)
    #
    # print('___---___')
    #
    # it = 1015808
    # print(dset.get_time(it))

    # grid = dset.get_grid(iteration=it)
    # print(grid.coordinates()[0]) ## gives 2d (or 3d) arrays of grid for each ref. level
    # print(grid.dim)
    # print(grid.mesh()[0]) ## gives 2d (or 3d) arrays of grid for each ref. level
    # rho = dset.get_grid_data(grid, iteration=it, variable="HYDROBASE::rho")

    # print(len(rho))
    # print(len(rho[0]))
    # print(rho[0])


    # d22 = LOAD_STORE_DATASETS("LS220_M13641364_M0_SR")
    # print(d22.get_dataset(1052672, 'xy', 'rho').metadata[0])
    # print(d22.get_dataset(1054720, 'xy', 'temperature').metadata[0])

    # ex = EXTRACT_STORE_DATA("LS220_M13641364_M0_SR")
    # print(ex.get_data(1052672, 'xy', 'rho'))

    # int_ = INTERPOLATE_STORE("LS220_M13641364_M0_SR")
    # int_.tmp()
    # print(int_.get_int(1052672, 'xy', 'rho'))
    # print(int_.get_int(1052672, 'xz', 'temperature'))
    print("Exit successful")
    # d22.load_h5_files('rho', 'xy')

    # exit(1)
    #
    # ''' DENSITY MODES D2 '''
    # d2 = SIM_2D_old('DD2_M13641364_M0_SR')
    # d2.save_densitymode()  # save the density decomposition