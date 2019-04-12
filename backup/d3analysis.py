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
from glob import glob
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
import gc

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

""" --- --- OTHERS --- --- """

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

    def __init__(self, grid_info, carpet_grid):

        self.grid_info = grid_info

        self.grid_type = grid_info["type"]

        self.carpet_grid = carpet_grid

        self.grid_v_ns = ["x_cyl", "y_cyl", "z_cyl",
                          "r_cyl", "phi_cyl",
                          "dr_cyl", "dphi_cyl", "dz_cyl"]

        print('-' * 25 + 'INITIALIZING CYLINDRICAL GRID' + '-' * 25)

        phi_cyl, r_cyl, z_cyl, \
        self.dphi_cyl_3d, self.dr_cyl_3d, self.dz_cyl_3d = self.get_phi_r_z_grid()

        self.r_cyl_3d, self.phi_cyl_3d, self.z_cyl_3d \
            = np.meshgrid(r_cyl, phi_cyl, z_cyl, indexing='ij')
        self.x_cyl_3d = self.r_cyl_3d * np.cos(self.phi_cyl_3d)
        self.y_cyl_3d = self.r_cyl_3d * np.sin(self.phi_cyl_3d)

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

    # generic methods to be present in all INTERPOLATION CLASSES
    def get_int_arr(self, arr_3d):

        # if not self.x_cyl_3d.shape == arr_3d.shape:
        #     raise ValueError("Passed for interpolation 3d array has wrong shape:\n"
        #                      "{} Expected {}".format(arr_3d.shape, self.x_cyl_3d.shape))
        xi = np.column_stack([self.x_cyl_3d.flatten(),
                              self.y_cyl_3d.flatten(),
                              self.z_cyl_3d.flatten()])
        F = Interpolator(self.carpet_grid, arr_3d, interp=1)
        res_arr_3d = F(xi).reshape(self.x_cyl_3d.shape)
        return res_arr_3d

    def get_int_grid(self, v_n):
        if v_n == "x_cyl":
            return self.x_cyl_3d
        elif v_n == "y_cyl":
            return self.y_cyl_3d
        elif v_n == "z_cyl":
            return self.z_cyl_3d
        elif v_n == "r_cyl":
            return self.r_cyl_3d
        elif v_n == "phi_cyl":
            return self.phi_cyl_3d
        elif v_n == "dr_cyl":
            return self.dr_cyl_3d
        elif v_n == "dphi_cyl":
            return self.dphi_cyl_3d
        elif v_n == "dz_cyl":
            return self.dz_cyl_3d
        else:
            raise NameError("v_n: {} not recogized in grid. Available:{}"
                            .format(v_n, self.grid_v_ns))

class FORMULAS:

    def __init__(self):
        pass

    @staticmethod
    def density(rho, w_lorentz, vol):
        return rho * w_lorentz * vol

    @staticmethod
    def vup(velx, vely, velz):
        return [velx, vely, velz]

    @staticmethod
    def metric(gxx, gxy, gxz, gyy, gyz, gzz):
        return [[gxx, gxy, gxz], [gxy, gyy, gyz], [gxz, gyz, gzz]]

    @staticmethod
    def enthalpy(eps, press, rho):
        return 1 + eps + (press / rho)

    @staticmethod
    def shift(betax, betay, betaz):
        return [betax, betay, betaz]

    @staticmethod
    def shvel(shift, vlow):
        shvel = np.zeros(shift[0].shape)
        for i in range(len(shift)):
            shvel += shift[i] * vlow[i]
        return shvel

    @staticmethod
    def u_0(w_lorentz, shvel, lapse):
        return w_lorentz * (shvel - lapse)

    @staticmethod
    def vlow(metric, vup):
        vlow = [np.zeros_like(vv) for vv in [vup[0], vup[1], vup[2]]]
        for i in range(3):  # for x, y
            for j in range(3):
                vlow[i][:] += metric[i][j][:] * vup[j][:]  # v_i = g_ij * v^j (lowering index) for x y
        return vlow

    @staticmethod
    def vphi(x, y, vlow):
        return -y * vlow[0] + x * vlow[1]

    @staticmethod
    def vr(x, y, vup):
        r = np.sqrt(x ** 2 + y ** 2)
        return (x / r) * vup[0] + (y / r) * vup[1]

    @staticmethod
    def theta(x, y, z):
        r = np.sqrt(x ** 2 + y ** 2 + z ** 2)
        return np.arccos(z/r)

    @staticmethod
    def dens_unb_geo(u_0, rho, w_lorentz, vol):

        c_geo = -u_0 - 1.0
        mask_geo = (c_geo > 0).astype(int)  # 1 or 0
        rho_unbnd_geo = rho * mask_geo  # if 1 -> same, if 0 -> masked
        dens_unbnd_geo = rho_unbnd_geo * w_lorentz * vol

        return dens_unbnd_geo

    @staticmethod
    def dens_unb_bern(enthalpy, u_0, rho, w_lorentz, vol):

        c_ber = -enthalpy * u_0 - 1.0
        mask_ber = (c_ber > 0).astype(int)
        rho_unbnd_bernoulli = rho * mask_ber
        density_unbnd_bernoulli = rho_unbnd_bernoulli * w_lorentz * vol

        return density_unbnd_bernoulli

    @staticmethod
    def dens_unb_garch(enthalpy, u_0, lapse, press, rho, w_lorentz, vol):

        c_ber = -enthalpy * u_0 - 1.0
        c_gar = c_ber - (lapse / w_lorentz) * (press / rho)
        mask_gar = (c_gar > 0).astype(int)
        rho_unbnd_garching = rho * mask_gar
        density_unbnd_garching = rho_unbnd_garching * w_lorentz * vol

        return density_unbnd_garching

    @staticmethod
    def ang_mom(rho, eps, press, w_lorentz, vol, vphi):
        return (rho * (1 + eps) + press) * w_lorentz * w_lorentz * vol * vphi

    @staticmethod
    def ang_mom_flux(ang_mom, lapse, vr):
        return ang_mom * lapse * vr

    # data manipulation methods
    @staticmethod
    def get_slice(x3d, y3d, z3d, data3d, slice='xy'):

        if slice == 'yz':
            ix0 = np.argmin(np.abs(x3d[:, 0, 0]))
            if abs(x3d[ix0, 0, 0]) < 1e-15:
                res = data3d[ix0, :, :]
            else:
                if x3d[ix0, 0, 0] > 0:
                    ix0 -= 1
                res = 0.5 * (data3d[ix0, :, :] + data3d[ix0 + 1, :, :])
        elif slice == 'xz':
            iy0 = np.argmin(np.abs(y3d[0, :, 0]))
            if abs(y3d[0, iy0, 0]) < 1e-15:
                res = data3d[:, iy0, :]
            else:
                if y3d[0, iy0, 0] > 0:
                    iy0 -= 1
                res = 0.5 * (data3d[:, iy0, :] + data3d[:, iy0 + 1, :])
        elif slice == 'xy':
            iz0 = np.argmin(np.abs(z3d[0, 0, :]))
            if abs(z3d[0, 0, iz0]) < 1e-15:
                res = data3d[:, :, iz0]
            else:
                if z3d[0, 0, iz0] > 0 and iz0 > 0:
                    iz0 -= 1
                res = 0.5 * (data3d[:, :, iz0] + data3d[:, :, iz0 + 1])
        else:
            raise ValueError("slice:{} not recognized. Use 'xy', 'yz' or 'xz' to get a slice")
        return res

""" --- --- INPUT DATA IS IN .dat.tar FORMAT --- --- """

class LOAD_STORE_DATASET:

    def __init__(self, sim):
        self.sim = sim
        self.gen_set = {'nlevels': 7,
                        'file_for_it': 'H.maximum.asc',
                        'iterations': 0,
                        'indir': Paths.gw170817 + sim + '/',
                        'outdir': Paths.ppr_sims + sim + '/res_3d/'}

        self.output_it_map = {}
        self.it_time = self.set_it_output_map()
        self.list_v_ns = ['lapse', 'betax', 'betay', 'betaz',
                          'gxx', 'gxy', 'gxz', 'gyy', 'gyz', 'gzz',
                          'rho', 'temperature',
                          'thc_M0_abs_energy', 'thc_M0_abs_nua', 'thc_M0_abs_nue',
                          'thc_M0_abs_number', 'thc_M0_eave_nua', 'thc_M0_eave_nue',
                          'thc_M0_eave_nux', 'thc_M0_E_nua', 'thc_M0_E_nue',
                          'thc_M0_E_nux', 'thc_M0_flux_fac.file', 'thc_M0_ndens_nua',
                          'thc_M0_ndens_nue', 'thc_M0_ndens_nux', 'thc_M0_N_nua',
                          'thc_M0_N_nue', 'thc_M0_N_nux',
                          # 'vel[0]', 'vel[1]', 'vel[2]',
                          'velx', 'vely', 'velz',
                          'vol', 'w_lorentz', 'Y_e']

        self.name_conversion_map = {'lapse': 'alp',
                                    'vol': 'volform',
                                    'velx': 'vel[0]',
                                    'vely': 'vel[1]',
                                    'velz': 'vel[2]'}


        # initialize data storage
        self.dataset_matrix = [[0
                                 for z in range(len(self.list_v_ns))]
                                 for s in range(len(self.output_it_map.keys()))]

    # --- index work

    def set_it_output_map(self):
        """
        Loads set of files that have '1:it 2:time ...' structure to get a map
        of what output-xxxx contains what iteration (and time)
        """

        folders = []

        files = glob(self.gen_set["indir"] + "output-*" + "/data/" + self.gen_set['file_for_it'])
        # print(len(files))
        # r=root, d=directories, f = files
        # for r, d, f in os.walk(self.gen_set["indir"]):
        #     for folder in d:
        #         folders.append(os.path.join(r, folder))



        print('-' * 25 + 'LOADING it list ({})'
              .format(self.gen_set['file_for_it']) + '-' * 25)
        print("\t loading from: {}".format(self.gen_set['indir']))
        # files = locate("data/" + self.gen_set['file_for_it'], root=self.gen_set["indir"], followlinks=True)
        # # files = [f for f in glob(self.gen_set["indir"] + 'data/' + self.gen_set['file_for_it'])]
        print("\t   files found: {}".format(len(files)))
        # # remove folders like 'collated'
        # selected = []
        # for file in files:
        #     if file.__contains__('output-'):
        #         selected.append(file)
        # for overall count of number of iterations and files
        it_time = np.zeros(2)
        for file in files:
            o_name = file.split('/')
            o_dir = ''
            for o_part in o_name:
                if o_part.__contains__('output-'):
                    o_dir = o_part
            if o_dir == '':
                raise NameError("Did not find output-xxxx in {}".format(o_name))
            it_time_i = np.loadtxt(file, usecols=(0,1))
            self.output_it_map[o_dir] = it_time_i
            it_time = np.vstack((it_time, it_time_i))
        it_time = np.delete(it_time, 0, 0)
        print('outputs:{} iterations:{} [{}->{}]'.format(len(files),
                                                         len(it_time[:,0]),
                                                         int(it_time[:,0].min()),
                                                         int(it_time[:,0].max())))
        if len(it_time[:,0]) != len(set(it_time[:,0])):
            Printcolor.yellow("Warning: repetitions found in the "
                              "loaded iterations")
            iterations = np.unique(it_time[:,0])
            timestpes  = np.unique(it_time[:,1])
            if not len(iterations) == len(timestpes):
                raise ValueError("Failed attmept to remove repetitions from "
                                 "\t it and time lists. Wrong lengths: {} {}"
                                 .format(len(iterations), len(timestpes)))
        else:
            Printcolor.blue("No repetitions found in loaded it list")
            iterations = np.unique(it_time[:,0])
            timestpes  = np.unique(it_time[:,1])

        print('-' * 30 + '------DONE-----' + '-' * 30)

        return np.vstack((iterations, timestpes)).T

    def check_v_n(self, v_n):
        if v_n not in self.list_v_ns:
            raise NameError("v_n:{} not in the v_n list\n{}"
                            .format(v_n, self.list_v_ns))

    def check_it(self, it):
        if int(it) not in self.it_time[:, 0]:
            idx = find_nearest_index(np.array(self.it_time[:,0], dtype=int), int(it))
            raise NameError("it:{} not in the it list\t"
                            "closes it:{}"
                            .format(it, int(self.it_time[idx, 0])))

    def check_output(self, o_dir):
        if o_dir not in self.output_it_map.keys():
            raise NameError("output-xxxx dir:{} not in the output_list\n{}"
                            .format(o_dir, self.output_it_map.keys()))

    def i_v_n(self, v_n):
        self.check_v_n(v_n)
        return int(self.list_v_ns.index(v_n))

    def i_output(self, o_dir):
        self.check_output(o_dir)
        return int(self.output_it_map.keys().index(o_dir))

    # --- loading data

    def load_dataset(self, o_dir, v_n):

        # flist = glob(dir + "/{}.file_*.h5".format(v_n))

        if v_n in self.name_conversion_map:
            v_n_ = self.name_conversion_map[v_n]
        else:
            v_n_ = v_n

        def get_number(file_):
            return int(str(file_.split('.file_')[-1]).split('.h5')[0])

        fname = v_n_ + '.file_*.h5'
        files = locate(fname, root=self.gen_set['indir'] + o_dir + '/data/', followlinks=True)
        files = sorted(files, key=get_number)
        if len(files) == 0:
            raise ValueError("For '{}' in {} found NO files \n searched:{}"
                             .format(v_n_, o_dir, fname))

        print("\t loading '{}' ({} with {} files) ".format(v_n, o_dir, len(files)))

        # carefully creating dataset, as some .h5 might be corrupted
        try:
            dset = h5.dataset(files)
        except IOError:
            cleared_files = []
            for file_ in files:
                try:
                    tmp = h5py.File(file_, "r")
                    tmp.close()
                except IOError:
                    Printcolor.red("Error! Corrupted file: {}".format(file_.split(self.sim)[-1]))
                    break
                cleared_files.append(file_)
            dset = h5.dataset(cleared_files)

        if not v_n_ in dset.contents.keys()[0]:
            raise NameError("Loaded dataset ({}) does not contain required v_n:{}"
                            .format(dset.contents.keys()[0], v_n_))

        self.dataset_matrix[self.i_output(o_dir)][self.i_v_n(v_n)] = copy.deepcopy(dset)
        dset.close_files()
        dset.close_files()
        dset.close_files()
        dset.close_files()




        # fname = v_n + '.' + plane + '.h5'
        # files = locate(fname, root=self.gen_set['indir'] + o_dir +'/', followlinks=True)
        # print("\t Loading: {} v_n:{} dataset"
        #       .format(o_dir, v_n))
        # if len(files) > 1:
        #     raise ValueError("More than 1 file ({}) found. \nFile:{} location:{}"
        #                      .format(len(files), fname, o_dir))
        # if len(files) == 0:
        #     raise ValueError("NO fils found. \nlocation:{}"
        #                      .format(fname, o_dir))
        # dset = h5.dataset(files)
        # self.dataset_matrix[self.i_output(o_dir)][self.i_plane(plane)][self.i_v_n(v_n)] = dset

    def is_dataset_loaded(self, o_dir, v_n):
        if isinstance(self.dataset_matrix[self.i_output(o_dir)][self.i_v_n(v_n)], int):
            self.load_dataset(o_dir, v_n)

    def it_to_output_dir(self, it):
        self.check_it(it)
        req_output_data_dir = []
        for output_data_dir in self.output_it_map.keys():
            if int(it) in np.array(self.output_it_map[output_data_dir], dtype=int)[:, 0]:
                req_output_data_dir.append(output_data_dir)

        if len(req_output_data_dir) > 1:
            raise ValueError("it:{} is found in multiple outputs:{}"
                             .format(it, req_output_data_dir))
        elif len(req_output_data_dir) == 0:
            raise ValueError("it:{} not found in a output_it_map:\n{}".format(it, self.output_it_map.keys()))
        else:
            return req_output_data_dir[0]

    def get_dataset(self, it, v_n):
        self.check_it(it)
        self.check_v_n(v_n)
        o_dir = self.it_to_output_dir(it)
        self.is_dataset_loaded(o_dir, v_n)
        dset = self.dataset_matrix[self.i_output(o_dir)][self.i_v_n(v_n)]
        # if not it in dset.iterations:
        #     it__ = int(dset.iterations[find_nearest_index(np.array(dset.iterations), it)])
        #     raise ValueError("Iteration it:{} (located in {}) \n"
        #                      "not in the dataset list. Closest:{} Full list:\n{}"
        #                      .format(it, o_dir, it__, dset.iterations))
        return dset

    def del_dataset(self, it, v_n):
        o_dir = self.it_to_output_dir(it)
        self.dataset_matrix[self.i_output(o_dir)][self.i_v_n(v_n)] = 0

class EXTRACT_STORE_DATA(LOAD_STORE_DATASET):

    def __init__(self, sim):

        LOAD_STORE_DATASET.__init__(self, sim)

        self.v_n_map = {
            'lapse':        "ADMBASE::alp",
            'betax':        "ADMBASE::betax",
            'betay':        "ADMBASE::betay",
            'betaz':        "ADMBASE::betaz",
            'gxx':          "ADMBASE::gxx",
            'gxy':          "ADMBASE::gxy",
            'gxz':          "ADMBASE::gxz",
            'gyy':          "ADMBASE::gyy",
            'gyz':          "ADMBASE::gyz",
            'gzz':          "ADMBASE::gzz",
            'rho':          "HYDROBASE::rho",
            'thc_M0_abs_energy': "THC_LEAKAGEM0::thc_M0_abs_energy",
            'thc_M0_abs_nua': "THC_LEAKAGEM0::thc_M0_abs_nua",
            'thc_M0_abs_nue': "THC_LEAKAGEM0::thc_M0_abs_nue",
            'thc_M0_abs_number':"THC_LEAKAGEM0::thc_M0_abs_number",
            'thc_M0_eave_nua': "THC_LEAKAGEM0::thc_M0_eave_nua",
            'thc_M0_eave_nue':"THC_LEAKAGEM0::thc_M0_eave_nue",
            'thc_M0_eave_nux': "THC_LEAKAGEM0::thc_M0_eave_nux",
            'thc_M0_E_nua':"THC_LEAKAGEM0::thc_M0_E_nua",
            'thc_M0_E_nue': "THC_LEAKAGEM0::thc_M0_E_nue",
            'thc_M0_E_nux': "THC_LEAKAGEM0::thc_M0_E_nux",
            'thc_M0_flux_fac':"THC_LEAKAGEM0::thc_M0_flux_fac",
            'thc_M0_ndens_nua':"THC_LEAKAGEM0::thc_M0_ndens_nua",
            'thc_M0_ndens_nue': "THC_LEAKAGEM0::thc_M0_ndens_nue",
            'thc_M0_ndens_nux': "THC_LEAKAGEM0::thc_M0_ndens_nux",
            'thc_M0_N_nua': "THC_LEAKAGEM0::thc_M0_N_nua",
            'thc_M0_N_nue': "THC_LEAKAGEM0::thc_M0_N_nue",
            'thc_M0_N_nux': "THC_LEAKAGEM0::thc_M0_N_nux",
            'velx':       "HYDROBASE::vel[0]",
            'vely':       "HYDROBASE::vel[1]",
            'velz':       "HYDROBASE::vel[2]",
            # 'vel[0]':       "HYDROBASE::vel[0]",
            # 'vel[1]':       "HYDROBASE::vel[1]",
            # 'vel[2]':       "HYDROBASE::vel[2]",
            # 'volform':      "THC_CORE::volform",
            'vol':          "THC_CORE::volform",
            'w_lorentz':    "HYDROBASE::w_lorentz",
            'Y_e':          "HYDROBASE::Y_e",
            'temperature':  "HYDROBASE::temperature"
        }

        self.eos_v_n_map = {
            'eps':     "internalEnergy",
            'press':   "pressure",
            'entropy': "entropy"
        }

        self.list_grid_v_ns = ["x", "y", "z", "delta"]

        self.gen_set["eos"] = Paths.SLy4_hydo

        # self.list_eos_v_ns = ['pres', 'eps']

        self.list_v_ns += self.eos_v_n_map.keys()

        self.data_matrix = [[0
                            for z in range(len(self.list_v_ns))]
                            for s in range(len(self.it_time[:,0]))]

        self.grid_matrix = [[0
                            for z in range(len(self.list_v_ns))]
                            for s in range(len(self.it_time[:,0]))]

        self.grid_matrix_v_n = [[0
                            for z in range(len(self.list_v_ns))]
                            for s in range(len(self.it_time[:,0]))]

    # --- index work

    def i_it(self, it):
        self.check_it(it)
        idx = list(self.it_time[:,0]).index(int(it))
        return idx

    def check_grid_v_n(self, v_n):
        if not v_n in self.list_grid_v_ns:
            raise NameError("v_n:{} not in the list of grind_v_ns \n {}"
                            .format(v_n, self.list_grid_v_ns))

    # --- GRID

    def extract_grid(self, it, v_n):
        print("\t extracting grid it:{} v_n:{}".format(it, v_n))
        dset = self.get_dataset(it, v_n)
        self.grid_matrix[self.i_it(it)][self.i_v_n(v_n)] = \
            dset.get_grid(iteration=it)

    def is_grid_extracted(self, it, v_n):

        if isinstance(self.grid_matrix[self.i_it(it)][self.i_v_n(v_n)], int):
            self.extract_grid(it, v_n)

    def get_grid(self, it, v_n):

        self.check_it(it)
        self.check_v_n(v_n)
        self.is_grid_extracted(it, v_n)

        return self.grid_matrix[self.i_it(it)][self.i_v_n(v_n)]

    # --- GRID v_ns

    def i_grid(self, v_n):
        self.check_grid_v_n(v_n)
        return int(self.list_grid_v_ns.index(v_n))

    def extract_grid_v_n(self, it):
        grid = None
        for v_n in self.list_v_ns:
            if not isinstance(self.grid_matrix[self.i_it(it)][self.i_v_n(v_n)], int):
                grid = self.get_grid(it, v_n)
        if grid == None:
            print("\tNo datasets were loaded. Loading 'rho' for grid")
            grid = self.get_grid(it, 'rho')

        x = []
        y = []
        z = []
        delta = []

        # ig grid is not a carpet_grid() it won't work
        print("\t Extracting grid_v_ns for it:{}".format(it))
        for rl in range(self.gen_set["nlevels"]):
            print("\t\trl:{}".format(rl)),
            x_, y_, z_ = grid.mesh()[rl]
            delta_ = grid[rl].delta
            x.append(x_)
            y.append(y_)
            z.append(z_)
            delta.append(delta_)
            print("{} d:{}".format(np.array(x_).shape, delta_[-1]))

        self.grid_matrix_v_n[self.i_it(it)][self.i_grid("x")] = x
        self.grid_matrix_v_n[self.i_it(it)][self.i_grid("y")] = y
        self.grid_matrix_v_n[self.i_it(it)][self.i_grid("z")] = z

        self.grid_matrix_v_n[self.i_it(it)][self.i_grid("delta")] = delta

    def is_grid_v_n_extracted(self, it, v_n):

        if isinstance(self.grid_matrix_v_n[self.i_it(it)][self.i_grid(v_n)], int):
            self.extract_grid_v_n(it)

    def get_grid_v_n(self, it, v_n):

        self.is_grid_v_n_extracted(it, v_n)
        return self.grid_matrix[self.i_it(it)][self.i_grid(v_n)]

    def del_grid(self, it, v_n):

        self.check_it(it)
        self.check_v_n(v_n)

        self.grid_matrix[self.i_it(it)][self.i_v_n(v_n)] = 0

    # --- DATA

    def extract_from_eos(self, it, v_n):

        from scivis import eostable
        print("\t extracting from eos it:{} v_n:{}".format(it, v_n))
        rho  = self.get_data(it, 'rho') # list, 7 entries
        temp = self.get_data(it, 'temperature')
        ye   = self.get_data(it, 'Y_e')

        # tst = np.ma.core.masked_array()

        # print(rho[0])
        # exit(1)

        eostable.init(self.gen_set["eos"])
        # print("\t: eos contains: {}".format(eostable.get_names()))

        res = []

        for rl in range(self.gen_set["nlevels"]):

            data = eostable.evaluate(self.eos_v_n_map[v_n], rho[rl], temp[rl], ye[rl])
            res.append(data)
            assert np.array(data).shape == rho[rl].shape

        self.data_matrix[self.i_it(it)][self.i_v_n(v_n)] = res

        # dset = self.get_dataset(it, v_n)
        # try:
        #     data = dset.get_grid_data(self.get_grid(it, v_n),
        #                               iteration=it,
        #                               variable=self.v_n_map[v_n])
        #     self.data_matrix[self.i_it(it)][self.i_v_n(v_n)] = data
        # except KeyError:
        #     raise KeyError("Wrong Key. Data not found. dset contains:{} attmeped:{} it:{}".format(dset.metadata[0],
        #                                                                                           self.v_n_map[v_n],

    def extract(self, it, v_n):

        if v_n in self.eos_v_n_map.keys():
            self.extract_from_eos(it, v_n)
        else:
            print("\t extracting it:{} v_n:{}".format(it, v_n))
            dset = self.get_dataset(it, v_n)
            try:
                data = dset.get_grid_data(self.get_grid(it, v_n),
                                          iteration=it,
                                          variable=self.v_n_map[v_n])
                self.data_matrix[self.i_it(it)][self.i_v_n(v_n)] = data
            except KeyError:
                raise KeyError("Wrong Key. Data not found. dset contains:{} attmeped:{} it:{}".format(dset.metadata[0],
                                                                                                      self.v_n_map[v_n],
                                                                                                      it))

    def is_data_extracted(self, it, v_n):

        if isinstance(self.data_matrix[self.i_it(it)][self.i_v_n(v_n)], int):
            self.extract(it, v_n)

    def get_data(self, it, v_n):
        self.check_it(it)
        self.check_v_n(v_n)

        self.is_data_extracted(it, v_n)

        return self.data_matrix[self.i_it(it)][self.i_v_n(v_n)]

    def del_data(self, it, v_n):

        self.check_it(it)
        self.check_v_n(v_n)

        self.data_matrix[self.i_it(it)][self.i_v_n(v_n)] = 0

""" --- --- INPUT DATA IS IN profile FORMAT --- --- """

class LOAD_STORE_PROFILE:

    def __init__(self):
        pass

class EXTRACT_STORE(LOAD_STORE_PROFILE):

    def __init__(self):
        LOAD_STORE_PROFILE.__init__(self)



""" --- --- COMPUTE ADDITIONAL QUANTITIES --- --- """

class COMPUTE_STORE:

    def __init__(self, data_object):

        self.compute_v_ns = ["density", "vup", "metric", "shift",
                          "enthalpy", "shvel", "u_0",
                          "vlow", "vphi", "vr",
                          "dens_unb_geo", "dens_unb_bern", "dens_unb_garch",
                          "ang_mom", "ang_mom_flux",
                          "theta"]

        self.grid_v_ns = data_object.list_grid_v_ns
        self.nlevels = data_object.gen_set["nlevels"]
        self.in_data_v_ns = data_object.list_v_ns
        self.in_data_cl = data_object
        self.all_iterations = data_object.it_time[:, 0]
        if len(self.all_iterations) == 0: raise ValueError("No iterations")

        self.data_matrix = [[0
                            for z in range(len(self.compute_v_ns))]
                            for s in range(len(self.all_iterations))]

    # --- INDEX

    def check_v_n(self, v_n):
        if not v_n in self.in_data_v_ns and \
                not v_n in self.compute_v_ns and \
                not v_n in self.grid_v_ns:
            raise NameError("v_n:{} is not in v_ns that can be computed:\n{}"
                            "\n neither is it in the initial data v_ns:\n{}"
                            "\n neither is it in the grid_v_ns:\n{}"
                            .format(v_n, self.in_data_v_ns, self.compute_v_ns, self.grid_v_ns))

    def i_it(self, it):
        return self.in_data_cl.i_it(it)

    def i_v_n(self, v_n):
        self.check_v_n(v_n)
        return int(self.compute_v_ns.index(v_n))

    def check_it(self, it):
        self.in_data_cl.check_it(it)

    # --- DATA

    def compute_for_rl(self, it, v_n, rl):

        if v_n == 'density':
            arr = FORMULAS.density(self.get_data(it, "rho")[rl],
                                   self.get_data(it, "w_lorentz")[rl],
                                   self.get_data(it, "vol")[rl])

        elif v_n == 'vup':
            arr = FORMULAS.vup(self.get_data(it, "velx")[rl],
                               self.get_data(it, "vely")[rl],
                               self.get_data(it, "velz")[rl])

        elif v_n == 'metric':  # gxx, gxy, gxz, gyy, gyz, gzz
            arr = FORMULAS.metric(self.get_data(it, "gxx")[rl],
                                  self.get_data(it, "gxy")[rl],
                                  self.get_data(it, "gxz")[rl],
                                  self.get_data(it, "gyy")[rl],
                                  self.get_data(it, "gyz")[rl],
                                  self.get_data(it, "gzz")[rl])

        elif v_n == 'shift':
            arr = FORMULAS.shift(self.get_data(it, "betax")[rl],
                                 self.get_data(it, "betay")[rl],
                                 self.get_data(it, "betaz")[rl])

        elif v_n == 'enthalpy':
            arr = FORMULAS.enthalpy(self.get_data(it, "eps")[rl],
                                    self.get_data(it, "press")[rl],
                                    self.get_data(it, "rho")[rl])

        elif v_n == 'shvel':
            arr = FORMULAS.shvel(self.get_data(it, "shift")[rl],
                                 self.get_data(it, "vlow")[rl])

        elif v_n == 'u_0':
            arr = FORMULAS.u_0(self.get_data(it, "w_lorentz")[rl],
                               self.get_data(it, "shvel")[rl],  # not input
                               self.get_data(it, "lapse")[rl])

        elif v_n == 'vlow':
            arr = FORMULAS.vlow(self.get_data(it, "metric")[rl],
                                self.get_data(it, "vup")[rl])

        elif v_n == 'vphi':
            arr = FORMULAS.vphi(self.get_data(it, "x")[rl],
                                self.get_data(it, "y")[rl],
                                self.get_data(it, "vlow")[rl])

        elif v_n == 'vr':
            arr = FORMULAS.vr(self.get_data(it, "x")[rl],
                              self.get_data(it, "y")[rl],
                              self.get_data(it, "vup")[rl])

        elif v_n == 'theta':
            arr = FORMULAS.theta(self.get_data(it, "x")[rl],
                                 self.get_data(it, "y")[rl],
                                 self.get_data(it, "z")[rl])

        elif v_n == 'ang_mom':
            arr = FORMULAS.ang_mom(self.get_data(it, "rho")[rl],
                                   self.get_data(it, "eps")[rl],
                                   self.get_data(it, "press")[rl],
                                   self.get_data(it, "w_lorentz")[rl],
                                   self.get_data(it, "vol")[rl],
                                   self.get_data(it, "vphi")[rl])

        elif v_n == 'ang_mom_flux':
            arr = FORMULAS.ang_mom_flux(self.get_data(it, "ang_mom")[rl],
                                        self.get_data(it, "lapse")[rl],
                                        self.get_data(it, "vr")[rl])

        elif v_n == 'dens_unb_geo':
            arr = FORMULAS.dens_unb_geo(self.get_data(it, "u_0")[rl],
                                        self.get_data(it, "rho")[rl],
                                        self.get_data(it, "w_lorentz")[rl],
                                        self.get_data(it, "vol")[rl])

        elif v_n == 'dens_unb_bern':
            arr = FORMULAS.dens_unb_bern(self.get_data(it, "enthalpy")[rl],
                                         self.get_data(it, "u_0")[rl],
                                         self.get_data(it, "rho")[rl],
                                         self.get_data(it, "w_lorentz")[rl],
                                         self.get_data(it, "vol")[rl])

        elif v_n == 'dens_unb_garch':
            arr = FORMULAS.dens_unb_garch(self.get_data(it, "enthalpy")[rl],
                                          self.get_data(it, "u_0")[rl],
                                          self.get_data(it, "lapse")[rl],
                                          self.get_data(it, "press")[rl],
                                          self.get_data(it, "rho")[rl],
                                          self.get_data(it, "w_lorentz")[rl],
                                          self.get_data(it, "vol")[rl])

        else:
            raise NameError("No method found for v_n:{} Add entry to 'compute_for_rl()'"
                            .format(v_n))

        return arr

    def compute(self, it, v_n):

        res = []

        for rl in range(self.nlevels):
            res.append(self.compute_for_rl(it, v_n, rl))

        self.data_matrix[self.i_it(it)][self.i_v_n(v_n)] = res

    def is_computed(self, it, v_n):

        if isinstance(self.data_matrix[self.i_it(it)][self.i_v_n(v_n)], int):
            self.compute(it, v_n)

    def get_data(self, it, v_n):
        self.check_it(it)
        self.check_v_n(v_n)
        if v_n in self.in_data_v_ns:
            return self.in_data_cl.get_data(it, v_n)
        elif v_n in self.in_data_cl.list_grid_v_ns:
            return self.in_data_cl.get_grid_v_n(it, v_n)
        elif v_n in self.compute_v_ns:
            self.is_computed(it, v_n)
            return self.data_matrix[self.i_it(it)][self.i_v_n(v_n)]



if __name__ == '__main__':

    # def get_number(file_):
    #     return int(str(file_.split('.file_')[-1]).split('.h5')[0])
    # path = "/data1/numrel/WhiskyTHC/Backup/2018/GW170817/SLy4_M13641364_M0_SR/output-0010/data/"
    # files = locate("gzz.file_*", root=path, followlinks=True)
    # print(files)
    # files = sorted(files, key=get_number)

    # file = "/data1/numrel/WhiskyTHC/Backup/2018/GW170817/SLy4_M13641364_M0_SR/output-0010/data/gzz.file_54.h5"
    # h5py.File(file, "r")
    # for file in files:
    #     print(file)
    #     loaded = h5py.File(file, "r")


    # dset = h5.dataset(files, maxnfiles=70)

    ''' --- DEBUGGING --- '''
    # ds_ = LOAD_STORE_DATASET("SLy4_M13641364_M0_SR")
    # ds_.load_dataset('output-0010', 'velx')
    # ds_.get_dataset(219904, 'gzz')
    # exit(1)

    ''''''
    # data_ = EXTRACT_STORE_DATA("SLy4_M13641364_M0_SR")
    # data_.get_data(245760, 'x')
    # data_.get_data(245760, 'y')
    # data_.extract_from_eos(245760, 'eps')
    # data_.extract_from_eos(245760, 'press')
    ''''''
    data_ = EXTRACT_STORE_DATA("SLy4_M13641364_M0_SR")
    comp = COMPUTE_STORE(data_)
    # comp.get_data(180224, 'vlow')
    print(comp.get_data(180224, 'dens_unb_bern'))