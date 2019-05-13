# ///////////////////////////////////////////////////////////////////////////////
# // disk: analysis of the disk in a THC run from the 3D carpet data
# // Copyright (C) 2019, Vsevolod Nedora <vsevolod.nedora@uni-jena.de>
# //
# // This program is free software: you can redistribute it and/or modify
# // it under the terms of the GNU General Public License as published by
# // the Free Software Foundation, either version 3 of the License, or
# // (at your option) any later version.
# //
# // This program is distributed in the hope that it will be useful,
# // but WITHOUT ANY WARRANTY; without even the implied warranty of
# // MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# // GNU General Public License for more details.
# //
# // You should have received a copy of the GNU General Public License
# // along with this program.  If not, see <http://www.gnu.org/licenses/>.
# ///////////////////////////////////////////////////////////////////////////////
# // This package contains utilities for
# // . Parsing 3D profile output of the carpet and hydro thorns
# // . Produce histograms-correlation maps and other properties
# //   of the disk for a cross analysis of multiple variables,
# //   including, density, angular momentum, its flux, density unbound.
# // . Produce total mass of the disk, interpolated 2D slices
# //
# // Usage
# // . Set the setup() according to required tasks to do
# // . python disk /path/to/profile.h5
# ///////////////////////////////////////////////////////////////////////////////

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

""" --- --- SETUP --- --- """

def setup():
    rho_const = 6.176269145886162e+17

    tasks_for_rl = {
        "mask": {'rm_rl': True, 'rho': [6.e4 / rho_const, 1.e13 / rho_const], 'lapse': [0.15, 1.]},  # rho in cgs

        # "task1::correlation": [
        #     {"v_n": "rho",  "edges": 10.0 ** np.linspace(4.0, 16.0, 500) / rho_const},  # not in CGS :^
        #     {"v_n": "temp", "edges": 10.0 ** np.linspace(-2, 2, 300)},
        #     {"v_n": "Ye",   "edges": np.linspace(0, 0.5, 300)}
        # ],

        # "task::mass": {}, # computes total mass

        # "task2::correlation": [
        #     {"v_n": "ang_mom_flux",  "edges": 10.0 ** np.linspace(-12., -7, 500)},  # not in CGS :^
        #     {"v_n": "dens_unb_bern", "edges": 10.0 ** np.linspace(-9., -7., 500)}
        # ]

        # "task3::correlation": [
        #     {"v_n": "ang_mom_flux", "edges": 10.0 ** np.linspace(-12., -5, 500)},  # not in CGS :^
        #     {"v_n": "Ye", "edges": np.linspace(0.01, 0.5, 500)}
        # ]

        "task3::correlation": [
            {"v_n": "density", "edges": 10.0 ** np.linspace(-12., -5, 500)},  # not in CGS :^
            {"v_n": "theta", "edges": np.linspace(0, 3.2, 500)}
        ]

    }

    tasks_for_int = {
        "grid": {'type': 'cyl', 'n_r': 150, 'n_phi': 150, 'n_z': 100},
        "save": {"grid_v_ns":["phi_cyl", "r_cyl", "z_cyl",
                              "dphi_cyl", "dr_cyl", "dz_cyl"],
                 "v_ns": ['ang_mom', 'ang_mom_flux', 'density', 'dens_unb_geo',
                          'dens_unb_bern','rho', 'temp', 'Ye']}
    }

    general_settings = {}
    general_settings['indataformat'] = 'profile' # profile output-xxxx
    general_settings['indatadir'] = './'
    general_settings['outdir'] = './postprocess/'
    general_settings['figdir'] = './postprocess/'
    general_settings['nlevels'] = 7

    return general_settings, tasks_for_rl, tasks_for_int

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

    def __init__(self, grid_info):

        self.grid_info = grid_info



        self.grid_type = grid_info["type"]

        # self.carpet_grid = carpet_grid

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

    def get_xi(self):
        return np.column_stack([self.x_cyl_3d.flatten(),
                                self.y_cyl_3d.flatten(),
                                self.z_cyl_3d.flatten()])

    def get_shape(self):
        return self.x_cyl_3d.shape

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
                            .format(v_n, self.list_int_grid_v_ns))

class FORMULAS:

    def __init__(self):
        pass

    @staticmethod
    def r(x, y):
        return np.sqrt(x ** 2 + y ** 2)# + z ** 2)

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
    def vr(x, y, r, vup):
        # r = np.sqrt(x ** 2 + y ** 2)
        return (x / r) * vup[0] + (y / r) * vup[1]

    @staticmethod
    def theta(r, z):
        # r = np.sqrt(x ** 2 + y ** 2 + z ** 2)
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

class SAVE_RESULT:

    def __init__(self):
        pass

    @staticmethod
    def correlation(it, v_ns, edges, corr, outdir):


        name = outdir + str(int(it)) + '_corr_'
        for v_n in v_ns:
            name += v_n
            if v_n != v_ns[-1]:
                name += '_'
        name += '.h5'
        print('-' * 30 + 'SAVING CORRELATION' + '-' * 30)
        print(' ' * 30 + '{}'.format(name) + ' ' * 30)
        outfile = h5py.File(name, "w")
        for v_n, edge in zip(v_ns, edges):
            outfile.create_dataset(v_n, data=edge)
        outfile.create_dataset("mass", data=corr)
        print('-' * 30 + '------DONE-----' + '-' * 30)
        print('\n')

        return name

""" --- --- DATA PROCESSING --- --- """

class LOAD_PROFILE:

    def __init__(self, fname):

        self.nlevels = 7
        self.profile = fname
        self.dfile = h5py.File(fname, "r")
        group_0 = self.dfile["reflevel={}".format(0)]
        self.time = group_0.attrs["time"] * 0.004925794970773136 * 1e-3 # [sec]
        self.iteration = group_0.attrs["iteration"]
        print("\t\t time: {}".format(self.time))
        print("\t\t iteration: {}".format(self.iteration))
        self.grid = self.read_carpet_grid(self.dfile)

        self.list_prof_v_ns = [
                             "rho", "w_lorentz", "vol",  # basic
                             "press", "eps", "lapse",    # basic + lapse
                             "velx", "vely", "velz",     # velocities
                             "gxx", "gxy", "gxz", "gyy", "gyz", "gzz",  # metric
                             "betax", "betay", "betaz",  # shift components
                             'temp', 'Ye']

        self.list_grid_v_ns = ["x", "y", "z", "delta"]

    def read_carpet_grid(self, dfile):
        import scidata.carpet.grid as grid
        L = []
        for il in range(self.nlevels):
            gname = "reflevel={}".format(il)
            group = dfile[gname]
            level = grid.basegrid()
            level.delta = np.array(group.attrs["delta"])
            level.dim = 3
            level.time = group.attrs["time"]
            # level.timestep = group.attrs["timestep"]
            level.directions = range(3)
            level.iorigin = np.array([0, 0, 0], dtype=np.int32)
            level.origin = np.array(group.attrs["extent"][0::2])
            level.n = np.array(group["rho"].shape, dtype=np.int32)
            level.rlevel = il
            L.append(level)
        return grid.grid(sorted(L, key=lambda x: x.rlevel))

    def get_prof_arr(self, rl, v_n):
        group = self.dfile["reflevel={}".format(rl)]
        try:
            arr = np.array(group[v_n])
        except:
            print('\nAvailable Parameters:')
            print(list(v_n_aval for v_n_aval in group))
            print('\n')
            raise ValueError('v_n:{} not in file:{}'.format(v_n, self.profile))
        return arr

    def get_prof_delta(self, rl):
        return self.grid[rl].delta

    def get_prof_x_y_z(self, rl):
        x, y, z = self.grid.mesh()[rl]
        return x, y, z

    def __delete__(self, instance):
        instance.dfile.close()

class COMPUTE_STORE(LOAD_PROFILE):

    def __init__(self, fname):
        LOAD_PROFILE.__init__(self, fname)

        rho_const = 6.176269145886162e+17
        self.task_setup = {'rm_rl': True, # REMOVE previouse ref. level from the next
                           'rho': [6.e4 / rho_const, 1.e13 / rho_const], # REMOVE atmo and NS
                           'lapse': [0.15, 1.]} # remove apparent horizon

        self.list_comp_v_ns = [
            "density", "vup", "metric", "shift",
            "enthalpy", "shvel", "u_0",
            "vlow", "vphi", "vr",
            "dens_unb_geo", "dens_unb_bern", "dens_unb_garch",
            "ang_mom", "ang_mom_flux",
            "theta", "r" # assumes cylindircal coordinates. r = x^2 + y^2
        ]

        self.list_all_v_ns = self.list_prof_v_ns + \
                             self.list_grid_v_ns + \
                             self.list_comp_v_ns

        self.data_matrix = [[np.zeros(0,)
                             for x in range(self.nlevels)]
                             for y in range(len(self.list_all_v_ns))]

    def check_v_n(self, v_n):
        if v_n not in self.list_all_v_ns:
            raise NameError("v_n:{} not in the v_n list \n{}"
                            .format(v_n, self.list_all_v_ns))

    def i_v_n(self, v_n):
        self.check_v_n(v_n)
        return int(self.list_all_v_ns.index(v_n))

    def set_data(self, rl, v_n, arr):
        self.data_matrix[self.i_v_n(v_n)][rl] = arr

    def extract_data(self, rl, v_n):
        data = self.get_prof_arr(rl, v_n)
        self.data_matrix[self.i_v_n(v_n)][rl] = data

    def extract_grid_data(self, rl, v_n):
        if v_n in ["x", "y", "z"]:
            x, y, z = self.get_prof_x_y_z(rl)
            self.data_matrix[self.i_v_n("x")][rl] = x
            self.data_matrix[self.i_v_n("y")][rl] = y
            self.data_matrix[self.i_v_n("z")][rl] = z
        elif v_n == "delta":
            delta = self.get_prof_delta(rl)
            self.data_matrix[self.i_v_n("x")][rl] = delta
        else:
            raise NameError("Grid variable {} not recognized".format(v_n))

    def compute_data(self, rl, v_n):

        if v_n == 'density':
            arr = FORMULAS.density(self.get_comp_data(rl, "rho"),
                                   self.get_comp_data(rl, "w_lorentz"),
                                   self.get_comp_data(rl, "vol"))

        elif v_n == 'vup':
            arr = FORMULAS.vup(self.get_comp_data(rl, "velx"),
                               self.get_comp_data(rl, "vely"),
                               self.get_comp_data(rl, "velz"))

        elif v_n == 'metric':  # gxx, gxy, gxz, gyy, gyz, gzz
            arr = FORMULAS.metric(self.get_comp_data(rl, "gxx"),
                                  self.get_comp_data(rl, "gxy"),
                                  self.get_comp_data(rl, "gxz"),
                                  self.get_comp_data(rl, "gyy"),
                                  self.get_comp_data(rl, "gyz"),
                                  self.get_comp_data(rl, "gzz"))

        elif v_n == 'shift':
            arr = FORMULAS.shift(self.get_comp_data(rl, "betax"),
                                 self.get_comp_data(rl, "betay"),
                                 self.get_comp_data(rl, "betaz"))

        elif v_n == 'enthalpy':
            arr = FORMULAS.enthalpy(self.get_comp_data(rl, "eps"),
                                    self.get_comp_data(rl, "press"),
                                    self.get_comp_data(rl, "rho"))

        elif v_n == 'shvel':
            arr = FORMULAS.shvel(self.get_comp_data(rl, "shift"),
                                 self.get_comp_data(rl, "vlow"))

        elif v_n == 'u_0':
            arr = FORMULAS.u_0(self.get_comp_data(rl, "w_lorentz"),
                               self.get_comp_data(rl, "shvel"),  # not input
                               self.get_comp_data(rl, "lapse"))

        elif v_n == 'vlow':
            arr = FORMULAS.vlow(self.get_comp_data(rl, "metric"),
                                self.get_comp_data(rl, "vup"))

        elif v_n == 'vphi':
            arr = FORMULAS.vphi(self.get_comp_data(rl, "x"),
                                self.get_comp_data(rl, "y"),
                                self.get_comp_data(rl, "vlow"))

        elif v_n == 'vr':
            arr = FORMULAS.vr(self.get_comp_data(rl, "x"),
                              self.get_comp_data(rl, "y"),
                              self.get_comp_data(rl, "r"),
                              self.get_comp_data(rl, "vup"))

        elif v_n == "r":
            arr = FORMULAS.r(self.get_comp_data(rl, "x"),
                             self.get_comp_data(rl, "y"))

        elif v_n == 'theta':
            arr = FORMULAS.theta(self.get_comp_data(rl, "r"),
                                 self.get_comp_data(rl, "z"))

        elif v_n == 'ang_mom':
            arr = FORMULAS.ang_mom(self.get_comp_data(rl, "rho"),
                                   self.get_comp_data(rl, "eps"),
                                   self.get_comp_data(rl, "press"),
                                   self.get_comp_data(rl, "w_lorentz"),
                                   self.get_comp_data(rl, "vol"),
                                   self.get_comp_data(rl, "vphi"))

        elif v_n == 'ang_mom_flux':
            arr = FORMULAS.ang_mom_flux(self.get_comp_data(rl, "ang_mom"),
                                        self.get_comp_data(rl, "lapse"),
                                        self.get_comp_data(rl, "vr"))

        elif v_n == 'dens_unb_geo':
            arr = FORMULAS.dens_unb_geo(self.get_comp_data(rl, "u_0"),
                                        self.get_comp_data(rl, "rho"),
                                        self.get_comp_data(rl, "w_lorentz"),
                                        self.get_comp_data(rl, "vol"))

        elif v_n == 'dens_unb_bern':
            arr = FORMULAS.dens_unb_bern(self.get_comp_data(rl, "enthalpy"),
                                         self.get_comp_data(rl, "u_0"),
                                         self.get_comp_data(rl, "rho"),
                                         self.get_comp_data(rl, "w_lorentz"),
                                         self.get_comp_data(rl, "vol"))

        elif v_n == 'dens_unb_garch':
            arr = FORMULAS.dens_unb_garch(self.get_comp_data(rl, "enthalpy"),
                                          self.get_comp_data(rl, "u_0"),
                                          self.get_comp_data(rl, "lapse"),
                                          self.get_comp_data(rl, "press"),
                                          self.get_comp_data(rl, "rho"),
                                          self.get_comp_data(rl, "w_lorentz"),
                                          self.get_comp_data(rl, "vol"))

        else:
            raise NameError("No method found for v_n:{} rl:{} Add entry to 'compute()'"
                            .format(v_n, rl))

        self.data_matrix[self.i_v_n(v_n)][rl] = arr

    def is_available(self, rl, v_n):
        self.check_v_n(v_n)
        data = self.data_matrix[self.i_v_n(v_n)][rl]
        if len(data) == 0:
            if v_n in self.list_prof_v_ns:
                self.extract_data(rl, v_n)
            elif v_n in self.list_grid_v_ns:
                self.extract_grid_data(rl, v_n)
            elif v_n in self.list_comp_v_ns:
                self.compute_data(rl, v_n)
            else:
                raise NameError("v_n is not recognized: '{}' [COMPUTE STORE]".format(v_n))

    def get_comp_data(self, rl, v_n):
        self.check_v_n(v_n)
        self.is_available(rl, v_n)

        return self.data_matrix[self.i_v_n(v_n)][rl]

    def __delete__(self, instance):
        instance.dfile.close()
        instance.data_matrix = [[np.zeros(0, )
                             for x in range(self.nlevels)]
                            for y in range(len(self.list_all_v_ns))]

class MASK_STORE(COMPUTE_STORE):

    def __init__(self, fname):
        COMPUTE_STORE.__init__(self, fname)

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
            if ii > 0 and self.task_setup["rm_rl"]:
                x_ = (x[ii][:, :, :] <= x[ii - 1][:, 0, 0].max()) & (
                        x[ii][:, :, :] >= x[ii - 1][:, 0, 0].min())
                y_ = (y[ii][:, :, :] <= y[ii - 1][0, :, 0].max()) & (
                        y[ii][:, :, :] >= y[ii - 1][0, :, 0].min())
                z_ = (z[ii][:, :, :] <= z[ii - 1][0, 0, :].max()) & (
                        z[ii][:, :, :] >= z[ii - 1][0, 0, :].min())
                mask = mask & np.invert((x_ & y_ & z_))

            for v_n in self.task_setup.keys()[1:]:
                self.check_v_n(v_n)
                if len(self.task_setup[v_n]) != 2:
                    raise NameError("Error. 2 values are required to set a limit. Give {} for {}"
                                     .format(self.task_setup[v_n], v_n))
                arr_1 = self.get_comp_data(rl, v_n)[3:-3, 3:-3, 3:-3]
                min_val = float(self.task_setup[v_n][0])
                max_val = float(self.task_setup[v_n][1])
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

class MAINMETHODS_STORE(MASK_STORE):

    def __init__(self, fname, sim):

        MASK_STORE.__init__(self, fname)

        self.sim = sim

        # "v_n": "temp", "edges": np.array()
        ''''''
        # "v_n": "temp", "points: number, "scale": "log", (and "min":number, "max":number)

        rho_const = 6.176269145886162e+17
        self.corr_task_dic_temp_ye = [
            # {"v_n": "rho",  "edges": 10.0 ** np.linspace(4.0, 16.0, 500) / rho_const},  # not in CGS :^
            {"v_n": "temp", "edges": 10.0 ** np.linspace(-2, 2, 300)},
            {"v_n": "Ye",   "edges": np.linspace(0, 0.5, 300)}
        ]

        self.corr_task_dic_rho_ye = [
            # {"v_n": "temp", "edges": 10.0 ** np.linspace(-2, 2, 300)},
            {"v_n": "rho",  "edges": 10.0 ** np.linspace(4.0, 13.0, 500) / rho_const},  # not in CGS :^
            {"v_n": "Ye",   "edges": np.linspace(0, 0.5, 300)}
        ]

        self.corr_task_dic_rho_theta = [
            {"v_n": "rho", "edges": 10.0 ** np.linspace(4.0, 13.0, 500) / rho_const},  # not in CGS :^
            {"v_n": "theta", "edges": np.linspace(0, 0.5*np.pi, 300)}
        ]

        self.corr_task_dic_rho_r = [
            {"v_n": "rho", "edges": 10.0 ** np.linspace(4.0, 13.0, 500) / rho_const},  # not in CGS :^
            {"v_n": "r", "edges": np.linspace(0, 100, 500)}
        ]

        self.corr_task_dic_rho_ang_mom = [
            {"v_n": "rho", "edges": 10.0 ** np.linspace(4.0, 13.0, 500) / rho_const},  # not in CGS :^
            {"v_n": "ang_mom", "points": 300, "scale": "log", "min":1e-9} # find min, max yourself
        ]

        self.corr_task_dic_rho_ang_mom_flux = [
            {"v_n": "rho", "edges": 10.0 ** np.linspace(4.0, 13.0, 500) / rho_const},  # not in CGS :^
            {"v_n": "ang_mom_flux", "points": 300, "scale": "log", "min":1e-12}
        ]

        self.corr_task_dic_rho_dens_unb_bern = [
            {"v_n": "rho", "edges": 10.0 ** np.linspace(4.0, 13.0, 500) / rho_const},  # not in CGS :^
            {"v_n": "dens_unb_bern", "edges": 10.0 ** np.linspace(-12., -6., 300)}
        ]

        self.corr_task_dic_ang_mom_flux_theta = [
            {"v_n": "ang_mom_flux", "points": 300, "scale": "log", "min":1e-12},  # not in CGS :^
            {"v_n": "theta", "edges": np.linspace(0, 0.5*np.pi, 500)}
        ]

        self.corr_task_dic_ang_mom_flux_dens_unb_bern = [
            {"v_n": "ang_mom_flux", "points": 500, "scale": "log", "min":1e-12},  # not in CGS :^
            {"v_n": "dens_unb_bern", "edges": 10.0 ** np.linspace(-12., -6., 500)}
        ]

        self.corr_task_dic_inv_ang_mom_flux_dens_unb_bern = [
            {"v_n": "inv_ang_mom_flux", "points": 500, "scale": "log", "min":1e-12},  # not in CGS :^
            {"v_n": "dens_unb_bern", "edges": 10.0 ** np.linspace(-12., -6., 500)}
        ]

    def get_total_mass(self, multiplier=2., save=False):
        mass = 0.
        for rl in range(self.nlevels):
            density = np.array(self.get_masked_data(rl, "density"))
            delta = self.get_prof_delta(rl)
            mass += float(multiplier * np.sum(density) * np.prod(delta))

        print("it:{} mass:{:3f}Msun".format(self.iteration, mass))

        if save:
            path = Paths.ppr_sims + self.sim + "/res_3d/"  + str(self.iteration) + '/'
            fname = "disk_mass.txt".format(self.iteration)

            if not os.path.exists(path):
                os.makedirs(path)

            np.savetxt(path + fname, np.array([mass]), fmt='%.5f')

        return mass

    def get_min_max(self, v_n):
        # self.check_v_n(v_n)
        min_, max_ = [], []
        for rl in range(self.nlevels):

            if v_n == 'inv_ang_mom_flux':
                v_n = 'ang_mom_flux'
                data = -1. * self.get_masked_data(rl, v_n)
            else:
                data = self.get_masked_data(rl, v_n)
            min_.append(data.min())
            max_.append(data.max())
        min_ = np.array(min_)
        max_ = np.array(max_)
        return min_.min(), max_.max()
            # print("rl:{} min:{} max:{}".format(rl, data.min(), data.max()))

    def get_edges(self, corr_task_dic):

        dic = dict(corr_task_dic)

        if "edges" in dic.keys():
            return dic["edges"]

        if "points" in dic.keys() and "scale" in dic.keys():
            min_, max_ = self.get_min_max(dic["v_n"])
            if "min" in dic.keys(): min_ = dic["min"]
            if "max" in dic.keys(): max_ = dic["max"]
            if dic["scale"] == "log":
                print("v_n: {} is in ({}->{}) range"
                      .format(dic["v_n"], min_, max_))
                if min_ <= 0: raise ValueError("for Logscale min cannot be < 0. "
                                               "found: {}".format(min_))
                if max_ <= 0:raise ValueError("for Logscale max cannot be < 0. "
                                               "found: {}".format(max_))
                edges = 10.0 ** np.linspace(np.log10(min_), np.log10(max_), dic["points"])

            elif dic["scale"] == "linear":
                edges = np.linspace(min_, max_, dic["points"])
            else:
                raise NameError("Unrecoginzed scale: {}".format(dic["scale"]))
            return edges

    def get_correlation(self, corr_task_dic, multiplier=2., save=False):

        v_ns = []
        edges = []
        for setup_dictionary in corr_task_dic:
            v_ns.append(setup_dictionary["v_n"])
            edges.append(self.get_edges(setup_dictionary))
        edges = tuple(edges)

        correlation = np.zeros([len(edge) - 1 for edge in edges])
        for rl in range(self.nlevels):
            data = []
            weights = self.get_masked_data(rl, "density").flatten() * \
                      np.prod(self.get_prof_delta(rl)) * multiplier
            for i_vn, v_n in enumerate(v_ns):

                if v_n == 'inv_ang_mom_flux':
                    v_n = 'ang_mom_flux'
                    data.append(-1. * self.get_masked_data(rl, v_n).flatten())
                else:
                    data.append(self.get_masked_data(rl, v_n).flatten())


            data = tuple(data)
            tmp, _ = np.histogramdd(data, bins=edges, weights=weights)
            correlation += tmp

        if save:
            path = Paths.ppr_sims + self.sim + "/res_3d/" + str(self.iteration) + '/'
            fname = "corr_".format(self.iteration)
            for v_n in v_ns:
                fname += v_n
                if v_n != v_ns[-1]:
                    fname += '_'
            fname += '.h5'

            if not os.path.exists(path):
                os.makedirs(path)

            outfile = h5py.File(path + fname, "w")
            for v_n, edge in zip(v_ns, edges):
                outfile.create_dataset(v_n, data=edge)
            outfile.create_dataset("mass", data=correlation)
            outfile.close()

        return correlation

    def __delete__(self, instance):
        instance.dfile.close()
        instance.data_matrix = [[np.zeros(0, )
                                 for x in range(self.nlevels)]
                                 for y in range(len(self.list_all_v_ns))]
        instance.mask_matrix = [np.ones(0, dtype=bool) for x in range(self.nlevels)]





class INTERPOLATE_STORE(MAINMETHODS_STORE):

    def __init__(self, fname, sim, grid_object):
        """
            fname - of the profile

            sim - name of the simulation (for directory searching)

            grid_object -
                object of the class with the interpolated grid. Must contain:

                list(list_grid_v_ns) that comtains the list of variable names of new grid,
                    for examply x_cyl ... z_cyl, r_cyl ... z_cyl, dr_cyl ... dz_cyl
                get_xi() function that returns array of the type
                    return np.column_stack([self.x_cyl_3d.flatten(),
                                self.y_cyl_3d.flatten(),
                                self.z_cyl_3d.flatten()])
                get_shape() function that returns the shape of the new grid such as
                    example: self.x_cyl_3d.shape
                get_int_grid(v_n) fucntion that returns the array of the new grid
                    for variable v_n. For ecample for v_n = "r_cyl"

        :param fname:
        :param sim:
        :param grid_object:
        """

        MAINMETHODS_STORE.__init__(self, fname, sim)

        self.new_grid = grid_object

        self.list_int_grid_v_ns = grid_object.list_int_grid_v_ns
        self.list_int_v_ns = self.list_prof_v_ns + \
                             self.list_comp_v_ns + \
                             self.list_grid_v_ns

        self.int_data_matrix = [np.zeros(0,) for y in range(len(self.list_int_v_ns))]


    def check_int_v_n(self, v_n):
        if v_n not in self.list_int_v_ns:
            raise NameError("v_n: '{}' not in the v_n list \n{}"
                            .format(v_n, self.list_int_v_ns))

    def i_int_v_n(self, v_n):
        self.check_int_v_n(v_n)
        return int(self.list_int_v_ns.index(v_n))

    def do_append_grid_var(self, v_n):
        self.int_data_matrix[self.i_int_v_n(v_n)] = \
            self.new_grid.get_int_grid(v_n)

    def do_interpolate(self, v_n):

        tmp = []
        for rl in range(self.nlevels):
            data = self.get_comp_data(rl, v_n)
            tmp.append(data)

        xi = self.new_grid.get_xi()
        shape = self.new_grid.get_shape()

        print("\t\tInterpolating: {}".format(v_n))
        F = Interpolator(self.grid, tmp, interp=1)
        arr = F(xi).reshape(shape)

        self.int_data_matrix[self.i_int_v_n(v_n)] = arr

    def is_data_interpolated(self, v_n):

        if len(self.int_data_matrix[self.i_int_v_n(v_n)]) == 0:
            if v_n in self.list_int_grid_v_ns:
                self.do_append_grid_var(v_n)
            else:
                self.do_interpolate(v_n)




    def get_int(self, v_n):
        self.check_int_v_n(v_n)
        self.is_data_interpolated(v_n)
        return self.int_data_matrix[self.i_int_v_n(v_n)]


class INTMETHODS_STORE(INTERPOLATE_STORE):

    def __init__(self, fname, sim, grid_object):

        INTERPOLATE_STORE.__init__(self, fname, sim, grid_object)

    def save_new_grid(self):

        grid_type = self.new_grid.grid_info['type']

        path = Paths.ppr_sims + self.sim + "/res_3d/" + str(self.iteration) + '/'
        outfile = h5py.File(path + grid_type + '_grid.h5', "w")

        if not os.path.exists(path):
            os.makedirs(path)

        # print("Saving grid...")
        for v_n in self.list_int_grid_v_ns:
            outfile.create_dataset(v_n, data=self.new_grid.get_int_grid(v_n))
        outfile.close()


    def save_int_v_n(self, v_n):

        path = Paths.ppr_sims + self.sim + "/res_3d/" + str(self.iteration) + '/'
        grid_type = self.new_grid.grid_info['type']
        outfile = h5py.File(path + grid_type + '_' + v_n + '.h5', "w")
        outfile.create_dataset(v_n, data=self.get_int(v_n))
        outfile.close()


""" --- --- LOADING & PLOTTING RESILTS --- --- """

class LOAD_RES_CORR:

    def __init__(self, sim):

        self.sim = sim
        fname = '*.h5'
        root = Paths.gw170817+ sim +'/profiles/3d/'
        files = locate(fname, root=root, followlinks=False)
        iterations = []

        if len(files) == 0:
            raise NameError("No iterations found in the root:{}".format(root))
        if os.path.isfile(root+"ittime.txt"):
            list_iterations, self.times = np.loadtxt(root+"ittime.txt", usecols=(0, 1), unpack=True)
            self.list_iterations = list(list_iterations)
        else:
            for file_ in files:
                iterations.append(int(str(file_.split('/')[-1]).split('.h5')[0]))
            self.times = interpoate_time_form_it(iterations, Paths.gw170817+sim+'/')
            self.list_iterations = iterations

        # self.list_corr_v_ns = ["temp_Ye", "rho_Ye", "rho_theta", "rho_r",
        #                        "rho_ang_mom", "rho_ang_mom_flux", "rho_dens_unb_bern",
        #                        "ang_mom_flux_dens_unb_bern", "ang_mom_flux_theta"]

        self.list_corr_v_ns = ["temp", "Ye", "rho", "theta", "r",
                               "ang_mom", "ang_mom_flux", "dens_unb_bern",
                               "inv_ang_mom_flux"
                               ]

        self.corr_matrix = [[np.zeros(0,)
                             for x in range(2 * len(self.list_corr_v_ns) + 2)] # Here 2 * () as for correlation 2 v_ns are aneeded
                             for y in range(len(self.list_iterations))]

    def check_it(self, it):
        if not it in self.list_iterations:
            raise NameError("it:{} not in the list of iterations\n{}"
                            .format(it, self.list_iterations))

    def check_v_n(self, v_n):
        if not v_n in self.list_corr_v_ns:
            raise NameError("v_n:{} not in list of corr_v_ns\n{}"
                            .format(v_n, self.list_corr_v_ns))

    def i_v_n(self, v_n_x, v_n_y):
        self.check_v_n(v_n_x)
        self.check_v_n(v_n_y)
        idx1 = int(self.list_corr_v_ns.index(v_n_x))
        idx2 = int(self.list_corr_v_ns.index(v_n_y))
        # shift = len(self.list_corr_v_ns)
        return int(idx1 + idx2)

    def i_it(self, it):
        self.check_it(it)
        return int(self.list_iterations.index(it))

    def get_corr_fpath(self, it, v_n):
        self.check_it(it)
        # self.check_v_n(v_n)
        fpath = Paths.ppr_sims + self.sim + "/res_3d/" + str(it) + "/corr_" + v_n + ".h5"
        if not os.path.isfile(fpath):
            raise IOError("Correlation file not found:\n{}".format(fpath))
        return fpath

    def load_corr_file_old(self, it, v_n):

        v_n = str(v_n)

        fpath = self.get_corr_fpath(it, v_n)

        dfile = h5py.File(fpath, "r")

        v_ns_in_data = []
        for v_n_ in dfile:
            v_ns_in_data.append(v_n_)

        if not "mass" in v_ns_in_data:
            raise NameError("mass is not found in file:{}".format(fpath))

        if len(v_ns_in_data) > 3:
            raise NameError("More than 3 datasets found in corr file: {}".format(fpath))

        v_ns_in_data.remove("mass")

        for v_n__ in v_ns_in_data:
            if not v_n__ in v_n:
                raise NameError("in_data_v_n: {} is not in corr name v_n: {}"
                                .format(v_n__, v_n))


        part1 = v_n.split(v_ns_in_data[0])
        part2 = v_n.split(v_ns_in_data[1])
        if v_ns_in_data[0] + '_' == part1[0]:
            v_n1 = v_ns_in_data[0]
            v_n2 = v_ns_in_data[1]
        elif '_' + v_ns_in_data[0] == part1[1]:
            v_n1 = v_ns_in_data[1]
            v_n2 = v_ns_in_data[0]
        elif v_ns_in_data[1] + '_' == part1[0]:
            v_n1 = v_ns_in_data[1]
            v_n2 = v_ns_in_data[0]
        elif '_' + v_ns_in_data[1] == part1[1]:
            v_n1 = v_ns_in_data[0]
            v_n2 = v_ns_in_data[1]
        else:
            print("v_n: {}".format(v_n))
            print("v_n_in_data: {}".format(v_ns_in_data))
            print("v_n.split({}): {}".format(v_ns_in_data[0], part1))
            print("v_n.split({}): {}".format(v_ns_in_data[1], part2))
            print("v_ns_in_data[0]: {}".format(v_ns_in_data[0]))
            print("v_ns_in_data[1]: {}".format(v_ns_in_data[1]))
            raise NameError("Get simpler for f*ck sake...")

        print("v_n1: {}".format(v_n1))
        print("v_n2: {}".format(v_n2))
        edge_x = np.array(dfile[v_n1])
        edge_y = np.array(dfile[v_n2])
        mass = np.array(dfile["mass"]).T

        arr_x = 0.5 * (edge_x[1:] + edge_x[:-1])
        arr_y = 0.5 * (edge_y[1:] + edge_y[:-1])

        result = combine(arr_x, arr_y, mass)

        self.corr_matrix[self.i_it(it)][self.i_v_n(v_n)] = result

    def load_corr_file(self, it, v_n_x, v_n_y):

        v_n_x = str(v_n_x)
        v_n_y = str(v_n_y)

        self.check_v_n(v_n_x)
        self.check_v_n(v_n_y)

        # check if the direct file exists or the inverse
        fpath_direct = Paths.ppr_sims + self.sim + "/res_3d/" + str(it) + "/corr_" + v_n_x + '_' + v_n_y + ".h5"
        fpath_inverse = Paths.ppr_sims + self.sim + "/res_3d/" + str(it) + "/corr_" + v_n_y + '_' + v_n_x + ".h5"
        if os.path.isfile(fpath_direct):
            fpath = fpath_direct
        elif os.path.isfile(fpath_inverse):
            fpath = fpath_inverse
        else:
            raise IOError("Correlation files not found:\n{}\n nor \n{}".format(fpath_direct, fpath_inverse))

        # check if the data inside is in the right format
        dfile = h5py.File(fpath, "r")
        v_ns_in_data = []
        for v_n_ in dfile:
            v_ns_in_data.append(v_n_)
        if not "mass" in v_ns_in_data:
            raise NameError("mass is not found in file:{}".format(fpath))
        if len(v_ns_in_data) > 3:
            raise NameError("More than 3 datasets found in corr file: {}".format(fpath))

        # extract edges and convert them into center of bins
        edge_x = np.array(dfile[v_n_x])
        edge_y = np.array(dfile[v_n_y])
        arr_x = 0.5 * (edge_x[1:] + edge_x[:-1])  # from edges to center of bins
        arr_y = 0.5 * (edge_y[1:] + edge_y[:-1])

        # extract mass (weights)
        if fpath == fpath_direct:
            mass = np.array(dfile["mass"]).T
        else:
            mass = np.array(dfile["mass"])

        # create a 2D table of the data (convenient format)
        result = combine(arr_x, arr_y, mass)
        self.corr_matrix[self.i_it(it)][self.i_v_n(v_n_x, v_n_y)] = result

        # fpath = Paths.ppr_sims + self.sim + "/res_3d/" + str(it) + "/corr_" + v_n_x + '_' + v_n_y + ".h5"
        #
        # if os.path.isfile(fpath):
        #     dfile = h5py.File(fpath, "r")
        #     v_ns_in_data = []
        #     for v_n_ in dfile:
        #         v_ns_in_data.append(v_n_)
        #     if not "mass" in v_ns_in_data:
        #         raise NameError("mass is not found in file:{}".format(fpath))
        #     if len(v_ns_in_data) > 3:
        #         raise NameError("More than 3 datasets found in corr file: {}".format(fpath))
        #     edge_x = np.array(dfile[v_n_x])
        #     edge_y = np.array(dfile[v_n_y])
        #     mass = np.array(dfile["mass"]).T
        #
        # if not os.path.isfile(fpath):
        #     print("Correlation file not found:\n{}".format(fpath))
        #     fpath_in = Paths.ppr_sims + self.sim + "/res_3d/" + str(it) + "/corr_" + v_n_y + '_' + v_n_x + ".h5"
        #     print("Loading inverse file:\n{}".format(fpath_in))
        #
        #     dfile = h5py.File(fpath_in, "r")
        #     v_ns_in_data = []
        #     for v_n_ in dfile:
        #         v_ns_in_data.append(v_n_)
        #     if not "mass" in v_ns_in_data:
        #         raise NameError("mass is not found in file:{}".format(fpath))
        #     if len(v_ns_in_data) > 3:
        #         raise NameError("More than 3 datasets found in corr file: {}".format(fpath))
        #     edge_x = np.array(dfile[v_n_x])
        #     edge_y = np.array(dfile[v_n_y])
        #     mass = np.array(dfile["mass"]).T
        #
        #     if not os.path.isfile(fpath_in):
        #         raise IOError("Correlation files not found:\n{}\n or \n{}".format(fpath, fpath_in))
        #
        #
        #
        # dfile = h5py.File(fpath, "r")
        #
        # v_ns_in_data = []
        # for v_n_ in dfile:
        #     v_ns_in_data.append(v_n_)
        #
        # if not "mass" in v_ns_in_data:
        #     raise NameError("mass is not found in file:{}".format(fpath))
        #
        # if len(v_ns_in_data) > 3:
        #     raise NameError("More than 3 datasets found in corr file: {}".format(fpath))
        #
        # v_ns_in_data.remove("mass")
        #
        # # for v_n__ in v_ns_in_data:
        # #     if not v_n__ in v_n:
        # #         raise NameError("in_data_v_n: {} is not in corr name v_n: {}"
        # #                         .format(v_n__, v_n))
        #
        #
        # # part1 = v_n.split(v_ns_in_data[0])
        # # part2 = v_n.split(v_ns_in_data[1])
        # # if v_ns_in_data[0] + '_' == part1[0]:
        # #     v_n1 = v_ns_in_data[0]
        # #     v_n2 = v_ns_in_data[1]
        # # elif '_' + v_ns_in_data[0] == part1[1]:
        # #     v_n1 = v_ns_in_data[1]
        # #     v_n2 = v_ns_in_data[0]
        # # elif v_ns_in_data[1] + '_' == part1[0]:
        # #     v_n1 = v_ns_in_data[1]
        # #     v_n2 = v_ns_in_data[0]
        # # elif '_' + v_ns_in_data[1] == part1[1]:
        # #     v_n1 = v_ns_in_data[0]
        # #     v_n2 = v_ns_in_data[1]
        # # else:
        # #     print("v_n: {}".format(v_n))
        # #     print("v_n_in_data: {}".format(v_ns_in_data))
        # #     print("v_n.split({}): {}".format(v_ns_in_data[0], part1))
        # #     print("v_n.split({}): {}".format(v_ns_in_data[1], part2))
        # #     print("v_ns_in_data[0]: {}".format(v_ns_in_data[0]))
        # #     print("v_ns_in_data[1]: {}".format(v_ns_in_data[1]))
        # #     raise NameError("Get simpler for f*ck sake...")
        # #
        # # print("v_n1: {}".format(v_n1))
        # # print("v_n2: {}".format(v_n2))
        # edge_x = np.array(dfile[v_n_x])
        # edge_y = np.array(dfile[v_n_y])
        # mass = np.array(dfile["mass"]).T
        #
        # arr_x = 0.5 * (edge_x[1:] + edge_x[:-1]) # from edges to center of bins
        # arr_y = 0.5 * (edge_y[1:] + edge_y[:-1])
        #
        # result = combine(arr_x, arr_y, mass)
        #
        # self.corr_matrix[self.i_it(it)][self.i_v_n(v_n_x, v_n_y)] = result

    def is_corr_loaded(self, it, v_n_x, v_n_y):

        if len(self.corr_matrix[self.i_it(it)]) < self.i_v_n(v_n_x, v_n_y):
            raise ValueError("{} < {}".format(len(self.corr_matrix[self.i_it(it)]), self.i_v_n(v_n_x, v_n_y)))

        corr = self.corr_matrix[self.i_it(it)][self.i_v_n(v_n_x, v_n_y)]
        if len(corr) == 0:
            self.load_corr_file(it, v_n_x, v_n_y)
        else:
            Printcolor.yellow("Warning. Rewriting loaded data: v_n_x:{} v_n_y:{}, it:{}"
                              .format(v_n_x, v_n_y, it))

    def get_res_corr(self, it, v_n_x, v_n_y):
        self.check_v_n(v_n_x)
        self.check_v_n(v_n_y)
        self.check_it(it)
        self.is_corr_loaded(it, v_n_x, v_n_y)
        return self.corr_matrix[self.i_it(it)][self.i_v_n(v_n_x, v_n_y)]

    def get_time(self, it):
        self.check_it(it)
        return self.times[self.list_iterations.index(it)]


class RES(LOAD_RES_CORR):

    def __init__(self, sim):
        LOAD_RES_CORR.__init__(self, sim)





class PLOT_TASK(RES):

    def __init__(self, sim):
        RES.__init__(self, sim)

    def plot_correlation(self, ax, dic):

        table = self.get_res_corr(dic["it"], dic["v_n_x"], dic["v_n_y"])

        table = np.array(table)
        # table[0, 1:] = table[0, 1:] * 6.176269145886162e+17

        x_arr = table[0,1:] #* 6.176269145886162e+17
        y_arr = table[1:,0]
        z_arr = table[1:,1:]

        # normalization
        z_arr = z_arr/np.sum(z_arr)

        # special treatment
        if dic["v_n_x"] == "theta": x_arr = 90 - (180 * x_arr / np.pi)
        if dic["v_n_y"] == "theta": y_arr = 90 - (180 * y_arr / np.pi)

        if dic["v_n_x"] == "rho": x_arr *= 6.176269145886162e+17
        if dic["v_n_y"] == "rho": y_arr *= 6.176269145886162e+17



        # limits
        if dic["xmin"] != None and dic["xmax"] != None:
            ax.set_xlim(dic["xmin"], dic["xmax"])

        if dic["ymin"] != None and dic["ymax"] != None:
            ax.set_ylim(dic["ymin"], dic["ymax"])

        if dic["vmin"] == None: dic["vmin"] = z_arr.min()
        if dic["vmax"] == None: dic["vmax"] = z_arr.max()

        if dic["norm"] == "norm":
            norm = Normalize(vmin=dic["vmin"], vmax=dic["vmax"])
        elif dic["norm"] == "log":
            norm = LogNorm(vmin=dic["vmin"], vmax=dic["vmax"])
        else:
            raise NameError("unrecognized norm: {} in task {}"
                            .format(dic["norm"], dic["v_n"]))

        ax.set_xlabel(dic["v_n_x"].replace('_', '\_'))
        ax.set_ylabel(dic["v_n_y"].replace('_', '\_'))

        if dic["xscale"] == 'log':
            ax.set_xscale("log")
        if dic["yscale"] == 'log':
            ax.set_yscale("log")

        im = ax.pcolormesh(x_arr, y_arr, z_arr, norm=norm, cmap=dic["cmap"])
        im.set_rasterized(True)

        return im


class PLOT_MANY_TASKS(PLOT_TASK):

    def __init__(self, sim):
        RES.__init__(self, sim)

        it = 1818738

        self.gen_set = {
            "figdir": Paths.ppr_sims + self.sim + "/res_3d/{}/".format(it),
            "figname": "inv_ang_mom_flux.png",
            # "figsize": (13.5, 3.5), # <->, |
            "figsize": (3.8, 3.5),  # <->, |
            "type": "cartesian",
            "subplots_adjust_h": 0.2,
            "subplots_adjust_w": 0.3
        }

        self.set_plot_dics = []

        corr_dic_temp_Ye = { # relies on the "get_res_corr(self, it, v_n): " method of data object
            'name': 'corr', 'position': (1, 1), 'title': 'time [ms]', 'cbar': 'right .05 .0',
            'it': 2237972, 'v_n_x': 'temp', 'v_n_y': 'Ye', 'v_n': 'mass',
            'xmin': 2., 'xmax': 15., 'ymin': 0., 'ymax': 0.2, 'vmin': 1e-4, 'vmax': None,
            'mask_below': None, 'mask_above': None, 'cmap': 'inferno_r', 'norm': 'log', 'todo': None
        }
        # self.set_plot_dics.append(corr_dic_temp_Ye)

        corr_dic_rho_ang_mom = { # relies on the "get_res_corr(self, it, v_n): " method of data object
            'name': 'corr', 'position': (1, 1), 'title': 'time [ms]', 'cbar':  'right .05 .0',
            'it': 761856, 'v_n_x': 'rho', 'v_n_y': 'ang_mom', 'v_n': 'mass',
            'xmin': None, 'xmax': None, 'ymin': None, 'ymax': None, 'vmin': 1e-8, 'vmax': None,
            'xscale':'log', 'yscale':'log',
            'mask_below': None, 'mask_above': None, 'cmap': 'inferno_r', 'norm': 'log', 'todo': None
        }
        # self.set_plot_dics.append(corr_dic_rho_ang_mom)

        corr_dic_rho_dens_unb_bern = { # relies on the "get_res_corr(self, it, v_n): " method of data object
            'name': 'corr', 'position': (1, 1), 'title': 'time [ms]', 'cbar':  None, #  'right .05 .0',
            'it': 1081344, 'v_n_x': 'rho', 'v_n_y': 'dens_unb_bern', 'v_n': 'mass',
            'xmin': None, 'xmax': None, 'ymin': None, 'ymax': 1e-3, 'vmin': 1e-7, 'vmax': None,
            'xscale':'log', 'yscale':'log',
            'mask_below': None, 'mask_above': None, 'cmap': 'inferno_r', 'norm': 'log', 'todo': None
        }
        # self.set_plot_dics.append(corr_dic_rho_dens_unb_bern)

        corr_dic_ang_mom_flux_theta = { # relies on the "get_res_corr(self, it, v_n): " method of data object
            'name': 'corr', 'position': (1, 3), 'title': 'time [ms]', 'cbar':  None, #'right .05 .0',
            'it': it, 'v_n_x': 'theta', 'v_n_y': 'ang_mom_flux', 'v_n': 'mass',
            'xmin': None, 'xmax': None, 'ymin': None, 'ymax': None, 'vmin': 1e-7, 'vmax': None,
            'xscale':'line', 'yscale':'log',
            'mask_below': None, 'mask_above': None, 'cmap': 'inferno_r', 'norm': 'log', 'todo': None
        }
        # self.set_plot_dics.append(corr_dic_ang_mom_flux_theta)

        corr_dic_ang_mom_flux_dens_unb_bern = { # relies on the "get_res_corr(self, it, v_n): " method of data object
            'name': 'corr', 'position': (1, 1), 'title': 'time [ms]', 'cbar': 'right .03 .0',
            'it': it, 'v_n_x': 'dens_unb_bern', 'v_n_y': 'ang_mom_flux', 'v_n': 'mass',
            'xmin': 1e-11, 'xmax': 1e-7, 'ymin': 1e-11, 'ymax': 1e-7, 'vmin': 1e-7, 'vmax': None,
            'xscale':'log', 'yscale':'log',
            'mask_below': None, 'mask_above': None, 'cmap': 'inferno_r', 'norm': 'log', 'todo': None
        }
        # self.set_plot_dics.append(corr_dic_ang_mom_flux_dens_unb_bern)

        corr_dic_inv_ang_mom_flux_dens_unb_bern = { # relies on the "get_res_corr(self, it, v_n): " method of data object
            'name': 'corr', 'position': (1, 1), 'title': 'time [ms]', 'cbar': 'right .03 .0',
            'it': it, 'v_n_x': 'dens_unb_bern', 'v_n_y': 'inv_ang_mom_flux', 'v_n': 'mass',
            'xmin': 1e-11, 'xmax': 1e-7, 'ymin': 1e-11, 'ymax': 1e-7, 'vmin': 1e-7, 'vmax': None,
            'xscale':'log', 'yscale':'log',
            'mask_below': None, 'mask_above': None, 'cmap': 'inferno_r', 'norm': 'log', 'todo': None
        }
        self.set_plot_dics.append(corr_dic_inv_ang_mom_flux_dens_unb_bern)

        corr_dic_rho_ang_mom_flux = { # relies on the "get_res_corr(self, it, v_n): " method of data object
            'name': 'corr', 'position': (1, 2), 'title': 'time [ms]', 'cbar':  'right .05 .0',
            'it': it, 'v_n_x': 'rho', 'v_n_y': 'ang_mom_flux', 'v_n': 'mass',
            'xmin': None, 'xmax': None, 'ymin': None, 'ymax': 1e-3, 'vmin': 1e-7, 'vmax': None,
            'xscale':'log', 'yscale':'log',
            'mask_below': None, 'mask_above': None, 'cmap': 'inferno_r', 'norm': 'log', 'todo': None
        }
        # self.set_plot_dics.append(corr_dic_rho_ang_mom_flux)

        corr_dic_rho_theta = { # relies on the "get_res_corr(self, it, v_n): " method of data object
            'name': 'corr', 'position': (1, 2), 'title': 'time [ms]', 'cbar':None, #   'right .05 .0',
            'it': it, 'v_n_x': 'theta', 'v_n_y': 'rho', 'v_n': 'mass',
            'xmin': None, 'xmax': None, 'ymin': None, 'ymax': None, 'vmin': 1e-7, 'vmax': None,
            'xscale': 'line', 'yscale': 'log',
            'mask_below': None, 'mask_above': None, 'cmap': 'inferno_r', 'norm': 'log', 'todo': None
        }
        # self.set_plot_dics.append(corr_dic_rho_theta)

        corr_dic_rho_r = { # relies on the "get_res_corr(self, it, v_n): " method of data object
            'name': 'corr', 'position': (1, 1), 'title': 'time [ms]', 'cbar': None,#  'right .05 .0',
            'it': it, 'v_n_x': 'r', 'v_n_y': 'rho', 'v_n': 'mass',
            'xmin': None, 'xmax': None, 'ymin': None, 'ymax': None, 'vmin': 1e-7, 'vmax': None,
            'xscale': 'line', 'yscale': 'log',
            'mask_below': None, 'mask_above': None, 'cmap': 'inferno_r', 'norm': 'log', 'todo': None
        }
        # self.set_plot_dics.append(corr_dic_rho_r)

        corr_dic_rho_Ye = { # relies on the "get_res_corr(self, it, v_n): " method of data object
            'name': 'corr', 'position': (1, 2), 'title': 'time [ms]', 'cbar': 'right .05 .0',
            'it': it, 'v_n_x': 'rho', 'v_n_y': 'Ye', 'v_n': 'mass',
            'xmin': None, 'xmax': None, 'ymin': 0., 'ymax': 0.4, 'vmin': 1e-8, 'vmax': None,
            'mask_below': None, 'mask_above': None, 'cmap': 'inferno_r', 'norm': 'log', 'todo': None
        }
        # self.set_plot_dics.append(corr_dic_rho_Ye)

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
                                                                      )#sharex=self.sbplot_matrix[n_col][0])
                    elif n_col > 0 and n_row == 0:
                        sbplot_matrix[n_col][n_row] = fig.add_subplot(self.n_rows, self.n_cols, i,
                                                                      )#sharey=self.sbplot_matrix[0][n_row])
                    else:
                        sbplot_matrix[n_col][n_row] = fig.add_subplot(self.n_rows, self.n_cols, i,
                                                                      #sharex=self.sbplot_matrix[n_col][0],
                                                                      )#sharey=self.sbplot_matrix[0][n_row])

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

    def plot_one_task(self, ax, dic):

        if dic["name"] == "corr":
            im = self.plot_correlation(ax, dic)
        else:
            raise NameError("name:{} is not recognized"
                            .format(dic["name"]))
        # self.time
        return im

    def set_plot_title(self, ax, plot_dic):
        if plot_dic["title"] != '' and plot_dic["title"] != None:

            if plot_dic["title"] == 'it':
                title = plot_dic["it"]
            elif plot_dic["title"] == 'time [s]' or \
                plot_dic["title"] == 'time':
                title = "%.3f" % self.get_time(plot_dic["it"]) + " [s]"
            elif plot_dic["title"] == 'time [ms]':
                title = "%.1f" % (self.get_time(plot_dic["it"]) * 1000) + " [ms]"
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
                    dic = dict(dic)
                    im = self.plot_one_task(ax, dic)
                    self.set_plot_title(ax, dic)
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
                cbar = plt.colorbar(im, cax=cax1, extend='both')#, format='%.1e')
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

""" --- --- LOADING & PLOTTING R-phi Phi-Z --- ---"""

class LOAD_INT_DATA:

    def __init__(self, sim):

        self.sim = sim

        # getting list of iterations ( for storage capacity reasons)
        iterations = []
        fname = '*.h5'
        root = Paths.gw170817 + sim + '/profiles/3d/'
        files = locate(fname, root=root, followlinks=False)
        if len(files) == 0:
            raise NameError("No iterations found in the root:{}".format(root))
        if os.path.isfile(root+"ittime.txt"):
            list_iterations, self.times = np.loadtxt(root+"ittime.txt", usecols=(0, 1), unpack=True)
            self.list_iterations = list(list_iterations)
        else:
            for file_ in files:
                iterations.append(int(str(file_.split('/')[-1]).split('.h5')[0]))
            self.times = interpoate_time_form_it(iterations, Paths.gw170817+sim+'/')
            self.list_iterations = iterations
        if len(self.list_iterations) == 0:
            raise IOError("No iterations found")

        self.grid_type = "cyl"
        self.list_grid_v_ns = ["x_cyl", "y_cyl", "z_cyl",
                          "r_cyl", "phi_cyl",
                          "dr_cyl", "dphi_cyl", "dz_cyl"]

        self.grid_data_matrix = [[np.zeros(0,)
                                 for x in range(len(self.list_grid_v_ns))]
                                 for y in range(len(self.list_iterations))]


        self.list_of_v_ns = ["ang_mom", "ang_mom_flux", "density", "dens_unb_geo",
                             "dens_unb_bern","rho", "temp", "Ye"]

        self.data_matrix = [[np.zeros(0,)
                                 for x in range(len(self.list_of_v_ns))]
                                 for y in range(len(self.list_iterations))]

    def check_grid_v_n(self, v_n):
        if not v_n in self.list_grid_v_ns:
            raise NameError("v_n:{} not in list of grid v_ns\n{}"
                            .format(v_n, self.list_grid_v_ns))

    def check_data_v_n(self, v_n):
        if not v_n in self.list_of_v_ns:
            raise NameError("v_n:{} not in list of data v_ns\n{}"
                            .format(v_n, self.list_of_v_ns))

    def check_it(self, it):
        if not it in self.list_iterations:
            raise NameError("it:{} not in the list of iterations\n{}"
                            .format(it, self.list_iterations))

    def i_data_v_n(self, v_n):
        self.check_data_v_n(v_n)
        return int(self.list_of_v_ns.index(v_n))

    def i_grid_v_n(self, v_n):
        self.check_grid_v_n(v_n)
        return int(self.list_grid_v_ns.index(v_n))

    def i_it(self, it):
        self.check_it(it)
        return int(self.list_iterations.index(it))

    def load_grid(self, it):

        path = Paths.ppr_sims + self.sim + "/res_3d/" + str(int(it)) + '/'
        fname = path + self.grid_type + '_grid.h5'
        grid_file = h5py.File(fname, "r")

        for v_n in self.list_grid_v_ns:
            if v_n not in grid_file:
                raise NameError("Loaded grid file {} does not have v_n:{} Expected only:\n{}"
                                .format(fname, v_n, self.list_grid_v_ns))

            grid_data = np.array(grid_file[v_n], dtype=np.float)
            self.grid_data_matrix[self.i_it(it)][self.i_grid_v_n(v_n)] = grid_data

    def load_data(self, it, v_n):

        path = Paths.ppr_sims + self.sim + "/res_3d/" + str(int(it)) + '/'
        fname = path + self.grid_type + '_' + v_n + ".h5"
        data_file = h5py.File(fname, "r")

        if len(data_file) > 1:
            raise IOError("More than one v_n is found in data_file: {}".format(fname))
        if len(data_file) == 0:
            raise IOError("No datasets found in data_file: {}".format(fname))

        for v_n_ in data_file:
            if v_n_ != v_n:
                raise NameError("required v_n:{} not the same as the one in datafile:{}"
                                .format(v_n, v_n_))

        data = np.array(data_file[v_n], dtype=np.float)

        self.data_matrix[self.i_it(it)][self.i_data_v_n(v_n)] = data

    def is_grid_loaded(self, it):

        if len(self.grid_data_matrix[self.i_it(it)][self.i_grid_v_n(self.list_grid_v_ns[0])]) == 0:
            self.load_grid(it)

    def is_data_loaded(self, it, v_n):

        if len(self.data_matrix[self.i_it(it)][self.i_data_v_n(v_n)]) == 0:
            self.load_data(it, v_n)


    def get_grid_data(self, it, v_n):
        self.check_it(it)
        self.check_grid_v_n(v_n)

        self.is_grid_loaded(it)

        return self.grid_data_matrix[self.i_it(it)][self.i_grid_v_n(v_n)]

    def get_data(self, it, v_n):
        self.check_it(it)
        self.check_data_v_n(v_n)

        self.is_data_loaded(it, v_n)

        return self.data_matrix[self.i_it(it)][self.i_data_v_n(v_n)]





""" --- --- TASK SPECIFIC FUNCTIONS --- --- --- """

def do_histogram_processing_of_iterations():

    sim = "DD2_M13641364_M0_SR"
    profs_loc = "/data/numrel/WhiskyTHC/Backup/2018/GW170817/DD2_M13641364_M0_SR/profiles/3d/"

    list_iterations, times = np.loadtxt(profs_loc+"ittime.txt", usecols=(0, 1), unpack=True)
    iterations = list_iterations[times>0.030]
    # print(iterations)
    # exit(1)
    processed = []

    for it in np.array(iterations, dtype=int):
        print("| processing iteration: {} ({} out {})|".format(it, len(processed), len(iterations)))
        o_methods = MAINMETHODS_STORE(profs_loc+"{}.h5".format(it), sim)
        #
        o_methods.get_total_mass(save=True)
        print(np.sum(o_methods.get_correlation(o_methods.corr_task_dic_rho_ye, save=True)))
        print(np.sum(o_methods.get_correlation(o_methods.corr_task_dic_temp_ye, save=True)))
        print(np.sum(o_methods.get_correlation(o_methods.corr_task_dic_rho_r, save=True)))
        print(np.sum(o_methods.get_correlation(o_methods.corr_task_dic_rho_theta, save=True)))
        print(np.sum(o_methods.get_correlation(o_methods.corr_task_dic_rho_ang_mom, save=True)))
        print(np.sum(o_methods.get_correlation(o_methods.corr_task_dic_rho_ang_mom_flux, save=True)))
        print(np.sum(o_methods.get_correlation(o_methods.corr_task_dic_rho_dens_unb_bern, save=True)))
        print(np.sum(o_methods.get_correlation(o_methods.corr_task_dic_ang_mom_flux_theta, save=True)))
        print(np.sum(o_methods.get_correlation(o_methods.corr_task_dic_ang_mom_flux_dens_unb_bern, save=True)))
        print(np.sum(o_methods.get_correlation(o_methods.corr_task_dic_inv_ang_mom_flux_dens_unb_bern, save=True)))
        o_methods.__delete__(o_methods)
        processed.append(it)
    exit(1)


if __name__ == "__main__":

    do_histogram_processing_of_iterations()




    # times, files = get_profiles("DD2_M13641364_M0_LR_R04",
    #                             time_list=[], it_list=[], n_more=0, ftype='.h5', time_units='s', add_path='profiles/3d/')
    # exit(1)


    # --- DEBUGGING ---
    # o_prof = LOAD_PROFILE("/data1/numrel/WhiskyTHC/Backup/2018/GW170817/SLy4_M13641364_M0_SR/profiles/434176.h5")

    # o_data = COMPUTE_STORE("/data1/numrel/WhiskyTHC/Backup/2018/GW170817/SLy4_M13641364_M0_SR/profiles/434176.h5")
    # print(o_data.get_prof_arr(0, "rho").shape); exit(1) # (262, 262, 134)
    # print(o_data.get_comp_data(2, "dens_unb_bern"))

    # o_mask = MASK_STORE("/data1/numrel/WhiskyTHC/Backup/2018/GW170817/SLy4_M13641364_M0_SR/profiles/434176.h5")
    # print(o_mask.get_masked_data(2, "rho"))
    ''' MAIN '''
    sim = "DD2_M13641364_M0_LK_SR_R04"
    it = 1111116

    # sim = "DD2_M13641364_M0_LK_SR_R04"; times: 0.084 0.073 0.062 0.051
    # sim = "DD2_M13641364_M0_LK_HR_R04"; times: 0.042 0.038 0.031
    # sim = "DD2_M13641364_M0_SR_R04";    times: 0.098, 0.094, 0.083


    #  ang_mom is in (-8.03540793958e-09->0.000251518721152)
    # DD2:  DD2_M13641364_M0_SR/profiles/3d/2025090.h5
    # SLy4: SLy4_M13641364_M0_SR/profiles/3d/761856.h5
    # LS220:LS220_M13641364_M0_SR/profiles/3d/1081344.h5
    #

    o_methods = MAINMETHODS_STORE("/data1/numrel/WhiskyTHC/Backup/2018/GW170817/"
                                  "{}/profiles/3d/{}.h5".format(sim, it), sim)
    #
    o_methods.get_total_mass(save=True)
    print(np.sum(o_methods.get_correlation(o_methods.corr_task_dic_rho_ye, save=True)))
    print(np.sum(o_methods.get_correlation(o_methods.corr_task_dic_temp_ye, save=True)))
    print(np.sum(o_methods.get_correlation(o_methods.corr_task_dic_rho_r, save=True)))
    print(np.sum(o_methods.get_correlation(o_methods.corr_task_dic_rho_theta, save=True)))
    print(np.sum(o_methods.get_correlation(o_methods.corr_task_dic_rho_ang_mom, save=True)))
    print(np.sum(o_methods.get_correlation(o_methods.corr_task_dic_rho_ang_mom_flux, save=True)))
    print(np.sum(o_methods.get_correlation(o_methods.corr_task_dic_rho_dens_unb_bern, save=True)))
    print(np.sum(o_methods.get_correlation(o_methods.corr_task_dic_ang_mom_flux_theta, save=True)))
    print(np.sum(o_methods.get_correlation(o_methods.corr_task_dic_ang_mom_flux_dens_unb_bern, save=True)))
    print(np.sum(o_methods.get_correlation(o_methods.corr_task_dic_inv_ang_mom_flux_dens_unb_bern, save=True)))
    exit(1)
    ''' PLOTTING '''
    # o_plot = PLOT_MANY_TASKS(sim)
    # o_plot.main()
    # exit(1)
    ''' TESTING '''
    _, _, task_for_int = setup()
    # int_ = INTERPOLATE_STORE("/data1/numrel/WhiskyTHC/Backup/2018/GW170817/"
    #                               "SLy4_M13641364_M0_SR/profiles/3d/761856.h5",
    #                          "SLy4_M13641364_M0_SR",
    #                          CYLINDRICAL_GRID(task_for_int["grid"]))
    # int_.get_int("x_cyl")
    # print(int_.get_int("rho"))

    # int_ = INTMETHODS_STORE("/data1/numrel/WhiskyTHC/Backup/2018/GW170817/"
    #                               "DD2_M13641364_M0_LK_SR_R04/profiles/3d/1818738.h5",
    #                          "DD2_M13641364_M0_LK_SR_R04",
    #                          CYLINDRICAL_GRID(task_for_int["grid"]))
    # int_.save_new_grid()
    # int_.save_int_v_n("ang_mom_flux")
    # int_.save_int_v_n("dens_unb_bern")

    # load interpolated data for plotting and processing
    load_int = LOAD_INT_DATA("DD2_M13641364_M0_LK_SR_R04")
    print(load_int.get_grid_data(1818738, "r_cyl"))




    ''''''

    # for rl in range(o_methods.nlevels):
    #     data = o_methods.get_masked_data(rl, "dens_unb_bern")
    #     print("rl:{} min:{} max:{}".format(rl, data.min(), data.max()))
    # exit(1)

    # o_res = LOAD_RES_CORR("DD2_M13641364_M0_SR")
    # print(np.sum(o_res.get_res_corr(2237972, "temp_Ye")[1:,1:]))# should yeid the disk mass


    # o_plot = PLOT_MANY_TASKS("SLy4_M13641364_M0_SR")
    # o_plot.main()

    # o_plot = PLOT_MANY()