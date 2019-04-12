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
        "grid": {'type': 'cylindrical', 'n_r': 150, 'n_phi': 150, 'n_z': -100},
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
            "theta"
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
                              self.get_comp_data(rl, "vup"))

        elif v_n == 'theta':
            arr = FORMULAS.theta(self.get_comp_data(rl, "x"),
                                 self.get_comp_data(rl, "y"),
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
                raise NameError("v_n is not recognized: {}".format(v_n))

    def get_comp_data(self, rl, v_n):
        self.check_v_n(v_n)
        self.is_available(rl, v_n)

        return self.data_matrix[self.i_v_n(v_n)][rl]


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


class MAINMETHODS_STORE(MASK_STORE):

    def __init__(self, fname, sim):

        MASK_STORE.__init__(self, fname)

        self.sim = sim

        rho_const = 6.176269145886162e+17
        self.task_correlation = [
            # {"v_n": "rho",  "edges": 10.0 ** np.linspace(4.0, 16.0, 500) / rho_const},  # not in CGS :^
            {"v_n": "temp", "edges": 10.0 ** np.linspace(-2, 2, 300)},
            {"v_n": "Ye",   "edges": np.linspace(0, 0.5, 300)}
        ]

    def get_total_mass(self, multiplier=2., save=False):
        mass = 0.
        for rl in range(self.nlevels):
            density = np.array(self.get_masked_data(rl, "density"))
            delta = self.get_prof_delta(rl)
            mass += float(multiplier * np.sum(density) * np.prod(delta))

        print("it:{} mass:{:3f}Msun".format(self.iteration, mass))

        if save:
            path = Paths.ppr_sims + self.sim + "/res_3d/"
            fname = "{}_disk_mass.txt".format(self.iteration)
            np.savetxt(path + fname, np.array([mass]), fmt='%.5f')

        return mass

    def get_correlation(self, multiplier=2., save=False):

        v_ns = []
        edges = []
        for setup_dictionary in self.task_correlation:
            v_ns.append(setup_dictionary["v_n"])
            edges.append(setup_dictionary["edges"])
        edges = tuple(edges)

        correlation = np.zeros([len(edge) - 1 for edge in edges])
        for rl in range(self.nlevels):
            data = []
            weights = self.get_masked_data(rl, "density").flatten() * \
                      np.prod(self.get_prof_delta(rl)) * multiplier
            for i_vn, v_n in enumerate(v_ns):
                data.append(self.get_masked_data(rl, v_n).flatten())

            data = tuple(data)
            tmp, _ = np.histogramdd(data, bins=edges, weights=weights)
            correlation += tmp

        if save:
            path = Paths.ppr_sims + self.sim + "/res_3d/"
            fname = "{}_corr_".format(self.iteration)
            for v_n in v_ns:
                fname += v_n
                if v_n != v_ns[-1]:
                    fname += '_'
            fname += '.h5'

            outfile = h5py.File(path + fname, "w")
            for v_n, edge in zip(v_ns, edges):
                outfile.create_dataset(v_n, data=edge)
            outfile.create_dataset("mass", data=correlation)
            outfile.close()

        return correlation


class INTERPOLATE_STORE(MAINMETHODS_STORE):

    def __init__(self, fname, sim):
        MAINMETHODS_STORE.__init__(self, fname, sim)

if __name__ == "__main__":

    # --- DEBUGGING ---
    # o_prof = LOAD_PROFILE("/data1/numrel/WhiskyTHC/Backup/2018/GW170817/SLy4_M13641364_M0_SR/profiles/434176.h5")

    # o_data = COMPUTE_STORE("/data1/numrel/WhiskyTHC/Backup/2018/GW170817/SLy4_M13641364_M0_SR/profiles/434176.h5")
    # print(o_data.get_prof_arr(0, "rho").shape); exit(1) # (262, 262, 134)
    # print(o_data.get_comp_data(2, "dens_unb_bern"))

    # o_mask = MASK_STORE("/data1/numrel/WhiskyTHC/Backup/2018/GW170817/SLy4_M13641364_M0_SR/profiles/434176.h5")
    # print(o_mask.get_masked_data(2, "rho"))

    o_methods = MAINMETHODS_STORE("/data1/numrel/WhiskyTHC/Backup/2018/GW170817/"
                                  "DD2_M13641364_M0_SR/profiles/3d_all/2237972.h5", "DD2_M13641364_M0_SR")
    o_methods.get_total_mass(save=True)
    o_methods.get_correlation(save=True)
    # corr = o_methods.get_correlation()
    # SAVE_RESULT.correlation()
    #
    #
    # print(o_methods.get_total_mass())