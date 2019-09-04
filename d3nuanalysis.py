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
import matplotlib
matplotlib.use("Agg")
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
import click
import h5py
import csv
import os
import gc
# from visit_utils import *
from scidata.utils import locate
import scidata.carpet.hdf5 as h5
import scidata.xgraph as xg
from scidata.carpet.interp import Interpolator
import scivis.data.carpet2d
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

from preanalysis import LOAD_ITTIME

class LOAD_NU_PROFILE(LOAD_ITTIME):

    def __init__(self, sim, symmetry=None):

        LOAD_ITTIME.__init__(self, sim)

        self.nuprof_name = "nu" # -> 12345nu.h5

        self.symmetry = symmetry

        # self.rmax = 500

        self.profpath = Paths.gw170817 + sim + '/' + "profiles/3d/"

        isnuprofs, itnuprofs, timenuprofs = \
            self.get_ittime("nuprofiles", "nuprof")
        if not isnuprofs:
            is3ddata, it3d, t3d = self.get_ittime("overall", d1d2d3prof="d3")
            if is3ddata:
                raise IOError("ittime.h5 says there are NO nuprofiles, while there IS 3D data for times:\n{}"
                              "\n Extract nuprofiles before proceeding"
                              .format(t3d))
            else:
                raise IOError("ittime.h5 says there are no profiles, and no 3D data found.")

        self.list_iterations = list(itnuprofs)
        self.list_times = timenuprofs

        self.list_nuprof_v_ns = ['abs_energy', 'abs_nua', 'abs_nue', 'abs_number', 'eave_nua', 'eave_nue',
                                 'eave_nux', 'E_nua', 'E_nue', 'E_nux', 'flux_fac', 'ndens_nua', 'ndens_nue',
                                 'ndens_nux','N_nua', 'N_nue', 'N_nux']

        self.list_nugrid_v_ns = ["x", "y", "z", "r", "theta", "phi"]

        self.nudfile_matrix = [0 for it in range(len(self.list_iterations))]

        # self.nugrid_matrix = [0 for it in range(len(self.list_iterations))]

        self.nuprof_arr_matrix = [[np.zeros(0, )
                                   for v_n in range(len(self.list_nuprof_v_ns))]
                                   for it in range(len(self.list_iterations))]

        self.nuprof_grid_params_matrix = [[-1.
                                   for v_n in range(len(self.list_nugrid_v_ns))]
                                   for it in range(len(self.list_iterations))]

    def check_nuprof_v_n(self, v_n):
        if not v_n in self.list_nuprof_v_ns:
            raise NameError("v_n:{} not in list of nuprofile v_ns:{}"
                            .format(v_n, self.list_nuprof_v_ns))

    def check_it(self, it):
        if not int(it) in self.list_iterations:
            raise NameError("it:{} not in list of iterations:{}"
                            .format(it, self.list_iterations))

    def i_nu_it(self, it):
        return int(self.list_iterations.index(it))


    # --- ---

    def load_nudfile(self, it):
        fpath = self.profpath + str(it) + self.nuprof_name + ".h5"
        if not os.path.isfile(fpath):
            raise IOError("Expected file:{} NOT found"
                          .format(fpath))
        dfile = h5py.File(fpath, "r")
        self.nudfile_matrix[self.i_nu_it(it)] = dfile

    def is_nudfile_loaded(self, it):
        if isinstance(self.nudfile_matrix[self.i_nu_it(it)], int):
            self.load_nudfile(it)

    def get_nuprofile_dfile(self, it):
        self.check_it(it)
        self.is_nudfile_loaded(it)
        return self.nudfile_matrix[self.i_nu_it(it)]

        # self.symmetry = symmetry
        # self.nlevels = 7
        # self.profile = fname
        # self.dfile = h5py.File(fname, "r")
        # group_0 = self.dfile["reflevel={}".format(0)]
        # self.time = group_0.attrs["time"] * 0.004925794970773136 * 1e-3 # [sec]
        # self.iteration = group_0.attrs["iteration"]
        # print("\t\t symmetry: {}".format(self.symmetry))
        # print("\t\t time: {}".format(self.time))
        # print("\t\t iteration: {}".format(self.iteration))
        # self.grid = self.read_carpet_grid(self.dfile)
        #
        # # print("grid: {}".format(self.grid))
        #
        #
        #
        # if self.symmetry == "pi" and not str(self.profile).__contains__("_PI"):
        #     raise NameError("profile {} does not seem to have a pi symmetry. Check"
        #                     .format(self.profile))

    def i_nu_v_n(self, v_n):
        return int(self.list_nuprof_v_ns.index(v_n))

    def extract_arr_from_nuprof(self, it, v_n):
        nudfile = self.get_nuprofile_dfile(it)
        arr = np.array(nudfile[v_n])
        self.nuprof_arr_matrix[self.i_nu_it(it)][self.i_nu_v_n(v_n)] = arr

    def is_nuprofarr_extracted(self, it, v_n):
        if len(self.nuprof_arr_matrix[self.i_nu_it(it)][self.i_nu_v_n(v_n)]) == 0:
            self.extract_arr_from_nuprof(it, v_n)

    def get_nuprof_arr(self, it, v_n):

        self.check_nuprof_v_n(v_n)

        self.is_nuprofarr_extracted(it, v_n)

        return self.nuprof_arr_matrix[self.i_nu_it(it)][self.i_nu_v_n(v_n)]

    # grid

    def get_nrad(self, it):
        nudfile = self.get_nuprofile_dfile(it)
        return int(nudfile.attrs["nrad"])

    def get_nphi(self, it):
        nudfile = self.get_nuprofile_dfile(it)
        return int(nudfile.attrs["nphi"])

    def get_ntheta(self, it):
        nudfile = self.get_nuprofile_dfile(it)
        return int(nudfile.attrs["ntheta"])

    def get_sph_grid(self, it, nextra=0):

        rad, phi, theta = np.mgrid[0:self.get_nrad(it) + nextra, 0:self.get_nphi(it) + nextra, \
                          0:self.get_ntheta(it) + nextra].astype(np.float32)

        return rad, phi, theta

    def get_x_y_z_grid(self, it, plane=None, dual=False, rmax=50):

        if dual:
            nextra = 1
            shift  = -0.5
        else:
            nextra = 0
            shift  = 0.0

        nrad, nphi, ntheta = self.get_nrad(it), self.get_nphi(it), self.get_ntheta(it)


        if plane is None:
            rad, phi, theta = np.mgrid[0:nrad+nextra, 0:nphi+nextra,\
                    0:ntheta+nextra].astype(np.float32)
            rad = (rad + shift) * rmax/(nrad - 1)
            phi = (phi + shift) * (2*pi)/(nphi - 1)
            theta = (theta + shift) * pi/(ntheta - 1)
            x = rad * np.cos(phi) * np.sin(theta)
            y = rad * np.sin(phi) * np.sin(theta)
            z = rad * np.cos(theta)
            return x, y, z
        if plane == "xy":
            rad, phi = np.mgrid[0:nrad+nextra,\
                    0:nphi+nextra].astype(np.float32)
            rad = (rad + shift) * rmax/(nrad - 1)
            phi = (phi + shift) * (2*pi)/(nphi - 1)
            x = rad * np.cos(phi)
            y = rad * np.sin(phi)
            return x, y
        if plane == "xz" or plane == "yz":
            rad, theta = np.mgrid[0:nrad+nextra,\
                    0:2*ntheta+nextra-1].astype(np.float32)
            rad = (rad  + shift) * rmax/(nrad - 1)
            theta = (theta + shift) * pi/(ntheta - 1)
            x = rad * np.sin(theta)
            z = rad * np.cos(theta)
            return x, z
        raise Exception("This is a bug in the code")

    # def check_nugrid_v_n(self, v_n):
    #     if not v_n in self.list_nugrid_v_ns:
    #         raise NameError("v_n:{} is not in the list of nugrid v_ns:{}"
    #                         .format(v_n, self.list_nugrid_v_ns))
    #
    # def is_grid_params_extracted(self, it, v_n):
    #     pass
    #
    # def get_sph_grid_params(self, it, v_n):
    #     self.check_nugrid_v_n(v_n)
    #     self.check_it(it)
    #     self.is_grid_params_extracted(it, v_n)


class MODIFY_LOADED_NU_DATA(LOAD_NU_PROFILE):

    def __init__(self, sim, symmetry=None):
        LOAD_NU_PROFILE.__init__(self, sim, symmetry)

    def get_nuprof_arr_sph(self, it, v_n):

        nrad = self.get_nrad(it)
        nphi = self.get_nphi(it)
        ntheta = self.get_ntheta(it)

        arr = self.get_nuprof_arr(it, v_n)
        reshaped_arr = arr.reshape((nrad, nphi, ntheta))

        return reshaped_arr

    def get_nuprof_arr_slice(self, it, plane, v_n):

        if not plane in ["xy", "xz", "yz"]:
            raise NameError("plane:{} is not recognized"
                            .format(plane))

        nrad = self.get_nrad(it)
        nphi = self.get_nphi(it)
        ntheta = self.get_ntheta(it)

        fnew = self.get_nuprof_arr_sph(it, v_n)

        if plane == "xy":
            out = np.empty((nrad, nphi), dtype=fnew.dtype)
            out[:] = np.NAN
            if 0 != ntheta % 2:
                out[:,:] = fnew[:,:,ntheta/2]
            else:
                itheta = ntheta/2
                out[:,:] = 0.5*(fnew[:,:,itheta-1] + fnew[:,:,itheta])
        elif plane == "xz":
            out = np.empty((nrad, 2*ntheta-1), dtype=fnew.dtype)
            out[:] = np.NAN
            out[:,:ntheta] = fnew[:,0,:]
            iphi = int(nphi/2)
            out[:,ntheta:] = 0.5*(fnew[:,iphi-1,-2::-1] + fnew[:,iphi,-2::-1])
        elif plane == "yz":
            out = np.empty((nrad, 2*ntheta-1), dtype=fnew.dtype)
            out[:] = np.NAN
            iphi1 = int(nphi/4)
            iphi2 = int(3*iphi1)
            out[:,:ntheta] = 0.5*(fnew[:,iphi1+1,:] + fnew[:,iphi1,:])
            out[:,ntheta:] = 0.5*(fnew[:,iphi2+1,-2::-1] + \
                    fnew[:,iphi2,-2::-1])
        else:
            raise Exception("This is a bug in the code")
        return np.ma.masked_invalid(out)



if __name__ == '__main__':
    nuprof = MODIFY_LOADED_NU_DATA("LS220_M13641364_M0_LK_SR_restart")
    # arr = nuprof.get_nuprof_arr(516096, "eave_nua")
    # print(arr.shape)
    # rarr = nuprof.get_nuprof_arr_sph(516096, "eave_nua")
    # print(rarr.shape)


    v_n = 'E_nua'

    data = nuprof.get_nuprof_arr_slice(516096, "xz", v_n)
    x, z = nuprof.get_x_y_z_grid(516096, plane="xz", rmax=512)

    from matplotlib import colors

    fig = plt.figure(figsize=(6.0, 3.0))
    ax = fig.add_subplot(111)
    ax.set_aspect(1.)

    min = data[(x>20)&(z>20)&(data>0.)].min()
    max = data[(x>20)&(z>20)].max()

    print("min:{}".format(min))
    print("max:{}".format(max))

    norm = colors.Normalize(vmin=min, vmax=max)
    im = ax.pcolormesh(x, z, data, norm=norm, cmap="inferno_r")

    # divider = make_axes_locatable(ax)
    # cax = divider.new_horizontal(size="5%", pad=0.7, pack_start=True)
    # fig.add_axes(cax)
    # fig.colorbar(im, cax=cax, orientation="vertical")
    fig.colorbar(im, ax=ax)
    # cbar = ax.cax.colorbar(im)
    # cbar = grid.cbar_axes[0].colorbar(im)

    plt.xlim(-70, 70)
    plt.ylim(0, 70)
    plt.title("${}$".format(v_n).replace("_", "\_"))
    plt.savefig('{}'.format(Paths.plots + "nu_test/nu_plot_{}_test.png".format(v_n)), bbox_inches='tight', dpi=128)
    # print("saved pi_test2.png")
    plt.close()

    # exit(1)
    for v_n in nuprof.list_nuprof_v_ns:
        data = nuprof.get_nuprof_arr_slice(516096, "xz", v_n)
        x, z = nuprof.get_x_y_z_grid(516096, plane="xz", rmax=512)

        from matplotlib import colors

        fig = plt.figure(figsize=(6.0, 3.0))
        ax = fig.add_subplot(111)
        ax.set_aspect(1.)

        min = data[(x > 20) & (z > 20) & (data > 0.)].min()
        max = data[(x > 20) & (z > 20)].max()

        print("min:{}".format(min))
        print("max:{}".format(max))

        if v_n.__contains__("E_") or v_n == "flux_fac" or v_n == "abs_nue" or v_n.__contains__("eave_"):
            norm = colors.Normalize(vmin=min, vmax=max)
        else:
            norm = colors.LogNorm(vmin=min, vmax=max)

        # norm = colors.LogNorm(vmin=data[(x > 20) & (z > 20) & (data > 0.)].min(), vmax=data[(x > 20) & (z > 20)].max())
        im = ax.pcolormesh(x, z, data, norm=norm, cmap="inferno_r")

        fig.colorbar(im, ax=ax)

        plt.xlim(-70, 70)
        plt.ylim(0, 70)
        plt.title("M0: {}".format(v_n).replace("_", "\_"))
        plt.savefig('{}'.format(Paths.plots + "nu_test/nu_plot_{}.png".format(v_n)), bbox_inches='tight', dpi=128)
        # print("saved pi_test2.png")
        plt.close()


    # print(plane)
    pass
