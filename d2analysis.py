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


class SIM_2D:

    def __init__(self, sim, outdir=''):

        # new grid parameters: polar
        self.n_phi = 150
        self.n_r = 140


        simdir = Paths.gw170817 + sim + '/'

        if outdir == '':
            self.outdir = Paths.ppr_sims + sim + '/'
        else:
            Printcolor.yellow("Setitng output directory to ./")
            self.outdir = './'

        files = locate("rho.xy.h5", root=simdir, followlinks=True)
        start_t = time.time()

        print("\t Loading [{}] files...".format(len(files))),
        self.dset = h5.dataset(files)
        print(" done! (%.2f sec)" % (time.time() - start_t))

    @staticmethod
    def make_stretched_grid(x0, x1, x2, nlin, nlog):
        assert x1 > 0
        assert x2 > 0
        x_lin_f = np.linspace(x0, x1, nlin)
        x_log_f = 10.0 ** np.linspace(log10(x1), log10(x2), nlog)
        return np.concatenate((x_lin_f, x_log_f))
    def get_2d_phi_r_grid(self):

        r_cyl_f = self.make_stretched_grid(0., 15., 512., self.n_r, self.n_phi)
        phi_cyl_f = np.linspace(0, 2 * np.pi, self.n_phi)

        r_cyl   = 0.5 * (r_cyl_f[1:] + r_cyl_f[:-1])
        phi_cyl = 0.5 * (phi_cyl_f[1:] + phi_cyl_f[:-1])

        dr_cyl      = np.diff(r_cyl_f)[:, np.newaxis]
        dphi_cyl    = np.diff(phi_cyl_f)[np.newaxis, :]

        return phi_cyl, r_cyl, dphi_cyl, dr_cyl

    def interp_onto_polar_grid(self, it, variable="HYDROBASE::rho", interp=0):

        from scidata.carpet.interp import Interpolator
        # assert x_cyl_2d.shape == arr_3d.shape

        phi_cyl, r_cyl, dphi_cyl_2d, dr_cyl_2d = self.get_2d_phi_r_grid()

        r_cyl_2d, phi_cyl_2d = np.meshgrid(r_cyl, phi_cyl, indexing='ij')
        x_cyl_2d = r_cyl_2d * np.cos(phi_cyl_2d)
        y_cyl_2d = r_cyl_2d * np.sin(phi_cyl_2d)

        xi = np.column_stack([x_cyl_2d.flatten(), y_cyl_2d.flatten()])

        grid = self.dset.get_grid(iteration=it)
        rho = self.dset.get_grid_data(grid, iteration=it, variable=variable)
        irho = Interpolator(grid, rho, interp=interp)
        res_arr_2d = irho(xi).reshape(x_cyl_2d.shape)

        return res_arr_2d, r_cyl_2d, phi_cyl_2d, dr_cyl_2d, dphi_cyl_2d

    def get_dens_decomp(self, dens_2d, phi_2d, dphi_2d, dr_2d, m=1):

        integ_over_phi = np.sum(dens_2d * np.exp(1j * m * phi_2d) * dphi_2d, axis=1)

        integ_over_phi_r = np.sum(integ_over_phi * dr_2d[:, 0])

        return integ_over_phi, integ_over_phi_r

    def save_densitymode(self):

        m0_integ_over_phi = []
        m0_integ_over_phi_r = []

        m1_integ_over_phi = []
        m1_integ_over_phi_r = []

        m2_integ_over_phi = []
        m2_integ_over_phi_r = []

        r_grid = []
        times = []
        for idx, it in enumerate(self.dset.iterations):
            print("Processing frame {}/{}...".format(idx, len(self.dset.iterations) - 1)),
            start_t = time.time()

            res_arr_2d, r_cyl_2d, phi_cyl_2d, dr_cyl_2d, dphi_cyl_2d = \
                self.interp_onto_polar_grid(it, variable="HYDROBASE::rho", interp=0)

            # --- m0

            m0_integ_over_phi_it, m0_integ_over_phi_r_it = \
                self.get_dens_decomp(res_arr_2d, phi_cyl_2d, dphi_cyl_2d, dr_cyl_2d, m=0)

            # append the absolute values (array - for every r) of the complex number
            m0_integ_over_phi.append(m0_integ_over_phi_it) # np.absolute(

            # appending the absolute value of the final value (summation is done for complex numbers)
            m0_integ_over_phi_r.append(m0_integ_over_phi_r_it) # np.absolute(

            # --- m1

            m1_integ_over_phi_it, m1_integ_over_phi_r_it = \
                self.get_dens_decomp(res_arr_2d, phi_cyl_2d, dphi_cyl_2d, dr_cyl_2d, m=1)

            # append the absolute values (array - for every r) of the complex number
            m1_integ_over_phi.append(m1_integ_over_phi_it) # np.absolute(

            # appending the absolute value of the final value (summation is done for complex numbers)
            m1_integ_over_phi_r.append(m1_integ_over_phi_r_it) # np.absolute(

            # --- m2

            m2_integ_over_phi_it, m2_integ_over_phi_r_it = \
                self.get_dens_decomp(res_arr_2d, phi_cyl_2d, dphi_cyl_2d, dr_cyl_2d, m=2)

            # append the absolute values (array - for every r) of the complex number
            m2_integ_over_phi.append(m2_integ_over_phi_it) # np.absolute(

            # appending the absolute value of the final value (summation is done for complex numbers)
            m2_integ_over_phi_r.append(m2_integ_over_phi_r_it) # np.absolute(

            times.append(self.dset.get_time(it))

            # saving the r_grid
            if it == self.dset.iterations[-1]:  r_grid = r_cyl_2d[:, 0]



            del res_arr_2d
            del r_cyl_2d
            del phi_cyl_2d
            del dr_cyl_2d
            del dphi_cyl_2d

            print("done! (%.2f sec)" % (time.time() - start_t))

        dfile = h5py.File(self.outdir + '/densitymode.h5', "w")
        # --- m0
        dfile.create_dataset("m0_integ_over_phi",
                             data=np.array(m0_integ_over_phi).reshape(len(self.dset.iterations), len(r_grid)))
        dfile.create_dataset("m0_integ_over_phi_r", data=np.array(m0_integ_over_phi_r))
        # --- m1
        dfile.create_dataset("m1_integ_over_phi",
                             data=np.array(m1_integ_over_phi).reshape(len(self.dset.iterations), len(r_grid)))
        dfile.create_dataset("m1_integ_over_phi_r", data=np.array(m1_integ_over_phi_r))
        # --- m2
        dfile.create_dataset("m2_integ_over_phi",
                             data=np.array(m2_integ_over_phi).reshape(len(self.dset.iterations), len(r_grid)))
        dfile.create_dataset("m2_integ_over_phi_r", data=np.array(m2_integ_over_phi_r))

        dfile.create_dataset("r_cyl", data=r_grid)
        dfile.create_dataset("iterations", data=self.dset.iterations)
        dfile.create_dataset("time", data=times)
        dfile.close()

def plot_dens_modes(sim, outname):

    fname = Paths.ppr_sims + sim + '/densitymode.h5'

    file = h5py.File(fname, 'r')
    print("Data File: {} [".format(fname)),
    for i in file:
        print(i),
    print(']')


    # COMPLEX NUMBERS WARNING
    # --- m0
    m0_integ_over_phi = np.array(file["m0_integ_over_phi"])      # 2D
    m0_integ_over_phi_r = np.array(file["m0_integ_over_phi_r"])  # 1D
    # --- m1
    m1_integ_over_phi = np.array(file["m1_integ_over_phi"])      # 2D
    m1_integ_over_phi_r = np.array(file["m1_integ_over_phi_r"])  # 1D
    # --- m2
    m2_integ_over_phi = np.array(file["m2_integ_over_phi"])      # 2D
    m2_integ_over_phi_r = np.array(file["m2_integ_over_phi_r"])  # 1D


    r_cyl = np.array(file["r_cyl"])  # 1D
    iterations = np.array(file["iterations"])  # 1D
    time = np.array(file["time"]) * time_constant * 1e-2 # s

    # --- m1

    m1_integ_over_phi_r_phase = np.angle(m1_integ_over_phi_r, deg=True)

    m1_integ_over_phi_r_normed = m1_integ_over_phi_r / m0_integ_over_phi_r
    m1_integ_over_phi_normed = m1_integ_over_phi / m0_integ_over_phi

    m1_integ_over_phi_r_normed_abs = np.absolute(m1_integ_over_phi_r_normed)
    m1_integ_over_phi_r_normed_phase = np.angle(m1_integ_over_phi_r_normed, deg=True)

    m1_integ_over_phi_normed_abs = np.absolute(m1_integ_over_phi_normed)
    m1_integ_over_phi_normed_phase = np.angle(m1_integ_over_phi_normed, deg=True)

    # --- m2

    m2_integ_over_phi_r_normed = m2_integ_over_phi_r / m0_integ_over_phi_r
    m2_integ_over_phi_r_normed_abs = np.absolute(m2_integ_over_phi_r_normed)

    # fig = plt.figure()
    # ax = fig.add_subplot(211)
    # ax.plot(time * 1e2, m1_integ_over_phi_r_normed_phase, '-', color='black', label=r'Phase $C_1 / C_0$')
    #
    # ax2 = fig.add_subplot(212)
    # ax2.plot(time * 1e2, m1_integ_over_phi_r_normed_abs, '-', color='black', label=r'Mag $C_1$')
    #
    #
    # plt.savefig('{}.png'.format('tst_phase'), bbox_inches='tight', dpi=128)
    # plt.legend()
    # plt.close()

    # exit(1)
    # --- plotting

    n_rows = 4
    n_cols = 1

    fig = plt.figure(figsize=(6.5, 3.5 * 3.6))  # figsize=(4.5, 2.5 * 3.6)  # (<->; v)
    axs = []
    for n in range(1, n_rows + 1):
        if n == 1:
            axs.append(fig.add_subplot(n_rows, n_cols, n))
        else:
            axs.append(fig.add_subplot(n_rows, n_cols, n, sharex=axs[n - 2]))  # sharex=axs[n - 2]))

    i = 0

    divider1 = make_axes_locatable(axs[i])
    cax = divider1.append_axes("right", size="5%", pad=0.05)
    cax.xaxis.set_ticklabels([])
    cax.get_xaxis().set_visible(False)
    cax.get_yaxis().set_visible(False)
    # fig.patch.set_visible(False)
    cax.axis('off')
    axs[i].plot(time * 1e2, m1_integ_over_phi_r_normed_abs, '-', color='black', label=r'$C_1/C_0$')
    axs[i].plot(time * 1e2, m2_integ_over_phi_r_normed_abs, '-', color='gray', label=r'$C_2/C_0$')
    axs[i].set_yscale('log')
    axs[i].set_ylabel(r'Mag $C_m = \int{\rho\cdot\exp(i m \phi) dr d\phi}$')
    axs[i].legend()

    i = 1

    divider1 = make_axes_locatable(axs[i])
    cax = divider1.append_axes("right", size="5%", pad=0.05)
    cax.xaxis.set_ticklabels([])
    cax.get_xaxis().set_visible(False)
    cax.get_yaxis().set_visible(False)
    # fig.patch.set_visible(False)
    cax.axis('off')
    axs[i].plot(time * 1e2, np.unwrap(m1_integ_over_phi_r_normed_phase), '-', color='black', label=r'$C_1/C_0$')
    # axs[i].plot(time * 1e2, m2_integ_over_phi_r / m0_integ_over_phi_r, '-', color='gray', label=r'$C_2/C_0$')
    # axs[i].set_yscale('log')
    axs[i].set_ylabel(r'Phase $C_m = \int{\rho\cdot\exp(i m \phi) dr d\phi}$')
    axs[i].legend()

    i = 2


    from matplotlib import colors

    norm = colors.LogNorm(vmin=1.e-2, vmax=1.)
    TIME, R = np.meshgrid(time, r_cyl)


    im = axs[i].pcolormesh(TIME*1e2, R, m1_integ_over_phi_normed_abs.T, cmap='Blues', norm=norm)#, vmin=0., vmax=0.0002)#, vmin=self.jmin, vmax=self.jmax)
    im.set_rasterized(True)

    divider1 = make_axes_locatable(axs[i])
    cax = divider1.append_axes("right", size="5%", pad=0.05)
    cbar = plt.colorbar(im, cax=cax, norm=norm)
    cbar.ax.set_title(r"$Mag: C_1(r) / C_0(r)$")
    cbar.ax.minorticks_off()

    axs[i].set_ylabel(r'$R$')
    axs[i].set_ylim(0, 50)
    # axs[i].set_xlabel('time [ms]')

    i = 3

    norm = colors.Normalize(vmin=-180, vmax=180)
    TIME, R = np.meshgrid(time, r_cyl)


    im = axs[i].pcolormesh(TIME*1e2, R, np.unwrap(m1_integ_over_phi_normed_phase.T), cmap='RdBu_r', norm=norm)#, vmin=0., vmax=0.0002)#, vmin=self.jmin, vmax=self.jmax)
    im.set_rasterized(True)

    divider1 = make_axes_locatable(axs[i])
    cax = divider1.append_axes("right", size="5%", pad=0.05)
    cbar = plt.colorbar(im, cax=cax, norm=norm)
    cbar.ax.set_title(r"$Phase: C_1(r) / C_0(r)$")
    cbar.ax.minorticks_off()

    axs[i].set_ylabel(r'$R$')
    axs[i].set_ylim(0, 50)
    axs[i].set_xlabel('time [ms]')


    # plt.savefig('{}.png'.format('tst_densmodes'), bbox_inches='tight', dpi=128)
    # plt.close()
    #
    #
    #
    #
    #
    #
    #
    #
    #
    # fig, ax = plt.subplots()
    # ax.plot(time*1e2, m1_integ_over_phi_r / m0_integ_over_phi_r, '-', color='black', label='m1/m0')
    # ax.plot(time*1e2, m2_integ_over_phi_r / m0_integ_over_phi_r, '-', color='gray', label='m2/m0')
    # ax.set_yscale('log')
    # plt.savefig('{}.png'.format('tst_dens_modes'), bbox_inches='tight', dpi=128)
    # plt.legend()
    # plt.close()
    #
    #
    # fig, ax = plt.subplots()
    # from matplotlib import colors
    # norm = colors.LogNorm(vmin=1.e-2, vmax=1.)
    #
    # TIME, R = np.meshgrid(time, r_cyl)
    # arr = m1_integ_over_phi / m0_integ_over_phi
    # im = ax.pcolormesh(TIME*1e2, R, arr.T, cmap='Blues', norm=norm)#, vmin=0., vmax=0.0002)#, vmin=self.jmin, vmax=self.jmax)
    # im.set_rasterized(True)
    #
    # cax1 = fig.add_axes([0.9, 0.1, 0.03, 0.8])
    # cbar = plt.colorbar(im, cax=cax1)
    # cbar.ax.set_title(r"$C_1(r) / C_0(r)$")
    #
    # ax.set_ylim(0, 50)

    plt.savefig('{}{}.png'.format(Paths.plots, outname), bbox_inches='tight', dpi=128)
    plt.close()

if __name__ == '__main__':

    ''' DENSITY MODES D2 '''
    d2 = SIM_2D('DD2_M13641364_M0_SR')
    d2.save_densitymode()  # save the density decomposition