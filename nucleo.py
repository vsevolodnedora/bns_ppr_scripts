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
from math import pi, log10, sqrt
import scipy.optimize as opt
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

from matplotlib.ticker import AutoMinorLocator, FixedLocator, NullFormatter, \
    MultipleLocator
from matplotlib.colors import LogNorm, Normalize
from matplotlib.colors import Normalize, LogNorm

from general import *
from lists import *
from filework import *
from units import time_constant, volume_constant, energy_constant


class SIM_NUCLEO:

    def __init__(self, sim, extension='_0'):

        self.sim = sim

        self.outflowdir = MakePath.outflow(sim, extension)
        self.path = self.outflowdir + 'yields.h5'

    # @staticmethod
    # def normalize_yields(As, Ys):
    #     Anrm = np.arange(As.max() + 1)
    #     Ynrm = np.zeros(int(As.max()) + 1)
    #     for i in range(Ynrm.shape[0]):  # changed xrange to range
    #         Ynrm[i] = Ys[As == i].sum()
    #     return Anrm, Ynrm

    def load_res(self):

        h5tbl = h5py.File(self.path, "r")
        ye_final = np.array(h5tbl["Y_final"])
        a_arr = np.array(h5tbl["A"])
        z_arr = np.array(h5tbl["Z"])
        h5tbl.close()

        return a_arr, ye_final, z_arr

    # def from_res(self, normalize=True):
    #
    #     ye_final, a_arr, z_arr = self.load_res()
    #
    #     if normalize == True:
    #         a_arr, ye_final = self.normalize_yields(a_arr, ye_final)
    #         norm = ye_final.sum()
    #         ye_final /= norm
    #         return ye_final, a_arr, z_arr

    # def from_res_norm_to_sol(self, element_a, a_sol, ye_sol):
    #
    #     if len(a_sol) == 0: raise ValueError('No solar arrays prvided')
    #
    #     ye_final, a_arr, z_arr = self.load_res()
    #     a_arr, ye_final = self.normalize_yields(a_arr, ye_final)
    #
    #     if element_a not in a_sol: raise ValueError('Element: a:{} not in solar A'.format(element_a, a_sol))
    #     if element_a not in a_arr: raise ValueError('Element: a:{} not in a_arr'.format(element_a, a_arr))
    #
    #     delta = np.float(ye_sol[np.where(a_sol == element_a)] / ye_final[np.where(a_arr == element_a)])
    #
    #
    #     print(a_sol[np.where(a_sol == element_a)])
    #     print(a_arr[np.where(a_sol == element_a)])
    #     print('delta:{}'.format(delta))
    #     ye_final *= delta
    #     # # print(a_sol)
    #     # # print(ye_final)
    #     # exit(1)
    #
    #     return ye_final, a_arr

    @staticmethod
    def sum_for_all_charge_z(As, Ys):
        '''Sums all Ys for a given A (for all Z)'''
        Anrm = np.arange(As.max() + 1)
        Ynrm = np.zeros(int(As.max()) + 1)
        for i in range(Ynrm.shape[0]):  # changed xrange to range
            Ynrm[i] = Ys[As == i].sum()
        return Anrm, Ynrm

    @staticmethod
    def normalize_yields_to_sum(a_arr, ye_final):
        ''' Normalizes to a sum of all A '''
        norm = ye_final.sum()
        ye_final /= norm
        return a_arr, ye_final

    def normalize_yields_to_a(self, element_a, a_arr, ye_final):

        a_sol, ye_sol = self.solar()

        if element_a not in a_sol: raise ValueError('Element: a:{} not in solar A'.format(element_a, a_sol))
        if element_a not in a_arr: raise ValueError('Element: a:{} not in a_arr'.format(element_a, a_arr))

        delta = np.float(ye_sol[np.where(a_sol == element_a)] / ye_final[np.where(a_arr == element_a)])

        # print(a_sol[np.where(a_sol == element_a)])
        # print(a_arr[np.where(a_sol == element_a)])
        # print('delta:{}'.format(delta))
        ye_final *= delta
        # print(a_sol)
        # print(ye_final)
        # exit(1)
        # print(a_arr)
        # print(ye_final)

        return a_arr, ye_final

    def load_from_outflow(self, normalization=None):

        a_arr, ye_final, z_arr = self.load_res() # ~8000 entries, 1D array, for Z and A
        a_arr, ye_final = self.sum_for_all_charge_z(a_arr, ye_final) # ~200 entries, 1D 1 arrays

        # normalization options
        if normalization == 'sum':
            return self.normalize_yields_to_sum(a_arr, ye_final)
        elif normalization == None:
            return a_arr, ye_final
        else:
            a_to_norm = int(normalization)
            if a_to_norm < a_arr.min():
                raise ValueError('Give normalization A:{} is < a_arr.min():{}'.format(a_to_norm, a_arr.min()))
            if a_to_norm > a_arr.max():
                raise ValueError('Give normalization A:{} is > a_arr.max():{}'.format(a_to_norm, a_arr.max()))
            return self.normalize_yields_to_a(a_to_norm, a_arr, ye_final)

    def solar(self):

        Asun, Ysun = np.loadtxt(Paths.skynet + Files.solar_r, unpack=True)
        Asun, Ysun = self.sum_for_all_charge_z(Asun, Ysun)
        Ysun /= np.sum(Ysun)

        return Asun, Ysun

class PLOT_NUCLEO:

    def __init__(self):

        # self.task_dic = {'sims':     ['DD2_M13641364_M0_SR']#, 'DD2_M13641364_M0_SR'],  #   corr_vn1
        #                  'criteria': ['geo']#, 'bern wind'],                  # /_b_w/corr_vn1_vn2
        #                  'det':      [0]#, 0],                                # /outflowed_0_b_w/corr_vn1_vn2
        #                  'norm':     ['195']#, 'sum'],
        #                  'labels':   ['']#, ''],
        #                  'colors':   ['blue']#, 'red'],
        #                  'yscale': 'log',
        #                  }

        self.task_dic = {'sims':     ['DD2_M13641364_M0_SR', 'SLy4_M13641364_M0_SR', 'LS220_M13641364_M0_SR', 'SFHo_M13641364_M0_SR'],#, 'DD2_M13641364_M0_SR'],  #   corr_vn1
                         'extension': ['_0', '_0', '_0', '_0'],#, 'bern wind'],                  # /_b_w/corr_vn1_vn2                              # /outflowed_0_b_w/corr_vn1_vn2
                         'norm':     ['195', '195', '195', '195'],#, 'sum'],
                         'labels':   ['DD2', 'SLy4', "LS220", 'SFHo'],#, ''],
                         'colors':   ['blue', 'green', 'orange', 'red'],#, 'red'],
                         'yscale': 'log',
                         }
        self.fig_name = 'yields'
        self.plot_solar = True

    def plot_from_dic(self):

        n_rows = 1
        n_cols = 1

        fig = plt.figure(figsize=(6.5, 3.6))  # figsize=(4.5, 2.5 * 3.6)  # (<->; v)

        axs = []
        for n in range(1, n_rows + 1):
            if n == 1:
                axs.append(fig.add_subplot(n_rows, n_cols, n))
            else:
                axs.append(fig.add_subplot(n_rows, n_cols, n, sharex=axs[n - 2]))  # sharex=axs[n - 2]))

        for i_sim, sim in enumerate(self.task_dic['sims']):

            cl_nuc = SIM_NUCLEO(sim, self.task_dic['extension'][i_sim])

            if self.plot_solar:
                a_sol, y_sol = cl_nuc.solar()
                axs[0].plot(a_sol, y_sol, '.', color='black')

            a_arr, ye_final = cl_nuc.load_from_outflow(self.task_dic['norm'][i_sim])
            axs[0].step(a_arr, ye_final, label=self.task_dic['labels'][i_sim], color=self.task_dic['colors'][i_sim])

            axs[0].legend(loc='best', numpoints=1)
            if self.task_dic['yscale'] == 'log':
                axs[0].set_yscale("log")
            axs[0].tick_params(labelsize=12)
            axs[0].set_ylabel("Relative final abundances", fontsize=12)


        plt.ylim(ymin=1e-5, ymax=2e-1)
        plt.xlim(xmin=50, xmax=210)
        # plt.ylabel("Relative final abundances")
        plt.xlabel("A", fontsize=12)
        # plt.yscale("log")
        plt.tick_params(axis='both', which='both', labelleft=True, labelright=False, tick1On=True, tick2On=True,
                        labelsize=12, direction='in')  # labeltop
        # plt.xticks(fontsize=12)
        # plt.yticks(fontsize=12)

        plt.minorticks_on()
        plt.savefig('{}{}.png'.format(Paths.plots, self.fig_name), bbox_inches='tight') # , dpi=128
        plt.close()

if __name__ == '__main__':

    ''' PLOT YIELDS '''
    pl_cl = PLOT_NUCLEO()
    pl_cl.plot_from_dic()
    exit(1)