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

    def __init__(self, sim, extension='_0', det=0):

        self.sim = sim
        self.det = det

        self.outflowdir = MakePath.outflow(sim, extension)
        self.path = MakePath.outflow(sim, extension)

        self.path = self.path_to_outflow_dir + 'yields.h5'

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

        Asun, Ysun = np.loadtxt(Paths., unpack=True)
        Asun, Ysun = self.sum_for_all_charge_z(Asun, Ysun)
        Ysun /= np.sum(Ysun)

        return Asun, Ysun