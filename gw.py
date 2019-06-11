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

class SIM_GW:

    fold = 'waveforms/'
    l2_m2_fname = 'waveform_l2_m2.dat' # created by convert.py
    tmerg_fname = 'tmerger.dat'
    tcoll_fname = 'tcoll.dat'

    def __init__(self, sim):

        self.sim = sim
        self.l2_m2_file = MakePath.waveforms(sim) + Files.l2_m2
        self.tmerg_file = MakePath.waveforms(sim) + Files.tmerg
        self.tmerg_file = MakePath.waveforms(sim) + SIM_GW.tcoll_fname

    def load_l2_m2(self):

        time, h_plus, h_cross = np.loadtxt(self.l2_m2_file, usecols=(0, 1, 2), unpack=True)
        # time = ut.conv_time(ut.cactus, ut.cgs, time)  # in s

        return time, h_plus, h_cross # time is in [s]

    def load_tmerg(self):

        try:
            tmerg = np.loadtxt(self.tmerg_file, unpack=True)
        except IOError:
            Printcolor.yellow("Warning: file {} not found".format(self.tmerg_file))
            return np.nan

        tmerg = np.float(tmerg) # in cactus units
        tmerg = ut.conv_time(ut.cactus, ut.cgs, tmerg)  # in s

        return tmerg # in s

    def load_tcoll(self):

        try:
            tcoll = np.loadtxt(self.tmerg_fname, unpack=True)
        except IOError:
            Printcolor.yellow("Warning: file {} not found".format(self.tmerg_fname))
            return np.nan

        tcoll = np.float(tcoll) # in cactus units
        tcoll = ut.conv_time(ut.cactus, ut.cgs, tcoll)  # in s

        return tcoll

if __name__ == '__main__':
    gw = SIM_GW('DD2_M13641364_M0_SR'); print("tmerg: {} [s]".format(gw.load_tmerg()))