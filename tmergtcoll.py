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
import scivis.units as ut  # for tmerg
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

strain_fname = "strain_l2_m2.dat"
outfile_tmerg = "tmerger2.dat"
outfile_tcoll = "tcoll.dat"
plot_name = "tmergtcoll.png"

for sim in ["LS220_M13641364_M0_LK_SR"]:#os.listdir(Paths.ppr_sims):
    print("Processing: {}".format(sim))

    if os.path.isdir(Paths.ppr_sims + sim + '/' + "waveforms/") and \
            os.path.isfile(Paths.ppr_sims + sim + '/' + "waveforms/" + strain_fname):# \
            #and sim == "SLy4_M13641364_M0_SR":

        w_dir = Paths.ppr_sims + sim + '/' + "waveforms/"
        strain_file = w_dir + strain_fname

        time, reh, imh = np.genfromtxt(strain_file, unpack=True, usecols=[0, 1, 2])
        h = reh + 1j * imh
        amp = np.abs(h)
        tmerg = time[np.argmax(amp)]

        print "tmerg: {:.4f}".format(tmerg)
        Printcolor.blue("\tSaving t merger:   {}".format(w_dir + outfile_tmerg))
        open(w_dir + outfile_tmerg, "w").write("{}\n".format(float(tmerg)))

        # reh = reh[time > tmerg]
        # amp = amp[time > tmerg]
        # time = time[time > tmerg]

        maxA = np.max(amp[time > tmerg])
        mins = np.where(amp[time > tmerg] < 0.01 * maxA)
        if len(mins) > 0:

            # TODO WROTE A WAY TO FIND A BEGINNING OF A LARGEST FLAT RGION < 0.05% OF MAXIMUM

            coll_part = amp[amp[time > tmerg] < 0.01 * maxA]
            t_coll_parts = time[amp[time > tmerg] < 0.01 * maxA]
            # diff_t_coll = np.diff(t_coll_parts)
            # t_begins = []
            # for i in range(len(t_coll_parts)-1):
            #     if t_coll_parts[i] != t_coll_parts[i+1] - t_coll_parts[i] != t_coll_parts[i+1]:
            #         t_begins.append(t_coll_parts[i])
            #
            # print(len(t_begins)); #exit(1)


            indx_collapse1 = np.min(np.where(amp[time > tmerg] < 0.01 * maxA))
            indx_collapse_last = np.max(np.where(amp[time > tmerg] < 0.01 * maxA))
            indx_collapse = indx_collapse1

            if indx_collapse_last <= indx_collapse1:
                raise ValueError("Index error 0<last")

            region = amp[indx_collapse1:indx_collapse_last]
            while max(region) > 0.01 * maxA and len(region) > 2:
                region = np.delete(region, 0, 0)
                indx_collapse1 += 1
                if len(region) == 3:
                    raise ValueError("Failed to fild indx_collapse1 iteratively "
                                                      "substracting the parts of > 0.01maxA ")


            tcoll = time[indx_collapse1]
            if tcoll < tmerg:
                Printcolor.yellow("\ttcoll < tmerg")
            print "tcoll: {:.4f}".format(tcoll)

            Printcolor.blue("\tSaving t collapse: {}".format(w_dir + outfile_tmerg))
            open(w_dir + outfile_tcoll, "w").write("{}\n".format(float(tcoll)))

        f = plt.figure()
        plt.plot(time, reh, c='k')
        plt.plot(time, amp, c='r', label="amplitude")
        plt.axvline(tmerg, ls='-.', c='b', label="tmerg")

        if len(mins) > 0:
            plt.axvline(tcoll, ls='-.', c='c', label="tcoll")

        # for t in t_coll_parts[::10]:
        plt.axvline(time[indx_collapse1], ls=':', c='gray')
        # plt.axvline(time[indx_collapse_last], ls=':', c='gray')

        plt.ylabel(r'strain [Re]', fontsize=12)
        plt.xlabel(r'time [s]', fontsize=12)
        plt.minorticks_on()
        plt.xticks(fontsize=12)
        plt.yticks(fontsize=12)
        plt.title('GW analysis', fontsize=20)
        plt.legend(loc='upper right', numpoints=1)
        plt.savefig(w_dir + plot_name, bbox_inches='tight', dpi=128)
        plt.close()

    else:
        Printcolor.red("Error: no {} in wafeforms/ found".format(strain_fname))
