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
import warnings
import os.path
import cPickle
import time
import copy
import click
import h5py
import csv
import os
import gc
warnings.filterwarnings("ignore",category=matplotlib.mplDeprecation)
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


from plotting_nexus import PLOT_MANY_TASKS
from d1analysis import COMPUTE_STORE_PAR, NORMALIZE_NUCLEO, ADD_METHODS_1D

# from scidata import units
# print(units.conv_length(units.cactus, units.cgs, 60))/100000.
# print(units.conv_length(units.cactus, units.cgs, 40))/100000.
# exit(1)


# --------------------------------
sim = "DD2_M15091235_M0_LK_HR"
simdir = "/data1/numrel/WhiskyTHC/Backup/2018/GW170817/" + sim +'/'
resdir = "/data01/numrel/vsevolod.nedora/postprocessed3/" + sim +'/'
scrdic = "/data01/numrel/vsevolod.nedora/scripts_server/"
reslist= ["*.done", "output-0003/parfile.par", "berntime.*", "WARNING*",
          "outflow*", "collated*", "waveforms"]


# -------------------------------

def get_it_time_cum_mass(simdir, asci_file="outflow_det_0.asc"):
    fpath = simdir + "output-*" + '/data/' + asci_file
    files = glob(fpath)
    print("\tcum_muss files found {}.".format(len(files)))
    if len(files) == 0:
        raise ValueError("No cum_muss files found for {} found in outputs"
                         .format(fpath))

    times, fluxes1, fluxes2 = np.zeros(0,), np.zeros(0,), np.zeros(0,)
    for file_ in files:
        time_, flux1, flux2 = np.loadtxt(file_, usecols=(1, 2, 5), unpack=True)
        times = np.hstack((times, time_))
        fluxes1 = np.hstack((fluxes1, flux1))
        fluxes2 = np.hstack((fluxes2, flux2))

    times = np.array(times)
    fluxes1 = np.array(fluxes1)
    fluxes2 = np.array(fluxes2)

    print(times.shape, fluxes1.shape, fluxes2.shape)

    times, fluxes1, fluxes2 = x_y_z_sort(times,
                                         fluxes1,
                                         fluxes2, sort_by_012=0)

    mass1 = mass2 = 0.0
    masses1, masses2 = [], []

    for i in range(1, times.shape[0]):
        mass1 += fluxes1[i-1] * (times[i] - times[i-1])
        mass2 += fluxes2[i-1] * (times[i] - times[i-1])
        masses1.append(mass1)
        masses2.append(mass2)

    return np.array(times[1:]), np.array(masses1), np.array(masses2)

def get_bern_time(simdir, fraction = 0.98, berntimefile = "berntime.txt", plot = "berntime.png"):

    times, masses1, masses2 = get_it_time_cum_mass(simdir)

    idx_ = find_nearest_index(masses2, fraction * masses2.max())

    if idx_ == len(times) - 1:
        raise ValueError("{}percent of ejecta is found at last iteration. Cannot use for Bernoulli"
                         .format(fraction * 100))

    if times[idx_] > 0.95 * times.max():
        Printcolor.yellow("\tWarning! Berntime is {}percent of the total time".format(times[idx_] * 100 / times.max()))

    Printcolor.blue("\t Bern. time is {:.3f}s (out of {:.3f}s total)"
                    .format(times[idx_] * 0.004925794970773136 * 1e-3,
                            times.max() * 0.004925794970773136 * 1e-3))

    berntimefile = simdir + berntimefile
    if berntimefile != None and berntimefile != '':
        if os.path.isfile(berntimefile):
            Printcolor.yellow("\tWarning. File:{} already exist. Rewriting".format(berntimefile))
            os.remove(berntimefile)
        else:
            Printcolor.blue("\tSaving {}".format(berntimefile))
        open(berntimefile, "w").write("{}\n".format(float(times[idx_])))

    if plot != None and plot != '':

        import matplotlib.pyplot as plt

        Printcolor.blue("\tPlotting {}".format(simdir + '{}.png'.format(plot)))
        fig = plt.figure()
        ax = fig.add_subplot(111)

        ax.axvline(x=times[idx_] * 0.004925794970773136 * 1e-3, color='black', lw='0.5', label='Bern time')
        ax.plot(times * 0.004925794970773136 * 1e-3, masses2 * 1e2, '-', color='black', label='Bernoulli wind')
        ax.plot(times * 0.004925794970773136 * 1e-3, masses1 * 1e2, '-', color='gray', label='Geodesic')
        ax.set_xlabel(r'time $[M_{\odot}$]', fontsize=12)
        ax.set_ylabel(r'$M_{ej}$ $[10^{-2}M_{\odot}]$', fontsize=12)
        ax.minorticks_on()
        ax.set_title('Outflow', fontsize=20)

        ax.legend(loc='upper right', numpoints=1)
        ax.tick_params(
            axis='both', which='both', labelleft=True,
            labelright=False, tick1On=True, tick2On=True,
            labelsize=12,
            direction='in',
            bottom=True, top=True, left=True, right=True
        )
        plt.savefig(simdir + plot, bbox_inches='tight', dpi=128)
        plt.close()


    return float(times[idx_])



''' --- --- --- D1 --- --- --- '''

def ejecta_profiles(d1class, plot_bern = True, fig_dir ='res_1d/', figname = 'ejecta_profile.png'):

    Printcolor.blue("\tInitializing ejecta profile plotting.")

    if not os.path.isdir(resdir):
        raise IOError("No simulation results frolder found: {}".format(resdir))
    if not os.path.isdir(resdir + fig_dir):
        Printcolor.yellow("\tCreating {} ".format(fig_dir))
        os.mkdir(resdir+fig_dir)

    o_plot = PLOT_MANY_TASKS()
    o_plot.gen_set["figdir"] = resdir + fig_dir
    o_plot.gen_set["type"] = "cartesian"
    o_plot.gen_set["figsize"] = (4.0, 3.5)  # <->, |]
    o_plot.gen_set["figname"] = figname
    o_plot.gen_set["sharex"] = False
    o_plot.gen_set["sharey"] = False
    o_plot.gen_set["subplots_adjust_h"] = 0.3
    o_plot.gen_set["subplots_adjust_w"] = 0.0
    o_plot.set_plot_dics = []

    # label1 = "Dyn"
    # label2 = "Wind"

    """ -------------- DD2 LK ---------------------"""

    Printcolor.blue("\tSetting ejecta profile plot parameters")
    # geo
    dic_ej_prof = {
        'task': 'ejprof', 'ptype': 'cartesian',
        'position': (1, 1),
        'data': d1class, 'criterion': '_0',
        'v_n_x': 't_tot_flux', 'v_n_y': 'mass_tot_flux',
        'color': 'red', 'ls': '--', 'lw': 1., 'ds': 'default','alpha':1.,
        'ymin': 0, 'ymax': None,
        'xlabel': r"time [ms]", 'ylabel': r'$M_{ej}$ $[10^{-2}M_{\odot}]$',
        'xunits': 'ms', '-t': d1class.get_par("tmerger_gw"), 'yunits': '1e-2Msun',
        'label': "Dyn.", 'yscale': None, 'title': {},
        'fancyticks': True, 'minorticks': True,
        'fontsize': 16,
        'labelsize': 16,
        'legend': {'loc': 'best', 'ncol': 1, 'fontsize': 14}
    }
    o_plot.set_plot_dics.append(dic_ej_prof)
    # bern
    if plot_bern:
        dic_ej_prof = {
            'task': 'ejprof', 'ptype': 'cartesian',
            'position': (1, 1),
            'data': d1class, 'criterion': '_0_b_w',
            'v_n_x': 't_tot_flux', 'v_n_y': 'mass_tot_flux',
            'color': 'red', 'ls': '-', 'lw': 1., 'ds': 'default','alpha':1.,
            # 'ymin': 1e-4, 'ymax': 1e0,
            'xlabel': r"time [ms]", 'ylabel': r'$M_{ej}$ $[10^{-2}M_{\odot}]$',
            'xunits': 'ms', '-t': d1class.get_par("tmerger_gw"), 'yunits': '1e-2Msun',
            'label': "Wind ", 'yscale': None, 'title': {},
            'fancyticks': True, 'minorticks': True,
            'fontsize': 16,
            'labelsize': 16,
            'legend': {'loc': 'best', 'ncol': 1, 'fontsize': 14}
        }
        o_plot.set_plot_dics.append(dic_ej_prof)
    Printcolor.blue("\tPlotting ejecta profile")
    o_plot.main()
    Printcolor.blue("\tPlotting ejecta profile is done.")
    print('\n')

def ejecta_properties(d1class, plot_bern = True, fig_dir ='res_1d/', figname = 'histograms.png'):

    Printcolor.blue("\tInitializing ejecta profile plotting.")

    if not os.path.isdir(resdir):
        raise IOError("No simulation results frolder found: {}".format(resdir))
    if not os.path.isdir(resdir + fig_dir):
        Printcolor.yellow("\tCreating {} ".format(fig_dir))
        os.mkdir(resdir+fig_dir)

    o_plot = PLOT_MANY_TASKS()
    o_plot.gen_set["figdir"] = resdir + fig_dir
    o_plot.gen_set["type"] = "cartesian"
    o_plot.gen_set["figsize"] = (9.0, 3.)  # <->, |]
    o_plot.gen_set["figname"] = figname
    o_plot.gen_set["sharex"] = False
    o_plot.gen_set["sharey"] = False
    o_plot.gen_set["subplots_adjust_h"] = 0.3
    o_plot.gen_set["subplots_adjust_w"] = 0.0
    o_plot.set_plot_dics = []

    ylabel = r'$M_{ej}$ [normed to 1]'
    label1 =  "Dyn."
    label2 =  "Wind"
    fontsize = 16
    labelsize = 16
    ''' --- --- sim 1 --- --- --- '''

    color = "red"
    Printcolor.blue("\tSetting histograms profile plot parameters")
    dic_hist_theta = {
        'task': 'hist1d', 'ptype': 'cartesian',
        'position': (1, 1),
        'data': d1class, 'norm': True, 'criterion': '_0',
        'v_n_x': 'hist_theta', 'v_n_y': 'hist_theta_m',
        'color': color, 'ls': ':', 'lw': 0.6, 'ds': 'steps', 'alpha':0.6,
        'ymin': 1e-4, 'ymax': 1e0,
        'xlabel': r"Angle from orbital plane", 'ylabel': ylabel,
        'label': label1, 'yscale': 'log',
        'fancyticks': True, 'minorticks': True,
        'fontsize': fontsize,
        'labelsize': labelsize,
        # 'legend': {'loc': 'best', 'ncol': 2, 'fontsize': 18}
    }
    o_plot.set_plot_dics.append(dic_hist_theta)
    dic_hist_theta_b = {
        'task': 'hist1d', 'ptype': 'cartesian',
        'position': (1, 1),
        'data': d1class, 'norm': True, 'criterion': '_0_b_w',
        'v_n_x': 'hist_theta', 'v_n_y': 'hist_theta_m',
        'color': color, 'ls': '-', 'lw': 1., 'ds': 'steps', 'alpha':1.,
        'ymin': 1e-4, 'ymax': 1e0,
        'xlabel': r"Angle from orbital plane", 'ylabel': ylabel,
        'label': label2,
        'yscale': 'log',
        'fancyticks': True, 'minorticks': True,
        'fontsize': fontsize,
        'labelsize': labelsize,
        # 'legend': {'loc': 'best', 'ncol': 2, 'fontsize': 18}
    }
    if plot_bern: o_plot.set_plot_dics.append(dic_hist_theta_b)
    dic_hist_vel_inf = {
        'task': 'hist1d', 'ptype': 'cartesian',
        'position': (1, 2),
        'data': d1class, 'norm': True, 'criterion': '_0',
        'v_n_x': 'hist_vel_inf', 'v_n_y': 'hist_vel_inf_m',
        'color': color, 'ls': ':', 'lw': 0.6, 'ds': 'steps', 'alpha':0.8,
        'xmin': 0., 'xmax': 1., 'ymin': 1e-4, 'ymax': 1e0,
        'xlabel': r"$\upsilon_{\infty}$ [c]", 'ylabel': ylabel,
        'label': label1, 'yscale': 'log',
        'fancyticks': True, 'minorticks': True,
        'sharey': True,
        'fontsize': fontsize,
        'labelsize': labelsize,
        'legend': {'loc': 'best', 'ncol': 1, 'fontsize': 14}
    }
    o_plot.set_plot_dics.append(dic_hist_vel_inf)
    dic_hist_vel_inf_b = {
        'task': 'hist1d', 'ptype': 'cartesian',
        'position': (1, 2),
        'data': d1class, 'norm': True, 'criterion': '_0_b_w',
        'v_n_x': 'hist_vel_inf', 'v_n_y': 'hist_vel_inf_m',
        'color': color, 'ls': '-', 'lw': 1., 'ds': 'steps', 'alpha':1.,
        'xmin': 0., 'xmax': 1., 'ymin': 1e-4, 'ymax': 1e0,
        'xlabel': r"$\upsilon_{\infty}$ [c]", 'ylabel': ylabel,
        'label': label2, 'yscale': 'log',
        'fancyticks': True, 'minorticks': True,
        'sharey': True,
        'fontsize': fontsize,
        'labelsize': labelsize,
        'legend': {'loc': 'best', 'ncol': 1, 'fontsize': 14}
    }
    if plot_bern: o_plot.set_plot_dics.append(dic_hist_vel_inf_b)
    dic_hist_ye = {
        'task': 'hist1d', 'ptype': 'cartesian',
        'position': (1, 3),
        'data': d1class, 'norm': True, 'criterion': '_0',
        'v_n_x': 'hist_ye', 'v_n_y': 'hist_ye_m',
        'color': color, 'ls': ':', 'lw': 0.6, 'ds': 'steps', 'alpha':0.8,
        'xmin': 0., 'xmax': 0.5, 'ymin': 1e-4, 'ymax': 1e0,
        'xlabel': r"$Y_e$", 'ylabel': ylabel,
        'label': None, 'yscale': 'log',
        'fancyticks': True, 'minorticks': True,
        'sharey': True,
        'fontsize': fontsize,
        'labelsize': labelsize,
        # 'legend': {'loc': 'best', 'ncol': 2, 'fontsize': 18}
    }
    o_plot.set_plot_dics.append(dic_hist_ye)
    dic_hist_ye_b = {
        'task': 'hist1d', 'ptype': 'cartesian',
        'position': (1, 3),
        'data': d1class, 'norm': True, 'criterion': '_0_b_w',
        'v_n_x': 'hist_ye', 'v_n_y': 'hist_ye_m',
        'color': color, 'ls': '-', 'lw': 1., 'ds': 'steps', 'alpha':1.,
        'xmin': 0., 'xmax': 0.5, 'ymin': 1e-4, 'ymax': 1e0,
        'xlabel': r"$Y_e$", 'ylabel': ylabel,
        'label': None, 'yscale': 'log',
        'fancyticks': True, 'minorticks': True,
        'sharey': True,
        'fontsize': fontsize,
        'labelsize': labelsize,
        # 'legend': {'loc': 'best', 'ncol': 2, 'fontsize': 18}
    }
    if plot_bern: o_plot.set_plot_dics.append(dic_hist_ye_b)


    Printcolor.blue("\tPlotting ejecta histograms")
    o_plot.main()
    Printcolor.blue("\tPlotting ejecta histograms is done.")
    # exit(1)

def ejecta_correlations_vinf(d1class, plot_bern = True, fig_dir ='res_1d/', figname = 'correlations_vinf.png'):

    Printcolor.blue("\tInitializing ejecta profile plotting.")

    if not os.path.isdir(resdir):
        raise IOError("No simulation results frolder found: {}".format(resdir))
    if not os.path.isdir(resdir + fig_dir):
        Printcolor.yellow("\tCreating {} ".format(fig_dir))
        os.mkdir(resdir + fig_dir)

    o_plot = PLOT_MANY_TASKS()
    o_plot.gen_set["figdir"] = resdir + fig_dir
    o_plot.gen_set["type"] = "cartesian"
    o_plot.gen_set["figsize"] = (9.0, 3.)  # <->, |]
    o_plot.gen_set["figname"] = figname
    o_plot.gen_set["sharex"] = False
    o_plot.gen_set["sharey"] = False
    o_plot.gen_set["subplots_adjust_h"] = 0.3
    o_plot.gen_set["subplots_adjust_w"] = 0.0
    o_plot.set_plot_dics = []

    v_n = r'$M_{ej}$ [normed to 1]'

    fontsize = 16
    labelsize = 16

    """ --- --- --- --- DD2 --- --- --- --- """

    Printcolor.blue("\tSetting correlations profile plot parameters")
    corr_vinf_theta = {  # relies on the "get_res_corr(self, it, v_n): " method of data object
        'task': 'outflow corr', 'dtype': 'corr', 'ptype': 'cartesian',
        'data': d1class, 'criterion': '_0',
        'position': (1, 1),
        'cbar': {},
        'v_n_x': 'vel_inf', 'v_n_y': 'theta',  'v_n': v_n,
        'xlabel': r'$\upsilon_{\infty}$ [c]', 'ylabel': r"Angle from orbital plane",
        'xmin': 0.05, 'xmax': 0.65, 'ymin': 0, 'ymax': 90, 'vmin': 1e-4, 'vmax': 1e-1,
        'xscale': None, 'yscale': None, 'normalize': True,
        'mask_below': 1e-4, 'mask_above': None, 'cmap': 'Reds', 'norm': 'log', 'todo': None,
        'fancyticks': True, 'minorticks': True,
        'fontsize': fontsize,
        'labelsize': labelsize,
        'title': {},
        'legend': {},
        'text':{'coords':(0.8, 0.8), 'text':"Dyn.", 'color':'black', 'fs':16}
    }
    o_plot.set_plot_dics.append(corr_vinf_theta)

    if plot_bern:
        corr_vinf_theta = {  # relies on the "get_res_corr(self, it, v_n): " method of data object
            'task': 'outflow corr', 'dtype': 'corr', 'ptype': 'cartesian',
            'data': d1class, 'criterion': '_0_b_w',
            'position': (1, 2),
            'cbar': {'location': 'right .03 .0', 'label': v_n,
                     'labelsize': 16, 'fontsize': 16}, # 'fmt': '%.1f',
            'v_n_x': 'vel_inf_bern', 'v_n_y': 'theta',  'v_n': v_n,
            'xlabel': r'$\upsilon_{\infty}$ [c]', 'ylabel': r"Angle from orbital plane",
            'xmin': 0.05, 'xmax': 0.65, 'ymin': 0, 'ymax': 90, 'vmin': 1e-4, 'vmax': 1e-1,
            'xscale': None, 'yscale': None, 'normalize': True,
            'mask_below': 1e-4, 'mask_above': None, 'cmap': 'Reds', 'norm': 'log', 'todo': None,
            'fancyticks': True, 'minorticks': True,
            'fontsize': fontsize,
            'labelsize': labelsize,
            'title': {},
            'legend': {},
            'text':{'coords':(0.8, 0.8), 'text':"Wind", 'color':'black', 'fs':16},
            'sharey': True
        }
        o_plot.set_plot_dics.append(corr_vinf_theta)



    coll_time_vert = {
        'task': 'vertline', 'dtype': '-', 'ptype': 'cartesian',
        'value': 0.25,
        'position': (1, 2),
        'ls': '-.', 'color': 'black', 'lw': 0.8,
    }
    o_plot.set_plot_dics.append(coll_time_vert)


    Printcolor.blue("\tPlotting ejecta correlations")
    o_plot.main()
    Printcolor.blue("\tPlotting ejecta correlations is done.")
    print('\n')

def ejecta_correlations_ye(d1class, plot_bern = True, fig_dir ='res_1d/', figname = 'correlations_ye.png'):

    Printcolor.blue("\tInitializing ejecta profile plotting.")

    if not os.path.isdir(resdir):
        raise IOError("No simulation results frolder found: {}".format(resdir))
    if not os.path.isdir(resdir + fig_dir):
        Printcolor.yellow("\tCreating {} ".format(fig_dir))
        os.mkdir(resdir + fig_dir)

    o_plot = PLOT_MANY_TASKS()
    o_plot.gen_set["figdir"] = resdir + fig_dir
    o_plot.gen_set["type"] = "cartesian"
    o_plot.gen_set["figsize"] = (9.0, 3.)  # <->, |]
    o_plot.gen_set["figname"] = figname
    o_plot.gen_set["sharex"] = False
    o_plot.gen_set["sharey"] = False
    o_plot.gen_set["subplots_adjust_h"] = 0.3
    o_plot.gen_set["subplots_adjust_w"] = 0.0
    o_plot.set_plot_dics = []

    v_n = r'$M_{ej}$ [normed to 1]'

    fontsize = 16
    labelsize = 16

    """ --- --- --- --- DD2 --- --- --- --- """

    Printcolor.blue("\tSetting correlations profile plot parameters")
    corr_vinf_theta = {  # relies on the "get_res_corr(self, it, v_n): " method of data object
        'task': 'outflow corr', 'dtype': 'corr', 'ptype': 'cartesian',
        'data': d1class, 'criterion': '_0',
        'position': (1, 1),
        'cbar': {}, # 'fmt': '%.1f',
        'v_n_x': 'ye', 'v_n_y': 'theta', 'v_n': 'mass',
        'xlabel': r'$Y_e$', 'ylabel': r"Angle from orbital plane",
        'xmin': 0.00, 'xmax': 0.5, 'ymin': 0, 'ymax': 90, 'vmin': 1e-4, 'vmax': 1e-1,
        'xscale': None, 'yscale': None, 'normalize': True,
        'mask_below': 1e-4, 'mask_above': None, 'cmap': 'Reds', 'norm': 'log', 'todo': None,
        'fancyticks': True, 'minorticks': True,
        'fontsize': fontsize,
        'labelsize': labelsize,
        'legend': {},  # {'loc': 'best', 'ncol': 2, 'fontsize': 18}
        'text':{'coords':(0.8, 0.8), 'text':"Dyn.", 'color':'black', 'fs':16}
    }
    o_plot.set_plot_dics.append(corr_vinf_theta)

    if plot_bern:
        corr_vinf_theta = {  # relies on the "get_res_corr(self, it, v_n): " method of data object
            'task': 'outflow corr', 'dtype': 'corr', 'ptype': 'cartesian',
            'data': d1class, 'criterion': '_0_b_w',
            'position': (1, 2),
            'cbar': {'location': 'right .03 .0', 'label': v_n,
                     'labelsize': 16, 'fontsize': 16},  # 'fmt': '%.1f',
            'v_n_x': 'ye', 'v_n_y': 'theta', 'v_n': 'mass',
            'xlabel': r'$Y_e$', 'ylabel': r"Angle from orbital plane",
            'xmin': 0.00, 'xmax': 0.5, 'ymin': 0, 'ymax': 90, 'vmin': 1e-4, 'vmax': 1e-1,
            'xscale': None, 'yscale': None, 'normalize': True,
            'mask_below': 1e-4, 'mask_above': None, 'cmap': 'Reds', 'norm': 'log', 'todo': None,
            'fancyticks': True, 'minorticks': True,
            'fontsize': fontsize,
            'labelsize': labelsize,
            'legend':  {}, # {'loc': 'best', 'ncol': 2, 'fontsize': 18},
            'text':{'coords':(0.8, 0.8), 'text':"Wind", 'color':'black', 'fs':16},
            'sharey': True
        }
        o_plot.set_plot_dics.append(corr_vinf_theta)


    coll_time_vert = {
        'task': 'vertline', 'dtype': '-', 'ptype': 'cartesian',
        'value': 0.25,
        'position': (1, 2),
        'ls': '-.', 'color': 'black', 'lw': 0.8,
    }
    o_plot.set_plot_dics.append(coll_time_vert)


    Printcolor.blue("\tPlotting ejecta correlations")
    o_plot.main()
    Printcolor.blue("\tPlotting ejecta correlations is done.")
    print('\n')

def plot_nucleo_yields(nuc_class, plot_bern = True, fig_dir ='res_1d/', figname = 'yields.png'):

    Printcolor.blue("\tInitializing ejecta profile plotting.")

    if not os.path.isdir(resdir):
        raise IOError("No simulation results frolder found: {}".format(resdir))
    if not os.path.isdir(resdir + fig_dir):
        Printcolor.yellow("\tCreating {} ".format(fig_dir))
        os.mkdir(resdir + fig_dir)

    o_plot = PLOT_MANY_TASKS()
    o_plot.gen_set["figdir"] = resdir + fig_dir
    o_plot.gen_set["type"] = "cartesian"
    o_plot.gen_set["figsize"] = (6.0, 3.)  # <->, |]
    o_plot.gen_set["figname"] = figname
    o_plot.gen_set["sharex"] = False
    o_plot.gen_set["sharey"] = False
    o_plot.gen_set["subplots_adjust_h"] = 0.3
    o_plot.gen_set["subplots_adjust_w"] = 0.0
    o_plot.set_plot_dics = []

    label1 =  "Dyn."
    label2 =  "Dyn.+Wind"
    fontsize = 16
    labelsize = 16
    ''' --- --- -- SIM 1 --- --- --- '''

    color='red'

    sim_nucleo = {
        'task': 'nucleo', 'ptype': 'cartesian',
        'position': (1, 1),
        'data': nuc_class, 'criterion': '_0', 'method': 'Asol=195',
        'v_n_x': 'A', 'v_n_y': 'Y_final',
        'color': color, 'ls': ':', 'lw': 0.6, 'ds': 'steps', 'alpha': 0.6,
        'ymin': 1e-5, 'ymax': 2e-1, 'xmin': 50, 'xmax': 210,
        'xlabel': r"A", 'ylabel': r'Relative final abundances',
        'label': label1, 'yscale': 'log',
        'fancyticks': True, 'minorticks': True,
        'fontsize': fontsize,
        'labelsize': labelsize,
    }
    o_plot.set_plot_dics.append(sim_nucleo)

    sim_nucleo_b = {
        'task': 'nucleo', 'ptype': 'cartesian',
        'position': (1, 1),
        'data': nuc_class, 'criterion': '_0 _0_b_w', 'method': 'Asol=195',
        'v_n_x': 'A', 'v_n_y': 'Y_final',
        'color': color, 'ls': '-', 'lw': 1.0, 'ds': 'steps', 'alpha': 1.0,
        'ymin': 1e-5, 'ymax': 2e-1, 'xmin': 50, 'xmax': 210,
        'xlabel': r"A", 'ylabel': r'Relative final abundances',
        'label': label2, 'yscale': 'log',
        'fancyticks': True, 'minorticks': True,
        'fontsize': fontsize,
        'labelsize': labelsize,
        'legend':  {}# {'loc': 'best', 'ncol': 2, 'fontsize': 18},
    }
    if plot_bern: o_plot.set_plot_dics.append(sim_nucleo_b)

    sol_yeilds = {
        'task': 'nucleo', 'ptype': 'cartesian',
        'position': (1, 1),
        'data': nuc_class, 'criterion': '_0', 'method': 'sum',
        'v_n_x': 'Asun', 'v_n_y': 'Ysun',
        'color': 'gray', 'marker': 'o', 'ms': 4, 'alpha':0.4,
        'ymin': 8e-5, 'ymax': 8e-1, 'xmin': 50, 'xmax': 210,
        'xlabel': r"A", 'ylabel': r'Relative final abundances',
        'label': 'solar', 'yscale': 'log',
        'fancyticks': True, 'minorticks': True,
        'fontsize': fontsize,
        'labelsize': labelsize,
        'legend': {'loc': 'best', 'ncol': 1, 'fontsize': 14},

    }
    o_plot.set_plot_dics.append(sol_yeilds)
    Printcolor.blue("\tPlotting ejecta nucleo yields")
    o_plot.main()
    Printcolor.blue("\tPlotting ejecta nucleo yeilds is done.")
    print('\n')

''' --- --- --- D2 --- --- --- '''

def movie_2d(d1class, d2class, v_n='rho', rl=3, fig_dir ='res_2d/', rewrite=False):



    rho_dic_xz = {'task': 'slice', 'dtype': '2d rl', 'ptype': 'cartesian',
        'data': d2class, 'it': 00000, 'plane': 'xz', 'rl':rl,
        'position': (1, 1),  # 'title': '[{:.1f} ms]'.format(time_),
        'cbar': {'location': 'right .03 .0', 'label': r'$\rho$ [geo]',  # 'fmt': '%.1e',
                 'labelsize': 14,
                 'fontsize': 14},
        'v_n_x': 'x', 'v_n_y': 'z', 'v_n': 'rho',
        'xmin': None, 'xmax': None, 'ymin': None, 'ymax': None, 'vmin': 1e-10, 'vmax': 1e-4,
        'fill_vmin': False,  # fills the x < vmin with vmin
        'xscale': None, 'yscale': None,
        'mask': None, 'cmap': 'inferno_r', 'norm': "log",
        'fancyticks': False,
        'title': {"text": r'$t-t_{merg}:$' + r'${:.1f}$'.format(0), 'fontsize': 14},
        'sharex': True,  # removes angular citkscitks
        'fontsize': 14,
        'labelsize': 14
        }
    rho_dic_xy = {'task': 'slice', 'dtype': '2d rl', 'ptype': 'cartesian',
        'data': d2class, 'it': 00000, 'plane': 'xy', 'rl':rl,
        'position': (2, 1),  # 'title': '[{:.1f} ms]'.format(time_),
        'cbar': {},
        'v_n_x': 'x', 'v_n_y': 'y', 'v_n': 'rho',
        'xmin': None, 'xmax': None, 'ymin': None, 'ymax': None, 'vmin': 1e-10, 'vmax': 1e-4,
        'fill_vmin': False,  # fills the x < vmin with vmin
        'xscale': None, 'yscale': None,
        'mask': None, 'cmap': 'inferno_r', 'norm': "log",
        'fancyticks': False,
        'title': {},
        'sharex': False,  # removes angular citkscitks
        'fontsize': 14,
        'labelsize': 14
        }

    ye_dic_xz = {'task': 'slice', 'dtype': '2d rl', 'ptype': 'cartesian',
        'data': d2class, 'it': 00000, 'plane': 'xz', 'rl':rl,
        'position': (1, 1),  # 'title': '[{:.1f} ms]'.format(time_),
        'cbar': {'location': 'right .03 .0', 'label': r'$Y_e$ [geo]',  # 'fmt': '%.1e',
                 'labelsize': 14,
                 'fontsize': 14},
        'v_n_x': 'x', 'v_n_y': 'z', 'v_n': 'Y_e',
        'xmin': None, 'xmax': None, 'ymin': None, 'ymax': None, 'vmin': 0.05, 'vmax': 0.45,
        'fill_vmin': False,  # fills the x < vmin with vmin
        'xscale': None, 'yscale': None,
        'mask': None, 'cmap': 'inferno_r', 'norm': "linear",
        'fancyticks': False,
        'title': {"text": r'$t-t_{merg}:$' + r'${:.1f}$'.format(0), 'fontsize': 14},
        'sharex': True,  # removes angular citkscitks
        'fontsize': 14,
        'labelsize': 14
        }
    ye_dic_xy = {'task': 'slice', 'dtype': '2d rl', 'ptype': 'cartesian',
        'data': d2class, 'it': 00000, 'plane': 'xy', 'rl':rl,
        'position': (2, 1),  # 'title': '[{:.1f} ms]'.format(time_),
        'cbar': {},
        'v_n_x': 'x', 'v_n_y': 'y', 'v_n': 'Y_e',
        'xmin': None, 'xmax': None, 'ymin': None, 'ymax': None, 'vmin': 0.05, 'vmax': 0.45,
        'fill_vmin': False,  # fills the x < vmin with vmin
        'xscale': None, 'yscale': None,
        'mask': None, 'cmap': 'inferno_r', 'norm': "linear",
        'fancyticks': False,
        'title': {},
        'sharex': False,  # removes angular citkscitks
        'fontsize': 14,
        'labelsize': 14
        }

    temp_dic_xz = {'task': 'slice', 'dtype': '2d rl', 'ptype': 'cartesian',
        'data': d2class, 'it': 00000, 'plane': 'xz', 'rl':rl,
        'position': (1, 1),  # 'title': '[{:.1f} ms]'.format(time_),
        'cbar': {'location': 'right .03 .0', 'label': r'$Temperature$ [geo]',  # 'fmt': '%.1e',
                 'labelsize': 14,
                 'fontsize': 14},
        'v_n_x': 'x', 'v_n_y': 'z', 'v_n': 'temperature',
        'xmin': None, 'xmax': None, 'ymin': None, 'ymax': None, 'vmin': 1e-2, 'vmax': 1e2,
        'fill_vmin': False,  # fills the x < vmin with vmin
        'xscale': None, 'yscale': None,
        'mask': None, 'cmap': 'inferno_r', 'norm': "log",
        'fancyticks': False,
        'title': {"text": r'$t-t_{merg}:$' + r'${:.1f}$'.format(0), 'fontsize': 14},
        'sharex': True,  # removes angular citkscitks
        'fontsize': 14,
        'labelsize': 14
        }
    temp_dic_xy = {'task': 'slice', 'dtype': '2d rl', 'ptype': 'cartesian',
        'data': d2class, 'it': 00000, 'plane': 'xy', 'rl':rl,
        'position': (2, 1),  # 'title': '[{:.1f} ms]'.format(time_),
        'cbar': {},
        'v_n_x': 'x', 'v_n_y': 'y', 'v_n': 'temperature',
        'xmin': None, 'xmax': None, 'ymin': None, 'ymax': None, 'vmin': 1e-2, 'vmax': 1e2,
        'fill_vmin': False,  # fills the x < vmin with vmin
        'xscale': None, 'yscale': None,
        'mask': None, 'cmap': 'inferno_r', 'norm': "log",
        'fancyticks': False,
        'title': {},
        'sharex': False,  # removes angular citkscitks
        'fontsize': 14,
        'labelsize': 14
        }

    s_dic_xz = {'task': 'slice', 'dtype': '2d rl', 'ptype': 'cartesian',
        'data': d2class, 'it': 00000, 'plane': 'xz', 'rl':rl,
        'position': (1, 1),  # 'title': '[{:.1f} ms]'.format(time_),
        'cbar': {'location': 'right .03 .0', 'label': r'$Entropy$ [geo]',  # 'fmt': '%.1e',
                 'labelsize': 14,
                 'fontsize': 14},
        'v_n_x': 'x', 'v_n_y': 'z', 'v_n': 'entropy',
        'xmin': None, 'xmax': None, 'ymin': None, 'ymax': None, 'vmin': 1e-1, 'vmax': 1e2,
        'fill_vmin': False,  # fills the x < vmin with vmin
        'xscale': None, 'yscale': None,
        'mask': None, 'cmap': 'inferno_r', 'norm': "log",
        'fancyticks': False,
        'title': {"text": r'$t-t_{merg}:$' + r'${:.1f}$'.format(0), 'fontsize': 14},
        'sharex': True,  # removes angular citkscitks
        'fontsize': 14,
        'labelsize': 14
        }
    s_dic_xy = {'task': 'slice', 'dtype': '2d rl', 'ptype': 'cartesian',
        'data': d2class, 'it': 00000, 'plane': 'xy', 'rl':rl,
        'position': (2, 1),  # 'title': '[{:.1f} ms]'.format(time_),
        'cbar': {},
        'v_n_x': 'x', 'v_n_y': 'y', 'v_n': 'entropy',
        'xmin': None, 'xmax': None, 'ymin': None, 'ymax': None, 'vmin': 1e-1, 'vmax': 1e2,
        'fill_vmin': False,  # fills the x < vmin with vmin
        'xscale': None, 'yscale': None,
        'mask': None, 'cmap': 'inferno_r', 'norm': "log",
        'fancyticks': False,
        'title': {},
        'sharex': False,  # removes angular citkscitks
        'fontsize': 14,
        'labelsize': 14
        }

    if v_n == 'rho':
        dic_xy = rho_dic_xy
        dic_xz = rho_dic_xz
    elif v_n == 'Y_e':
        dic_xy = ye_dic_xy
        dic_xz = ye_dic_xz
    elif v_n == 'temperature':
        dic_xy = temp_dic_xy
        dic_xz = temp_dic_xz
    elif v_n == 'entropy':
        dic_xy = s_dic_xy
        dic_xz = s_dic_xz
    else:
        raise NameError("v_n:{} not recogmized".format(v_n))

    """ --- --- --- """

    Printcolor.blue("\tInitializing profile for plotting the movie.")

    if not os.path.isdir(resdir):
        raise IOError("No simulation results frolder found: {}".format(resdir))
    movie_dir = v_n + '_movie/'
    moviename = v_n + ".mp4"
    if not os.path.isdir(resdir + fig_dir):
        Printcolor.blue("\tCreating {} ".format(fig_dir))
        os.mkdir(resdir + fig_dir)

    if not os.path.isdir(resdir + fig_dir + movie_dir):
        Printcolor.blue("\tCreating {} ".format(fig_dir + movie_dir))
        os.mkdir(resdir + fig_dir + movie_dir)

    outfpath = resdir + fig_dir + movie_dir

    # d2class.set_all_it_times_from_outputs(rho_dic_xy["plane"], rho_dic_xy["v_n"])
    d2class.load_all(dic_xy["plane"], dic_xy["v_n"])
    d2class.load_all(dic_xz["plane"], dic_xz["v_n"])
    d2class.set_use_new_output_if_duplicated = True

    o_plot = PLOT_MANY_TASKS()
    o_plot.gen_set["figdir"] = outfpath
    o_plot.gen_set["type"] = "cartesian"
    o_plot.gen_set["figsize"] = (4.2, 8.0)  # <->, |] # to match hists with (8.5, 2.7)
    o_plot.gen_set["figname"] = "0000000.png"
    o_plot.gen_set["sharex"] = False
    o_plot.gen_set["sharey"] = False
    o_plot.gen_set["subplots_adjust_h"] = 0.0
    o_plot.gen_set["subplots_adjust_w"] = 0.2
    o_plot.set_plot_dics = []

    i = 1
    tmerg = d1class.get_par("tmerger_gw")
    for it, t in zip(d2class.iterations, d2class.times):  # zip([346112],[0.020]):# #

        figname = "{0:07d}.png".format(int(it))
        tr = (t - tmerg) * 1e3  # ms
        # print("\nProcessing v_n:{} it:{} ".format(v_n, it) + "{:.2f}ms".format(tr))

        if os.path.isfile(outfpath+figname) and not rewrite:
            Printcolor.blue("Skipping v_n:{} it:{} [{}/{}] t:"
                            .format(v_n, it, i, len(d2class.iterations)) + "{:.2f}ms".format(tr))
            continue
        else:
            Printcolor.blue("\nProcessing v_n:{} it:{} [{}/{}] t:"
                            .format(v_n, it, i, len(d2class.iterations)) + "{:.2f}ms".format(tr))
            dic_xz["it"] = int(it)
            dic_xz["title"]["text"] = r'$t-t_{merg}:$' + r'${:.2f}ms$'.format(float(tr))
            o_plot.gen_set["figname"] = figname  # 7 digit output
            o_plot.set_plot_dics.append(dic_xz)

            dic_xy["it"] = int(it)
            # rho_dic_xy["title"]["text"] = r'$t-t_{merg}:$' + r'${:.2f}ms$'.format(float(tr))
            o_plot.gen_set["figname"] = figname  # 7 digit output
            o_plot.set_plot_dics.append(dic_xy)

            o_plot.main()
            o_plot.set_plot_dics = []
        i = i+1

    if os.path.isfile(outfpath + moviename):
        Printcolor.blue("Rewriting movie: {}".format(outfpath + moviename))
        os.remove(outfpath + moviename)
    else:
        Printcolor.blue("Saving movie: {}".format(outfpath + moviename))
    os.system("ffmpeg -framerate 4 -pattern_type glob -i '{}*.png' {}"
              .format(outfpath, outfpath + moviename))

''' --- --- --- D3 --- --- --- '''




# class D2MOVIES:
#
#     def __init__(self, sim, o_int_data, o_d1_data):
#         self.sim = sim
#         self
#
#         rho_dic = {
#             'task': '2d movie', 'dtype': 'int', 'ptype': 'polar',
#             'data': o_int_data, 'it': 00000, 'plane': 'xy',
#             'position': (1, 1),  # 'title': '[{:.1f} ms]'.format(time_),
#             'cbar': {'location': 'right .03 .0', 'label': r'$\rho$ [geo]',  # 'fmt': '%.1e',
#                      'labelsize': 14,
#                      'fontsize': 14},
#             'v_n_x': 'phi_cyl', 'v_n_y': 'r_cyl', 'v_n': 'rho', 'mod': 'fill_phi',
#             'phimin': None, 'phimax': None, 'rmin': 0, 'rmax': 50, 'vmin': 1e-10, 'vmax': 1e-4,
#             'fill_vmin': False,  # fills the x < vmin with vmin
#             'xscale': None, 'yscale': None,
#             'mask': None, 'cmap': 'inferno_r', 'norm': "log",
#             'fancyticks': False,
#             'title': {"text": r'$t-t_{merg}:$' + r'${:.1f}$'.format(0), 'fontsize': 14},
#             'sharex': True,  # removes angular citkscitks
#         }

#
#
#
#
#
#
# def movie_2d_rho(d1class, d2class, fig_dir ='res_2d/', movie_dir = 'rho_xz/'):
#
#     plane = 'xy'
#     v_n = 'rho'
#     moviename = "rho_xy.mp4"
#
#     Printcolor.blue("\tInitializing profile for plotting the movie.")
#
#     if not os.path.isdir(resdir):
#         raise IOError("No simulation results frolder found: {}".format(resdir))
#
#     if not os.path.isdir(resdir + fig_dir):
#         Printcolor.blue("\tCreating {} ".format(fig_dir))
#         os.mkdir(resdir + fig_dir)
#
#     if not os.path.isdir(resdir + fig_dir + movie_dir):
#         Printcolor.blue("\tCreating {} ".format(fig_dir + movie_dir))
#         os.mkdir(resdir + fig_dir + movie_dir)
#
#     outfpath = resdir + fig_dir + movie_dir
#
#     o_plot = PLOT_MANY_TASKS()
#     o_plot.gen_set["figdir"] = outfpath
#     o_plot.gen_set["type"] = "cartesian"
#     o_plot.gen_set["figsize"] = (4.2, 3.6)  # <->, |] # to match hists with (8.5, 2.7)
#     o_plot.gen_set["figname"] = "0000000.png"
#     o_plot.gen_set["sharex"] = False
#     o_plot.gen_set["sharey"] = False
#     o_plot.gen_set["subplots_adjust_h"] = 0.3
#     o_plot.gen_set["subplots_adjust_w"] = 0.2
#     o_plot.set_plot_dics = []
#
#     rho_dic_xy = {
#         'task': '2d movie xy', 'dtype': 'int', 'ptype': 'polar',
#         'data': d2class, 'it': 00000, 'plane': 'xy',
#         'position': (1, 1),  # 'title': '[{:.1f} ms]'.format(time_),
#         'cbar': {'location': 'right .03 .0', 'label': r'$\rho$ [geo]',  # 'fmt': '%.1e',
#                  'labelsize': 14,
#                  'fontsize': 14},
#         'v_n_x': 'phi_cyl', 'v_n_y': 'r_cyl', 'v_n': 'rho', 'mod': 'fill_phi',
#         'phimin': None, 'phimax': None, 'rmin': 0, 'rmax': 50, 'vmin': 1e-10, 'vmax': 1e-4,
#         'fill_vmin': False,  # fills the x < vmin with vmin
#         'xscale': None, 'yscale': None,
#         'mask': None, 'cmap': 'inferno_r', 'norm': "log",
#         'fancyticks': False,
#         'title': {"text": r'$t-t_{merg}:$' + r'${:.1f}$'.format(0), 'fontsize': 14},
#         'sharex': True,  # removes angular citkscitks
#     }
#
#     rho_dic_xz = {
#         'task': '2d movie xz', 'dtype': 'int', 'ptype': 'cartesian',
#         'data': d2class, 'it': 00000, 'plane': 'xz',
#         'position': (1, 1),  # 'title': '[{:.1f} ms]'.format(time_),
#         'cbar': {'location': 'right .03 .0', 'label': r'$\rho$ [geo]',  # 'fmt': '%.1e',
#                  'labelsize': 14,
#                  'fontsize': 14},
#         'v_n_x': 'phi_cyl', 'v_n_y': 'z_cyl', 'v_n': 'rho', #'mod': 'fill_phi',
#         'xmin': None, 'xmax': None, 'ymin': 0.1, 'ymax': 100, 'vmin': 1e-10, 'vmax': 1e-4,
#         'fill_vmin': False,  # fills the x < vmin with vmin
#         'xscale': None, 'yscale': None,
#         'mask': None, 'cmap': 'inferno_r', 'norm': "log",
#         'fancyticks': False,
#         'title': {"text": r'$t-t_{merg}:$' + r'${:.1f}$'.format(0), 'fontsize': 14},
#         # 'sharex': True,  # removes angular citkscitks
#         'centery': True,
#         'fontsize': 14,
#         'labelsize': 14
#     }
#
#     tmerg = d1class.get_par("tmerger_gw")
#     for it, t in zip([346112],[0.020]):# np.array([346112],dtype=int):#zip(d2class.iterations, d2class.times):  #
#         # --- xy --- #
#         # print('\n')
#         # print("it:{}".format(it))
#         # tr = (t - tmerg) * 1e3  # ms
#         # print("t-t_{merg}:" + "{:.2f}ms".format(tr))
#         # rho_dic_xy["it"] = int(it)
#         # rho_dic_xy["title"]["text"] = r'$t-t_{merg}:$' + r'${:.2f}ms$'.format(float(tr))
#         # o_plot.gen_set["figname"] = "{0:07d}.png".format(int(it))  # 7 digit output
#         # o_plot.set_plot_dics.append(rho_dic_xy)
#         # o_plot.main()
#         # o_plot.set_plot_dics = []
#         # --- xz --- #
#         print('\n')
#         print("it:{}".format(it))
#         tr = (t - tmerg) * 1e3  # ms
#         print("t-t_{merg}:" + "{:.2f}ms".format(tr))
#         rho_dic_xz["it"] = int(it)
#         rho_dic_xz["title"]["text"] = r'$t-t_{merg}:$' + r'${:.2f}ms$'.format(float(tr))
#         o_plot.gen_set["figname"] = "{0:07d}.png".format(int(it))  # 7 digit output
#         o_plot.set_plot_dics.append(rho_dic_xz)
#         o_plot.main()
#         o_plot.set_plot_dics = []
#
#
#         # exit(1)
#         # try:
#         #     plot_dic["title"]["text"] = r'$t-t_{merg}:$'+r'${:.1f}ms$'.format(tr)
#         #     o_plot.gen_set["figname"] = "{0:07d}.png".format(int(it)) # 7 digit output
#         #     o_plot.main()
#         # except ValueError:
#         #     Printcolor.yellow("Warning. it:{} failed with ValueError".format(it))
#
#     if os.path.isfile(outfpath + moviename):
#         Printcolor.blue("Rewriting movie: {}".format(outfpath + moviename))
#         os.remove(outfpath + moviename)
#     else:
#         Printcolor.blue("Saving movie: {}".format(outfpath + moviename))
#     os.system("ffmpeg -framerate 4 -pattern_type glob -i '{}*.png' {}"
#               .format(outfpath, outfpath + moviename))
#
#
# def movie_2d_rho__(d1class, d2class, fig_dir ='res_2d/', movie_dir = 'rho_xy/'):
#
#     plane = 'xy'
#     v_n = 'rho'
#     moviename = "rho_xy.mp4"
#
#     Printcolor.blue("\tInitializing profile for plotting the movie.")
#
#     if not os.path.isdir(resdir):
#         raise IOError("No simulation results frolder found: {}".format(resdir))
#
#     if not os.path.isdir(resdir + fig_dir):
#         Printcolor.blue("\tCreating {} ".format(fig_dir))
#         os.mkdir(resdir + fig_dir)
#
#     if not os.path.isdir(resdir + fig_dir + movie_dir):
#         Printcolor.blue("\tCreating {} ".format(fig_dir + movie_dir))
#         os.mkdir(resdir + fig_dir + movie_dir)
#
#     outfpath = resdir + fig_dir + movie_dir
#
#     o_plot = PLOT_MANY_TASKS()
#     o_plot.gen_set["figdir"] = outfpath
#     o_plot.gen_set["type"] = "polar"
#     o_plot.gen_set["figsize"] = (4.2, 3.6)  # <->, |] # to match hists with (8.5, 2.7)
#     o_plot.gen_set["figname"] = "0000000.png"
#     o_plot.gen_set["sharex"] = False
#     o_plot.gen_set["sharey"] = False
#     o_plot.gen_set["subplots_adjust_h"] = 0.3
#     o_plot.gen_set["subplots_adjust_w"] = 0.2
#     o_plot.set_plot_dics = []
#
#     # load data for all iteration and get a sorted list of iterations
#     d2class.load_all(plane, v_n)
#     d2class.iterations = d2class.get_all_iterations_times(plane, v_n)
#     d2class.iterations = np.array(d2class.iterations, dtype=int)
#     d2class.iterations.sort(axis=0)
#     d2class.set_use_new_output_if_duplicated = True
#
#     plot_dic = {
#         'task': '2d movie', 'dtype': 'int', 'ptype': 'polar',
#         'data': d2class, 'it': 00000, 'plane': plane,
#         'position': (1, 1),# 'title': '[{:.1f} ms]'.format(time_),
#         'cbar': {'location': 'right .03 .0', 'label': r'$\rho$ [geo]', # 'fmt': '%.1e',
#                  'labelsize': 14,
#                  'fontsize': 14},
#         'v_n_x': 'phi_cyl', 'v_n_y': 'r_cyl', 'v_n': v_n, 'mod': 'fill_phi',
#         'phimin': None, 'phimax': None, 'rmin': 0, 'rmax': 50, 'vmin': 1e-10, 'vmax': 1e-4,
#         'fill_vmin': False,  # fills the x < vmin with vmin
#         'xscale': None, 'yscale': None,
#         'mask': None, 'cmap': 'inferno_r', 'norm': "log",
#         'fancyticks': False,
#         'title': {"text":r'$t-t_{merg}:$'+r'${:.1f}$'.format(0), 'fontsize':14},
#         'sharex': True,  # removes angular citkscitks
#     }
#
#     print(d2class.iterations)
#     tmerg = d1class.get_par("tmerger_gw")
#     for it in  np.sort(d2class.iterations): # np.array([346112],dtype=int):#
#         print('\n')
#         print("it:{}".format(it))
#         t = d2class.get_time(it)
#         tr = (t - tmerg) * 1e3  # ms
#         print("t-t_{merg}:" + "{:.2f}ms".format(tr))
#         plot_dic["it"] = int(it)
#         plot_dic["title"]["text"] = r'$t-t_{merg}:$' + r'${:.2f}ms$'.format(float(tr))
#         o_plot.gen_set["figname"] = "{0:07d}.png".format(int(it))  # 7 digit output
#         o_plot.set_plot_dics.append(plot_dic)
#         o_plot.main()
#         o_plot.set_plot_dics = []
#         # exit(1)
#         # try:
#         #     plot_dic["title"]["text"] = r'$t-t_{merg}:$'+r'${:.1f}ms$'.format(tr)
#         #     o_plot.gen_set["figname"] = "{0:07d}.png".format(int(it)) # 7 digit output
#         #     o_plot.main()
#         # except ValueError:
#         #     Printcolor.yellow("Warning. it:{} failed with ValueError".format(it))
#
#     if os.path.isfile(outfpath + moviename):
#         Printcolor.blue("Rewriting movie: {}".format(outfpath + moviename))
#         os.remove(outfpath + moviename)
#     else:
#         Printcolor.blue("Saving movie: {}".format(outfpath + moviename))
#     os.system("ffmpeg -framerate 4 -pattern_type glob -i '{}*.png' {}"
#               .format(outfpath, outfpath + moviename))

if __name__ == '__main__':

    ''' --- --- --- PREPROCESSING --- --- --- '''

    # get_bern_time(simdir)
    # os.system(scrdic + 'analyze.sh' + ' ' + simdir)
    #
    # for item in reslist:
    #     fname = simdir + item
    #     os.system("cp -r " + fname + ' ' + resdir)
    # exit(0)

    ''' --- --- --- ANALYSIS 1D --- --- --- '''

    d1_class = ADD_METHODS_1D(sim)
    # ejecta_profiles(d1_class)
    # ejecta_properties(d1_class)
    # ejecta_correlations_vinf(d1_class)
    # ejecta_correlations_ye(d1_class)
    # nuc_class = NORMALIZE_NUCLEO(sim)
    # plot_nucleo_yields(nuc_class)

    ''' --- --- --- ANALYSIS 2D --- --- --- '''

    from d2analysis import INTERPOLATE_STORE_2D, COMPUTE_STORE
    # d2_class = INTERPOLATE_STORE_2D(sim)
    # d2_class.new_grid_cl.save_grid(sim, "xy")
    # d2_class.new_grid_cl.save_grid(sim, "xz")
    # d2_class.save_all("xy", "rho")
    # d2_class.save_all("xz", "rho")

    # exit(0)

    d2_class = COMPUTE_STORE(sim)

    movie_2d(d1_class, d2_class, rl=3, v_n="rho", rewrite=True)
    movie_2d(d1_class, d2_class, rl=3, v_n="Y_e", rewrite=True)
    movie_2d(d1_class, d2_class, rl=3, v_n="temperature", rewrite=True)
    movie_2d(d1_class, d2_class, rl=3, v_n="entropy", rewrite=True)