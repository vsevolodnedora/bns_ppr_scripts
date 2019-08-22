from __future__ import division
from sys import path
path.append('modules/')
from general import *
from lists import *
from filework import *
# import matplotlib
# from matplotlib import pyplot
import matplotlib.pyplot as plt
plt.rc('text', usetex=True)
plt.rc('font', family='serif')

import copy

from plotting_nexus import PLOT_MANY_TASKS
from d3analysis import LOAD_INT_DATA, LOAD_DENSITY_MODES, LOAD_RES_CORR, ADD_METHODS_FOR_INT_DATA
from d1analysis import COMPUTE_STORE_PAR, NORMALIZE_NUCLEO, ADD_METHODS_1D
from mkn_interface import EXTRACT_LIGHTCURVE, COMBINE_LIGHTCURVES

"""================================================ EJ POFS ========================================================="""

def ejecta_profiles_dd2_ls22_qnot1():

    o_plot = PLOT_MANY_TASKS()
    o_plot.gen_set["figdir"] = Paths.plots
    o_plot.gen_set["type"] = "cartesian"
    o_plot.gen_set["figsize"] = (9.0, 2.7)  # <->, |]
    o_plot.gen_set["figname"] = "ejecta_dd2_dd2_ls220.png"
    o_plot.gen_set["sharex"] = False
    o_plot.gen_set["sharey"] = False
    o_plot.gen_set["subplots_adjust_h"] = 0.3
    o_plot.gen_set["subplots_adjust_w"] = 0.0
    o_plot.set_plot_dics = []

    # label1 = "Dyn"
    # label2 = "Wind"
    label1v = "Dyn."
    label2v = "Wind"
    fontsize = 12
    labelsize = 12



    """ -------------- DD2 LK ---------------------"""

    sim = "LS220_M14691268_M0_LK_SR"
    o_data2 = COMPUTE_STORE_PAR(sim)
    # geo
    dic_ej_prof = {
        'task': 'ejprof', 'ptype': 'cartesian',
        'position': (1, 1),
        'data': o_data2, 'criterion': '_0',
        'v_n_x': 't_tot_flux', 'v_n_y': 'mass_tot_flux',
        'color': 'red', 'ls': '--', 'lw': 1., 'ds': 'default','alpha':1.,
        'ymin': 0, 'ymax': 1.4,
        'xlabel': r"time [ms]", 'ylabel': r'$M_{ej}$ $[10^{-2}M_{\odot}]$',
        'xunits': 'ms', '-t': o_data2.get_par("tmerger_gw"), 'yunits': '1e-2Msun',
        'label': "LS220 14691268 "+label1v, 'yscale': None, 'title': {},
        'fancyticks': True, 'minorticks': True,
        'fontsize': fontsize,
        'labelsize': labelsize
        # 'legend': True
    }
    o_plot.set_plot_dics.append(dic_ej_prof)
    # bern
    dic_ej_prof = {
        'task': 'ejprof', 'ptype': 'cartesian',
        'position': (1, 1),
        'data': o_data2, 'criterion': '_0_b_w',
        'v_n_x': 't_tot_flux', 'v_n_y': 'mass_tot_flux',
        'color': 'red', 'ls': '-', 'lw': 1., 'ds': 'default','alpha':1.,
        # 'ymin': 1e-4, 'ymax': 1e0,
        'xlabel': r"time [ms]", 'ylabel': r'$M_{ej}$ $[10^{-2}M_{\odot}]$',
        'xunits': 'ms', '-t': o_data2.get_par("tmerger_gw"), 'yunits': '1e-2Msun',
        'label': "LS220 14691268 "+label2v, 'yscale': None, 'title': {},
        'fancyticks': True, 'minorticks': True,
        'fontsize': fontsize,
        'labelsize': labelsize
        # 'legend':  True
    }
    o_plot.set_plot_dics.append(dic_ej_prof)

    """ -------------- DD2 LK ---------------------"""

    sim = "LS220_M13641364_M0_LK_SR"
    o_data2 = COMPUTE_STORE_PAR(sim)
    # geo
    dic_ej_prof = {
        'task': 'ejprof', 'ptype': 'cartesian',
        'position': (1, 1),
        'data': o_data2, 'criterion': '_0',
        'v_n_x': 't_tot_flux', 'v_n_y': 'mass_tot_flux',
        'color': 'purple', 'ls': '--', 'lw': 1., 'ds': 'default','alpha':1.,
        'ymin': 0, 'ymax': 1.4,
        'xlabel': r"time [ms]", 'ylabel': r'$M_{ej}$ $[10^{-2}M_{\odot}]$',
        'xunits': 'ms', '-t': o_data2.get_par("tmerger_gw"), 'yunits': '1e-2Msun',
        'label': "LS220 13641364 "+label1v, 'yscale': None, 'title': {},
        'fancyticks': True, 'minorticks': True,
        'fontsize': fontsize,
        'labelsize': labelsize
        # 'legend': True
    }
    o_plot.set_plot_dics.append(dic_ej_prof)
    # bern
    dic_ej_prof = {
        'task': 'ejprof', 'ptype': 'cartesian',
        'position': (1, 1),
        'data': o_data2, 'criterion': '_0_b_w',
        'v_n_x': 't_tot_flux', 'v_n_y': 'mass_tot_flux',
        'color': 'purple', 'ls': '-', 'lw': 1., 'ds': 'default','alpha':1.,
        # 'ymin': 1e-4, 'ymax': 1e0,
        'xlabel': r"time [ms]", 'ylabel': r'$M_{ej}$ $[10^{-2}M_{\odot}]$',
        'xunits': 'ms', '-t': o_data2.get_par("tmerger_gw"), 'yunits': '1e-2Msun',
        'label': "LS220 13641364 "+label2v, 'yscale': None, 'title': {},
        'fancyticks': True, 'minorticks': True,
        'fontsize': fontsize,
        'labelsize': labelsize
        # 'legend':  True
    }
    o_plot.set_plot_dics.append(dic_ej_prof)

    """ -------------- DD2 PI ----------------------"""

    sim = "DD2_M13641364_M0_LK_SR_R04"
    o_data2 = COMPUTE_STORE_PAR(sim)
    # geo
    dic_ej_prof = {
        'task': 'ejprof', 'ptype': 'cartesian',
        'position': (1, 1),
        'data': o_data2, 'criterion': '_0',
        'v_n_x': 't_tot_flux', 'v_n_y': 'mass_tot_flux',
        'color': 'green', 'ls': '--', 'lw': 1., 'ds': 'default','alpha':1.,
        'ymin': 0, 'ymax': 1.4,
        'xlabel': r"time [ms]", 'ylabel': r'$M_{ej}$ $[10^{-2}M_{\odot}]$',
        'xunits': 'ms', '-t': o_data2.get_par("tmerger_gw"), 'yunits': '1e-2Msun',
        'label': "DD2 13641364 "+label1v, 'yscale': None, 'title': {},
        'fontsize': fontsize,
        'labelsize': labelsize,
        'fancyticks': True, 'minorticks': True
    }
    o_plot.set_plot_dics.append(dic_ej_prof)
    # bern
    dic_ej_prof = {
        'task': 'ejprof', 'ptype': 'cartesian',
        'position': (1, 1),
        'data': o_data2, 'criterion': '_0_b_w',
        'v_n_x': 't_tot_flux', 'v_n_y': 'mass_tot_flux',
        'color': 'green', 'ls': '-', 'lw': 1., 'ds': 'default','alpha':1.,
        'xmin':0., 'xmax':100, 'ymin': 0, 'ymax': 1.4,
        'xlabel': r"$t-t_{merg}$ [ms]", 'ylabel': r'$M_{ej}$ $[10^{-2}M_{\odot}]$',
        'xunits': 'ms', '-t': o_data2.get_par("tmerger_gw"), 'yunits': '1e-2Msun',
        'label': "DD2 13641364 "+label2v, 'yscale': None, 'title': {},
        'fancyticks': True, 'minorticks': True,
        # 'sharey': True,
        'fontsize': fontsize,
        'labelsize': labelsize,
        'legend':{'loc': 'center right', 'ncol': 2, 'fontsize':fontsize,
                  'bbox_to_anchor':(1.0, 0.3)}
    }
    o_plot.set_plot_dics.append(dic_ej_prof)
    # coll time


    """ -------------- DD2 PI ----------------------"""

    sim = "DD2_M15091235_M0_LK_SR"
    o_data2 = COMPUTE_STORE_PAR(sim)
    # geo
    dic_ej_prof = {
        'task': 'ejprof', 'ptype': 'cartesian',
        'position': (1, 1),
        'data': o_data2, 'criterion': '_0',
        'v_n_x': 't_tot_flux', 'v_n_y': 'mass_tot_flux',
        'color': 'blue', 'ls': '--', 'lw': 1., 'ds': 'default','alpha':1.,
        'ymin': 0, 'ymax': 1.4,
        'xlabel': r"time [ms]", 'ylabel': r'$M_{ej}$ $[10^{-2}M_{\odot}]$',
        'xunits': 'ms', '-t': o_data2.get_par("tmerger_gw"), 'yunits': '1e-2Msun',
        'label': "DD2 15091235 "+label1v, 'yscale': None, 'title': {},
        'fontsize': fontsize,
        'labelsize': labelsize,
        'fancyticks': True, 'minorticks': True
    }
    o_plot.set_plot_dics.append(dic_ej_prof)
    # bern
    dic_ej_prof = {
        'task': 'ejprof', 'ptype': 'cartesian',
        'position': (1, 1),
        'data': o_data2, 'criterion': '_0_b_w',
        'v_n_x': 't_tot_flux', 'v_n_y': 'mass_tot_flux',
        'color': 'blue', 'ls': '-', 'lw': 1., 'ds': 'default','alpha':1.,
        'xmin':0., 'xmax':90, 'ymin': 0, 'ymax': 2.80,
        'xlabel': r"$t-t_{merg}$ [ms]", 'ylabel': r'$M_{ej}$ $[10^{-2}M_{\odot}]$',
        'xunits': 'ms', '-t': o_data2.get_par("tmerger_gw"), 'yunits': '1e-2Msun',
        'label': "DD2 15091235 "+label2v, 'yscale': None, 'title': {},
        'fancyticks': True, 'minorticks': True,
        # 'sharey': True,
        'fontsize': fontsize,
        'labelsize': labelsize,
        'legend':{'loc': 'center right', 'ncol': 2, 'fontsize':11,
                  'bbox_to_anchor':(0.65, 0.7)} # (1.0, 1.2)
    }
    o_plot.set_plot_dics.append(dic_ej_prof)
    # coll time


    o_plot.main()
    exit(0)
# ejecta_profiles_dd2_ls22_qnot1()

def plot_hists_dd2_ls22_qnot1():

    o_plot = PLOT_MANY_TASKS()
    o_plot.gen_set["figdir"] = Paths.plots
    o_plot.gen_set["type"] = "cartesian"
    o_plot.gen_set["figsize"] = (9.0, 2.7)  # <->, |]
    o_plot.gen_set["figname"] = "hist_dd2_dd2_ls220_ls220.png"
    o_plot.gen_set["sharex"] = False
    o_plot.gen_set["sharey"] = False
    o_plot.gen_set["subplots_adjust_h"] = 0.3
    o_plot.gen_set["subplots_adjust_w"] = 0.0

    label1 =  "Dyn."
    label2 =  "Wind"
    fontsize = 12
    labelsize = 12
    ''' --- --- sim 1 --- --- --- '''

    sim = "LS220_M14691268_M0_LK_SR"
    color = "red"
    o_data = COMPUTE_STORE_PAR(sim)

    dic_hist_theta = {
        'task': 'hist1d', 'ptype': 'cartesian',
        'position': (1, 1),
        'data': o_data, 'norm': True, 'criterion': '_0',
        'v_n_x': 'hist_theta', 'v_n_y': 'hist_theta_m',
        'color': color, 'ls': ':', 'lw': 0.6, 'ds': 'steps', 'alpha':0.6,
        'ymin': 1e-4, 'ymax': 1e0,
        'xlabel': r"Angle from orbital plane", 'ylabel': 'mass',
        'label': None,
        'yscale': 'log', #"LS220 M14691268 "+label1
        'fancyticks': True, 'minorticks': True,
        'fontsize': fontsize,
        'labelsize': labelsize,
        'legend': {}#{'loc': 'best', 'ncol': 2, 'fontsize': 18}
    }
    o_plot.set_plot_dics.append(dic_hist_theta)
    dic_hist_theta_b = {
        'task': 'hist1d', 'ptype': 'cartesian',
        'position': (1, 1),
        'data': o_data, 'norm': True, 'criterion': '_0_b_w',
        'v_n_x': 'hist_theta', 'v_n_y': 'hist_theta_m',
        'color': color, 'ls': '-', 'lw': 1., 'ds': 'steps', 'alpha':1.,
        'ymin': 1e-4, 'ymax': 1e0,
        'xlabel': r"Angle from orbital plane", 'ylabel': 'mass',
        'label': None, #"LS220 M14691268 "+label2,
        'yscale': 'log',
        'fancyticks': True, 'minorticks': True,
        'fontsize': fontsize,
        'labelsize': labelsize,
        'legend': {}#{'loc': 'best', 'ncol': 2, 'fontsize': 18}
    }
    o_plot.set_plot_dics.append(dic_hist_theta_b)

    dic_hist_vel_inf = {
        'task': 'hist1d', 'ptype': 'cartesian',
        'position': (1, 2),
        'data': o_data, 'norm': True, 'criterion': '_0',
        'v_n_x': 'hist_vel_inf', 'v_n_y': 'hist_vel_inf_m',
        'color': color, 'ls': ':', 'lw': 0.6, 'ds': 'steps', 'alpha':0.8,
        'xmin': 0., 'xmax': 1., 'ymin': 1e-4, 'ymax': 1e0,
        'xlabel': r"$\upsilon_{\infty}$ [c]", 'ylabel': 'mass',
        'label': None, #"LS220 M14691268 "+label2,
        'yscale': 'log',
        'fancyticks': True, 'minorticks': True,
        'sharey': True,
        'fontsize': fontsize,
        'labelsize': labelsize,
        'legend': {}#{'loc': 'best', 'ncol': 2, 'fontsize': 18}
    }
    o_plot.set_plot_dics.append(dic_hist_vel_inf)
    dic_hist_vel_inf_b = {
        'task': 'hist1d', 'ptype': 'cartesian',
        'position': (1, 2),
        'data': o_data, 'norm': True, 'criterion': '_0_b_w',
        'v_n_x': 'hist_vel_inf', 'v_n_y': 'hist_vel_inf_m',
        'color': color, 'ls': '-', 'lw': 1., 'ds': 'steps', 'alpha':1.,
        'xmin': 0., 'xmax': 1., 'ymin': 1e-4, 'ymax': 1e0,
        'xlabel': r"$\upsilon_{\infty}$ [c]", 'ylabel': 'mass',
        'label': None,
        'yscale': 'log',
        'fancyticks': True, 'minorticks': True,
        'sharey': True,
        'fontsize': fontsize,
        'labelsize': labelsize,
        'legend': {}# {'loc': 'best', 'ncol': 2, 'fontsize': 18}
    }
    o_plot.set_plot_dics.append(dic_hist_vel_inf_b)

    dic_hist_ye = {
        'task': 'hist1d', 'ptype': 'cartesian',
        'position': (1, 3),
        'data': o_data, 'norm': True, 'criterion': '_0',
        'v_n_x': 'hist_ye', 'v_n_y': 'hist_ye_m',
        'color': color, 'ls': ':', 'lw': 0.6, 'ds': 'steps', 'alpha':0.8,
        'xmin': 0., 'xmax': 0.5, 'ymin': 1e-4, 'ymax': 1e0,
        'xlabel': r"$Y_e$", 'ylabel': 'mass',
        'label': None,
        'yscale': 'log',
        'fancyticks': True, 'minorticks': True,
        'sharey': True,
        'fontsize': fontsize,
        'labelsize': labelsize,
        'legend': {}#{'loc': 'best', 'ncol': 2, 'fontsize': 18}
    }
    o_plot.set_plot_dics.append(dic_hist_ye)
    dic_hist_ye_b = {
        'task': 'hist1d', 'ptype': 'cartesian',
        'position': (1, 3),
        'data': o_data, 'norm': True, 'criterion': '_0_b_w',
        'v_n_x': 'hist_ye', 'v_n_y': 'hist_ye_m',
        'color': color, 'ls': '-', 'lw': 1., 'ds': 'steps', 'alpha':1.,
        'xmin': 0., 'xmax': 0.5, 'ymin': 1e-4, 'ymax': 1e0,
        'xlabel': r"$Y_e$", 'ylabel': 'mass',
        'label': None,
        'yscale': 'log',
        'fancyticks': True, 'minorticks': True,
        'sharey': True,
        'fontsize': fontsize,
        'labelsize': labelsize,
        'legend': {}#{'loc': 'best', 'ncol': 2, 'fontsize': 18}
    }
    o_plot.set_plot_dics.append(dic_hist_ye_b)

    ''' --- --- sim 1 --- --- --- '''

    sim = "LS220_M13641364_M0_LK_SR"
    color = "purple"
    o_data = COMPUTE_STORE_PAR(sim)

    dic_hist_theta = {
        'task': 'hist1d', 'ptype': 'cartesian',
        'position': (1, 1),
        'data': o_data, 'norm': True, 'criterion': '_0',
        'v_n_x': 'hist_theta', 'v_n_y': 'hist_theta_m',
        'color': color, 'ls': ':', 'lw': 0.6, 'ds': 'steps', 'alpha':0.6,
        'ymin': 1e-4, 'ymax': 1e0,
        'xlabel': r"Angle from orbital plane", 'ylabel': 'mass',
        'label': None,
        'yscale': 'log', #"LS220 M14691268 "+label1
        'fancyticks': True, 'minorticks': True,
        'fontsize': fontsize,
        'labelsize': labelsize,
        'legend': {}#{'loc': 'best', 'ncol': 2, 'fontsize': 18}
    }
    o_plot.set_plot_dics.append(dic_hist_theta)
    dic_hist_theta_b = {
        'task': 'hist1d', 'ptype': 'cartesian',
        'position': (1, 1),
        'data': o_data, 'norm': True, 'criterion': '_0_b_w',
        'v_n_x': 'hist_theta', 'v_n_y': 'hist_theta_m',
        'color': color, 'ls': '-', 'lw': 1., 'ds': 'steps', 'alpha':1.,
        'ymin': 1e-4, 'ymax': 1e0,
        'xlabel': r"Angle from orbital plane", 'ylabel': 'mass',
        'label': None, #"LS220 M14691268 "+label2,
        'yscale': 'log',
        'fancyticks': True, 'minorticks': True,
        'fontsize': fontsize,
        'labelsize': labelsize,
        'legend': {}#{'loc': 'best', 'ncol': 2, 'fontsize': 18}
    }
    o_plot.set_plot_dics.append(dic_hist_theta_b)

    dic_hist_vel_inf = {
        'task': 'hist1d', 'ptype': 'cartesian',
        'position': (1, 2),
        'data': o_data, 'norm': True, 'criterion': '_0',
        'v_n_x': 'hist_vel_inf', 'v_n_y': 'hist_vel_inf_m',
        'color': color, 'ls': ':', 'lw': 0.6, 'ds': 'steps', 'alpha':0.8,
        'xmin': 0., 'xmax': 1., 'ymin': 1e-4, 'ymax': 1e0,
        'xlabel': r"$\upsilon_{\infty}$ [c]", 'ylabel': 'mass',
        'label': None,#"LS220 M14691268 "+label2,
        'yscale': 'log',
        'fancyticks': True, 'minorticks': True,
        'sharey': True,
        'fontsize': fontsize,
        'labelsize': labelsize,
        'legend': {}#{'loc': 'best', 'ncol': 2, 'fontsize': 18}
    }
    o_plot.set_plot_dics.append(dic_hist_vel_inf)
    dic_hist_vel_inf_b = {
        'task': 'hist1d', 'ptype': 'cartesian',
        'position': (1, 2),
        'data': o_data, 'norm': True, 'criterion': '_0_b_w',
        'v_n_x': 'hist_vel_inf', 'v_n_y': 'hist_vel_inf_m',
        'color': color, 'ls': '-', 'lw': 1., 'ds': 'steps', 'alpha':1.,
        'xmin': 0., 'xmax': 1., 'ymin': 1e-4, 'ymax': 1e0,
        'xlabel': r"$\upsilon_{\infty}$ [c]", 'ylabel': 'mass',
        'label': None,
        'yscale': 'log',
        'fancyticks': True, 'minorticks': True,
        'sharey': True,
        'fontsize': fontsize,
        'labelsize': labelsize,
        'legend': {}# {'loc': 'best', 'ncol': 2, 'fontsize': 18}
    }
    o_plot.set_plot_dics.append(dic_hist_vel_inf_b)

    dic_hist_ye = {
        'task': 'hist1d', 'ptype': 'cartesian',
        'position': (1, 3),
        'data': o_data, 'norm': True, 'criterion': '_0',
        'v_n_x': 'hist_ye', 'v_n_y': 'hist_ye_m',
        'color': color, 'ls': ':', 'lw': 0.6, 'ds': 'steps', 'alpha':0.8,
        'xmin': 0., 'xmax': 0.5, 'ymin': 1e-4, 'ymax': 1e0,
        'xlabel': r"$Y_e$", 'ylabel': 'mass',
        'label': None,
        'yscale': 'log',
        'fancyticks': True, 'minorticks': True,
        'sharey': True,
        'fontsize': fontsize,
        'labelsize': labelsize,
        'legend': {}#{'loc': 'best', 'ncol': 2, 'fontsize': 18}
    }
    o_plot.set_plot_dics.append(dic_hist_ye)
    dic_hist_ye_b = {
        'task': 'hist1d', 'ptype': 'cartesian',
        'position': (1, 3),
        'data': o_data, 'norm': True, 'criterion': '_0_b_w',
        'v_n_x': 'hist_ye', 'v_n_y': 'hist_ye_m',
        'color': color, 'ls': '-', 'lw': 1., 'ds': 'steps', 'alpha':1.,
        'xmin': 0., 'xmax': 0.5, 'ymin': 1e-4, 'ymax': 1e0,
        'xlabel': r"$Y_e$", 'ylabel': 'mass',
        'label': None,
        'yscale': 'log',
        'fancyticks': True, 'minorticks': True,
        'sharey': True,
        'fontsize': fontsize,
        'labelsize': labelsize,
        'legend': {}#{'loc': 'best', 'ncol': 2, 'fontsize': 18}
    }
    o_plot.set_plot_dics.append(dic_hist_ye_b)

    ''' --- --- sim 1 --- --- --- '''

    sim = "DD2_M13641364_M0_LK_SR_R04"
    color = "green"
    o_data = COMPUTE_STORE_PAR(sim)

    dic_hist_theta = {
        'task': 'hist1d', 'ptype': 'cartesian',
        'position': (1, 1),
        'data': o_data, 'norm': True, 'criterion': '_0',
        'v_n_x': 'hist_theta', 'v_n_y': 'hist_theta_m',
        'color': color, 'ls': ':', 'lw': 0.6, 'ds': 'steps', 'alpha':0.6,
        'ymin': 1e-4, 'ymax': 1e0,
        'xlabel': r"Angle from orbital plane", 'ylabel': 'mass',
        'label': None,
        'yscale': 'log', #"LS220 M14691268 "+label1
        'fancyticks': True, 'minorticks': True,
        'fontsize': fontsize,
        'labelsize': labelsize,
        'legend': {}#{'loc': 'best', 'ncol': 2, 'fontsize': 18}
    }
    o_plot.set_plot_dics.append(dic_hist_theta)
    dic_hist_theta_b = {
        'task': 'hist1d', 'ptype': 'cartesian',
        'position': (1, 1),
        'data': o_data, 'norm': True, 'criterion': '_0_b_w',
        'v_n_x': 'hist_theta', 'v_n_y': 'hist_theta_m',
        'color': color, 'ls': '-', 'lw': 1., 'ds': 'steps', 'alpha':1.,
        'ymin': 1e-4, 'ymax': 1e0,
        'xlabel': r"Angle from orbital plane", 'ylabel': 'mass',
        'label': None, #"LS220 M14691268 "+label2,
        'yscale': 'log',
        'fancyticks': True, 'minorticks': True,
        'fontsize': fontsize,
        'labelsize': labelsize,
        'legend': {}#{'loc': 'best', 'ncol': 2, 'fontsize': 18}
    }
    o_plot.set_plot_dics.append(dic_hist_theta_b)

    dic_hist_vel_inf = {
        'task': 'hist1d', 'ptype': 'cartesian',
        'position': (1, 2),
        'data': o_data, 'norm': True, 'criterion': '_0',
        'v_n_x': 'hist_vel_inf', 'v_n_y': 'hist_vel_inf_m',
        'color': color, 'ls': ':', 'lw': 0.6, 'ds': 'steps', 'alpha':0.8,
        'xmin': 0., 'xmax': 1., 'ymin': 1e-4, 'ymax': 1e0,
        'xlabel': r"$\upsilon_{\infty}$ [c]", 'ylabel': 'mass',
        'label': None,#"LS220 M14691268 "+label2,
        'yscale': 'log',
        'fancyticks': True, 'minorticks': True,
        'sharey': True,
        'fontsize': fontsize,
        'labelsize': labelsize,
        'legend': {}#{'loc': 'best', 'ncol': 2, 'fontsize': 18}
    }
    o_plot.set_plot_dics.append(dic_hist_vel_inf)
    dic_hist_vel_inf_b = {
        'task': 'hist1d', 'ptype': 'cartesian',
        'position': (1, 2),
        'data': o_data, 'norm': True, 'criterion': '_0_b_w',
        'v_n_x': 'hist_vel_inf', 'v_n_y': 'hist_vel_inf_m',
        'color': color, 'ls': '-', 'lw': 1., 'ds': 'steps', 'alpha':1.,
        'xmin': 0., 'xmax': 1., 'ymin': 1e-4, 'ymax': 1e0,
        'xlabel': r"$\upsilon_{\infty}$ [c]", 'ylabel': 'mass',
        'label': None,
        'yscale': 'log',
        'fancyticks': True, 'minorticks': True,
        'sharey': True,
        'fontsize': fontsize,
        'labelsize': labelsize,
        'legend': {}# {'loc': 'best', 'ncol': 2, 'fontsize': 18}
    }
    o_plot.set_plot_dics.append(dic_hist_vel_inf_b)

    dic_hist_ye = {
        'task': 'hist1d', 'ptype': 'cartesian',
        'position': (1, 3),
        'data': o_data, 'norm': True, 'criterion': '_0',
        'v_n_x': 'hist_ye', 'v_n_y': 'hist_ye_m',
        'color': color, 'ls': ':', 'lw': 0.6, 'ds': 'steps', 'alpha':0.8,
        'xmin': 0., 'xmax': 0.5, 'ymin': 1e-4, 'ymax': 1e0,
        'xlabel': r"$Y_e$", 'ylabel': 'mass',
        'label': None,
        'yscale': 'log',
        'fancyticks': True, 'minorticks': True,
        'sharey': True,
        'fontsize': fontsize,
        'labelsize': labelsize,
        'legend': {}#{'loc': 'best', 'ncol': 2, 'fontsize': 18}
    }
    o_plot.set_plot_dics.append(dic_hist_ye)
    dic_hist_ye_b = {
        'task': 'hist1d', 'ptype': 'cartesian',
        'position': (1, 3),
        'data': o_data, 'norm': True, 'criterion': '_0_b_w',
        'v_n_x': 'hist_ye', 'v_n_y': 'hist_ye_m',
        'color': color, 'ls': '-', 'lw': 1., 'ds': 'steps', 'alpha':1.,
        'xmin': 0., 'xmax': 0.5, 'ymin': 1e-4, 'ymax': 1e0,
        'xlabel': r"$Y_e$", 'ylabel': 'mass',
        'label': None,
        'yscale': 'log',
        'fancyticks': True, 'minorticks': True,
        'sharey': True,
        'fontsize': fontsize,
        'labelsize': labelsize,
        'legend': {}#{'loc': 'best', 'ncol': 2, 'fontsize': 18}
    }
    o_plot.set_plot_dics.append(dic_hist_ye_b)

    ''' --- --- sim 2 --- --- --- '''

    sim = "DD2_M15091235_M0_LK_SR"
    color = "blue"
    o_data = COMPUTE_STORE_PAR(sim)

    dic_hist_theta = {
        'task': 'hist1d', 'ptype': 'cartesian',
        'position': (1, 1),
        'data': o_data, 'norm': True, 'criterion': '_0',
        'v_n_x': 'hist_theta', 'v_n_y': 'hist_theta_m',
        'color': color, 'ls': ':', 'lw': 0.6, 'ds': 'steps', 'alpha':0.6,
        'ymin': 1e-4, 'ymax': 1e0,
        'xlabel': r"Angle from orbital plane", 'ylabel': 'mass',
        'label': None,
        'yscale': 'log', #"LS220 M14691268 "+label1
        'fancyticks': True, 'minorticks': True,
        'fontsize': fontsize,
        'labelsize': labelsize,
        'legend': {}#{'loc': 'best', 'ncol': 2, 'fontsize': 18}
    }
    o_plot.set_plot_dics.append(dic_hist_theta)
    dic_hist_theta_b = {
        'task': 'hist1d', 'ptype': 'cartesian',
        'position': (1, 1),
        'data': o_data, 'norm': True, 'criterion': '_0_b_w',
        'v_n_x': 'hist_theta', 'v_n_y': 'hist_theta_m',
        'color': color, 'ls': '-', 'lw': 1., 'ds': 'steps', 'alpha':1.,
        'ymin': 1e-4, 'ymax': 1e0,
        'xlabel': r"Angle from orbital plane", 'ylabel': 'mass',
        'label': None, #"LS220 M14691268 "+label2,
        'yscale': 'log',
        'fancyticks': True, 'minorticks': True,
        'fontsize': fontsize,
        'labelsize': labelsize,
        'legend': {}#{'loc': 'best', 'ncol': 2, 'fontsize': 18}
    }
    o_plot.set_plot_dics.append(dic_hist_theta_b)

    dic_hist_vel_inf = {
        'task': 'hist1d', 'ptype': 'cartesian',
        'position': (1, 2),
        'data': o_data, 'norm': True, 'criterion': '_0',
        'v_n_x': 'hist_vel_inf', 'v_n_y': 'hist_vel_inf_m',
        'color': color, 'ls': ':', 'lw': 0.6, 'ds': 'steps', 'alpha':0.8,
        'xmin': 0., 'xmax': 1., 'ymin': 1e-4, 'ymax': 1e0,
        'xlabel': r"$\upsilon_{\infty}$ [c]", 'ylabel': 'mass',
        'label': None,#"LS220 M14691268 "+label2,
        'yscale': 'log',
        'fancyticks': True, 'minorticks': True,
        'sharey': True,
        'fontsize': fontsize,
        'labelsize': labelsize,
        'legend': {}#{'loc': 'best', 'ncol': 2, 'fontsize': 18}
    }
    o_plot.set_plot_dics.append(dic_hist_vel_inf)
    dic_hist_vel_inf_b = {
        'task': 'hist1d', 'ptype': 'cartesian',
        'position': (1, 2),
        'data': o_data, 'norm': True, 'criterion': '_0_b_w',
        'v_n_x': 'hist_vel_inf', 'v_n_y': 'hist_vel_inf_m',
        'color': color, 'ls': '-', 'lw': 1., 'ds': 'steps', 'alpha':1.,
        'xmin': 0., 'xmax': 1., 'ymin': 1e-4, 'ymax': 1e0,
        'xlabel': r"$\upsilon_{\infty}$ [c]", 'ylabel': 'mass',
        'label': None,
        'yscale': 'log',
        'fancyticks': True, 'minorticks': True,
        'sharey': True,
        'fontsize': fontsize,
        'labelsize': labelsize,
        'legend': {}# {'loc': 'best', 'ncol': 2, 'fontsize': 18}
    }
    o_plot.set_plot_dics.append(dic_hist_vel_inf_b)

    dic_hist_ye = {
        'task': 'hist1d', 'ptype': 'cartesian',
        'position': (1, 3),
        'data': o_data, 'norm': True, 'criterion': '_0',
        'v_n_x': 'hist_ye', 'v_n_y': 'hist_ye_m',
        'color': color, 'ls': ':', 'lw': 0.6, 'ds': 'steps', 'alpha':0.8,
        'xmin': 0., 'xmax': 0.5, 'ymin': 1e-4, 'ymax': 1e0,
        'xlabel': r"$Y_e$", 'ylabel': 'mass',
        'label': None,
        'yscale': 'log',
        'fancyticks': True, 'minorticks': True,
        'sharey': True,
        'fontsize': fontsize,
        'labelsize': labelsize,
        'legend': {}#{'loc': 'best', 'ncol': 2, 'fontsize': 18}
    }
    o_plot.set_plot_dics.append(dic_hist_ye)
    dic_hist_ye_b = {
        'task': 'hist1d', 'ptype': 'cartesian',
        'position': (1, 3),
        'data': o_data, 'norm': True, 'criterion': '_0_b_w',
        'v_n_x': 'hist_ye', 'v_n_y': 'hist_ye_m',
        'color': color, 'ls': '-', 'lw': 1., 'ds': 'steps', 'alpha':1.,
        'xmin': 0., 'xmax': 0.5, 'ymin': 1e-4, 'ymax': 1e0,
        'xlabel': r"$Y_e$", 'ylabel': 'mass',
        'label': None,
        'yscale': 'log',
        'fancyticks': True, 'minorticks': True,
        'sharey': True,
        'fontsize': fontsize,
        'labelsize': labelsize,
        'legend': {}#{'loc': 'best', 'ncol': 2, 'fontsize': 18}
    }
    o_plot.set_plot_dics.append(dic_hist_ye_b)

    o_plot.main()
    exit(1)
# plot_hists_dd2_ls22_qnot1()

def plot_m1_dd2_ls220_qnot1():

    o_plot = PLOT_MANY_TASKS()
    o_plot.gen_set["figdir"] = Paths.plots
    o_plot.gen_set["type"] = "cartesian"
    o_plot.gen_set["figsize"] = (9.0, 2.7)  # <->, |]
    o_plot.gen_set["figname"] = "dm_dd2_dd2_ls220.png"
    o_plot.gen_set["sharex"] = True
    o_plot.gen_set["sharey"] = False
    o_plot.set_plot_dics = []

    fname  =  "density_modes_lap15.h5"
    fontsize = 12
    labelsize = 12

    " --- --- --- --- --- --- ---- --- --- --- --- --- -----"
    # ------------------------------------------------------#
    " --- --- --- --- --- SIM0 --- --- --- ---  --- --- --- "

    sim = "LS220_M14691268_M0_LK_SR"
    o_1d_data = COMPUTE_STORE_PAR(sim)

    o_dm_data = LOAD_DENSITY_MODES(sim)
    o_dm_data.gen_set['fname'] = Paths.ppr_sims + sim + '/res_3d/' + fname
    color='red'
    densmode_m1 = {
        'task': 'line', 'dtype': 'dm', 'ptype': 'cartesian',
        'data': o_dm_data,
        'position': (1, 1),
        'v_n_x': 'times', 'v_n_y': 'int_phi_r abs',
        'xunits': 'ms', '-t': o_1d_data.get_par("tmerger_gw"),
        'mode': 1,  'norm_to_m': 0,
        'ls': '-', 'color': color, 'lw': 1., 'ds': 'default','alpha':1.,
        'ylabel': r'$C_m/C_0$ Magnitude', 'xlabel': r'time [ms]',
        'label': None, #r'DD2 $m=1$',
        'xmin': None, 'xmax': None, 'ymin': None, 'ymax': None,
        'xscale': None, 'yscale': 'log', 'legend': {},
        'fancyticks': True, 'minorticks': True,
        'fontsize': fontsize,
        'labelsize': labelsize,
    }
    densmode_m2 = {
        'task': 'line', 'dtype': 'dm', 'ptype': 'cartesian',
        'data': o_dm_data,
        'position': (1, 1),
        'v_n_x': 'times', 'v_n_y': 'int_phi_r abs',
        'xunits': 'ms', '-t': o_1d_data.get_par("tmerger_gw"),
        'mode': 2, 'norm_to_m': 0,
        'ls': ':', 'color': color, 'lw': 0.8, 'ds': 'default', 'alpha': 1.,
        'ylabel': r'$C_m/C_0$ Magnitude', 'xlabel': r'time [ms]',
        'label': None, #r'DD2 $m=2$',
        'xmin': None, 'xmax': None, 'ymin': None, 'ymax': None,
        'xscale': None, 'yscale': 'log', 'legend': {},
        'fancyticks': True, 'minorticks': True,
        'fontsize': fontsize,
        'labelsize': labelsize,
    }
    # o_plot.set_plot_dics.append(densmode_m0)
    o_plot.set_plot_dics.append(densmode_m1)
    o_plot.set_plot_dics.append(densmode_m2)

    " --- --- --- --- --- --- ---- --- --- --- --- --- -----"
    # ------------------------------------------------------#
    " --- --- --- --- --- SIM0 --- --- --- ---  --- --- --- "

    sim = "DD2_M13641364_M0_LK_SR_R04"
    o_1d_data = COMPUTE_STORE_PAR(sim)

    o_dm_data = LOAD_DENSITY_MODES(sim)
    o_dm_data.gen_set['fname'] = Paths.ppr_sims + sim + '/res_3d/' + fname
    color = 'green'
    densmode_m1 = {
        'task': 'line', 'dtype': 'dm', 'ptype': 'cartesian',
        'data': o_dm_data,
        'position': (1, 1),
        'v_n_x': 'times', 'v_n_y': 'int_phi_r abs',
        'xunits': 'ms', '-t': o_1d_data.get_par("tmerger_gw"),
        'mode': 1, 'norm_to_m': 0,
        'ls': '-', 'color': color, 'lw': 1., 'ds': 'default', 'alpha': 1.,
        'ylabel': r'$C_m/C_0$ Magnitude', 'xlabel': r'time [ms]',
        'label': None,  # r'DD2 $m=1$',
        'xmin': None, 'xmax': None, 'ymin': None, 'ymax': None,
        'xscale': None, 'yscale': 'log', 'legend': {},
        'fancyticks': True, 'minorticks': True,
        'fontsize': fontsize,
        'labelsize': labelsize,
    }
    densmode_m2 = {
        'task': 'line', 'dtype': 'dm', 'ptype': 'cartesian',
        'data': o_dm_data,
        'position': (1, 1),
        'v_n_x': 'times', 'v_n_y': 'int_phi_r abs',
        'xunits': 'ms', '-t': o_1d_data.get_par("tmerger_gw"),
        'mode': 2, 'norm_to_m': 0,
        'ls': ':', 'color': color, 'lw': 0.8, 'ds': 'default', 'alpha': 1.,
        'ylabel': r'$C_m/C_0$ Magnitude', 'xlabel': r'time [ms]',
        'label': None,  # r'DD2 $m=2$',
        'xmin': None, 'xmax': None, 'ymin': None, 'ymax': None,
        'xscale': None, 'yscale': 'log', 'legend': {},
        'fancyticks': True, 'minorticks': True,
        'fontsize': fontsize,
        'labelsize': labelsize,
    }
    # o_plot.set_plot_dics.append(densmode_m0)
    o_plot.set_plot_dics.append(densmode_m1)
    o_plot.set_plot_dics.append(densmode_m2)

    " --- --- --- --- --- --- ---- --- --- --- --- --- -----"
    # ------------------------------------------------------#
    " --- --- --- --- --- SIM0 --- --- --- ---  --- --- --- "

    sim = "DD2_M15091235_M0_LK_SR"
    o_1d_data = COMPUTE_STORE_PAR(sim)

    o_dm_data = LOAD_DENSITY_MODES(sim)
    o_dm_data.gen_set['fname'] = Paths.ppr_sims + sim + '/res_3d/' + fname
    color = 'blue'
    densmode_m1 = {
        'task': 'line', 'dtype': 'dm', 'ptype': 'cartesian',
        'data': o_dm_data,
        'position': (1, 1),
        'v_n_x': 'times', 'v_n_y': 'int_phi_r abs',
        'xunits': 'ms', '-t': o_1d_data.get_par("tmerger_gw"),
        'mode': 1, 'norm_to_m': 0,
        'ls': '-', 'color': color, 'lw': 1., 'ds': 'default', 'alpha': 1.,
        'ylabel': r'$C_m/C_0$ Magnitude', 'xlabel': r'time [ms]',
        'label': None,  # r'DD2 $m=1$',
        'xmin': None, 'xmax': None, 'ymin': None, 'ymax': None,
        'xscale': None, 'yscale': 'log', 'legend': {},
        'fancyticks': True, 'minorticks': True,
        'fontsize': fontsize,
        'labelsize': labelsize,
    }
    densmode_m2 = {
        'task': 'line', 'dtype': 'dm', 'ptype': 'cartesian',
        'data': o_dm_data,
        'position': (1, 1),
        'v_n_x': 'times', 'v_n_y': 'int_phi_r abs',
        'xunits': 'ms', '-t': o_1d_data.get_par("tmerger_gw"),
        'mode': 2, 'norm_to_m': 0,
        'ls': ':', 'color': color, 'lw': 0.8, 'ds': 'default', 'alpha': 1.,
        'ylabel': r'$C_m/C_0$ Magnitude', 'xlabel': r'time [ms]',
        'label': None,  # r'DD2 $m=2$',
        'xmin': None, 'xmax': None, 'ymin': 1e-3, 'ymax': 1e-1,
        'xscale': None, 'yscale': 'log', 'legend': {},
        'fancyticks': True, 'minorticks': True,
        'fontsize': fontsize,
        'labelsize': labelsize,
    }
    # o_plot.set_plot_dics.append(densmode_m0)
    o_plot.set_plot_dics.append(densmode_m1)
    o_plot.set_plot_dics.append(densmode_m2)

    o_plot.main()
    exit(1)
# plot_m1_dd2_ls220_qnot1()


"""================================================= NUCLEO ========================================================="""

def plot_nucleo_mult_sims_ls220():
    o_plot = PLOT_MANY_TASKS()
    o_plot.gen_set["figdir"] = Paths.plots + "all_ls220_SR/"
    o_plot.gen_set["type"] = "cartesian"
    o_plot.gen_set["figsize"] = (9.0, 6.0)  # <->, |]
    o_plot.gen_set["figname"] = "nucleo_mult_ls220_wind_sims.png"
    o_plot.gen_set["sharex"] = False
    o_plot.gen_set["sharey"] = False
    o_plot.gen_set["subplots_adjust_h"] = 0.3
    o_plot.gen_set["subplots_adjust_w"] = 0.0
    o_plot.set_plot_dics = []

    sims_dd2 = ["DD2_M13641364_M0_SR_R04",
                "DD2_M13641364_M0_LK_SR_R04",
                "DD2_M13641364_M0_SR",
                "DD2_M14971245_M0_SR",
                "DD2_M15091235_M0_LK_SR"]
    labels_dd2 = ['136136', '136136 LK', "136136 R03", "150125", "151124 LK"]

    sims = [
        "LS220_M13641364_M0_SR", # FULL 3D!
        "LS220_M14001330_M0_SR",
        "LS220_M14351298_M0_SR",
        "LS220_M14691268_M0_SR",

        "LS220_M13641364_M0_LK_SR",
        "LS220_M14691268_M0_LK_SR",
    ]
    labels = ['136136', '140133', "144130", "147127", "136136 LK", "147127 LK"]

    colors = ['red', 'orange', "pink", "magenta", "lightblue", "blue"]
    lss = ['--', '--', "--", '--', "--", "--"]
    dyns = [False, False, False, False, False, False] # [False, False, False, False, False]
    berns = [True, True, True, True, True, True] # [True, True, True, True, True]
    alphas = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0]
    lws = [0.8, 0.8, 0.8, 0.8, 0.8, 1.0]

    dyn_lbl =  " D."
    wind_lbl =  " DW."
    fontsize = 12
    labelsize = 12

    for sim, label, color, alpha, ls, if_bern, if_dyn in \
            zip(sims, labels, colors, alphas, lss, berns, dyns):
        o_data = NORMALIZE_NUCLEO(sim)

        if if_dyn:
            sim_nucleo = {
                'task': 'nucleo', 'ptype': 'cartesian',
                'position': (1, 1),
                'data': o_data, 'criterion': '_0', 'method': 'Asol=195',
                'v_n_x': 'A', 'v_n_y': 'Y_final',
                'color': color, 'ls': ls, 'lw': 1.0, 'ds': 'steps', 'alpha': 1.0,
                'ymin': 1e-5, 'ymax': 2e-1, 'xmin': 50, 'xmax': 210,
                'xlabel': r"A", 'ylabel': r'Relative final abundances',
                'label': label + dyn_lbl, 'yscale': 'log',
                'fancyticks': True, 'minorticks': True,
                'legend': {},
                'fontsize': 14,
                'labelsize': 14,
            }
            o_plot.set_plot_dics.append(sim_nucleo)
        if if_bern:
            sim_nucleo_b = {
                'task': 'nucleo', 'ptype': 'cartesian',
                'position': (1, 1),
                'data': o_data, 'criterion': '_0 _0_b_w', 'method': 'Asol=195',
                'v_n_x': 'A', 'v_n_y': 'Y_final',
                'color': color, 'ls': '-', 'lw': 1.0, 'ds': 'steps', 'alpha': 1.0,
                'ymin': 1e-5, 'ymax': 2e-1, 'xmin': 50, 'xmax': 210,
                'xlabel': r"A", 'ylabel': r'Relative final abundances',
                'label': label + wind_lbl, 'yscale': 'log',
                'fancyticks': True, 'minorticks': True,
                'legend': {},
                'fontsize': 14,
                'labelsize': 14,
            }
            o_plot.set_plot_dics.append(sim_nucleo_b)

        if sim == sims[-1]:
            # sol_yeilds['legend'] = {'loc': 'best', 'ncol': 2, 'fontsize': 14}
            sol_yeilds = {
                'task': 'nucleo', 'ptype': 'cartesian',
                'position': (1, 1),
                'data': o_data, 'criterion': '_0', 'method': 'sum',
                'v_n_x': 'Asun', 'v_n_y': 'Ysun',
                'color': 'gray', 'marker': 'o', 'ms': 4, 'alpha': alpha,
                'ymin': 8e-5, 'ymax': 8e-1, 'xmin': 50, 'xmax': 210,
                'xlabel': r"Mass number, A", 'ylabel': r'Relative final abundances',
                'label': 'solar', 'yscale': 'log',
                'fancyticks': True, 'minorticks': True,
                'legend':  {'loc': 'upper right', 'ncol': 3, 'fontsize': 14},
                'fontsize': 14,
                'labelsize': 14,
            }
            o_plot.set_plot_dics.append(sol_yeilds)

    o_plot.main()
    exit(1)
# plot_nucleo_mult_sims_ls220()
def plot_nucleo_mult_sims_dd2():
    o_plot = PLOT_MANY_TASKS()
    o_plot.gen_set["figdir"] = Paths.plots + "all_dd2_SR/"
    o_plot.gen_set["type"] = "cartesian"
    o_plot.gen_set["figsize"] = (9.0, 6.0)  # <->, |]
    o_plot.gen_set["figname"] = "nucleo_mult_dd2_dyn_sims.png"
    o_plot.gen_set["sharex"] = False
    o_plot.gen_set["sharey"] = False
    o_plot.gen_set["subplots_adjust_h"] = 0.3
    o_plot.gen_set["subplots_adjust_w"] = 0.0
    o_plot.set_plot_dics = []

    sims = ["DD2_M13641364_M0_SR_R04",
                "DD2_M13641364_M0_LK_SR_R04",
                "DD2_M13641364_M0_SR",
                "DD2_M14971245_M0_SR",
                "DD2_M15091235_M0_LK_SR"]
    labels = ['136136', '136136 LK', "136136 R03", "150125", "151124 LK"]
    colors = ['red', 'orange', "pink", "lightblue", "blue"]
    lss = ['--', '--', "--", '--', "--"]
    dyns = [True, True, True, True, True] # [False, False, False, False, False]
    berns = [False, False, False, False, False] # [True, True, True, True, True]
    alphas = [1.0, 1.0, 1.0, 1.0, 1.0]
    lws = [0.8, 0.8, 0.8, 0.8, 1.0]

    dyn_lbl =  " D."
    wind_lbl =  " DW."
    fontsize = 12
    labelsize = 12

    for sim, label, color, alpha, ls, if_bern, if_dyn in \
            zip(sims, labels, colors, alphas, lss, berns, dyns):
        o_data = NORMALIZE_NUCLEO(sim)

        if if_dyn:
            sim_nucleo = {
                'task': 'nucleo', 'ptype': 'cartesian',
                'position': (1, 1),
                'data': o_data, 'criterion': '_0', 'method': 'Asol=195',
                'v_n_x': 'A', 'v_n_y': 'Y_final',
                'color': color, 'ls': ls, 'lw': 1.0, 'ds': 'steps', 'alpha': 1.0,
                'ymin': 1e-5, 'ymax': 2e-1, 'xmin': 50, 'xmax': 210,
                'xlabel': r"A", 'ylabel': r'Relative final abundances',
                'label': label + dyn_lbl, 'yscale': 'log',
                'fancyticks': True, 'minorticks': True,
                'legend': {},
                'fontsize': 14,
                'labelsize': 14,
            }
            o_plot.set_plot_dics.append(sim_nucleo)
        if if_bern:
            sim_nucleo_b = {
                'task': 'nucleo', 'ptype': 'cartesian',
                'position': (1, 1),
                'data': o_data, 'criterion': '_0 _0_b_w', 'method': 'Asol=195',
                'v_n_x': 'A', 'v_n_y': 'Y_final',
                'color': color, 'ls': '-', 'lw': 1.0, 'ds': 'steps', 'alpha': 1.0,
                'ymin': 1e-5, 'ymax': 2e-1, 'xmin': 50, 'xmax': 210,
                'xlabel': r"A", 'ylabel': r'Relative final abundances',
                'label': label + wind_lbl, 'yscale': 'log',
                'fancyticks': True, 'minorticks': True,
                'legend': {},
                'fontsize': 14,
                'labelsize': 14,
            }
            o_plot.set_plot_dics.append(sim_nucleo_b)

        if sim == sims[-1]:
            # sol_yeilds['legend'] = {'loc': 'best', 'ncol': 2, 'fontsize': 14}
            sol_yeilds = {
                'task': 'nucleo', 'ptype': 'cartesian',
                'position': (1, 1),
                'data': o_data, 'criterion': '_0', 'method': 'sum',
                'v_n_x': 'Asun', 'v_n_y': 'Ysun',
                'color': 'gray', 'marker': 'o', 'ms': 4, 'alpha': alpha,
                'ymin': 8e-5, 'ymax': 8e-1, 'xmin': 50, 'xmax': 210,
                'xlabel': r"Mass number, A", 'ylabel': r'Relative final abundances',
                'label': 'solar', 'yscale': 'log',
                'fancyticks': True, 'minorticks': True,
                'legend':  {'loc': 'upper right', 'ncol': 3, 'fontsize': 14},
                'fontsize': 14,
                'labelsize': 14,
            }
            o_plot.set_plot_dics.append(sol_yeilds)

    o_plot.main()
    exit(1)
# plot_nucleo_mult_sims_dd2()
def plot_nucleo_mult_sims_sfho():
    o_plot = PLOT_MANY_TASKS()
    o_plot.gen_set["figdir"] = Paths.plots + "all_sfho_SR/"
    o_plot.gen_set["type"] = "cartesian"
    o_plot.gen_set["figsize"] = (9.0, 6.0)  # <->, |]
    o_plot.gen_set["figname"] = "nucleo_mult_sfho_wind_sims.png"
    o_plot.gen_set["sharex"] = False
    o_plot.gen_set["sharey"] = False
    o_plot.gen_set["subplots_adjust_h"] = 0.3
    o_plot.gen_set["subplots_adjust_w"] = 0.0
    o_plot.set_plot_dics = []

    sims = [
        "SFHo_M13641364_M0_SR", # 21ms
        "SFHo_M13641364_M0_LK_SR_2019pizza",  # 55ms
        "SFHo_M14521283_M0_LK_SR",# 30ms
        "SFHo_M14521283_M0_LK_SR_2019pizza",# 26ms
        "SFHo_M14521283_M0_SR", # 33ms
    ]
    labels = ['136136', '136136 LK 19', "145128 LK", "145128 LK 19", "145128"]
    colors = ['red', 'orange', "pink", "lightblue", "blue"]
    lss = ['--', '--', "--", '--', "--"]
    dyns =[False, False, False, False, False] # [False, False, False, False, False]
    berns =[True, True, True, True, True]# [True, True, True, True, True]
    alphas = [1.0, 1.0, 1.0, 1.0, 1.0]
    lws = [0.8, 0.8, 0.8, 0.8, 1.0]

    dyn_lbl =  " D."
    wind_lbl =  " DW."
    fontsize = 12
    labelsize = 12

    for sim, label, color, alpha, ls, if_bern, if_dyn in \
            zip(sims, labels, colors, alphas, lss, berns, dyns):
        o_data = NORMALIZE_NUCLEO(sim)

        if if_dyn:
            sim_nucleo = {
                'task': 'nucleo', 'ptype': 'cartesian',
                'position': (1, 1),
                'data': o_data, 'criterion': '_0', 'method': 'Asol=195',
                'v_n_x': 'A', 'v_n_y': 'Y_final',
                'color': color, 'ls': ls, 'lw': 1.0, 'ds': 'steps', 'alpha': 1.0,
                'ymin': 1e-5, 'ymax': 2e-1, 'xmin': 50, 'xmax': 210,
                'xlabel': r"A", 'ylabel': r'Relative final abundances',
                'label': label + dyn_lbl, 'yscale': 'log',
                'fancyticks': True, 'minorticks': True,
                'legend': {},
                'fontsize': 14,
                'labelsize': 14,
            }
            o_plot.set_plot_dics.append(sim_nucleo)
        if if_bern:
            sim_nucleo_b = {
                'task': 'nucleo', 'ptype': 'cartesian',
                'position': (1, 1),
                'data': o_data, 'criterion': '_0 _0_b_w', 'method': 'Asol=195',
                'v_n_x': 'A', 'v_n_y': 'Y_final',
                'color': color, 'ls': '-', 'lw': 1.0, 'ds': 'steps', 'alpha': 1.0,
                'ymin': 1e-5, 'ymax': 2e-1, 'xmin': 50, 'xmax': 210,
                'xlabel': r"A", 'ylabel': r'Relative final abundances',
                'label': label + wind_lbl, 'yscale': 'log',
                'fancyticks': True, 'minorticks': True,
                'legend': {},
                'fontsize': 14,
                'labelsize': 14,
            }
            o_plot.set_plot_dics.append(sim_nucleo_b)

        if sim == sims[-1]:
            # sol_yeilds['legend'] = {'loc': 'best', 'ncol': 2, 'fontsize': 14}
            sol_yeilds = {
                'task': 'nucleo', 'ptype': 'cartesian',
                'position': (1, 1),
                'data': o_data, 'criterion': '_0', 'method': 'sum',
                'v_n_x': 'Asun', 'v_n_y': 'Ysun',
                'color': 'gray', 'marker': 'o', 'ms': 4, 'alpha': alpha,
                'ymin': 8e-5, 'ymax': 8e-1, 'xmin': 50, 'xmax': 210,
                'xlabel': r"Mass number, A", 'ylabel': r'Relative final abundances',
                'label': 'solar', 'yscale': 'log',
                'fancyticks': True, 'minorticks': True,
                'legend':  {'loc': 'upper right', 'ncol': 3, 'fontsize': 14},
                'fontsize': 14,
                'labelsize': 14,
            }
            o_plot.set_plot_dics.append(sol_yeilds)

    o_plot.main()
    exit(1)
# plot_nucleo_mult_sims_sfho()
def plot_nucleo_mult_sims_sly4():
    o_plot = PLOT_MANY_TASKS()
    o_plot.gen_set["figdir"] = Paths.plots + "all_sly4_SR/"
    o_plot.gen_set["type"] = "cartesian"
    o_plot.gen_set["figsize"] = (9.0, 6.0)  # <->, |]
    o_plot.gen_set["figname"] = "nucleo_mult_sly4_dyn_sims.png"
    o_plot.gen_set["sharex"] = False
    o_plot.gen_set["sharey"] = False
    o_plot.gen_set["subplots_adjust_h"] = 0.3
    o_plot.gen_set["subplots_adjust_w"] = 0.0
    o_plot.set_plot_dics = []

    sims = [
        "SLy4_M13641364_M0_SR", # 21ms
        "SLy4_M13641364_M0_LK_SR",  # 55ms
        "SLy4_M14521283_M0_SR",# 30ms
    ]
    labels = ['136136', "136136 LK", "145128"]
    colors = ['red', "lightblue", "blue"]
    lss =    ['--', '--', "--"]
    dyns =   [True, True, True] # [False, False, False, False, False]
    berns =  [False, False, False]# [True, True, True, True, True]
    alphas = [1.0, 1.0, 1.0]
    lws =    [0.8, 0.8, 0.8]

    dyn_lbl =  " D."
    wind_lbl =  " DW."
    fontsize = 12
    labelsize = 12

    for sim, label, color, alpha, ls, if_bern, if_dyn in \
            zip(sims, labels, colors, alphas, lss, berns, dyns):
        o_data = NORMALIZE_NUCLEO(sim)

        if if_dyn:
            sim_nucleo = {
                'task': 'nucleo', 'ptype': 'cartesian',
                'position': (1, 1),
                'data': o_data, 'criterion': '_0', 'method': 'Asol=195',
                'v_n_x': 'A', 'v_n_y': 'Y_final',
                'color': color, 'ls': ls, 'lw': 1.0, 'ds': 'steps', 'alpha': 1.0,
                'ymin': 1e-5, 'ymax': 2e-1, 'xmin': 50, 'xmax': 210,
                'xlabel': r"A", 'ylabel': r'Relative final abundances',
                'label': label + dyn_lbl, 'yscale': 'log',
                'fancyticks': True, 'minorticks': True,
                'legend': {},
                'fontsize': 14,
                'labelsize': 14,
            }
            o_plot.set_plot_dics.append(sim_nucleo)
        if if_bern:
            sim_nucleo_b = {
                'task': 'nucleo', 'ptype': 'cartesian',
                'position': (1, 1),
                'data': o_data, 'criterion': '_0 _0_b_w', 'method': 'Asol=195',
                'v_n_x': 'A', 'v_n_y': 'Y_final',
                'color': color, 'ls': '-', 'lw': 1.0, 'ds': 'steps', 'alpha': 1.0,
                'ymin': 1e-5, 'ymax': 2e-1, 'xmin': 50, 'xmax': 210,
                'xlabel': r"A", 'ylabel': r'Relative final abundances',
                'label': label + wind_lbl, 'yscale': 'log',
                'fancyticks': True, 'minorticks': True,
                'legend': {},
                'fontsize': 14,
                'labelsize': 14,
            }
            o_plot.set_plot_dics.append(sim_nucleo_b)

        if sim == sims[-1]:
            # sol_yeilds['legend'] = {'loc': 'best', 'ncol': 2, 'fontsize': 14}
            sol_yeilds = {
                'task': 'nucleo', 'ptype': 'cartesian',
                'position': (1, 1),
                'data': o_data, 'criterion': '_0', 'method': 'sum',
                'v_n_x': 'Asun', 'v_n_y': 'Ysun',
                'color': 'gray', 'marker': 'o', 'ms': 4, 'alpha': alpha,
                'ymin': 8e-5, 'ymax': 8e-1, 'xmin': 50, 'xmax': 210,
                'xlabel': r"Mass number, A", 'ylabel': r'Relative final abundances',
                'label': 'solar', 'yscale': 'log',
                'fancyticks': True, 'minorticks': True,
                'legend':  {'loc': 'upper right', 'ncol': 3, 'fontsize': 14},
                'fontsize': 14,
                'labelsize': 14,
            }
            o_plot.set_plot_dics.append(sol_yeilds)

    o_plot.main()
    exit(1)
# plot_nucleo_mult_sims_sly4()

def ejecta_profs_mult_sims_ls220():

    o_plot = PLOT_MANY_TASKS()
    o_plot.gen_set["figdir"] = Paths.plots + 'all_ls220_SR/'
    o_plot.gen_set["type"] = "cartesian"
    o_plot.gen_set["figsize"] = (9.0, 3.7)  # <->, |]
    o_plot.gen_set["figname"] = "ejecta_ls220_sims.png"
    o_plot.gen_set["sharex"] = False
    o_plot.gen_set["sharey"] = False
    o_plot.gen_set["subplots_adjust_h"] = 0.3
    o_plot.gen_set["subplots_adjust_w"] = 0.0
    o_plot.set_plot_dics = []

    # label1 = "Dyn"
    # label2 = "Wind"

    fontsize = 14
    labelsize = 14

    sims = [
        "LS220_M13641364_M0_SR", # FULL 3D!
        "LS220_M14001330_M0_SR",
        "LS220_M14351298_M0_SR",
        "LS220_M14691268_M0_SR",

        "LS220_M13641364_M0_LK_SR",
        "LS220_M14691268_M0_LK_SR",
    ]
    labels = ['136136', '140133', "144130", "147127", "136136 LK", "147127 LK"]
    colors = ['red', 'orange', "pink", "magenta", "lightblue", "blue"]
    lss = ['--', '--', "--", "--", "--", "--"]
    dyns = [True, True, True, True, True, True]  # [False, False, False, False, False]
    berns = [True, True, True, True, True, True]  # [True, True, True, True, True]
    alphas = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0]
    lws = [0.8, 0.8, 0.8, 0.8, 0.8, 0.8]

    dyn_lbl = " D."
    wind_lbl = " DW."

    for sim, label, color, alpha, ls, if_bern, if_dyn in \
            zip(sims, labels, colors, alphas, lss, berns, dyns):

        o_data2 = COMPUTE_STORE_PAR(sim)

        if if_dyn:
            dic_ej_prof = {
                'task': 'ejprof', 'ptype': 'cartesian',
                'position': (1, 1),
                'data': o_data2, 'criterion': '_0',
                'v_n_x': 't_tot_flux', 'v_n_y': 'mass_tot_flux',
                'color': color, 'ls': ':', 'lw': 1., 'ds': 'default', 'alpha': 1.,
                'xmin': 0., 'xmax': 80, 'ymin': 0, 'ymax': 1.6,#2.80,
                'xlabel': r"$t-t_{merg}$ [ms]", 'ylabel': r'$M_{ej}$ $[10^{-2}M_{\odot}]$',
                'xunits': 'ms', '-t': o_data2.get_par("tmerger_gw"), 'yunits': '1e-2Msun',
                'label': label + dyn_lbl, 'yscale': None, 'title': {},
                'fancyticks': True, 'minorticks': True,
                # 'sharey': True,
                'fontsize': fontsize,
                'labelsize': labelsize,
                'legend': {'loc': 'lower right', 'ncol': 2, 'fontsize': 12}#,
                           # 'bbox_to_anchor': (0.65, 0.7)}  # (1.0, 1.2)
            }
            o_plot.set_plot_dics.append(dic_ej_prof)

        if if_bern:
            dic_ej_prof = {
                'task': 'ejprof', 'ptype': 'cartesian',
                'position': (1, 1),
                'data': o_data2, 'criterion': '_0_b_w',
                'v_n_x': 't_tot_flux', 'v_n_y': 'mass_tot_flux',
                'color': color, 'ls': '-', 'lw': 1., 'ds': 'default', 'alpha': 1.,
                'xmin': 0., 'xmax': 80, 'ymin': 0, 'ymax': 2.80,
                'xlabel': r"$t-t_{merg}$ [ms]", 'ylabel': r'$M_{ej}$ $[10^{-2}M_{\odot}]$',
                'xunits': 'ms', '-t': o_data2.get_par("tmerger_gw"), 'yunits': '1e-2Msun',
                'label': label + wind_lbl, 'yscale': None, 'title': {},
                'fancyticks': True, 'minorticks': True,
                # 'sharey': True,
                'fontsize': fontsize,
                'labelsize': labelsize,
                'legend': {'loc': 'lower right', 'ncol': 2, 'fontsize': 12}#,
                           # 'bbox_to_anchor': (0.65, 0.7)}  # (1.0, 1.2)
                           # 'bbox_to_anchor': (0.65, 0.7)}  # (1.0, 1.2)
            }
            o_plot.set_plot_dics.append(dic_ej_prof)

    o_plot.main()
    exit(1)
# ejecta_profs_mult_sims_ls220()
def ejecta_profs_mult_sims_dd2():

    o_plot = PLOT_MANY_TASKS()
    o_plot.gen_set["figdir"] = Paths.plots + 'all_dd2_SR/'
    o_plot.gen_set["type"] = "cartesian"
    o_plot.gen_set["figsize"] = (9.0, 3.7)  # <->, |]
    o_plot.gen_set["figname"] = "ejecta_dd2_sims.png"
    o_plot.gen_set["sharex"] = False
    o_plot.gen_set["sharey"] = False
    o_plot.gen_set["subplots_adjust_h"] = 0.3
    o_plot.gen_set["subplots_adjust_w"] = 0.0
    o_plot.set_plot_dics = []

    # label1 = "Dyn"
    # label2 = "Wind"

    fontsize = 14
    labelsize = 14

    sims = ["DD2_M13641364_M0_SR_R04",
            "DD2_M13641364_M0_LK_SR_R04",
            "DD2_M13641364_M0_SR",
            "DD2_M14971245_M0_SR",
            "DD2_M15091235_M0_LK_SR"]
    labels = ['136136', '136136 LK', "136136 R03", "150125", "151124 LK"]
    colors = ['red', 'orange', "pink", "lightblue", "blue"]
    lss = ['--', '--', "--", "--", "--"]
    dyns = [True, True, True, True, True]  # [False, False, False, False, False]
    berns = [True, True, True, True, True]  # [True, True, True, True, True]
    alphas = [1.0, 1.0, 1.0, 1.0, 1.0]
    lws = [0.8, 0.8, 0.8, 0.8, 0.8]

    dyn_lbl = " D."
    wind_lbl = " DW."

    for sim, label, color, alpha, ls, if_bern, if_dyn in \
            zip(sims, labels, colors, alphas, lss, berns, dyns):

        o_data2 = COMPUTE_STORE_PAR(sim)

        if if_dyn:
            dic_ej_prof = {
                'task': 'ejprof', 'ptype': 'cartesian',
                'position': (1, 1),
                'data': o_data2, 'criterion': '_0',
                'v_n_x': 't_tot_flux', 'v_n_y': 'mass_tot_flux',
                'color': color, 'ls': ':', 'lw': 1., 'ds': 'default', 'alpha': 1.,
                'xmin': 0., 'xmax': 90, 'ymin': 0, 'ymax': 1.6,#2.80,
                'xlabel': r"$t-t_{merg}$ [ms]", 'ylabel': r'$M_{ej}$ $[10^{-2}M_{\odot}]$',
                'xunits': 'ms', '-t': o_data2.get_par("tmerger_gw"), 'yunits': '1e-2Msun',
                'label': label + dyn_lbl, 'yscale': None, 'title': {},
                'fancyticks': True, 'minorticks': True,
                # 'sharey': True,
                'fontsize': fontsize,
                'labelsize': labelsize,
                'legend': {'loc': 'upper left', 'ncol': 2, 'fontsize': 12}#,
                           # 'bbox_to_anchor': (0.65, 0.7)}  # (1.0, 1.2)
            }
            o_plot.set_plot_dics.append(dic_ej_prof)

        if if_bern:
            dic_ej_prof = {
                'task': 'ejprof', 'ptype': 'cartesian',
                'position': (1, 1),
                'data': o_data2, 'criterion': '_0_b_w',
                'v_n_x': 't_tot_flux', 'v_n_y': 'mass_tot_flux',
                'color': color, 'ls': '-', 'lw': 1., 'ds': 'default', 'alpha': 1.,
                'xmin': 0., 'xmax': 110, 'ymin': 0, 'ymax': 1.6,##2.80,
                'xlabel': r"$t-t_{merg}$ [ms]", 'ylabel': r'$M_{ej}$ $[10^{-2}M_{\odot}]$',
                'xunits': 'ms', '-t': o_data2.get_par("tmerger_gw"), 'yunits': '1e-2Msun',
                'label': label + wind_lbl, 'yscale': None, 'title': {},
                'fancyticks': True, 'minorticks': True,
                # 'sharey': True,
                'fontsize': fontsize,
                'labelsize': labelsize,
                'legend': {'loc': 'upper left', 'ncol': 2, 'fontsize': 12}#,
                           # 'bbox_to_anchor': (0.65, 0.7)}  # (1.0, 1.2)
            }
            o_plot.set_plot_dics.append(dic_ej_prof)

    o_plot.main()
    exit(1)
# ejecta_profs_mult_sims_dd2()
def ejecta_profs_mult_sims_sfho():

    o_plot = PLOT_MANY_TASKS()
    o_plot.gen_set["figdir"] = Paths.plots + 'all_sfho_SR/'
    o_plot.gen_set["type"] = "cartesian"
    o_plot.gen_set["figsize"] = (9.0, 3.7)  # <->, |]
    o_plot.gen_set["figname"] = "ejecta_sfho_sims.png"
    o_plot.gen_set["sharex"] = False
    o_plot.gen_set["sharey"] = False
    o_plot.gen_set["subplots_adjust_h"] = 0.3
    o_plot.gen_set["subplots_adjust_w"] = 0.0
    o_plot.set_plot_dics = []

    # label1 = "Dyn"
    # label2 = "Wind"

    fontsize = 14
    labelsize = 14

    sims = [
        "SFHo_M13641364_M0_SR",  # 21ms
        "SFHo_M13641364_M0_LK_SR_2019pizza",  # 55ms
        "SFHo_M14521283_M0_LK_SR",  # 30ms
        "SFHo_M14521283_M0_LK_SR_2019pizza",  # 26ms
        "SFHo_M14521283_M0_SR",  # 33ms
    ]
    labels = ['136136', '136136 LK 19', "145128 LK", "145128 LK 19", "145128"]
    colors = ['red', 'orange', "pink", "lightblue", "blue"]
    lss = ['--', '--', "--", "--", "--"]
    dyns = [True, True, True, True, True]  # [False, False, False, False, False]
    berns = [True, True, True, True, True]  # [True, True, True, True, True]
    alphas = [1.0, 1.0, 1.0, 1.0, 1.0]
    lws = [0.8, 0.8, 0.8, 0.8, 0.8]

    dyn_lbl = " D."
    wind_lbl = " DW."

    for sim, label, color, alpha, ls, if_bern, if_dyn in \
            zip(sims, labels, colors, alphas, lss, berns, dyns):

        o_data2 = COMPUTE_STORE_PAR(sim)

        if if_dyn:
            dic_ej_prof = {
                'task': 'ejprof', 'ptype': 'cartesian',
                'position': (1, 1),
                'data': o_data2, 'criterion': '_0',
                'v_n_x': 't_tot_flux', 'v_n_y': 'mass_tot_flux',
                'color': color, 'ls': ':', 'lw': 1., 'ds': 'default', 'alpha': 1.,
                'xmin': 0., 'xmax': 90, 'ymin': 0, 'ymax': 1.6,#2.80,
                'xlabel': r"$t-t_{merg}$ [ms]", 'ylabel': r'$M_{ej}$ $[10^{-2}M_{\odot}]$',
                'xunits': 'ms', '-t': o_data2.get_par("tmerger_gw"), 'yunits': '1e-2Msun',
                'label': label + dyn_lbl, 'yscale': None, 'title': {},
                'fancyticks': True, 'minorticks': True,
                # 'sharey': True,
                'fontsize': fontsize,
                'labelsize': labelsize,
                'legend': {'loc': 'upper left', 'ncol': 2, 'fontsize': 12}#,
                           # 'bbox_to_anchor': (0.65, 0.7)}  # (1.0, 1.2)
            }
            o_plot.set_plot_dics.append(dic_ej_prof)

        if if_bern:
            dic_ej_prof = {
                'task': 'ejprof', 'ptype': 'cartesian',
                'position': (1, 1),
                'data': o_data2, 'criterion': '_0_b_w',
                'v_n_x': 't_tot_flux', 'v_n_y': 'mass_tot_flux',
                'color': color, 'ls': '-', 'lw': 1., 'ds': 'default', 'alpha': 1.,
                'xmin': 0., 'xmax': 30, 'ymin': 0, 'ymax': 1.6,##2.80,
                'xlabel': r"$t-t_{merg}$ [ms]", 'ylabel': r'$M_{ej}$ $[10^{-2}M_{\odot}]$',
                'xunits': 'ms', '-t': o_data2.get_par("tmerger_gw"), 'yunits': '1e-2Msun',
                'label': label + wind_lbl, 'yscale': None, 'title': {},
                'fancyticks': True, 'minorticks': True,
                # 'sharey': True,
                'fontsize': fontsize,
                'labelsize': labelsize,
                'legend': {'loc': 'upper left', 'ncol': 2, 'fontsize': 12}#,
                           # 'bbox_to_anchor': (0.65, 0.7)}  # (1.0, 1.2)
            }
            o_plot.set_plot_dics.append(dic_ej_prof)

    o_plot.main()
    exit(1)
# ejecta_profs_mult_sims_sfho()
def ejecta_profs_mult_sims_sly4():

    o_plot = PLOT_MANY_TASKS()
    o_plot.gen_set["figdir"] = Paths.plots + 'all_sly4_SR/'
    o_plot.gen_set["type"] = "cartesian"
    o_plot.gen_set["figsize"] = (9.0, 3.7)  # <->, |]
    o_plot.gen_set["figname"] = "ejecta_sly4_sims.png"
    o_plot.gen_set["sharex"] = False
    o_plot.gen_set["sharey"] = False
    o_plot.gen_set["subplots_adjust_h"] = 0.3
    o_plot.gen_set["subplots_adjust_w"] = 0.0
    o_plot.set_plot_dics = []

    # label1 = "Dyn"
    # label2 = "Wind"

    fontsize = 14
    labelsize = 14

    sims = [
        "SLy4_M13641364_M0_SR",  # 21ms
        "SLy4_M13641364_M0_LK_SR",  # 55ms
        "SLy4_M14521283_M0_SR",  # 30ms
    ]
    labels = ['136136', "136136 LK", "145128"]
    colors = ['red', "lightblue", "blue"]
    lss = ['--', '--', "--"]
    dyns = [True, True, True]  # [False, False, False, False, False]
    berns = [True, True, True]  # [True, True, True, True, True]
    alphas = [1.0, 1.0, 1.0]
    lws = [0.8, 0.8, 0.8]

    dyn_lbl = " D."
    wind_lbl = " DW."

    for sim, label, color, alpha, ls, if_bern, if_dyn in \
            zip(sims, labels, colors, alphas, lss, berns, dyns):

        o_data2 = COMPUTE_STORE_PAR(sim)

        if if_dyn:
            dic_ej_prof = {
                'task': 'ejprof', 'ptype': 'cartesian',
                'position': (1, 1),
                'data': o_data2, 'criterion': '_0',
                'v_n_x': 't_tot_flux', 'v_n_y': 'mass_tot_flux',
                'color': color, 'ls': ':', 'lw': 1., 'ds': 'default', 'alpha': 1.,
                'xmin': 0., 'xmax': 90, 'ymin': 0, 'ymax': 1.6,#2.80,
                'xlabel': r"$t-t_{merg}$ [ms]", 'ylabel': r'$M_{ej}$ $[10^{-2}M_{\odot}]$',
                'xunits': 'ms', '-t': o_data2.get_par("tmerger_gw"), 'yunits': '1e-2Msun',
                'label': label + dyn_lbl, 'yscale': None, 'title': {},
                'fancyticks': True, 'minorticks': True,
                # 'sharey': True,
                'fontsize': fontsize,
                'labelsize': labelsize,
                'legend': {'loc': 'upper left', 'ncol': 2, 'fontsize': 12}#,
                           # 'bbox_to_anchor': (0.65, 0.7)}  # (1.0, 1.2)
            }
            o_plot.set_plot_dics.append(dic_ej_prof)

        if if_bern:
            dic_ej_prof = {
                'task': 'ejprof', 'ptype': 'cartesian',
                'position': (1, 1),
                'data': o_data2, 'criterion': '_0_b_w',
                'v_n_x': 't_tot_flux', 'v_n_y': 'mass_tot_flux',
                'color': color, 'ls': '-', 'lw': 1., 'ds': 'default', 'alpha': 1.,
                'xmin': 0., 'xmax': 30, 'ymin': 0, 'ymax': 1.6,##2.80,
                'xlabel': r"$t-t_{merg}$ [ms]", 'ylabel': r'$M_{ej}$ $[10^{-2}M_{\odot}]$',
                'xunits': 'ms', '-t': o_data2.get_par("tmerger_gw"), 'yunits': '1e-2Msun',
                'label': label + wind_lbl, 'yscale': None, 'title': {},
                'fancyticks': True, 'minorticks': True,
                # 'sharey': True,
                'fontsize': fontsize,
                'labelsize': labelsize,
                'legend': {'loc': 'upper left', 'ncol': 2, 'fontsize': 12}#,
                           # 'bbox_to_anchor': (0.65, 0.7)}  # (1.0, 1.2)
            }
            o_plot.set_plot_dics.append(dic_ej_prof)

    o_plot.main()
    exit(1)
# ejecta_profs_mult_sims_sly4()

def plot_hists_mult_sims_ls220():
    o_plot = PLOT_MANY_TASKS()
    o_plot.gen_set["figdir"] = Paths.plots + "all_ls220_SR/"
    o_plot.gen_set["type"] = "cartesian"
    o_plot.gen_set["figsize"] = (9.0, 3.7)  # <->, |]
    o_plot.gen_set["figname"] = "hists_mult_ls220_wind_sims.png"
    o_plot.gen_set["sharex"] = False
    o_plot.gen_set["sharey"] = False
    o_plot.gen_set["subplots_adjust_h"] = 0.3
    o_plot.gen_set["subplots_adjust_w"] = 0.0

    sims = [
        "LS220_M13641364_M0_SR",  # FULL 3D!
        "LS220_M14001330_M0_SR",
        "LS220_M14351298_M0_SR",
        "LS220_M14691268_M0_SR",

        "LS220_M13641364_M0_LK_SR",
        "LS220_M14691268_M0_LK_SR",
    ]
    labels = ['136136', '140133', "144130", "147127", "136136 LK", "147127 LK"]
    colors = ['red', 'orange', "pink", "magenta", "lightblue", "blue"]
    lss = ['--', '--', "--", "--", "--", "--"]
    dyns = [False, False, False, False, False, False]  # [False, False, False, False, False, False]
    berns = [True, True, True, True, True, True]  # [True, True, True, True, True, True]
    alphas = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0]
    lws = [0.8, 0.8, 0.8, 0.8, 0.8, 0.8]

    dyn_lbl =  " D."
    wind_lbl =  " DW."
    fontsize = 14
    labelsize = 14

    for sim, label, color, alpha, ls, if_bern, if_dyn in \
            zip(sims, labels, colors, alphas, lss, berns, dyns):

        o_data = COMPUTE_STORE_PAR(sim)
        lw = 1.0
        alpha = 1.0
        if if_dyn:
            dic_hist_theta = {
                'task': 'hist1d', 'ptype': 'cartesian',
                'position': (1, 1),
                'data': o_data, 'norm': True, 'criterion': '_0',
                'v_n_x': 'hist_theta', 'v_n_y': 'hist_theta_m',
                'color': color, 'ls': ':', 'lw':lw, 'ds': 'steps', 'alpha': alpha,
                'ymin': 1e-4, 'ymax': 1e0,
                'xlabel': r"Angle from orbital plane", 'ylabel': r'$M_{ej}$ [normed]',
                'label': label + dyn_lbl,
                'yscale': 'log',  # "LS220 M14691268 "+label1
                'fancyticks': True, 'minorticks': True,
                'fontsize': fontsize,
                'labelsize': labelsize,
                'legend': {}  # {'loc': 'best', 'ncol': 2, 'fontsize': 18}
            }

            dic_hist_vel_inf = {
                'task': 'hist1d', 'ptype': 'cartesian',
                'position': (1, 2),
                'data': o_data, 'norm': True, 'criterion': '_0',
                'v_n_x': 'hist_vel_inf', 'v_n_y': 'hist_vel_inf_m',
                'color': color, 'ls': ':', 'lw': lw, 'ds': 'steps', 'alpha': alpha,
                'xmin': 0., 'xmax': 1., 'ymin': 1e-4, 'ymax': 1e0,
                'xlabel': r"$\upsilon_{\infty}$ [c]", 'ylabel': r'$M_{ej}$ [normed]',
                'label': label + dyn_lbl,  # "LS220 M14691268 "+label2,
                'yscale': 'log',
                'fancyticks': True, 'minorticks': True,
                'sharey': True,
                'fontsize': fontsize,
                'labelsize': labelsize,
                'legend': {}  # {'loc': 'best', 'ncol': 2, 'fontsize': 18}
            }

            dic_hist_ye = {
                'task': 'hist1d', 'ptype': 'cartesian',
                'position': (1, 3),
                'data': o_data, 'norm': True, 'criterion': '_0',
                'v_n_x': 'hist_ye', 'v_n_y': 'hist_ye_m',
                'color': color, 'ls': ':', 'lw': lw, 'ds': 'steps', 'alpha': alpha,
                'xmin': 0., 'xmax': 0.5, 'ymin': 1e-4, 'ymax': 1e0,
                'xlabel': r"$Y_e$", 'ylabel': r'$M_{ej}$ [normed]',
                'label': label + dyn_lbl,
                'yscale': 'log',
                'fancyticks': True, 'minorticks': True,
                'sharey': True,
                'fontsize': fontsize,
                'labelsize': labelsize,
                'legend': {}  # {'loc': 'best', 'ncol': 2, 'fontsize': 18}
            }

            if sim == sims[-1]:
                dic_hist_theta['legend'] =  {'loc': 'lower right', 'ncol': 3, 'fontsize': 12,
                                             'bbox_to_anchor': (2.20, 1.02)}

            o_plot.set_plot_dics.append(dic_hist_theta)
            o_plot.set_plot_dics.append(dic_hist_vel_inf)
            o_plot.set_plot_dics.append(dic_hist_ye)

        if if_bern:
            dic_hist_theta_b = {
                'task': 'hist1d', 'ptype': 'cartesian',
                'position': (1, 1),
                'data': o_data, 'norm': True, 'criterion': '_0_b_w',
                'v_n_x': 'hist_theta', 'v_n_y': 'hist_theta_m',
                'color': color, 'ls': '-', 'lw': lw, 'ds': 'steps', 'alpha': alpha,
                'ymin': 1e-4, 'ymax': 1e0,
                'xlabel': r"Angle from orbital plane", 'ylabel': r'$M_{ej}$ [normed]',
                'label': label + wind_lbl,  # "LS220 M14691268 "+label2,
                'yscale': 'log',
                'fancyticks': True, 'minorticks': True,
                'fontsize': fontsize,
                'labelsize': labelsize,
                'legend': {}  # {'loc': 'best', 'ncol': 2, 'fontsize': 18}
            }

            dic_hist_vel_inf_b = {
                'task': 'hist1d', 'ptype': 'cartesian',
                'position': (1, 2),
                'data': o_data, 'norm': True, 'criterion': '_0_b_w',
                'v_n_x': 'hist_vel_inf', 'v_n_y': 'hist_vel_inf_m',
                'color': color, 'ls': '-', 'lw': lw, 'ds': 'steps', 'alpha': alpha,
                'xmin': 0., 'xmax': 1., 'ymin': 1e-4, 'ymax': 1e0,
                'xlabel': r"$\upsilon_{\infty}$ [c]", 'ylabel': r'$M_{ej}$ [normed]',
                'label': label + wind_lbl,
                'yscale': 'log',
                'fancyticks': True, 'minorticks': True,
                'sharey': True,
                'fontsize': fontsize,
                'labelsize': labelsize,
                'legend': {}  # {'loc': 'best', 'ncol': 2, 'fontsize': 18}
            }

            dic_hist_ye_b = {
                'task': 'hist1d', 'ptype': 'cartesian',
                'position': (1, 3),
                'data': o_data, 'norm': True, 'criterion': '_0_b_w',
                'v_n_x': 'hist_ye', 'v_n_y': 'hist_ye_m',
                'color': color, 'ls': '-', 'lw': lw, 'ds': 'steps', 'alpha': alpha,
                'xmin': 0., 'xmax': 0.5, 'ymin': 1e-4, 'ymax': 1e0,
                'xlabel': r"$Y_e$", 'ylabel': r'$M_{ej}$ [normed]',
                'label': label + wind_lbl,
                'yscale': 'log',
                'fancyticks': True, 'minorticks': True,
                'sharey': True,
                'fontsize': fontsize,
                'labelsize': labelsize,
                'legend': {}  # {'loc': 'best', 'ncol': 2, 'fontsize': 18}
            }

            if sim == sims[-1]:
                dic_hist_theta_b['legend'] = {'loc': 'lower right', 'ncol': 3, 'fontsize': 12,
                                           'bbox_to_anchor': (2.40, 1.02)}

            o_plot.set_plot_dics.append(dic_hist_theta_b)
            o_plot.set_plot_dics.append(dic_hist_vel_inf_b)
            o_plot.set_plot_dics.append(dic_hist_ye_b)

    o_plot.main()
    exit(1)
# plot_hists_mult_sims_ls220()
def plot_hists_mult_sims_dd2():
    o_plot = PLOT_MANY_TASKS()
    o_plot.gen_set["figdir"] = Paths.plots + "all_dd2_SR/"
    o_plot.gen_set["type"] = "cartesian"
    o_plot.gen_set["figsize"] = (9.0, 3.7)  # <->, |]
    o_plot.gen_set["figname"] = "hists_mult_dd2_dyn_sims.png"
    o_plot.gen_set["sharex"] = False
    o_plot.gen_set["sharey"] = False
    o_plot.gen_set["subplots_adjust_h"] = 0.3
    o_plot.gen_set["subplots_adjust_w"] = 0.0

    sims = ["DD2_M13641364_M0_SR_R04",
            "DD2_M13641364_M0_LK_SR_R04",
            "DD2_M13641364_M0_SR",
            "DD2_M14971245_M0_SR",
            "DD2_M15091235_M0_LK_SR"]
    labels = ['136136', '136136 LK', "136136 R03", "150125", "151124 LK"]
    colors = ['red', 'orange', "pink", "lightblue", "blue"]
    lss = ['--', '--', "--", "--", "--"]
    dyns =  [True, True, True, True, True]  # [False, False, False, False, False]
    berns = [False, False, False, False, False]  # [True, True, True, True, True]
    alphas = [1.0, 1.0, 1.0, 1.0, 1.0]
    lws = [0.8, 0.8, 0.8, 0.8, 0.8]

    dyn_lbl =  " D."
    wind_lbl =  " DW."
    fontsize = 14
    labelsize = 14

    for sim, label, color, alpha, ls, if_bern, if_dyn in \
            zip(sims, labels, colors, alphas, lss, berns, dyns):

        o_data = COMPUTE_STORE_PAR(sim)
        lw = 1.0
        alpha = 1.0
        if if_dyn:
            dic_hist_theta = {
                'task': 'hist1d', 'ptype': 'cartesian',
                'position': (1, 1),
                'data': o_data, 'norm': True, 'criterion': '_0',
                'v_n_x': 'hist_theta', 'v_n_y': 'hist_theta_m',
                'color': color, 'ls': ':', 'lw':lw, 'ds': 'steps', 'alpha': alpha,
                'ymin': 1e-4, 'ymax': 1e0,
                'xlabel': r"Angle from orbital plane", 'ylabel': r'$M_{ej}$ [normed]',
                'label': label + dyn_lbl,
                'yscale': 'log',  # "LS220 M14691268 "+label1
                'fancyticks': True, 'minorticks': True,
                'fontsize': fontsize,
                'labelsize': labelsize,
                'legend': {}  # {'loc': 'best', 'ncol': 2, 'fontsize': 18}
            }

            dic_hist_vel_inf = {
                'task': 'hist1d', 'ptype': 'cartesian',
                'position': (1, 2),
                'data': o_data, 'norm': True, 'criterion': '_0',
                'v_n_x': 'hist_vel_inf', 'v_n_y': 'hist_vel_inf_m',
                'color': color, 'ls': ':', 'lw': lw, 'ds': 'steps', 'alpha': alpha,
                'xmin': 0., 'xmax': 1., 'ymin': 1e-4, 'ymax': 1e0,
                'xlabel': r"$\upsilon_{\infty}$ [c]", 'ylabel': r'$M_{ej}$ [normed]',
                'label': label + dyn_lbl,  # "LS220 M14691268 "+label2,
                'yscale': 'log',
                'fancyticks': True, 'minorticks': True,
                'sharey': True,
                'fontsize': fontsize,
                'labelsize': labelsize,
                'legend': {}  # {'loc': 'best', 'ncol': 2, 'fontsize': 18}
            }

            dic_hist_ye = {
                'task': 'hist1d', 'ptype': 'cartesian',
                'position': (1, 3),
                'data': o_data, 'norm': True, 'criterion': '_0',
                'v_n_x': 'hist_ye', 'v_n_y': 'hist_ye_m',
                'color': color, 'ls': ':', 'lw': lw, 'ds': 'steps', 'alpha': alpha,
                'xmin': 0., 'xmax': 0.5, 'ymin': 1e-4, 'ymax': 1e0,
                'xlabel': r"$Y_e$", 'ylabel': r'$M_{ej}$ [normed]',
                'label': label + dyn_lbl,
                'yscale': 'log',
                'fancyticks': True, 'minorticks': True,
                'sharey': True,
                'fontsize': fontsize,
                'labelsize': labelsize,
                'legend': {}  # {'loc': 'best', 'ncol': 2, 'fontsize': 18}
            }

            if sim == sims[-1]:
                dic_hist_theta['legend'] =  {'loc': 'lower right', 'ncol': 3, 'fontsize': 12,
                                             'bbox_to_anchor': (2.20, 1.02)}

            o_plot.set_plot_dics.append(dic_hist_theta)
            o_plot.set_plot_dics.append(dic_hist_vel_inf)
            o_plot.set_plot_dics.append(dic_hist_ye)

        if if_bern:
            dic_hist_theta_b = {
                'task': 'hist1d', 'ptype': 'cartesian',
                'position': (1, 1),
                'data': o_data, 'norm': True, 'criterion': '_0_b_w',
                'v_n_x': 'hist_theta', 'v_n_y': 'hist_theta_m',
                'color': color, 'ls': '-', 'lw': lw, 'ds': 'steps', 'alpha': alpha,
                'ymin': 1e-4, 'ymax': 1e0,
                'xlabel': r"Angle from orbital plane", 'ylabel': r'$M_{ej}$ [normed]',
                'label': label + wind_lbl,  # "LS220 M14691268 "+label2,
                'yscale': 'log',
                'fancyticks': True, 'minorticks': True,
                'fontsize': fontsize,
                'labelsize': labelsize,
                'legend': {}  # {'loc': 'best', 'ncol': 2, 'fontsize': 18}
            }

            dic_hist_vel_inf_b = {
                'task': 'hist1d', 'ptype': 'cartesian',
                'position': (1, 2),
                'data': o_data, 'norm': True, 'criterion': '_0_b_w',
                'v_n_x': 'hist_vel_inf', 'v_n_y': 'hist_vel_inf_m',
                'color': color, 'ls': '-', 'lw': lw, 'ds': 'steps', 'alpha': alpha,
                'xmin': 0., 'xmax': 1., 'ymin': 1e-4, 'ymax': 1e0,
                'xlabel': r"$\upsilon_{\infty}$ [c]", 'ylabel': r'$M_{ej}$ [normed]',
                'label': label + wind_lbl,
                'yscale': 'log',
                'fancyticks': True, 'minorticks': True,
                'sharey': True,
                'fontsize': fontsize,
                'labelsize': labelsize,
                'legend': {}  # {'loc': 'best', 'ncol': 2, 'fontsize': 18}
            }

            dic_hist_ye_b = {
                'task': 'hist1d', 'ptype': 'cartesian',
                'position': (1, 3),
                'data': o_data, 'norm': True, 'criterion': '_0_b_w',
                'v_n_x': 'hist_ye', 'v_n_y': 'hist_ye_m',
                'color': color, 'ls': '-', 'lw': lw, 'ds': 'steps', 'alpha': alpha,
                'xmin': 0., 'xmax': 0.5, 'ymin': 1e-4, 'ymax': 1e0,
                'xlabel': r"$Y_e$", 'ylabel': r'$M_{ej}$ [normed]',
                'label': label + wind_lbl,
                'yscale': 'log',
                'fancyticks': True, 'minorticks': True,
                'sharey': True,
                'fontsize': fontsize,
                'labelsize': labelsize,
                'legend': {}  # {'loc': 'best', 'ncol': 2, 'fontsize': 18}
            }

            if sim == sims[-1]:
                dic_hist_theta_b['legend'] = {'loc': 'lower right', 'ncol': 3, 'fontsize': 12,
                                           'bbox_to_anchor': (2.40, 1.02)}

            o_plot.set_plot_dics.append(dic_hist_theta_b)
            o_plot.set_plot_dics.append(dic_hist_vel_inf_b)
            o_plot.set_plot_dics.append(dic_hist_ye_b)

    o_plot.main()
    exit(1)
# plot_hists_mult_sims_dd2()
def plot_hists_mult_sims_sfho():
    o_plot = PLOT_MANY_TASKS()
    o_plot.gen_set["figdir"] = Paths.plots + "all_sfho_SR/"
    o_plot.gen_set["type"] = "cartesian"
    o_plot.gen_set["figsize"] = (9.0, 3.7)  # <->, |]
    o_plot.gen_set["figname"] = "hists_mult_sfho_wind_sims.png"
    o_plot.gen_set["sharex"] = False
    o_plot.gen_set["sharey"] = False
    o_plot.gen_set["subplots_adjust_h"] = 0.3
    o_plot.gen_set["subplots_adjust_w"] = 0.0

    sims = [
        "SFHo_M13641364_M0_SR",  # 21ms
        "SFHo_M13641364_M0_LK_SR_2019pizza",  # 55ms
        "SFHo_M14521283_M0_LK_SR",  # 30ms
        "SFHo_M14521283_M0_LK_SR_2019pizza",  # 26ms
        "SFHo_M14521283_M0_SR",  # 33ms
    ]
    labels = ['136136', '136136 LK 19', "145128 LK", "145128 LK 19", "145128"]
    colors = ['red', 'orange', "pink", "lightblue", "blue"]
    lss = ['--', '--', "--", "--", "--"]
    dyns =  [False, False, False, False, False]  # [False, False, False, False, False]
    berns = [True, True, True, True, True]  # [True, True, True, True, True]
    alphas = [1.0, 1.0, 1.0, 1.0, 1.0]
    lws = [0.8, 0.8, 0.8, 0.8, 0.8]

    dyn_lbl =  " D."
    wind_lbl =  " DW."
    fontsize = 14
    labelsize = 14

    for sim, label, color, alpha, ls, if_bern, if_dyn in \
            zip(sims, labels, colors, alphas, lss, berns, dyns):

        o_data = COMPUTE_STORE_PAR(sim)
        lw = 1.0
        alpha = 1.0
        if if_dyn:
            dic_hist_theta = {
                'task': 'hist1d', 'ptype': 'cartesian',
                'position': (1, 1),
                'data': o_data, 'norm': True, 'criterion': '_0',
                'v_n_x': 'hist_theta', 'v_n_y': 'hist_theta_m',
                'color': color, 'ls': ':', 'lw':lw, 'ds': 'steps', 'alpha': alpha,
                'ymin': 1e-4, 'ymax': 1e0,
                'xlabel': r"Angle from orbital plane", 'ylabel': r'$M_{ej}$ [normed]',
                'label': label + dyn_lbl,
                'yscale': 'log',  # "LS220 M14691268 "+label1
                'fancyticks': True, 'minorticks': True,
                'fontsize': fontsize,
                'labelsize': labelsize,
                'legend': {}  # {'loc': 'best', 'ncol': 2, 'fontsize': 18}
            }

            dic_hist_vel_inf = {
                'task': 'hist1d', 'ptype': 'cartesian',
                'position': (1, 2),
                'data': o_data, 'norm': True, 'criterion': '_0',
                'v_n_x': 'hist_vel_inf', 'v_n_y': 'hist_vel_inf_m',
                'color': color, 'ls': ':', 'lw': lw, 'ds': 'steps', 'alpha': alpha,
                'xmin': 0., 'xmax': 1., 'ymin': 1e-4, 'ymax': 1e0,
                'xlabel': r"$\upsilon_{\infty}$ [c]", 'ylabel': r'$M_{ej}$ [normed]',
                'label': label + dyn_lbl,  # "LS220 M14691268 "+label2,
                'yscale': 'log',
                'fancyticks': True, 'minorticks': True,
                'sharey': True,
                'fontsize': fontsize,
                'labelsize': labelsize,
                'legend': {}  # {'loc': 'best', 'ncol': 2, 'fontsize': 18}
            }

            dic_hist_ye = {
                'task': 'hist1d', 'ptype': 'cartesian',
                'position': (1, 3),
                'data': o_data, 'norm': True, 'criterion': '_0',
                'v_n_x': 'hist_ye', 'v_n_y': 'hist_ye_m',
                'color': color, 'ls': ':', 'lw': lw, 'ds': 'steps', 'alpha': alpha,
                'xmin': 0., 'xmax': 0.5, 'ymin': 1e-4, 'ymax': 1e0,
                'xlabel': r"$Y_e$", 'ylabel': r'$M_{ej}$ [normed]',
                'label': label + dyn_lbl,
                'yscale': 'log',
                'fancyticks': True, 'minorticks': True,
                'sharey': True,
                'fontsize': fontsize,
                'labelsize': labelsize,
                'legend': {}  # {'loc': 'best', 'ncol': 2, 'fontsize': 18}
            }

            if sim == sims[-1]:
                dic_hist_theta['legend'] =  {'loc': 'lower right', 'ncol': 3, 'fontsize': 12,
                                             'bbox_to_anchor': (2.20, 1.02)}

            o_plot.set_plot_dics.append(dic_hist_theta)
            o_plot.set_plot_dics.append(dic_hist_vel_inf)
            o_plot.set_plot_dics.append(dic_hist_ye)

        if if_bern:
            dic_hist_theta_b = {
                'task': 'hist1d', 'ptype': 'cartesian',
                'position': (1, 1),
                'data': o_data, 'norm': True, 'criterion': '_0_b_w',
                'v_n_x': 'hist_theta', 'v_n_y': 'hist_theta_m',
                'color': color, 'ls': '-', 'lw': lw, 'ds': 'steps', 'alpha': alpha,
                'ymin': 1e-4, 'ymax': 1e0,
                'xlabel': r"Angle from orbital plane", 'ylabel': r'$M_{ej}$ [normed]',
                'label': label + wind_lbl,  # "LS220 M14691268 "+label2,
                'yscale': 'log',
                'fancyticks': True, 'minorticks': True,
                'fontsize': fontsize,
                'labelsize': labelsize,
                'legend': {}  # {'loc': 'best', 'ncol': 2, 'fontsize': 18}
            }

            dic_hist_vel_inf_b = {
                'task': 'hist1d', 'ptype': 'cartesian',
                'position': (1, 2),
                'data': o_data, 'norm': True, 'criterion': '_0_b_w',
                'v_n_x': 'hist_vel_inf', 'v_n_y': 'hist_vel_inf_m',
                'color': color, 'ls': '-', 'lw': lw, 'ds': 'steps', 'alpha': alpha,
                'xmin': 0., 'xmax': 1., 'ymin': 1e-4, 'ymax': 1e0,
                'xlabel': r"$\upsilon_{\infty}$ [c]", 'ylabel': r'$M_{ej}$ [normed]',
                'label': label + wind_lbl,
                'yscale': 'log',
                'fancyticks': True, 'minorticks': True,
                'sharey': True,
                'fontsize': fontsize,
                'labelsize': labelsize,
                'legend': {}  # {'loc': 'best', 'ncol': 2, 'fontsize': 18}
            }

            dic_hist_ye_b = {
                'task': 'hist1d', 'ptype': 'cartesian',
                'position': (1, 3),
                'data': o_data, 'norm': True, 'criterion': '_0_b_w',
                'v_n_x': 'hist_ye', 'v_n_y': 'hist_ye_m',
                'color': color, 'ls': '-', 'lw': lw, 'ds': 'steps', 'alpha': alpha,
                'xmin': 0., 'xmax': 0.5, 'ymin': 1e-4, 'ymax': 1e0,
                'xlabel': r"$Y_e$", 'ylabel': r'$M_{ej}$ [normed]',
                'label': label + wind_lbl,
                'yscale': 'log',
                'fancyticks': True, 'minorticks': True,
                'sharey': True,
                'fontsize': fontsize,
                'labelsize': labelsize,
                'legend': {}  # {'loc': 'best', 'ncol': 2, 'fontsize': 18}
            }

            if sim == sims[-1]:
                dic_hist_theta_b['legend'] = {'loc': 'lower right', 'ncol': 3, 'fontsize': 12,
                                           'bbox_to_anchor': (2.40, 1.02)}

            o_plot.set_plot_dics.append(dic_hist_theta_b)
            o_plot.set_plot_dics.append(dic_hist_vel_inf_b)
            o_plot.set_plot_dics.append(dic_hist_ye_b)

    o_plot.main()
    exit(1)
# plot_hists_mult_sims_sfho()
def plot_hists_mult_sims_sly4():
    o_plot = PLOT_MANY_TASKS()
    o_plot.gen_set["figdir"] = Paths.plots + "all_sly4_SR/"
    o_plot.gen_set["type"] = "cartesian"
    o_plot.gen_set["figsize"] = (9.0, 3.7)  # <->, |]
    o_plot.gen_set["figname"] = "hists_mult_sly4_dyn_sims.png"
    o_plot.gen_set["sharex"] = False
    o_plot.gen_set["sharey"] = False
    o_plot.gen_set["subplots_adjust_h"] = 0.3
    o_plot.gen_set["subplots_adjust_w"] = 0.0

    sims = [
        "SLy4_M13641364_M0_SR",  # 21ms
        "SLy4_M13641364_M0_LK_SR",  # 55ms
        "SLy4_M14521283_M0_SR",  # 30ms
    ]
    labels = ['136136', "136136 LK", "145128"]
    colors = ['red', "lightblue", "blue"]
    lss = ['--', '--', "--"]
    dyns = [True, True, True]  # [False, False, False, False, False]
    berns = [False, False, False]  # [True, True, True, True, True]
    alphas = [1.0, 1.0, 1.0]
    lws = [0.8, 0.8, 0.8]

    dyn_lbl =  " D."
    wind_lbl =  " DW."
    fontsize = 14
    labelsize = 14

    for sim, label, color, alpha, ls, if_bern, if_dyn in \
            zip(sims, labels, colors, alphas, lss, berns, dyns):

        o_data = COMPUTE_STORE_PAR(sim)
        lw = 1.0
        alpha = 1.0
        if if_dyn:
            dic_hist_theta = {
                'task': 'hist1d', 'ptype': 'cartesian',
                'position': (1, 1),
                'data': o_data, 'norm': True, 'criterion': '_0',
                'v_n_x': 'hist_theta', 'v_n_y': 'hist_theta_m',
                'color': color, 'ls': ':', 'lw':lw, 'ds': 'steps', 'alpha': alpha,
                'ymin': 1e-4, 'ymax': 1e0,
                'xlabel': r"Angle from orbital plane", 'ylabel': r'$M_{ej}$ [normed]',
                'label': label + dyn_lbl,
                'yscale': 'log',  # "LS220 M14691268 "+label1
                'fancyticks': True, 'minorticks': True,
                'fontsize': fontsize,
                'labelsize': labelsize,
                'legend': {}  # {'loc': 'best', 'ncol': 2, 'fontsize': 18}
            }

            dic_hist_vel_inf = {
                'task': 'hist1d', 'ptype': 'cartesian',
                'position': (1, 2),
                'data': o_data, 'norm': True, 'criterion': '_0',
                'v_n_x': 'hist_vel_inf', 'v_n_y': 'hist_vel_inf_m',
                'color': color, 'ls': ':', 'lw': lw, 'ds': 'steps', 'alpha': alpha,
                'xmin': 0., 'xmax': 1., 'ymin': 1e-4, 'ymax': 1e0,
                'xlabel': r"$\upsilon_{\infty}$ [c]", 'ylabel': r'$M_{ej}$ [normed]',
                'label': label + dyn_lbl,  # "LS220 M14691268 "+label2,
                'yscale': 'log',
                'fancyticks': True, 'minorticks': True,
                'sharey': True,
                'fontsize': fontsize,
                'labelsize': labelsize,
                'legend': {}  # {'loc': 'best', 'ncol': 2, 'fontsize': 18}
            }

            dic_hist_ye = {
                'task': 'hist1d', 'ptype': 'cartesian',
                'position': (1, 3),
                'data': o_data, 'norm': True, 'criterion': '_0',
                'v_n_x': 'hist_ye', 'v_n_y': 'hist_ye_m',
                'color': color, 'ls': ':', 'lw': lw, 'ds': 'steps', 'alpha': alpha,
                'xmin': 0., 'xmax': 0.5, 'ymin': 1e-4, 'ymax': 1e0,
                'xlabel': r"$Y_e$", 'ylabel': r'$M_{ej}$ [normed]',
                'label': label + dyn_lbl,
                'yscale': 'log',
                'fancyticks': True, 'minorticks': True,
                'sharey': True,
                'fontsize': fontsize,
                'labelsize': labelsize,
                'legend': {}  # {'loc': 'best', 'ncol': 2, 'fontsize': 18}
            }

            if sim == sims[-1]:
                dic_hist_theta['legend'] =  {'loc': 'lower right', 'ncol': 3, 'fontsize': 12,
                                             'bbox_to_anchor': (2.20, 1.02)}

            o_plot.set_plot_dics.append(dic_hist_theta)
            o_plot.set_plot_dics.append(dic_hist_vel_inf)
            o_plot.set_plot_dics.append(dic_hist_ye)

        if if_bern:
            dic_hist_theta_b = {
                'task': 'hist1d', 'ptype': 'cartesian',
                'position': (1, 1),
                'data': o_data, 'norm': True, 'criterion': '_0_b_w',
                'v_n_x': 'hist_theta', 'v_n_y': 'hist_theta_m',
                'color': color, 'ls': '-', 'lw': lw, 'ds': 'steps', 'alpha': alpha,
                'ymin': 1e-4, 'ymax': 1e0,
                'xlabel': r"Angle from orbital plane", 'ylabel': r'$M_{ej}$ [normed]',
                'label': label + wind_lbl,  # "LS220 M14691268 "+label2,
                'yscale': 'log',
                'fancyticks': True, 'minorticks': True,
                'fontsize': fontsize,
                'labelsize': labelsize,
                'legend': {}  # {'loc': 'best', 'ncol': 2, 'fontsize': 18}
            }

            dic_hist_vel_inf_b = {
                'task': 'hist1d', 'ptype': 'cartesian',
                'position': (1, 2),
                'data': o_data, 'norm': True, 'criterion': '_0_b_w',
                'v_n_x': 'hist_vel_inf', 'v_n_y': 'hist_vel_inf_m',
                'color': color, 'ls': '-', 'lw': lw, 'ds': 'steps', 'alpha': alpha,
                'xmin': 0., 'xmax': 1., 'ymin': 1e-4, 'ymax': 1e0,
                'xlabel': r"$\upsilon_{\infty}$ [c]", 'ylabel': r'$M_{ej}$ [normed]',
                'label': label + wind_lbl,
                'yscale': 'log',
                'fancyticks': True, 'minorticks': True,
                'sharey': True,
                'fontsize': fontsize,
                'labelsize': labelsize,
                'legend': {}  # {'loc': 'best', 'ncol': 2, 'fontsize': 18}
            }

            dic_hist_ye_b = {
                'task': 'hist1d', 'ptype': 'cartesian',
                'position': (1, 3),
                'data': o_data, 'norm': True, 'criterion': '_0_b_w',
                'v_n_x': 'hist_ye', 'v_n_y': 'hist_ye_m',
                'color': color, 'ls': '-', 'lw': lw, 'ds': 'steps', 'alpha': alpha,
                'xmin': 0., 'xmax': 0.5, 'ymin': 1e-4, 'ymax': 1e0,
                'xlabel': r"$Y_e$", 'ylabel': r'$M_{ej}$ [normed]',
                'label': label + wind_lbl,
                'yscale': 'log',
                'fancyticks': True, 'minorticks': True,
                'sharey': True,
                'fontsize': fontsize,
                'labelsize': labelsize,
                'legend': {}  # {'loc': 'best', 'ncol': 2, 'fontsize': 18}
            }

            if sim == sims[-1]:
                dic_hist_theta_b['legend'] = {'loc': 'lower right', 'ncol': 3, 'fontsize': 12,
                                           'bbox_to_anchor': (2.40, 1.02)}

            o_plot.set_plot_dics.append(dic_hist_theta_b)
            o_plot.set_plot_dics.append(dic_hist_vel_inf_b)
            o_plot.set_plot_dics.append(dic_hist_ye_b)

    o_plot.main()
    exit(1)
# plot_hists_mult_sims_sly4()