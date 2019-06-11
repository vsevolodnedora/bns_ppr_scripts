from __future__ import division
from sys import path
path.append('modules/')
from general import *
from lists import *
from filework import *
import matplotlib
from matplotlib import pyplot
import matplotlib.pyplot as plt
plt.rc('text', usetex=True)
plt.rc('font', family='serif')

import copy

from plotting_nexus import PLOT_MANY_TASKS
from d3analysis import LOAD_INT_DATA, LOAD_DENSITY_MODES, LOAD_RES_CORR, ADD_METHODS_FOR_INT_DATA
from d1analysis import COMPUTE_STORE_PAR, NORMALIZE_NUCLEO
from mkn_interface import EXTRACT_LIGHTCURVE


""" -- --- FINALS --- --- """

def plot_jflux_m1_for_3sims():


    o_plot = PLOT_MANY_TASKS()
    o_plot.gen_set["figdir"] = Paths.plots + '/combine_test/'
    o_plot.gen_set["type"] = "polar"
    o_plot.gen_set["figsize"] = (9.0, 2.7)  # <->, |] # to match hists with (8.5, 2.7)
    o_plot.gen_set["figname"] = "jflux_final.pdf"
    o_plot.gen_set["sharex"] = False
    o_plot.gen_set["sharey"] = False
    o_plot.gen_set["subplots_adjust_h"] = 0.3
    o_plot.gen_set["subplots_adjust_w"] = 0.3
    o_plot.set_plot_dics = []


    ''' --- SIM 1 --- '''

    sim = "DD2_M13641364_M0_SR"

    o_int_data = ADD_METHODS_FOR_INT_DATA(sim)
    o_int_data.flag_force_unique_grid = True
    # o_int_data.it_for_unique_grid = 2254356 # grid is loaded only for this iteration and assumed to be constant


    int_ang_mom_flux_dic = {
        'task': '2d projection', 'dtype': 'int', 'ptype': 'polar',
        'data': o_int_data, 'it': 751278,
        'position': (1, 1), 'title': 'DD2', 'cbar': None,#'right .05 .0',
        'v_n_x': 'phi_cyl', 'v_n_y': 'r_cyl', 'v_n': 'ang_mom_flux',  'mod': 'integ_over_z',
        'phimin': None, 'phimax': None, 'rmin': 0, 'rmax': 50, 'vmin': -1e-6, 'vmax': 1e-6,
        'fill_vmin': False, # fills the x < vmin with vmin
        'xscale': None, 'yscale': None,
        'mask': None, 'cmap': 'RdBu_r', 'norm': "linear",
        'fancyticks': False
    }
    int_ang_mom_flux_dic_rev = {
        'task': '2d projection', 'dtype': 'int', 'ptype': 'polar',
        'data': o_int_data, 'it': 751278,
        'position': (1, 1), 'title': 'DD2', 'cbar': None,#'right .05 .0','cbar fmt':'%.1e',
        'v_n_x': 'phi_cyl', 'v_n_y': 'r_cyl', 'v_n': 'ang_mom_flux', 'mod': 'integ_over_z',
        'phimin': None, 'phimax': None, 'rmin': 0, 'rmax': 50, 'vmin': -5e-6, 'vmax': 5e-6,
        'fill_vmin': False, # fills the x < vmin with vmin
        'xscale': None, 'yscale': None,
        'mask': None, 'cmap': 'RdBu_r', 'norm': "linear",
        'fancyticks': False
    }

    o_dm_data = LOAD_DENSITY_MODES(sim)
    densmode = {
        'task': '2d projection', 'dtype': 'dm', 'ptype': 'polar',
        'data': o_dm_data, 'it': 751278,
        'name': 'densmode', 'position': (1, 1),
        'v_n_x': 'phi_cyl', 'v_n_y': 'r_cyl', 'v_n': 'int_phi',
        'mode': 1, #'int': 'spline',
        'rmax': 50, 'ls': '-', 'color': 'black',
    }
    densmode2 = copy.deepcopy(densmode)
    densmode2["v_n"] = 'int_phi_r'
    densmode2["ls"] = '--'

    o_plot.set_plot_dics.append(int_ang_mom_flux_dic)
    # o_plot.set_plot_dics.append(int_ang_mom_flux_dic_rev)
    o_plot.set_plot_dics.append(densmode)
    o_plot.set_plot_dics.append(densmode2)
    ''' --- SIM 2 --- '''

    sim = "LS220_M13641364_M0_SR"

    o_int_data = ADD_METHODS_FOR_INT_DATA(sim)
    o_int_data.flag_force_unique_grid = True
    # o_int_data.it_for_unique_grid = 2254356  # grid is loaded only for this iteration and assumed to be constant

    int_ang_mom_flux_dic = {
        'task': '2d projection', 'dtype': 'int', 'ptype': 'polar',
        'data': o_int_data, 'it': 753664,
        'position': (1, 2), 'title': 'LS220', 'cbar': None, #'right .05 .0',  'cbar fmt':'%.1e'
        'v_n_x': 'phi_cyl', 'v_n_y': 'r_cyl', 'v_n': 'ang_mom_flux', 'mod': 'integ_over_z',
        'phimin': None, 'phimax': None, 'rmin': 0, 'rmax': 50, 'vmin': -1e-6, 'vmax': 1e-6,
        'fill_vmin': False,  # fills the x < vmin with vmin
        'xscale': None, 'yscale': None,
        'mask': None, 'cmap': 'RdBu_r', 'norm': "linear",
        'fancyticks': False
    }
    int_ang_mom_flux_dic_rev = {
        'task': '2d projection', 'dtype': 'int', 'ptype': 'polar',
        'data': o_int_data, 'it': 753664,
        'position': (1, 2), 'title': 'LS220', 'cbar': None,#'right .05 .0',
        'v_n_x': 'phi_cyl', 'v_n_y': 'r_cyl', 'v_n': 'ang_mom_flux',
        'phimin': None, 'phimax': None, 'rmin': 0, 'rmax': 50, 'vmin': 1e-9, 'vmax': 1e-5,
        'fill_vmin': False,  # fills the x < vmin with vmin
        'xscale': None, 'yscale': None,
        'mask': 'positive', 'cmap': 'Blues', 'norm': "log", 'todo': None,
        'fancyticks': False
    }

    o_dm_data = LOAD_DENSITY_MODES(sim)

    densmode = {
        'task': '2d projection', 'dtype': 'dm', 'ptype': 'polar',
        'data': o_dm_data, 'it': 753664,
        'name': 'densmode', 'position': (1, 2),
        'v_n_x': 'phi_cyl', 'v_n_y': 'r_cyl', 'v_n': 'int_phi',
        'mode': 1,
        'rmax': 50, 'ls': '-', 'color': 'black',
    }
    densmode2 = copy.deepcopy(densmode)
    densmode2["v_n"] = 'int_phi_r'
    densmode2["ls"] = '--'
    # o_plot.set_plot_dics.append(int_ang_mom_flux_dic_rev)
    o_plot.set_plot_dics.append(int_ang_mom_flux_dic)
    # o_plot.set_plot_dics.append(int_ang_mom_flux_dic_rev)
    o_plot.set_plot_dics.append(densmode)
    o_plot.set_plot_dics.append(densmode2)

    ''' --- SIM 3 --- '''

    sim = "SLy4_M13641364_M0_SR"

    o_int_data = ADD_METHODS_FOR_INT_DATA(sim)
    o_int_data.flag_force_unique_grid = True
    # o_int_data.it_for_unique_grid = 2254356  # grid is loaded only for this iteration and assumed to be constant

    int_ang_mom_flux_dic = {
        'task': '2d projection', 'dtype': 'int', 'ptype': 'polar',
        'data': o_int_data, 'it': 761856,
        'position': (1, 3), 'title': 'SLy4', 'cbar': 'right .03 .0','cbar fmt':'%.1e',
        'v_n_x': 'phi_cyl', 'v_n_y': 'r_cyl', 'v_n': 'ang_mom_flux', 'mod': 'integ_over_z', 'cbar label':r"$\int(\dot{J}dz)$",
        'phimin': None, 'phimax': None, 'rmin': 0, 'rmax': 50, 'vmin': -1e-6, 'vmax': 1e-6,
        'fill_vmin': True,  # fills the x < vmin with vmin
        'xscale': None, 'yscale': None,
         'mask': None, 'cmap': 'RdBu_r', 'norm': "linear",
        'fancyticks': False
    }
    int_ang_mom_flux_dic_rev = {
        'task': '2d projection', 'dtype': 'int', 'ptype': 'polar',
        'data': o_int_data, 'it': 761856,
        'position': (1, 3), 'title': 'SLy4', 'cbar': None,
        'v_n_x': 'phi_cyl', 'v_n_y': 'r_cyl', 'v_n': 'ang_mom_flux',
        'phimin': None, 'phimax': None, 'rmin': 0, 'rmax': 50, 'vmin': 1e-9, 'vmax': 1e-5,
        'fill_vmin': True,  # fills the x < vmin with vmin
        'xscale': None, 'yscale': None,
        'mask': 'positive', 'cmap': 'Blues', 'norm': "log", 'todo': None,
        'fancyticks': False
    }

    o_dm_data = LOAD_DENSITY_MODES(sim)

    densmode = {
        'task': '2d projection', 'dtype': 'dm', 'ptype': 'polar',
        'data': o_dm_data, 'it': 761856,
        'name': 'densmode', 'position': (1, 3),
        'v_n_x': 'phi_cyl', 'v_n_y': 'r_cyl', 'v_n': 'int_phi',
        'mode': 1,
        'rmax': 50, 'ls': '-', 'color': 'black',
    }
    densmode2 = copy.deepcopy(densmode)
    densmode2["v_n"] = 'int_phi_r'
    densmode2["ls"] = '--'

    o_plot.set_plot_dics.append(int_ang_mom_flux_dic)
    # o_plot.set_plot_dics.append(int_ang_mom_flux_dic_rev)
    o_plot.set_plot_dics.append(densmode)
    o_plot.set_plot_dics.append(densmode2)


    o_plot.main()
# plot_jflux_m1_for_3sims()


def plot_ej_profiles_with_viscosity_3sim():

    o_plot = PLOT_MANY_TASKS()
    o_plot.gen_set["figdir"] = Paths.plots + '/combine_test/'
    o_plot.gen_set["type"] = "cartesian"
    o_plot.gen_set["figsize"] = (9.0, 2.7)  # <->, |]
    o_plot.gen_set["figname"] = "ej_profs_final.pdf"
    o_plot.gen_set["sharex"] = False
    o_plot.gen_set["sharey"] = False
    o_plot.gen_set["subplots_adjust_h"] = 0.3
    o_plot.gen_set["subplots_adjust_w"] = 0.0
    o_plot.set_plot_dics = []

    ''' --- --- sim 1 --- --- --- '''

    sim = "DD2_M13641364_M0_SR"
    o_data1 = COMPUTE_STORE_PAR(sim)


    sim = "DD2_M13641364_M0_LK_SR_R04"
    o_data2 = COMPUTE_STORE_PAR(sim)


    dic_band = {
        'task': 'ejband', 'ptype': 'cartesian',
        'position': (1, 1),
        'data1': o_data1, 'data2': o_data2,
        'criterion1': '_0', 'criterion2': '_0',
        'v_n_x': 't_tot_flux', 'v_n_y': 'mass_tot_flux',
        'color': 'green', 'alpha':0.4,
        # 'ymin': 1e-4, 'ymax': 1e0,
        'xlabel': r"time [ms]", 'ylabel': r'$M_{ej}$ $[10^{-2}M_{\odot}]$',
        'xunits': 'ms', 'yunits': '1e-2Msun',
        'label': 'Geodesic', 'yscale': None, 'title': 'DD2',
        'fancyticks': True, 'minorticks': True
    }
    o_plot.set_plot_dics.append(dic_band)
    dic_ej_prof = {
        'task': 'ejprof', 'ptype': 'cartesian',
        'position': (1, 1),
        'data': o_data1, 'criterion': '_0',
        'v_n_x': 't_tot_flux', 'v_n_y': 'mass_tot_flux',
        'color': 'green', 'ls': '--', 'lw': 1.,
        # 'ymin': 1e-4, 'ymax': 1e0,
        'xlabel': r"time [ms]", 'ylabel': r'$M_{ej}$ $[10^{-2}M_{\odot}]$',
        'xunits': 'ms', '-t': o_data1.get_par("tmerger_gw"), 'yunits': '1e-2Msun',
        'label': None, 'yscale': None, 'title': None,
        'fancyticks': True, 'minorticks': True
    }
    o_plot.set_plot_dics.append(dic_ej_prof)
    dic_ej_prof = {
        'task': 'ejprof', 'ptype': 'cartesian',
        'position': (1, 1),
        'data': o_data2, 'criterion': '_0',
        'v_n_x': 't_tot_flux', 'v_n_y': 'mass_tot_flux',
        'color': 'lightgreen', 'ls': '--', 'lw': 1.,
        'ymin': 0, 'ymax': 1.4,
        'xlabel': r"time [ms]", 'ylabel': r'$M_{ej}$ $[10^{-2}M_{\odot}]$',
        'xunits': 'ms', '-t': o_data2.get_par("tmerger_gw"), 'yunits': '1e-2Msun',
        'label': None, 'yscale': None, 'title': None,
        'fancyticks': True, 'minorticks': True,
        'legend': True
    }
    o_plot.set_plot_dics.append(dic_ej_prof)


    """ --------------------BERNULLI------------------- """

    dic_band = {
        'task': 'ejband', 'ptype': 'cartesian',
        'position': (1, 1),
        'data1': o_data1, 'data2': o_data2,
        'criterion1': '_0_b_w', 'criterion2': '_0_b_w',
        'v_n_x': 't_tot_flux', 'v_n_y': 'mass_tot_flux',
        'color': 'green', 'alpha':0.6,
        # 'ymin': 1e-4, 'ymax': 1e0,
        'xlabel': r"time [ms]", 'ylabel': r'$M_{ej}$ $[10^{-2}M_{\odot}]$',
        'xunits': 'ms', 'yunits': '1e-2Msun',
        'label': 'Bernoulli', 'yscale': None, 'title': 'DD2',
        'fancyticks': True, 'minorticks': True
    }
    o_plot.set_plot_dics.append(dic_band)
    dic_ej_prof = {
        'task': 'ejprof', 'ptype': 'cartesian',
        'position': (1, 1),
        'data': o_data1, 'criterion': '_0_b_w',
        'v_n_x': 't_tot_flux', 'v_n_y': 'mass_tot_flux',
        'color': 'green', 'ls': '-', 'lw': 1.,
        # 'ymin': 1e-4, 'ymax': 1e0,
        'xlabel': r"time [ms]", 'ylabel': r'$M_{ej}$ $[10^{-2}M_{\odot}]$',
        'xunits': 'ms', '-t': o_data1.get_par("tmerger_gw"), 'yunits': '1e-2Msun',
        'label': 'no viscosity', 'yscale': None, 'title': None,
        'fancyticks': True, 'minorticks': True,
        'legend': True
    }
    o_plot.set_plot_dics.append(dic_ej_prof)
    dic_ej_prof = {
        'task': 'ejprof', 'ptype': 'cartesian',
        'position': (1, 1),
        'data': o_data2, 'criterion': '_0_b_w',
        'v_n_x': 't_tot_flux', 'v_n_y': 'mass_tot_flux',
        'color': 'lightgreen', 'ls': '-', 'lw': 1.,
        # 'ymin': 1e-4, 'ymax': 1e0,
        'xlabel': r"time [ms]", 'ylabel': r'$M_{ej}$ $[10^{-2}M_{\odot}]$',
        'xunits': 'ms', '-t': o_data2.get_par("tmerger_gw"), 'yunits': '1e-2Msun',
        'label': 'viscosity', 'yscale': None, 'title': None,
        'fancyticks': True, 'minorticks': True,
        'legend':  True
    }
    o_plot.set_plot_dics.append(dic_ej_prof)

    ''' --- --- ---- ----- ---- - '''
    #
    #
    #
    ''' --- --- sim 2 --- --- --- '''

    sim = "LS220_M13641364_M0_SR"
    o_data1 = COMPUTE_STORE_PAR(sim)

    sim = "LS220_M13641364_M0_LK_SR"
    o_data2 = COMPUTE_STORE_PAR(sim)

    dic_band = {
        'task': 'ejband', 'ptype': 'cartesian',
        'position': (1, 2),
        'data1': o_data1, 'data2': o_data2,
        'criterion1': '_0', 'criterion2': '_0',
        'v_n_x': 't_tot_flux', 'v_n_y': 'mass_tot_flux',
        'color': 'blue', 'alpha': 0.4,
        # 'ymin': 1e-4, 'ymax': 1e0,
        'xlabel': r"time [ms]", 'ylabel': r'$M_{ej}$ $[10^{-2}M_{\odot}]$',
        'xunits': 'ms', 'yunits': '1e-2Msun',
        'label': 'Geodesic', 'yscale': None, 'title': 'LS220',
        'fancyticks': True, 'minorticks': True
    }
    o_plot.set_plot_dics.append(dic_band)
    dic_ej_prof = {
        'task': 'ejprof', 'ptype': 'cartesian',
        'position': (1, 2),
        'data': o_data1, 'criterion': '_0',
        'v_n_x': 't_tot_flux', 'v_n_y': 'mass_tot_flux',
        'color': 'blue', 'ls': '--', 'lw': 1.,
        # 'ymin': 1e-4, 'ymax': 1e0,
        'xlabel': r"time [ms]", 'ylabel': r'$M_{ej}$ $[10^{-2}M_{\odot}]$',
        'xunits': 'ms', '-t': o_data1.get_par("tmerger_gw"), 'yunits': '1e-2Msun',
        'label': None, 'yscale': None, 'title': None,
        'fancyticks': True, 'minorticks': True
    }
    o_plot.set_plot_dics.append(dic_ej_prof)
    dic_ej_prof = {
        'task': 'ejprof', 'ptype': 'cartesian',
        'position': (1, 2),
        'data': o_data2, 'criterion': '_0',
        'v_n_x': 't_tot_flux', 'v_n_y': 'mass_tot_flux',
        'color': 'lightblue', 'ls': '--', 'lw': 1.,
        'ymin': 0, 'ymax': 1.4,
        'xlabel': r"time [ms]", 'ylabel': r'$M_{ej}$ $[10^{-2}M_{\odot}]$',
        'xunits': 'ms', '-t': o_data2.get_par("tmerger_gw"), 'yunits': '1e-2Msun',
        'label': None, 'yscale': None, 'title': None,
        'fancyticks': True, 'minorticks': True
    }
    o_plot.set_plot_dics.append(dic_ej_prof)

    """ --------------------BERNULLI------------------- """

    dic_band = {
        'task': 'ejband', 'ptype': 'cartesian',
        'position': (1, 2),
        'data1': o_data1, 'data2': o_data2,
        'criterion1': '_0_b_w', 'criterion2': '_0_b_w',
        'v_n_x': 't_tot_flux', 'v_n_y': 'mass_tot_flux',
        'color': 'blue', 'alpha': 0.6,
        # 'ymin': 1e-4, 'ymax': 1e0,
        'xlabel': r"time [ms]", 'ylabel': r'$M_{ej}$ $[10^{-2}M_{\odot}]$',
        'xunits': 'ms', 'yunits': '1e-2Msun',
        'label': 'Bernoulli', 'yscale': None, 'title': 'LS220',
        'fancyticks': True, 'minorticks': True
    }
    o_plot.set_plot_dics.append(dic_band)
    dic_ej_prof = {
        'task': 'ejprof', 'ptype': 'cartesian',
        'position': (1, 2),
        'data': o_data1, 'criterion': '_0_b_w',
        'v_n_x': 't_tot_flux', 'v_n_y': 'mass_tot_flux',
        'color': 'blue', 'ls': '-', 'lw': 1.,
        # 'ymin': 1e-4, 'ymax': 1e0,
        'xlabel': r"time [ms]", 'ylabel': r'$M_{ej}$ $[10^{-2}M_{\odot}]$',
        'xunits': 'ms', '-t': o_data1.get_par("tmerger_gw"), 'yunits': '1e-2Msun',
        'label': 'no viscosity', 'yscale': None, 'title': None,
        'fancyticks': True, 'minorticks': True
    }
    o_plot.set_plot_dics.append(dic_ej_prof)
    dic_ej_prof = {
        'task': 'ejprof', 'ptype': 'cartesian',
        'position': (1, 2),
        'data': o_data2, 'criterion': '_0_b_w',
        'v_n_x': 't_tot_flux', 'v_n_y': 'mass_tot_flux',
        'color': 'lightblue', 'ls': '-', 'lw': 1.,
        'ymin': 0, 'ymax': 1.4,
        'xlabel': r"time [ms]", 'ylabel': r'$M_{ej}$ $[10^{-2}M_{\odot}]$',
        'xunits': 'ms', '-t': o_data2.get_par("tmerger_gw"), 'yunits': '1e-2Msun',
        'label': 'viscosity', 'yscale': None, 'title': None,
        'fancyticks': True, 'minorticks': True,
        'sharey': True,
        'legend':True
    }
    o_plot.set_plot_dics.append(dic_ej_prof)

    coll_time_vert = {
        'task': 'vertline', 'dtype': '-', 'ptype': 'cartesian',
        'value': (o_data1.get_par("tcoll_gw") - o_data1.get_par("tmerger_gw"))*1e3,
        'position': (1, 2),
        'ls': '-.', 'color': 'blue', 'lw': 0.5,
    }
    o_plot.set_plot_dics.append(coll_time_vert)
    coll_time_vert = {
        'task': 'vertline', 'dtype': '-', 'ptype': 'cartesian',
        'value': (o_data2.get_par("tcoll_gw") - o_data2.get_par("tmerger_gw"))*1e3,
        'position': (1, 2),
        'ls': '-.', 'color': 'lightblue', 'lw': 0.5,
    }
    o_plot.set_plot_dics.append(coll_time_vert)

    ''' --- --- ---- ----- ---- - '''
    #
    #
    #
    ''' --- --- sim 2 --- --- --- '''

    sim = "SLy4_M13641364_M0_SR"
    o_data1 = COMPUTE_STORE_PAR(sim)

    sim = "SLy4_M13641364_M0_LK_SR"
    # o_data2 = COMPUTE_STORE_PAR(sim)

    # dic_band = {
    #     'task': 'ejband', 'ptype': 'cartesian',
    #     'position': (1, 2),
    #     'data1': o_data1, 'data2': o_data2,
    #     'criterion1': '_0', 'criterion2': '_0',
    #     'v_n_x': 't_tot_flux', 'v_n_y': 'mass_tot_flux',
    #     'color': 'blue', 'alpha': 0.4,
    #     # 'ymin': 1e-4, 'ymax': 1e0,
    #     'xlabel': r"time [ms]", 'ylabel': r'$M_{ej}$ $[10^{-2}M_{\odot}]$',
    #     'xunits': 'ms',
    #     'label': None, 'yscale': None, 'title': 'LS220',
    #     'fancyticks': True, 'minorticks': True
    # }
    # o_plot.set_plot_dics.append(dic_band)
    dic_ej_prof = {
        'task': 'ejprof', 'ptype': 'cartesian',
        'position': (1, 3),
        'data': o_data1, 'criterion': '_0',
        'v_n_x': 't_tot_flux', 'v_n_y': 'mass_tot_flux',
        'color': 'red', 'ls': '--', 'lw': 1.,
        # 'ymin': 1e-4, 'ymax': 1e0,
        'xlabel': r"time [ms]", 'ylabel': r'$M_{ej}$ $[10^{-2}M_{\odot}]$',
        'xunits': 'ms', '-t': o_data1.get_par("tmerger_gw"), 'yunits': '1e-2Msun',
        'label': 'Geodesic', 'yscale': None, 'title': "SLy4",
        'fancyticks': True, 'minorticks': True,
        'sharey': True
    }
    o_plot.set_plot_dics.append(dic_ej_prof)
    # dic_ej_prof = {
    #     'task': 'ejprof', 'ptype': 'cartesian',
    #     'position': (1, 2),
    #     'data': o_data2, 'criterion': '_0',
    #     'v_n_x': 't_tot_flux', 'v_n_y': 'mass_tot_flux',
    #     'color': 'lightblue', 'ls': '--', 'lw': 1.,
    #     # 'ymin': 1e-4, 'ymax': 1e0,
    #     'xlabel': r"time [ms]", 'ylabel': r'$M_{ej}$ $[10^{-2}M_{\odot}]$',
    #     'xunits': 'ms', '-t': o_data2.get_par("tmerger_gw"),
    #     'label': None, 'yscale': None, 'title': None,
    #     'fancyticks': True, 'minorticks': True
    # }
    # o_plot.set_plot_dics.append(dic_ej_prof)

    """ --------------------BERNULLI------------------- """

    # dic_band = {
    #     'task': 'ejband', 'ptype': 'cartesian',
    #     'position': (1, 2),
    #     'data1': o_data1, 'data2': o_data2,
    #     'criterion1': '_0_b_w', 'criterion2': '_0_b_w',
    #     'v_n_x': 't_tot_flux', 'v_n_y': 'mass_tot_flux',
    #     'color': 'blue', 'alpha': 0.6,
    #     # 'ymin': 1e-4, 'ymax': 1e0,
    #     'xlabel': r"time [ms]", 'ylabel': r'$M_{ej}$ $[10^{-2}M_{\odot}]$',
    #     'xunits': 'ms',
    #     'label': None, 'yscale': None, 'title': 'DD2',
    #     'fancyticks': True, 'minorticks': True
    # }
    # o_plot.set_plot_dics.append(dic_band)
    dic_ej_prof = {
        'task': 'ejprof', 'ptype': 'cartesian',
        'position': (1, 3),
        'data': o_data1, 'criterion': '_0_b_w',
        'v_n_x': 't_tot_flux', 'v_n_y': 'mass_tot_flux',
        'color': 'red', 'ls': '-', 'lw': 1.,
        'ymin': 0, 'ymax': 1.4,
        'xlabel': r"time [ms]", 'ylabel': r'$M_{ej}$ $[10^{-2}M_{\odot}]$',
        'xunits': 'ms', '-t': o_data1.get_par("tmerger_gw"), 'yunits': '1e-2Msun',
        'label': 'Bernoulli', 'yscale': None, 'title': None,
        'fancyticks': True, 'minorticks': True,
        'sharey':True,
        'legend': True
    }
    o_plot.set_plot_dics.append(dic_ej_prof)
    # dic_ej_prof = {
    #     'task': 'ejprof', 'ptype': 'cartesian',
    #     'position': (1, 2),
    #     'data': o_data2, 'criterion': '_0_b_w',
    #     'v_n_x': 't_tot_flux', 'v_n_y': 'mass_tot_flux',
    #     'color': 'lightblue', 'ls': '-', 'lw': 1.,
    #     # 'ymin': 1e-4, 'ymax': 1e0,
    #     'xlabel': r"time [ms]", 'ylabel': r'$M_{ej}$ $[10^{-2}M_{\odot}]$',
    #     'xunits': 'ms', '-t': o_data2.get_par("tmerger_gw"),
    #     'label': None, 'yscale': None, 'title': None,
    #     'fancyticks': True, 'minorticks': True
    # }
    # o_plot.set_plot_dics.append(dic_ej_prof)

    coll_time_vert = {
        'task': 'vertline', 'dtype': '-', 'ptype': 'cartesian',
        'value': (o_data1.get_par("tcoll_gw") - o_data1.get_par("tmerger_gw")) * 1e3,
        'position': (1, 3),
        'ls': '-.', 'color': 'red', 'lw': 0.5,

    }
    o_plot.set_plot_dics.append(coll_time_vert)
    # coll_time_vert = {
    #     'task': 'vertline', 'dtype': '-', 'ptype': 'cartesian',
    #     'value': (o_data2.get_par("tcoll_gw") - o_data2.get_par("tmerger_gw")) * 1e3,
    #     'position': (1, 2),
    #     'ls': '-.', 'color': 'lightblue', 'lw': 0.5,
    # }
    # o_plot.set_plot_dics.append(coll_time_vert)

    o_plot.main()
# plot_ej_profiles_with_viscosity_3sim()

def plot_3hists_for_3sims_in_one_row():

    o_plot = PLOT_MANY_TASKS()
    o_plot.gen_set["figdir"] = Paths.plots + '/combine_test/'
    o_plot.gen_set["type"] = "cartesian"
    o_plot.gen_set["figsize"] = (9.0, 2.7)  # <->, |]
    o_plot.gen_set["figname"] = "histograms_final.pdf"
    o_plot.gen_set["sharex"] = False
    o_plot.gen_set["sharey"] = False
    o_plot.gen_set["subplots_adjust_h"] = 0.3
    o_plot.gen_set["subplots_adjust_w"] = 0.0

    ''' --- --- sim 1 --- --- --- '''

    sim = "DD2_M13641364_M0_SR"
    o_data = COMPUTE_STORE_PAR(sim)

    dic_hist_theta = {
        'task': 'hist1d', 'ptype': 'cartesian',
        'position': (1, 1),
        'data': o_data, 'norm': True, 'criterion': '_0',
        'v_n_x': 'hist_theta', 'v_n_y': 'hist_theta_m',
        'color': 'green', 'ls': ':', 'lw': 0.6, 'ds': 'steps', 'alpha':0.6,
        'ymin': 1e-4, 'ymax': 1e0,
        'xlabel': r"Angle from orbital plane", 'ylabel': 'mass',
        'label': None, 'yscale': 'log',
        'fancyticks': True, 'minorticks': True
    }
    dic_hist_theta_b = {
        'task': 'hist1d', 'ptype': 'cartesian',
        'position': (1, 1),
        'data': o_data, 'norm': True, 'criterion': '_0_b_w',
        'v_n_x': 'hist_theta', 'v_n_y': 'hist_theta_m',
        'color': 'green', 'ls': '-', 'lw': 1., 'ds': 'steps',
        'ymin': 1e-4, 'ymax': 1e0,
        'xlabel': r"Angle from orbital plane", 'ylabel': 'mass',
        'label': "DD2", 'legend':True,
        'yscale': 'log',
        'fancyticks': True, 'minorticks': True,
    }

    dic_hist_vel_inf = {
        'task': 'hist1d', 'ptype': 'cartesian',
        'position': (1, 2),
        'data': o_data, 'norm': True, 'criterion': '_0',
        'v_n_x': 'hist_vel_inf', 'v_n_y': 'hist_vel_inf_m',
        'color': 'green', 'ls': ':', 'lw': 0.6, 'ds': 'steps', 'alpha':0.6,
        'xmin': 0., 'xmax': 1., 'ymin': 1e-4, 'ymax': 1e0,
        'xlabel': r"$\upsilon_{\infty}$ [c]", 'ylabel': 'mass',
        'label': None, 'yscale': 'log',
        'fancyticks': True, 'minorticks': True,
        'sharey': True
    }
    dic_hist_vel_inf_b = {
        'task': 'hist1d', 'ptype': 'cartesian',
        'position': (1, 2),
        'data': o_data, 'norm': True, 'criterion': '_0_b_w',
        'v_n_x': 'hist_vel_inf', 'v_n_y': 'hist_vel_inf_m',
        'color': 'green', 'ls': '-', 'lw': 1., 'ds': 'steps',
        'xmin': 0., 'xmax': 1., 'ymin': 1e-4, 'ymax': 1e0,
        'xlabel': r"$\upsilon_{\infty}$ [c]", 'ylabel': 'mass',
        'label': None, 'yscale': 'log',
        'fancyticks': True, 'minorticks': True,
        'sharey': True
    }

    dic_hist_ye = {
        'task': 'hist1d', 'ptype': 'cartesian',
        'position': (1, 3),
        'data': o_data, 'norm': True, 'criterion': '_0',
        'v_n_x': 'hist_ye', 'v_n_y': 'hist_ye_m',
        'color': 'green', 'ls': ':', 'lw': 0.6, 'ds': 'steps', 'alpha':0.6,
        'xmin': 0., 'xmax': 0.5, 'ymin': 1e-4, 'ymax': 1e0,
        'xlabel': r"$Y_e$", 'ylabel': 'mass',
        'label': None, 'yscale': 'log',
        'fancyticks': True, 'minorticks': True,
        'sharey': True
    }
    dic_hist_ye_b = {
        'task': 'hist1d', 'ptype': 'cartesian',
        'position': (1, 3),
        'data': o_data, 'norm': True, 'criterion': '_0_b_w',
        'v_n_x': 'hist_ye', 'v_n_y': 'hist_ye_m',
        'color': 'green', 'ls': '-', 'lw': 1., 'ds': 'steps',
        'xmin': 0., 'xmax': 0.5, 'ymin': 1e-4, 'ymax': 1e0,
        'xlabel': r"$Y_e$", 'ylabel': 'mass',
        'label': None, 'yscale': 'log',
        'fancyticks': True, 'minorticks': True,
        'sharey': True
    }

    ''' --- --- sim 2 --- --- --- '''

    sim = "LS220_M13641364_M0_SR"
    o_data = COMPUTE_STORE_PAR(sim)

    dic_hist_theta_LS220 = {
        'task': 'hist1d', 'ptype': 'cartesian',
        'position': (1, 1),
        'data': o_data, 'norm': True, 'criterion': '_0',
        'v_n_x': 'hist_theta', 'v_n_y': 'hist_theta_m',
        'color': 'blue', 'ls': ':', 'lw': 0.6, 'ds': 'steps', 'alpha': 0.6,
        'ymin': 1e-4, 'ymax': 1e0,
        'xlabel': r"Angle from orbital plane", 'ylabel': 'mass',
        'label': None, 'yscale': 'log',
        'fancyticks': True, 'minorticks': True
    }
    dic_hist_theta_b_LS220 = {
        'task': 'hist1d', 'ptype': 'cartesian',
        'position': (1, 1),
        'data': o_data, 'norm': True, 'criterion': '_0_b_w',
        'v_n_x': 'hist_theta', 'v_n_y': 'hist_theta_m',
        'color': 'blue', 'ls': '-', 'lw': 1., 'ds': 'steps',
        'ymin': 1e-4, 'ymax': 1e0,
        'xlabel': r"Angle from orbital plane", 'ylabel': 'mass',
        'label': "LS220", 'yscale': 'log',
        'fancyticks': True, 'minorticks': True,
    }

    dic_hist_vel_inf_LS220 = {
        'task': 'hist1d', 'ptype': 'cartesian',
        'position': (1, 2),
        'data': o_data, 'norm': True, 'criterion': '_0',
        'v_n_x': 'hist_vel_inf', 'v_n_y': 'hist_vel_inf_m',
        'color': 'blue', 'ls': ':', 'lw': 0.6, 'ds': 'steps', 'alpha': 0.6,
        'xmin': 0., 'xmax': 1., 'ymin': 1e-4, 'ymax': 1e0,
        'xlabel': r"$\upsilon_{\infty}$ [c]", 'ylabel': 'mass',
        'label': "Geodesic", 'yscale': 'log',
        'fancyticks': True, 'minorticks': True,
        'sharey': True
    }
    dic_hist_vel_inf_b_LS220= {
        'task': 'hist1d', 'ptype': 'cartesian',
        'position': (1, 2),
        'data': o_data, 'norm': True, 'criterion': '_0_b_w',
        'v_n_x': 'hist_vel_inf', 'v_n_y': 'hist_vel_inf_m',
        'color': 'blue', 'ls': '-', 'lw': 1., 'ds': 'steps',
        'xmin': 0., 'xmax': 1., 'ymin': 1e-4, 'ymax': 1e0,
        'xlabel': r"$\upsilon_{\infty}$ [c]", 'ylabel': 'mass',
        'label': "Bernoulli", 'legend':True,
        'yscale': 'log',
        'fancyticks': True, 'minorticks': True,
        'sharey': True
    }

    dic_hist_ye_LS220 = {
        'task': 'hist1d', 'ptype': 'cartesian',
        'position': (1, 3),
        'data': o_data, 'norm': True, 'criterion': '_0',
        'v_n_x': 'hist_ye', 'v_n_y': 'hist_ye_m',
        'color': 'blue', 'ls': ':', 'lw': 0.6, 'ds': 'steps', 'alpha': 0.6,
        'xmin': 0., 'xmax': 0.5, 'ymin': 1e-4, 'ymax': 1e0,
        'xlabel': r"$Y_e$", 'ylabel': 'mass',
        'label': None, 'yscale': 'log',
        'fancyticks': True, 'minorticks': True,
        'sharey': True
    }
    dic_hist_ye_b_LS220 = {
        'task': 'hist1d', 'ptype': 'cartesian',
        'position': (1, 3),
        'data': o_data, 'norm': True, 'criterion': '_0_b_w',
        'v_n_x': 'hist_ye', 'v_n_y': 'hist_ye_m',
        'color': 'blue', 'ls': '-', 'lw': 1., 'ds': 'steps',
        'xmin': 0., 'xmax': 0.5, 'ymin': 1e-4, 'ymax': 1e0,
        'xlabel': r"$Y_e$", 'ylabel': 'mass',
        'label': None, 'yscale': 'log',
        'fancyticks': True, 'minorticks': True,
        'sharey': True
    }

    ''' --- --- sim 3 --- --- --- '''

    sim = "SLy4_M13641364_M0_SR"
    o_data = COMPUTE_STORE_PAR(sim)

    dic_hist_theta_DD2 = {
        'task': 'hist1d', 'ptype': 'cartesian',
        'position': (1, 1),
        'data': o_data, 'norm': True, 'criterion': '_0',
        'v_n_x': 'hist_theta', 'v_n_y': 'hist_theta_m',
        'color': 'red', 'ls': ':', 'lw': 0.6, 'ds': 'steps', 'alpha': 0.6,
        'ymin': 1e-4, 'ymax': 1e0,
        'xlabel': r"Angle from orbital plane", 'ylabel': 'mass',
        'label': None, 'yscale': 'log',
        'fancyticks': True, 'minorticks': True
    }
    dic_hist_theta_b_DD2 = {
        'task': 'hist1d', 'ptype': 'cartesian',
        'position': (1, 1),
        'data': o_data, 'norm': True, 'criterion': '_0_b_w',
        'v_n_x': 'hist_theta', 'v_n_y': 'hist_theta_m',
        'color': 'red', 'ls': '-', 'lw': 1., 'ds': 'steps',
        'ymin': 1e-4, 'ymax': 1e0,
        'xlabel': r"Angle from orbital plane", 'ylabel': 'mass',
        'label': "SLy4", 'yscale': 'log',
        'fancyticks': True, 'minorticks': True,
        'legend':True
    }

    dic_hist_vel_inf_DD2 = {
        'task': 'hist1d', 'ptype': 'cartesian',
        'position': (1, 2),
        'data': o_data, 'norm': True, 'criterion': '_0',
        'v_n_x': 'hist_vel_inf', 'v_n_y': 'hist_vel_inf_m',
        'color': 'red', 'ls': ':', 'lw': 0.6, 'ds': 'steps', 'alpha': 0.6,
        'xmin': 0., 'xmax': 1., 'ymin': 1e-4, 'ymax': 1e0,
        'xlabel': r"$\upsilon_{\infty}$ [c]", 'ylabel': 'mass',
        'label': None, 'yscale': 'log',
        'fancyticks': True, 'minorticks': True,
        'sharey': True
    }
    dic_hist_vel_inf_b_DD2= {
        'task': 'hist1d', 'ptype': 'cartesian',
        'position': (1, 2),
        'data': o_data, 'norm': True, 'criterion': '_0_b_w',
        'v_n_x': 'hist_vel_inf', 'v_n_y': 'hist_vel_inf_m',
        'color': 'red', 'ls': '-', 'lw': 1., 'ds': 'steps',
        'xmin': 0., 'xmax': 1., 'ymin': 1e-4, 'ymax': 1e0,
        'xlabel': r"$\upsilon_{\infty}$ [c]", 'ylabel': 'mass',
        'label': None, 'yscale': 'log',
        'fancyticks': True, 'minorticks': True,
        'sharey': True
    }

    dic_hist_ye_DD2 = {
        'task': 'hist1d', 'ptype': 'cartesian',
        'position': (1, 3),
        'data': o_data, 'norm': True, 'criterion': '_0',
        'v_n_x': 'hist_ye', 'v_n_y': 'hist_ye_m',
        'color': 'red', 'ls': ':', 'lw': 0.6, 'ds': 'steps', 'alpha': 0.6,
        'xmin': 0., 'xmax': 0.5, 'ymin': 1e-4, 'ymax': 1e0,
        'xlabel': r"$Y_e$", 'ylabel': 'mass',
        'label': None, 'yscale': 'log',
        'fancyticks': True, 'minorticks': True,
        'sharey': True
    }
    dic_hist_ye_b_DD2 = {
        'task': 'hist1d', 'ptype': 'cartesian',
        'position': (1, 3),
        'data': o_data, 'norm': True, 'criterion': '_0_b_w',
        'v_n_x': 'hist_ye', 'v_n_y': 'hist_ye_m',
        'color': 'red', 'ls': '-', 'lw': 1., 'ds': 'steps',
        'xmin': 0., 'xmax': 0.5, 'ymin': 1e-4, 'ymax': 1e0,
        'xlabel': r"$Y_e$", 'ylabel': 'mass',
        'label': None, 'legend':True,
        'yscale': 'log',
        'fancyticks': True, 'minorticks': True,
        'sharey': True
    }


    o_plot.set_plot_dics = [
        # dic_ej_prof, dic_ej_prof_b,
        dic_hist_theta, dic_hist_theta_b,
        dic_hist_vel_inf, dic_hist_vel_inf_b,
        dic_hist_ye, dic_hist_ye_b,

        dic_hist_theta_LS220, dic_hist_theta_b_LS220,
        dic_hist_vel_inf_LS220, dic_hist_vel_inf_b_LS220,
        dic_hist_ye_LS220, dic_hist_ye_b_LS220,
        #
        dic_hist_theta_DD2, dic_hist_theta_b_DD2,
        dic_hist_vel_inf_DD2, dic_hist_vel_inf_b_DD2,
        dic_hist_ye_DD2, dic_hist_ye_b_DD2

    ]
    o_plot.main()
# plot_3hists_for_3sims_in_one_row()

def plot_m1_in_time_3sims_tcoll():
    o_plot = PLOT_MANY_TASKS()
    o_plot.gen_set["figdir"] = Paths.plots + '/combine_test/'
    o_plot.gen_set["type"] = "cartesian"
    o_plot.gen_set["figsize"] = (9.0, 2.7)  # <->, |]
    o_plot.gen_set["figname"] = "dens_modes_final.pdf"
    o_plot.gen_set["sharex"] = True
    o_plot.gen_set["sharey"] = False
    o_plot.set_plot_dics = []

    " --- --- --- --- --- SIM0 --- --- --- ---  --- --- --- "

    sim = "DD2_M13641364_M0_SR"
    o_1d_data = COMPUTE_STORE_PAR(sim)

    o_dm_data = LOAD_DENSITY_MODES(sim)
    o_dm_data.gen_set['fname'] = Paths.ppr_sims + sim + '/res_3d/' + "density_modes_lap15.h5"
    # densmode_m0 = {
    #     'task': 'line', 'dtype': 'dm', 'ptype': 'cartesian',
    #     'data': o_dm_data,
    #     'position': (1, 1),
    #     'v_n_x': 'times', 'v_n_y': 'int_phi_r abs',
    #     'xunits': 'ms', '-t': o_1d_data.get_par("tmerger_gw"),
    #     'mode': 0, 'norm_to_m': 0,
    #     'ls': '-.', 'color': 'red', 'lw': 0.4,
    #     'label': None, 'ylabel': r'$C_m/C_0$ Magnitude', 'xlabel': r'time [ms]',
    #     'xmin': None, 'xmax': None, 'ymin': None, 'ymax': None,
    #     'xscale': None, 'yscale': 'log', 'legend': True,
    #     'fancyticks': True, 'minorticks': True
    # }
    densmode_m1 = {
        'task': 'line', 'dtype': 'dm', 'ptype': 'cartesian',
        'data': o_dm_data,
        'position': (1, 1),
        'v_n_x': 'times', 'v_n_y': 'int_phi_r abs',
        'xunits': 'ms', '-t': o_1d_data.get_par("tmerger_gw"),
        'mode': 1, 'norm_to_m': 0,
        'ls': '-', 'color': 'green', 'lw': 1.,
        'label': 'DD2 m=1', 'ylabel': r'$C_m/C_0$ Magnitude', 'xlabel': r'time [ms]',
        'xmin': None, 'xmax': None, 'ymin': None, 'ymax': None,
        'xscale': None, 'yscale': 'log', 'legend': False,
        'fancyticks': True, 'minorticks': True,
    }
    densmode_m2 = {
        'task': 'line', 'dtype': 'dm', 'ptype': 'cartesian',
        'data': o_dm_data,
        'position': (1, 1),
        'v_n_x': 'times', 'v_n_y': 'int_phi_r abs',
        'xunits': 'ms', '-t': o_1d_data.get_par("tmerger_gw"),
        'mode': 2, 'norm_to_m': 0,
        'ls': ':', 'color': 'green', 'lw': 0.6,
        'label': 'DD2 m=2', 'ylabel': r'$C_m/C_0$ Magnitude', 'xlabel': r'time [ms]',
        'xmin': None, 'xmax': None, 'ymin': None, 'ymax': None,
        'xscale': None, 'yscale': 'log', 'legend': False,
        'fancyticks': True, 'minorticks': True
    }
    # o_plot.set_plot_dics.append(densmode_m0)
    o_plot.set_plot_dics.append(densmode_m1)
    o_plot.set_plot_dics.append(densmode_m2)

    " --- --- --- --- --- --- ---- --- --- --- --- --- -----"
    # ------------------------------------------------------#
    " --- --- --- --- --- SIM0 --- --- --- ---  --- --- --- "

    sim = "LS220_M13641364_M0_SR"
    o_1d_data = COMPUTE_STORE_PAR(sim)

    o_dm_data = LOAD_DENSITY_MODES(sim)
    o_dm_data.gen_set['fname'] = Paths.ppr_sims + sim + '/res_3d/' + "density_modes_lap15.h5"
    # densmode_m0 = {
    #     'task': 'line', 'dtype': 'dm', 'ptype': 'cartesian',
    #     'data': o_dm_data,
    #     'position': (1, 1),
    #     'v_n_x': 'times', 'v_n_y': 'int_phi_r abs',
    #     'xunits': 'ms', '-t': o_1d_data.get_par("tmerger_gw"),
    #     'mode': 0, 'norm_to_m': 0,
    #     'ls': '-.', 'color': 'blue', 'lw': 0.4,
    #     'label': None, 'ylabel': r'$C_m/C_0$ Magnitude', 'xlabel': r'time [ms]',
    #     'xmin': None, 'xmax': None, 'ymin': None, 'ymax': None,
    #     'xscale': None, 'yscale': 'log', 'legend': True,
    #     'fancyticks': True, 'minorticks': True
    # }
    densmode_m1 = {
        'task': 'line', 'dtype': 'dm', 'ptype': 'cartesian',
        'data': o_dm_data,
        'position': (1, 1),
        'v_n_x': 'times', 'v_n_y': 'int_phi_r abs',
        'xunits': 'ms', '-t': o_1d_data.get_par("tmerger_gw"),
        'mode': 1, 'norm_to_m': 0,
        'ls': '-', 'color': 'blue', 'lw': 1.,
        'label': 'LS220', 'ylabel': r'$C_m/C_0$ Magnitude', 'xlabel': r'time [ms]',
        'xmin': None, 'xmax': None, 'ymin': None, 'ymax': None,
        'xscale': None, 'yscale': 'log', 'legend': False,
        'fancyticks': True, 'minorticks': True
    }
    densmode_m2 = {
        'task': 'line', 'dtype': 'dm', 'ptype': 'cartesian',
        'data': o_dm_data,
        'position': (1, 1),
        'v_n_x': 'times', 'v_n_y': 'int_phi_r abs',
        'xunits': 'ms', '-t': o_1d_data.get_par("tmerger_gw"),
        'mode': 2, 'norm_to_m': 0,
        'ls': ':', 'color': 'blue', 'lw': 0.6,
        'label': None, 'ylabel': r'$C_m/C_0$ Magnitude', 'xlabel': r'time [ms]',
        'xmin': None, 'xmax': None, 'ymin': None, 'ymax': None,
        'xscale': None, 'yscale': 'log', 'legend': False,
        'fancyticks': True, 'minorticks': True
    }
    coll_time_vert = {
        'task': 'vertline', 'dtype': '-', 'ptype': 'cartesian',
        'value': (o_1d_data.get_par("tcoll_gw") - o_1d_data.get_par("tmerger_gw")) * 1e3,
        'position': (1, 1),
        'ls': '-.', 'color': 'blue', 'lw': 0.5,
    }
    o_plot.set_plot_dics.append(coll_time_vert)
    # o_plot.set_plot_dics.append(densmode_m0)
    o_plot.set_plot_dics.append(densmode_m1)
    o_plot.set_plot_dics.append(densmode_m2)

    " --- --- --- --- --- --- ---- --- --- --- --- --- -----"
    # ------------------------------------------------------#
    " --- --- --- --- --- SIM0 --- --- --- ---  --- --- --- "

    sim = "SLy4_M13641364_M0_SR"
    o_1d_data = COMPUTE_STORE_PAR(sim)

    o_dm_data = LOAD_DENSITY_MODES(sim)
    o_dm_data.gen_set['fname'] = Paths.ppr_sims + sim + '/res_3d/' + "density_modes_lap15.h5"
    # densmode_m0 = {
    #     'task': 'line', 'dtype': 'dm', 'ptype': 'cartesian',
    #     'data': o_dm_data,
    #     'position': (1, 1),
    #     'v_n_x': 'times', 'v_n_y': 'int_phi_r abs',
    #     'xunits': 'ms', '-t': o_1d_data.get_par("tmerger_gw"),
    #     'mode': 0, 'norm_to_m': 0,
    #     'ls': '-.', 'color': 'red', 'lw': 0.4,
    #     'label': None, 'ylabel': r'$C_m/C_0$ Magnitude', 'xlabel': r'time [ms]',
    #     'xmin': None, 'xmax': None, 'ymin': None, 'ymax': None,
    #     'xscale': None, 'yscale': 'log', 'legend': True,
    #     'fancyticks': True, 'minorticks': True
    # }
    densmode_m1 = {
        'task': 'line', 'dtype': 'dm', 'ptype': 'cartesian',
        'data': o_dm_data,
        'position': (1, 1),
        'v_n_x': 'times', 'v_n_y': 'int_phi_r abs',
        'xunits': 'ms', '-t': o_1d_data.get_par("tmerger_gw"),
        'mode': 1, 'norm_to_m': 0,
        'ls': '-', 'color': 'red', 'lw': 1.,
        'label': 'SLy4', 'ylabel': r'$C_m/C_0$ Magnitude', 'xlabel': r'time [ms]',
        'xmin': None, 'xmax': None, 'ymin': None, 'ymax': None,
        'xscale': None, 'yscale': 'log', 'legend': False,
        'fancyticks': True, 'minorticks': True
    }
    densmode_m2 = {
        'task': 'line', 'dtype': 'dm', 'ptype': 'cartesian',
        'data': o_dm_data,
        'position': (1, 1),
        'v_n_x': 'times', 'v_n_y': 'int_phi_r abs',
        'xunits': 'ms', '-t': o_1d_data.get_par("tmerger_gw"),
        'mode': 2, 'norm_to_m': 0,
        'ls': ':', 'color': 'red', 'lw': 0.6,
        'label': None, 'ylabel': r'$C_m/C_0$ Magnitude', 'xlabel': r'time [ms]',
        'xmin': None, 'xmax': None, 'ymin': None, 'ymax': None,
        'xscale': None, 'yscale': 'log', 'legend': True,
        'fancyticks': True, 'minorticks': True
    }
    coll_time_vert = {
        'task': 'vertline', 'dtype': '-', 'ptype': 'cartesian',
        'value': (o_1d_data.get_par("tcoll_gw") - o_1d_data.get_par("tmerger_gw")) * 1e3,
        'position': (1, 1),
        'ls': '-.', 'color': 'red', 'lw': 0.5,
    }
    o_plot.set_plot_dics.append(coll_time_vert)
    # o_plot.set_plot_dics.append(densmode_m0)
    o_plot.set_plot_dics.append(densmode_m1)
    o_plot.set_plot_dics.append(densmode_m2)

    o_plot.main()
# plot_m1_in_time_3sims_tcoll()

def plot_mkn_3sims():


    o_plot = PLOT_MANY_TASKS()
    o_plot.gen_set["figdir"] = Paths.plots + '/combine_test/'
    o_plot.gen_set["type"] = "cartesian"
    o_plot.gen_set["figsize"] = (9.0, 2.7)  # <->, |]
    o_plot.gen_set["figname"] = "mkn_final.pdf"
    o_plot.gen_set["sharex"] = False
    o_plot.gen_set["sharey"] = False
    o_plot.gen_set["invert_y"] = False
    o_plot.gen_set["subplots_adjust_h"] = 0.3
    o_plot.gen_set["subplots_adjust_w"] = 0.0
    o_plot.set_plot_dics = []


    """-------------------------------------BAND1--------------------------------"""

    """ --- --- ---- --- --- --- SIM 1 --- -- --- --- --- --- """

    sim = "DD2_M13641364_M0_SR"
    o_lc_data = EXTRACT_LIGHTCURVE(sim)
    sim_lc = {
        'task': 'mkn model', "ptype": "cartesian",
        'position': (1, 1),
        'data': o_lc_data, 'band': 'g', 'obs': True, 'fname': 'mkn_model2.h5',
        'v_n_x': 'time', 'v_n_y': 'mag',
        'color': 'green', 'alpha': 0.5,
        'ymin': 25, 'ymax': 15, 'xmin': 3e-1, 'xmax': 3e1,
        'xlabel': r"time [days]", 'ylabel': r"AB magnitude at 40 Mpc",
        'label': 'DD2', 'xscale': 'log',
        'fancyticks': True, 'minorticks': True,
    }
    o_plot.set_plot_dics.append(sim_lc)

    """ --- --- ---- --- --- --- SIM 2 --- -- --- --- --- --- """

    sim = "LS220_M13641364_M0_SR"
    o_lc_data = EXTRACT_LIGHTCURVE(sim)
    sim_lc = {
        'task': 'mkn model', "ptype": "cartesian",
        'position': (1, 1),
        'data': o_lc_data, 'band': 'g', 'obs': True, 'fname': 'mkn_model2.h5',
        'v_n_x': 'time', 'v_n_y': 'mag',
        'color': 'blue', 'alpha': 0.5,
        'ymin': 25, 'ymax': 15,  'xmin': 3e-1, 'xmax': 3e1,
        'xlabel': r"time [days]", 'ylabel': r"AB magnitude at 40 Mpc",
        'label': 'LS220', 'xscale': 'log',
        'fancyticks': True, 'minorticks': True,
    }
    o_plot.set_plot_dics.append(sim_lc)

    """ --- --- ---- --- --- --- SIM 3 --- -- --- --- --- --- """

    sim = "SLy4_M13641364_M0_SR"
    o_lc_data = EXTRACT_LIGHTCURVE(sim)
    sim_lc = {
        'task': 'mkn model', "ptype": "cartesian",
        'position': (1, 1),
        'data': o_lc_data, 'band': 'g', 'obs': True, 'fname': 'mkn_model2.h5',
        'v_n_x': 'time', 'v_n_y': 'mag',
        'color': 'red', 'alpha': 0.5,
        'ymin': 25, 'ymax': 15,  'xmin': 3e-1, 'xmax': 3e1,
        'xlabel': r"time [days]", 'ylabel': r"AB magnitude at 40 Mpc",
        'label': 'SLy4', 'xscale': 'log',
        'fancyticks': True, 'minorticks': True,
    }
    o_plot.set_plot_dics.append(sim_lc)

    sim_lc_obs = {
        'task': 'mkn obs', "ptype": "cartesian",
        'position': (1, 1),
        'data': o_lc_data, 'band': 'g', 'obs': True,
        'v_n_x': 'time', 'v_n_y': 'mag',
        'color': 'gray', 'marker': 'o', 'ms': 5., 'alpha': 0.6,
        'ymin': 25, 'ymax': 15,  'xmin': 3e-1, 'xmax': 3e1,
        'xlabel': r"time [days]", 'ylabel': r"AB magnitude at 40 Mpc",
        'label': None, 'xscale': 'log',
        'fancyticks': True, 'minorticks': True,
        'legend': True, 'title': 'g band',
    }
    o_plot.set_plot_dics.append(sim_lc_obs)

    """-------------------------------------BAND2--------------------------------"""

    """ --- --- ---- --- --- --- SIM 1 --- -- --- --- --- --- """

    sim = "DD2_M13641364_M0_SR"
    o_lc_data = EXTRACT_LIGHTCURVE(sim)
    sim_lc = {
        'task': 'mkn model', "ptype": "cartesian",
        'position': (1, 2),
        'data': o_lc_data, 'band': 'z', 'obs': True, 'fname': 'mkn_model2.h5',
        'v_n_x': 'time', 'v_n_y': 'mag',
        'color': 'green', 'alpha': 0.5,
        'ymin': 25, 'ymax': 15,  'xmin': 3e-1, 'xmax': 3e1,
        'xlabel': r"time [days]", 'ylabel': r"AB magnitude at 40 Mpc",
        'label': None, 'xscale': 'log',
        'fancyticks': True, 'minorticks': True,
    }
    o_plot.set_plot_dics.append(sim_lc)

    """ --- --- ---- --- --- --- SIM 2 --- -- --- --- --- --- """

    sim = "LS220_M13641364_M0_SR"
    o_lc_data = EXTRACT_LIGHTCURVE(sim)
    sim_lc = {
        'task': 'mkn model', "ptype": "cartesian",
        'position': (1, 2),
        'data': o_lc_data, 'band': 'z', 'obs': True, 'fname': 'mkn_model2.h5',
        'v_n_x': 'time', 'v_n_y': 'mag',
        'color': 'blue', 'alpha': 0.5,
        'ymin': 25, 'ymax': 15, 'xmin': 3e-1, 'xmax': 3e1,
        'xlabel': r"time [days]", 'ylabel': r"AB magnitude at 40 Mpc",
        'label': None, 'xscale': 'log',
        'fancyticks': True, 'minorticks': True,
    }
    o_plot.set_plot_dics.append(sim_lc)

    """ --- --- ---- --- --- --- SIM 3 --- -- --- --- --- --- """

    sim = "SLy4_M13641364_M0_SR"
    o_lc_data = EXTRACT_LIGHTCURVE(sim)
    sim_lc = {
        'task': 'mkn model', "ptype": "cartesian",
        'position': (1, 2),
        'data': o_lc_data, 'band': 'z', 'obs': True, 'fname': 'mkn_model2.h5',
        'v_n_x': 'time', 'v_n_y': 'mag',
        'color': 'red', 'alpha': 0.5,
        'ymin': 25, 'ymax': 15, 'xmin': 3e-1, 'xmax': 3e1,
        'xlabel': r"time [days]", 'ylabel': r"AB magnitude at 40 Mpc",
        'label': None, 'xscale': 'log',
        'fancyticks': True, 'minorticks': True,
        'sharey': True
    }
    o_plot.set_plot_dics.append(sim_lc)

    sim_lc_obs = {
        'task': 'mkn obs', "ptype": "cartesian",
        'position': (1, 2),
        'data': o_lc_data, 'band': 'z', 'obs': True,
        'v_n_x': 'time', 'v_n_y': 'mag',
        'color': 'gray', 'marker': 'o', 'ms': 5., 'alpha': 0.6,
        'ymin': 25, 'ymax': 15, 'xmin': 3e-1, 'xmax': 3e1,
        'xlabel': r"time [days]", 'ylabel': r"AB magnitude at 40 Mpc",
        'label': None, 'xscale': 'log',
        'fancyticks': True, 'minorticks': True,
        'legend': True, 'title': 'z band',
        'sharey': True
    }
    o_plot.set_plot_dics.append(sim_lc_obs)

    """-------------------------------------BAND3--------------------------------"""

    """ --- --- ---- --- --- --- SIM 1 --- -- --- --- --- --- """

    sim = "DD2_M13641364_M0_SR"
    o_lc_data = EXTRACT_LIGHTCURVE(sim)
    sim_lc = {
        'task': 'mkn model', "ptype": "cartesian",
        'position': (1, 3),
        'data': o_lc_data, 'band': 'Ks', 'obs': True, 'fname': 'mkn_model2.h5',
        'v_n_x': 'time', 'v_n_y': 'mag',
        'color': 'green', 'alpha': 0.5,
        'ymin': 25, 'ymax': 15, 'xmin': 3e-1, 'xmax': 3e1,
        'xlabel': r"time [days]", 'ylabel': r"AB magnitude at 40 Mpc",
        'label': None, 'xscale': 'log',
        'fancyticks': True, 'minorticks': True,
    }
    o_plot.set_plot_dics.append(sim_lc)

    """ --- --- ---- --- --- --- SIM 2 --- -- --- --- --- --- """

    sim = "LS220_M13641364_M0_SR"
    o_lc_data = EXTRACT_LIGHTCURVE(sim)
    sim_lc = {
        'task': 'mkn model', "ptype": "cartesian",
        'position': (1, 3),
        'data': o_lc_data, 'band': 'Ks', 'obs': True, 'fname':'mkn_model2.h5',
        'v_n_x': 'time', 'v_n_y': 'mag',
        'color': 'blue', 'alpha': 0.5,
        'ymin': 25, 'ymax': 15, 'xmin': 3e-1, 'xmax': 3e1,
        'xlabel': r"time [days]", 'ylabel': r"AB magnitude at 40 Mpc",
        'label': None, 'xscale': 'log',
        'fancyticks': True, 'minorticks': True,
    }
    o_plot.set_plot_dics.append(sim_lc)

    """ --- --- ---- --- --- --- SIM 3 --- -- --- --- --- --- """

    sim = "SLy4_M13641364_M0_SR"
    o_lc_data = EXTRACT_LIGHTCURVE(sim)
    sim_lc = {
        'task': 'mkn model', "ptype": "cartesian",
        'position': (1, 3),
        'data': o_lc_data, 'band': 'Ks', 'obs': True, 'fname': 'mkn_model2.h5',
        'v_n_x': 'time', 'v_n_y': 'mag',
        'color': 'red', 'alpha': 0.5,
        'ymin': 25, 'ymax': 15, 'xmin': 3e-1, 'xmax': 3e1,
        'xlabel': r"time [days]", 'ylabel': r"AB magnitude at 40 Mpc",
        'label': None, 'xscale': 'log',
        'fancyticks': True, 'minorticks': True,
        'sharey': True
    }
    o_plot.set_plot_dics.append(sim_lc)

    sim_lc_obs = {
        'task': 'mkn obs', "ptype": "cartesian",
        'position': (1, 3),
        'data': o_lc_data, 'band': 'Ks', 'obs': True,
        'v_n_x': 'time', 'v_n_y': 'mag',
        'color': 'gray', 'marker': 'o', 'ms': 5., 'alpha': 0.6,
        'ymin': 25, 'ymax': 15, 'xmin': 3e-1, 'xmax': 3e1,
        'xlabel': r"time [days]", 'ylabel': r"AB magnitude at 40 Mpc",
        'label': None, 'xscale': 'log',
        'fancyticks': True, 'minorticks': True,
        'legend': True, 'title': 'Ks band',
        'sharey': True
    }
    o_plot.set_plot_dics.append(sim_lc_obs)

    o_plot.main()
    exit(1)
# plot_mkn_3sims()

def plot_nucleo_3sim_one_plot_combined_yields():
    o_plot = PLOT_MANY_TASKS()
    o_plot.gen_set["figdir"] = Paths.plots + '/combine_test/'
    o_plot.gen_set["type"] = "cartesian"
    o_plot.gen_set["figsize"] = (4.6, 3.6)  # <->, |]
    o_plot.gen_set["figname"] = "nucleo_final.pdf"
    o_plot.gen_set["sharex"] = False
    o_plot.gen_set["sharey"] = False
    o_plot.gen_set["subplots_adjust_h"] = 0.3
    o_plot.gen_set["subplots_adjust_w"] = 0.0
    o_plot.set_plot_dics = []

    ''' --- --- -- SIM 1 --- --- --- '''

    sim = "DD2_M13641364_M0_SR"
    o_data = NORMALIZE_NUCLEO(sim)

    sim_nucleo = {
        'task': 'nucleo', 'ptype': 'cartesian',
        'position': (1, 1),
        'data': o_data, 'criterion': '_0', 'method': 'Asol=195',
        'v_n_x': 'A', 'v_n_y': 'Y_final',
        'color': 'green', 'ls': ':', 'lw': 0.6, 'ds': 'steps', 'alpha': 0.6,
        'ymin': 1e-5, 'ymax': 2e-1, 'xmin': 50, 'xmax': 210,
        'xlabel': r"A", 'ylabel': r'Relative final abundances',
        'label': "DD2 Geodesic", 'yscale': 'log',
        'fancyticks': True, 'minorticks': True
    }
    o_plot.set_plot_dics.append(sim_nucleo)

    sim_nucleo = {
        'task': 'nucleo', 'ptype': 'cartesian',
        'position': (1, 1),
        'data': o_data, 'criterion': '_0_b_w', 'method': 'A=195Asun',
        'v_n_x': 'A', 'v_n_y': 'Y_final',
        'color': 'green', 'ls': '--', 'lw': 0.6, 'ds': 'steps', 'alpha': 0.6,
        'ymin': 1e-5, 'ymax': 2e-1, 'xmin': 50, 'xmax': 210,
        'xlabel': r"A", 'ylabel': r'Relative final abundances',
        'label': "DD2", 'yscale': 'log',
        'fancyticks': True, 'minorticks': True
    }
    # o_plot.set_plot_dics.append(sim_nucleo)

    sim_nucleo_b = {
        'task': 'nucleo', 'ptype': 'cartesian',
        'position': (1, 1),
        'data': o_data, 'criterion': '_0 _0_b_w', 'method': 'Asol=195',
        'v_n_x': 'A', 'v_n_y': 'Y_final',
        'color': 'green', 'ls': '-', 'lw': 1.0, 'ds': 'steps',
        'ymin': 1e-5, 'ymax': 2e-1, 'xmin': 50, 'xmax': 210,
        'xlabel': r"A", 'ylabel': r'Relative final abundances',
        'label': "DD2 combined", 'yscale': 'log',
        'fancyticks': True, 'minorticks': True,
        'legend': True
    }
    o_plot.set_plot_dics.append(sim_nucleo_b)

    sol_yeilds = {
        'task': 'nucleo', 'ptype': 'cartesian',
        'position': (1, 1),
        'data': o_data, 'criterion': '_0', 'method': 'sum',
        'v_n_x': 'Asun', 'v_n_y': 'Ysun',
        'color': 'gray', 'marker': 'o', 'ms': 4, 'alpha':0.4,
        'ymin': 1e-5, 'ymax': 2e-1, 'xmin': 50, 'xmax': 210,
        'xlabel': r"A", 'ylabel': r'Relative final abundances',
        'label': None, 'yscale': 'log',
        'fancyticks': True, 'minorticks': True
    }
    o_plot.set_plot_dics.append(sol_yeilds)


    ''' --- --- -- SIM 2 --- --- --- '''

    sim = "LS220_M13641364_M0_SR"
    o_data = NORMALIZE_NUCLEO(sim)

    sim_nucleo = {
        'task': 'nucleo', 'ptype': 'cartesian',
        'position': (1, 1),
        'data': o_data, 'criterion': '_0', 'method': 'Asol=195',
        'v_n_x': 'A', 'v_n_y': 'Y_final',
        'color': 'blue', 'ls': ':', 'lw': 0.6, 'ds': 'steps', 'alpha': 0.6,
        'ymin': 1e-5, 'ymax': 2e-1, 'xmin': 50, 'xmax': 210,
        'xlabel': r"A", 'ylabel': r'Relative final abundances',
        'label': None, 'yscale': 'log',
        'fancyticks': True, 'minorticks': True
    }

    sim_nucleo_b = {
        'task': 'nucleo', 'ptype': 'cartesian',
        'position': (1, 1),
        'data': o_data, 'criterion': '_0 _0_b_w', 'method': 'Asol=195',
        'v_n_x': 'A', 'v_n_y': 'Y_final',
        'color': 'blue', 'ls': '-', 'lw': 1.0, 'ds': 'steps',
        'ymin': 1e-5, 'ymax': 2e-1, 'xmin': 50, 'xmax': 210,
        'xlabel': r"A", 'ylabel': r'Relative final abundances',
        'label': "LS220", 'yscale': 'log',
        'fancyticks': True, 'minorticks': True
    }

    o_plot.set_plot_dics.append(sim_nucleo)
    o_plot.set_plot_dics.append(sim_nucleo_b)

    ''' --- --- -- SIM 3 --- --- --- '''

    sim = "SLy4_M13641364_M0_SR"
    o_data = NORMALIZE_NUCLEO(sim)

    sim_nucleo = {
        'task': 'nucleo', 'ptype': 'cartesian',
        'position': (1, 1),
        'data': o_data, 'criterion': '_0', 'method': 'Asol=195',
        'v_n_x': 'A', 'v_n_y': 'Y_final',
        'color': 'red', 'ls': ':', 'lw': 0.6, 'ds': 'steps', 'alpha': 0.6,
        'ymin': 1e-5, 'ymax': 2e-1, 'xmin': 50, 'xmax': 210,
        'xlabel': r"A", 'ylabel': r'Relative final abundances',
        'label': None, 'yscale': 'log',
        'fancyticks': True, 'minorticks': True
    }

    sim_nucleo_b = {
        'task': 'nucleo', 'ptype': 'cartesian',
        'position': (1, 1),
        'data': o_data, 'criterion': '_0 _0_b_w', 'method': 'Asol=195',
        'v_n_x': 'A', 'v_n_y': 'Y_final',
        'color': 'red', 'ls': '-', 'lw': 1.0, 'ds': 'steps',
        'ymin': 1e-5, 'ymax': 2e-1, 'xmin': 50, 'xmax': 210,
        'xlabel': r"A", 'ylabel': r'Relative final abundances',
        'label': "SLy4", 'yscale': 'log',
        'fancyticks': True, 'minorticks': True,
        'legend': True
    }

    o_plot.set_plot_dics.append(sim_nucleo)
    o_plot.set_plot_dics.append(sim_nucleo_b)

    o_plot.main()
# plot_nucleo_3sim_one_plot_combined_yields()


""" --- --- MAPS --- --- --- """

def plot_cartesian_edensity_jflux_for_1sim():


    o_plot = PLOT_MANY_TASKS()
    o_plot.gen_set["figdir"] = Paths.plots + '/combine_test/'
    o_plot.gen_set["type"] = "polar"
    o_plot.gen_set["figsize"] = (6.0, 2.7)  # <->, |] # to match hists with (8.5, 2.7)
    o_plot.gen_set["figname"] = "jflux_final.png"
    o_plot.gen_set["sharex"] = False
    o_plot.gen_set["sharey"] = False
    o_plot.gen_set["subplots_adjust_h"] = 0.3
    o_plot.gen_set["subplots_adjust_w"] = 0.3
    o_plot.set_plot_dics = []

    sim = "DD2_M13641364_M0_SR"
    o_int_data = ADD_METHODS_FOR_INT_DATA(sim)
    o_int_data.flag_force_unique_grid = True

    int_ang_mom_flux_dic = {
        'task': '2d projection', 'dtype': 'int', 'ptype': 'polar',
        'data': o_int_data, 'it': 751278,
        'position': (1, 1), 'title': r'$r\times\int \dot{J} dz $ [geo]', 'cbar':  'left .15 .0', 'cbar fmt': '%.1e',
        'cbar label': '',
        'v_n_x': 'phi_cyl', 'v_n_y': 'r_cyl', 'v_n': 'ang_mom_flux', 'mod': 'integ_over_z fill_phi *r',
        'phimin': None, 'phimax': None, 'rmin': 0, 'rmax': 50, 'vmin': -5e-5, 'vmax': 5e-5,
        'fill_vmin': False,  # fills the x < vmin with vmin
        'xscale': None, 'yscale': None,
        'mask': None, 'cmap': 'RdBu_r', 'norm': "linear",
        'fancyticks': False,
        'sharex': True # removes angular citks
    }
    o_plot.set_plot_dics.append(int_ang_mom_flux_dic)

    int_density_dic = {
        'task': '2d projection', 'dtype': 'int', 'ptype': 'polar',
        'data': o_int_data, 'it': 751278,
        'position': (1, 2), 'title': r'$(D - \sum D / N_D)(r)$, $D = \int W\rho\sqrt{g} dz $ [geo]', 'cbar': 'right .05 .0', 'cbar fmt': '%.1e',
        'cbar label': '',
        'v_n_x': 'phi_cyl', 'v_n_y': 'r_cyl', 'v_n': 'density', 'mod': 'integ_over_z fill_phi -ave(r)',
        'phimin': None, 'phimax': None, 'rmin': 0, 'rmax': 50, 'vmin': -5e-6, 'vmax': 5e-6,
        'fill_vmin': False,  # fills the x < vmin with vmin
        'xscale': None, 'yscale': None,
        'mask': None, 'cmap': 'RdBu_r', 'norm': "linear",
        'fancyticks': False,
        'sharex': True  # removes angular citks
    }
    o_plot.set_plot_dics.append(int_density_dic)

    o_plot.main()
    exit(0)

    ''' --- SIM 1 --- '''

    sim = "DD2_M13641364_M0_SR"

    o_int_data = ADD_METHODS_FOR_INT_DATA(sim)
    o_int_data.flag_force_unique_grid = True
    # o_int_data.it_for_unique_grid = 2254356 # grid is loaded only for this iteration and assumed to be constant


    int_ang_mom_flux_dic = {
        'task': '2d projection', 'dtype': 'int', 'ptype': 'polar',
        'data': o_int_data, 'it': 751278,
        'position': (1, 1), 'title': 'DD2', 'cbar': None,#'right .05 .0',
        'v_n_x': 'phi_cyl', 'v_n_y': 'r_cyl', 'v_n': 'ang_mom_flux',  'mod': 'integ_over_z',
        'phimin': None, 'phimax': None, 'rmin': 0, 'rmax': 50, 'vmin': -1e-6, 'vmax': 1e-6,
        'fill_vmin': False, # fills the x < vmin with vmin
        'xscale': None, 'yscale': None,
        'mask': None, 'cmap': 'RdBu_r', 'norm': "linear",
        'fancyticks': False
    }
    int_ang_mom_flux_dic_rev = {
        'task': '2d projection', 'dtype': 'int', 'ptype': 'polar',
        'data': o_int_data, 'it': 751278,
        'position': (1, 1), 'title': 'DD2', 'cbar': None,#'right .05 .0','cbar fmt':'%.1e',
        'v_n_x': 'phi_cyl', 'v_n_y': 'r_cyl', 'v_n': 'ang_mom_flux', 'mod': 'integ_over_z',
        'phimin': None, 'phimax': None, 'rmin': 0, 'rmax': 50, 'vmin': -5e-6, 'vmax': 5e-6,
        'fill_vmin': False, # fills the x < vmin with vmin
        'xscale': None, 'yscale': None,
        'mask': None, 'cmap': 'RdBu_r', 'norm': "linear",
        'fancyticks': False
    }

    o_dm_data = LOAD_DENSITY_MODES(sim)
    densmode = {
        'task': '2d projection', 'dtype': 'dm', 'ptype': 'polar',
        'data': o_dm_data, 'it': 751278,
        'name': 'densmode', 'position': (1, 1),
        'v_n_x': 'phi_cyl', 'v_n_y': 'r_cyl', 'v_n': 'int_phi',
        'mode': 1, #'int': 'spline',
        'rmax': 50, 'ls': '-', 'color': 'black',
    }
    densmode2 = copy.deepcopy(densmode)
    densmode2["v_n"] = 'int_phi_r'
    densmode2["ls"] = '--'

    o_plot.set_plot_dics.append(int_ang_mom_flux_dic)
    # o_plot.set_plot_dics.append(int_ang_mom_flux_dic_rev)
    o_plot.set_plot_dics.append(densmode)
    o_plot.set_plot_dics.append(densmode2)
    ''' --- SIM 2 --- '''

    sim = "LS220_M13641364_M0_SR"

    o_int_data = ADD_METHODS_FOR_INT_DATA(sim)
    o_int_data.flag_force_unique_grid = True
    # o_int_data.it_for_unique_grid = 2254356  # grid is loaded only for this iteration and assumed to be constant

    int_ang_mom_flux_dic = {
        'task': '2d projection', 'dtype': 'int', 'ptype': 'polar',
        'data': o_int_data, 'it': 753664,
        'position': (1, 2), 'title': 'LS220', 'cbar': None, #'right .05 .0',  'cbar fmt':'%.1e'
        'v_n_x': 'phi_cyl', 'v_n_y': 'r_cyl', 'v_n': 'ang_mom_flux', 'mod': 'integ_over_z',
        'phimin': None, 'phimax': None, 'rmin': 0, 'rmax': 50, 'vmin': -1e-6, 'vmax': 1e-6,
        'fill_vmin': False,  # fills the x < vmin with vmin
        'xscale': None, 'yscale': None,
        'mask': None, 'cmap': 'RdBu_r', 'norm': "linear",
        'fancyticks': False
    }
    int_ang_mom_flux_dic_rev = {
        'task': '2d projection', 'dtype': 'int', 'ptype': 'polar',
        'data': o_int_data, 'it': 753664,
        'position': (1, 2), 'title': 'LS220', 'cbar': None,#'right .05 .0',
        'v_n_x': 'phi_cyl', 'v_n_y': 'r_cyl', 'v_n': 'ang_mom_flux',
        'phimin': None, 'phimax': None, 'rmin': 0, 'rmax': 50, 'vmin': 1e-9, 'vmax': 1e-5,
        'fill_vmin': False,  # fills the x < vmin with vmin
        'xscale': None, 'yscale': None,
        'mask': 'positive', 'cmap': 'Blues', 'norm': "log", 'todo': None,
        'fancyticks': False
    }

    o_dm_data = LOAD_DENSITY_MODES(sim)

    densmode = {
        'task': '2d projection', 'dtype': 'dm', 'ptype': 'polar',
        'data': o_dm_data, 'it': 753664,
        'name': 'densmode', 'position': (1, 2),
        'v_n_x': 'phi_cyl', 'v_n_y': 'r_cyl', 'v_n': 'int_phi',
        'mode': 1,
        'rmax': 50, 'ls': '-', 'color': 'black',
    }
    densmode2 = copy.deepcopy(densmode)
    densmode2["v_n"] = 'int_phi_r'
    densmode2["ls"] = '--'
    # o_plot.set_plot_dics.append(int_ang_mom_flux_dic_rev)
    o_plot.set_plot_dics.append(int_ang_mom_flux_dic)
    # o_plot.set_plot_dics.append(int_ang_mom_flux_dic_rev)
    o_plot.set_plot_dics.append(densmode)
    o_plot.set_plot_dics.append(densmode2)

    ''' --- SIM 3 --- '''

    sim = "SLy4_M13641364_M0_SR"

    o_int_data = ADD_METHODS_FOR_INT_DATA(sim)
    o_int_data.flag_force_unique_grid = True
    # o_int_data.it_for_unique_grid = 2254356  # grid is loaded only for this iteration and assumed to be constant

    int_ang_mom_flux_dic = {
        'task': '2d projection', 'dtype': 'int', 'ptype': 'polar',
        'data': o_int_data, 'it': 761856,
        'position': (1, 3), 'title': 'SLy4', 'cbar': 'right .03 .0','cbar fmt':'%.1e',
        'v_n_x': 'phi_cyl', 'v_n_y': 'r_cyl', 'v_n': 'ang_mom_flux', 'mod': 'integ_over_z', 'cbar label':r"$\int(\dot{J}dz)$",
        'phimin': None, 'phimax': None, 'rmin': 0, 'rmax': 50, 'vmin': -1e-6, 'vmax': 1e-6,
        'fill_vmin': True,  # fills the x < vmin with vmin
        'xscale': None, 'yscale': None,
         'mask': None, 'cmap': 'RdBu_r', 'norm': "linear",
        'fancyticks': False
    }
    int_ang_mom_flux_dic_rev = {
        'task': '2d projection', 'dtype': 'int', 'ptype': 'polar',
        'data': o_int_data, 'it': 761856,
        'position': (1, 3), 'title': 'SLy4', 'cbar': None,
        'v_n_x': 'phi_cyl', 'v_n_y': 'r_cyl', 'v_n': 'ang_mom_flux',
        'phimin': None, 'phimax': None, 'rmin': 0, 'rmax': 50, 'vmin': 1e-9, 'vmax': 1e-5,
        'fill_vmin': True,  # fills the x < vmin with vmin
        'xscale': None, 'yscale': None,
        'mask': 'positive', 'cmap': 'Blues', 'norm': "log", 'todo': None,
        'fancyticks': False
    }

    o_dm_data = LOAD_DENSITY_MODES(sim)

    densmode = {
        'task': '2d projection', 'dtype': 'dm', 'ptype': 'polar',
        'data': o_dm_data, 'it': 761856,
        'name': 'densmode', 'position': (1, 3),
        'v_n_x': 'phi_cyl', 'v_n_y': 'r_cyl', 'v_n': 'int_phi',
        'mode': 1,
        'rmax': 50, 'ls': '-', 'color': 'black',
    }
    densmode2 = copy.deepcopy(densmode)
    densmode2["v_n"] = 'int_phi_r'
    densmode2["ls"] = '--'

    o_plot.set_plot_dics.append(int_ang_mom_flux_dic)
    # o_plot.set_plot_dics.append(int_ang_mom_flux_dic_rev)
    o_plot.set_plot_dics.append(densmode)
    o_plot.set_plot_dics.append(densmode2)


    # o_plot.main()
plot_cartesian_edensity_jflux_for_1sim()


""" --- COLORFILLED PLOTS"""

def plot_one_jflux():

    sim = "DD2_M13641364_M0_SR"

    o_int_data = LOAD_INT_DATA(sim)
    o_int_data.flag_force_unique_grid = True
    o_int_data.it_for_unique_grid = 2254356 # grid is loaded only for this iteration and assumed to be constant

    o_plot = PLOT_MANY_TASKS()
    o_plot.gen_set["figdir"] = Paths.plots + '/combine_test/'
    o_plot.gen_set["type"] = "polar"
    o_plot.gen_set["figsize"] = (3.5, 3.5)  # <->, |]
    o_plot.gen_set["figname"] = "{0:07d}.png".format(int(669810))

    int_ang_mom_flux_dic = {
        'task': '2d projection', 'dtype': 'int', 'data': o_int_data, 'it': 669810,
        'position': (1, 1), 'title': 'time [ms]', 'cbar': 'right .05 .0',
        'ptype':'polar', 'v_n_x': 'phi_cyl', 'v_n_y': 'r_cyl', 'v_n': 'ang_mom_flux',
        'phimin': None, 'phimax': None, 'rmin': 0, 'rmax': 50, 'vmin': 1e-8, 'vmax': 1e-5,
        'fill_vmin': True, # fills the x < vmin with vmin
        'xscale': None, 'yscale': None,
        'mask_below': None, 'mask_above': None, 'cmap': 'inferno', 'norm': "log", 'todo': None,
        'fancyticks': False
    }
    o_plot.set_plot_dics = [int_ang_mom_flux_dic]
    o_plot.main()
# plot_one_jflux()

def plot_one_jflux_m1():

    sim = "DD2_M13641364_M0_SR"

    o_int_data = LOAD_INT_DATA(sim)
    o_int_data.flag_force_unique_grid = True
    o_int_data.it_for_unique_grid = 2254356 # grid is loaded only for this iteration and assumed to be constant

    o_plot = PLOT_MANY_TASKS()
    o_plot.gen_set["figdir"] = Paths.plots + '/combine_test/'
    o_plot.gen_set["type"] = "polar"
    o_plot.gen_set["figsize"] = (3.5, 3.5)  # <->, |]
    o_plot.gen_set["figname"] = "{0:07d}.png".format(int(669810))

    int_ang_mom_flux_dic = {
        'task': '2d projection', 'dtype': 'int', 'ptype': 'polar',
        'data': o_int_data, 'it': 669810,
        'position': (1, 1), 'title': 'time [ms]', 'cbar': 'right .05 .0',
        'v_n_x': 'phi_cyl', 'v_n_y': 'r_cyl', 'v_n': 'ang_mom_flux',
        'phimin': None, 'phimax': None, 'rmin': 0, 'rmax': 50, 'vmin': 1e-8, 'vmax': 1e-5,
        'fill_vmin': True, # fills the x < vmin with vmin
        'xscale': None, 'yscale': None,
        'mask_below': None, 'mask_above': None, 'cmap': 'inferno', 'norm': "log", 'todo': None,
        'fancyticks': False
    }

    o_dm_data = LOAD_DENSITY_MODES(sim)

    densmode = {
        'task': '2d projection', 'dtype': 'dm', 'ptype': 'polar',
        'data': o_dm_data, 'it': 669810,
        'name': 'densmode', 'position': (1, 1),
        'v_n_x': 'phi_cyl', 'v_n_y': 'r_cyl', 'v_n': 'int_phi',
        'mode': 1,
        'rmax': 50, 'ls': '.-', 'color': 'white',
    }
    densmode2 = densmode.copy()
    densmode2["v_n"] = 'int_phi_r'
    densmode2["ls"] = '--'

    o_plot.set_plot_dics = [int_ang_mom_flux_dic, densmode, densmode2]
    o_plot.main()
# plot_one_jflux_m1()


def plot_test_velocity_jflux():


    o_plot = PLOT_MANY_TASKS()
    o_plot.gen_set["figdir"] = Paths.plots + '/combine_test/'
    o_plot.gen_set["type"] = "polar"
    o_plot.gen_set["figsize"] = (6.9, 2.9)  # <->, |] # to match hists with (8.5, 2.7)
    o_plot.gen_set["figname"] = "test_velocity_jflux.png"
    o_plot.gen_set["sharex"] = False
    o_plot.gen_set["sharey"] = False
    o_plot.gen_set["subplots_adjust_h"] = 0.3
    o_plot.gen_set["subplots_adjust_w"] = 0.3
    o_plot.set_plot_dics = []


    sim = "LS220_M13641364_M0_SR"

    o_int_data = LOAD_INT_DATA(sim)
    o_int_data.flag_force_unique_grid = True
    # o_int_data.it_for_unique_grid = 2254356  # grid is loaded only for this iteration and assumed to be constant

    int_ang_mom_flux_dic = {
        'task': '2d projection', 'dtype': 'int', 'ptype': 'polar',
        'data': o_int_data, 'it': 753664,
        'position': (1, 1), 'title': None, 'cbar': 'left .15 .0',
        'v_n_x': 'phi_cyl', 'v_n_y': 'r_cyl', 'v_n': 'ang_mom_flux',
        'phimin': None, 'phimax': None, 'rmin': 0, 'rmax': 50, 'vmin': 1e-9, 'vmax': 1e-5,
        'fill_vmin': False,  # fills the x < vmin with vmin
        'xscale': None, 'yscale': None,
        'cmap': 'Reds', 'norm': "log", 'todo': None,
        'fancyticks': False
    }
    int_ang_mom_flux_dic_rev = {
        'task': '2d projection', 'dtype': 'int', 'ptype': 'polar',
        'data': o_int_data, 'it': 753664,
        'position': (1, 1), 'title': 'time [ms]', 'cbar': None,#'right .05 .0',
        'v_n_x': 'phi_cyl', 'v_n_y': 'r_cyl', 'v_n': 'ang_mom_flux',
        'phimin': None, 'phimax': None, 'rmin': 0, 'rmax': 50, 'vmin': 1e-9, 'vmax': 1e-5,
        'fill_vmin': False,  # fills the x < vmin with vmin
        'xscale': None, 'yscale': None,
        'mask': 'positive', 'cmap': 'Blues', 'norm': "log", 'todo': None,
        'fancyticks': False
    }

    o_plot.set_plot_dics.append(int_ang_mom_flux_dic_rev)
    o_plot.set_plot_dics.append(int_ang_mom_flux_dic)
    # o_plot.set_plot_dics.append(int_ang_mom_flux_dic_rev)

    int_vel_r_dic = {
        'task': '2d projection', 'dtype': 'int', 'ptype': 'polar',
        'data': o_int_data, 'it': 753664,
        'position': (1, 2), 'title': None, 'cbar': 'right .05 .0',
        'v_n_x': 'phi_cyl', 'v_n_y': 'r_cyl', 'v_n': 'vr',
        'phimin': None, 'phimax': None, 'rmin': 0, 'rmax': 50, 'vmin': 1e-5, 'vmax': 1e-1,
        'fill_vmin': False,  # fills the x < vmin with vmin
        'xscale': None, 'yscale': None,
        'mask': 'negative', 'cmap': 'Reds', 'norm': "log", 'todo': None,
        'fancyticks': False
    }
    int_vel_r_dic_rev = {
        'task': '2d projection', 'dtype': 'int', 'ptype': 'polar',
        'data': o_int_data, 'it': 753664,
        'position': (1, 2), 'title': None, 'cbar': 'right .05 .0',
        'v_n_x': 'phi_cyl', 'v_n_y': 'r_cyl', 'v_n': 'vr',
        'phimin': None, 'phimax': None, 'rmin': 0, 'rmax': 50, 'vmin': 1e-5, 'vmax': 1e-1,
        'fill_vmin': False,  # fills the x < vmin with vmin
        'xscale': None, 'yscale': None,
        'mask': 'positive', 'cmap': 'Blues', 'norm': "log", 'todo': None,
        'fancyticks': False
    }
    o_plot.set_plot_dics.append(int_vel_r_dic_rev)
    o_plot.set_plot_dics.append(int_vel_r_dic)

    o_plot.main()
# plot_test_velocity_jflux()

def plot_dung_jflux_corr():

    sim = "DD2_M13641364_M0_SR"

    o_corr_data = LOAD_RES_CORR(sim)

    o_plot = PLOT_MANY_TASKS()
    o_plot.gen_set["figdir"] = Paths.plots + '/combine_test/'
    o_plot.gen_set["type"] = "cartesian"
    o_plot.gen_set["figsize"] = (3.5, 3.5)  # <->, |]
    o_plot.gen_set["figname"] = "corr_{0:07d}.png".format(int(669810))

    corr_dic_ang_mom_flux_dens_unb_bern = {  # relies on the "get_res_corr(self, it, v_n): " method of data object
        'task': '2d colormesh', 'dtype': 'corr', 'ptype': 'cartesian',
        'data': o_corr_data, 'it': 669810,
        'position': (1, 1),
        'v_n_x': 'dens_unb_bern', 'v_n_y': 'ang_mom_flux', 'v_n': 'mass',
        'xmin': 1e-11, 'xmax': 1e-7, 'ymin': 1e-11, 'ymax': 1e-7, 'vmin': 1e-7, 'vmax': None,
        'xscale': 'log', 'yscale': 'log',
        'mask_below': None, 'mask_above': None, 'cmap': 'inferno_r', 'norm': 'log', 'todo': None
    }

    o_plot.set_plot_dics = [corr_dic_ang_mom_flux_dens_unb_bern]
    o_plot.main()
# plot_dung_jflux_hist()

def plot_dung_jflux_corr_3sim():

    o_plot = PLOT_MANY_TASKS()
    o_plot.gen_set["figdir"] = Paths.plots + '/combine_test/'
    o_plot.gen_set["type"] = "cartesian"
    o_plot.gen_set["figsize"] = (8.5, 2.7)  # <->, |]
    o_plot.gen_set["figname"] = "corr_Dunb_Jflux_3sims.png"
    o_plot.gen_set["sharex"] = False
    o_plot.gen_set["sharey"] = False
    o_plot.gen_set["subplots_adjust_h"] = 0.3
    o_plot.gen_set["subplots_adjust_w"] = 0.4
    ''' --- --- sim 1 --- --- --- '''

    sim = "DD2_M13641364_M0_SR"
    o_corr_data = LOAD_RES_CORR(sim)

    corr_dic_ang_mom_flux_dens_unb_bern = {  # relies on the "get_res_corr(self, it, v_n): " method of data object
        'task': '2d colormesh', 'dtype': 'corr', 'ptype': 'cartesian',
        'data': o_corr_data, 'it': 751278,
        'position': (1, 1),
        'v_n_x': 'dens_unb_bern', 'v_n_y': 'ang_mom_flux', 'v_n': 'mass',
        'xmin': 1e-11, 'xmax': 1e-7, 'ymin': 1e-11, 'ymax': 1e-7, 'vmin': 1e-7, 'vmax': None,
        'xscale': 'log', 'yscale': 'log',
        'mask_below': None, 'mask_above': None, 'cmap': 'inferno_r', 'norm': 'log', 'todo': None,
        'title': o_corr_data.sim.replace('_', '\_')
    }

    ''' --- --- sim 2 --- --- --- '''

    sim = "LS220_M13641364_M0_SR"
    o_corr_data = LOAD_RES_CORR(sim)

    corr_dic_ang_mom_flux_dens_unb_bern2 = {  # relies on the "get_res_corr(self, it, v_n): " method of data object
        'task': '2d colormesh', 'dtype': 'corr', 'ptype': 'cartesian',
        'data': o_corr_data, 'it': 753664,
        'position': (1, 2),
        'v_n_x': 'dens_unb_bern', 'v_n_y': 'ang_mom_flux', 'v_n': 'mass',
        'xmin': 1e-11, 'xmax': 1e-7, 'ymin': 1e-11, 'ymax': 1e-7, 'vmin': 1e-7, 'vmax': None,
        'xscale': 'log', 'yscale': 'log',
        'mask_below': None, 'mask_above': None, 'cmap': 'inferno_r', 'norm': 'log', 'todo': None,
        'title': o_corr_data.sim.replace('_', '\_')
    }

    ''' --- --- sim 3 --- --- --- '''

    sim = "SLy4_M13641364_M0_SR"
    o_corr_data = LOAD_RES_CORR(sim)

    corr_dic_ang_mom_flux_dens_unb_bern3 = {  # relies on the "get_res_corr(self, it, v_n): " method of data object
        'task': '2d colormesh', 'dtype': 'corr', 'ptype': 'cartesian',
        'data': o_corr_data, 'it': 761856,
        'position': (1, 3),
        'v_n_x': 'dens_unb_bern', 'v_n_y': 'ang_mom_flux', 'v_n': 'mass', 'cbar': 'right .03 .0',
        'xmin': 1e-11, 'xmax': 1e-7, 'ymin': 1e-11, 'ymax': 1e-7, 'vmin': 1e-7, 'vmax': None,
        'xscale': 'log', 'yscale': 'log',
        'mask_below': None, 'mask_above': None, 'cmap': 'inferno_r', 'norm': 'log', 'todo': None,
        'title': o_corr_data.sim.replace('_', '\_')
    }

    o_plot.set_plot_dics = [corr_dic_ang_mom_flux_dens_unb_bern,
                            corr_dic_ang_mom_flux_dens_unb_bern2,
                            corr_dic_ang_mom_flux_dens_unb_bern3]
    o_plot.main()
# plot_dung_jflux_corr_3sim()



""" --- LINE PLOTS ---- """

def plot_m1_in_time():

    sim = "DD2_M13641364_M0_SR"

    o_dm_data = LOAD_DENSITY_MODES(sim)

    o_plot = PLOT_MANY_TASKS()
    o_plot.gen_set["figdir"] = Paths.plots + '/combine_test/'
    o_plot.gen_set["type"] = "cartesian"
    o_plot.gen_set["figsize"] = (3.5, 3.5)  # <->, |]
    o_plot.gen_set["figname"] = "dens_modes.png"

    densmode_m1 = {
        'task': 'line', 'dtype': 'dm', 'ptype': 'cartesian',
        'data': o_dm_data,
        'position': (1, 1),
        'v_n_x': 'times', 'v_n_y': 'int_phi_r abs',
        'x_units': 'ms',
        'mode': 1,
        'ls': '-', 'color': 'black', 'lw': 1.,
        'label': 'sim', 'ylabel': r'$C_m/C_0$ Magnitude', 'xlabel': r'time [ms]',
        'xmin': None, 'xmax': None, 'ymin': None, 'ymax': None,
        'xscale': None, 'yscale': 'log', 'legend': True,
        'fancyticks': True
    }

    densmode_m2 = copy.deepcopy(densmode_m1)
    densmode_m2['mode'] = 2
    densmode_m2['lw'] = 0.5
    densmode_m2['color'] = 'gray'

    # densmode2 = densmode.copy()
    # densmode2["v_n"] = 'int_phi_r'
    # densmode2["ls"] = '--'

    o_plot.set_plot_dics = [densmode_m1, densmode_m2]
    o_plot.main()
# plot_m1_in_time()


def plot_m1_in_time_3sims_tcoll_tsub():


    o_plot = PLOT_MANY_TASKS()
    o_plot.gen_set["figdir"] = Paths.plots + '/combine_test/'
    o_plot.gen_set["type"] = "cartesian"
    o_plot.gen_set["figsize"] = (9.0, 3.7)  # <->, |]
    o_plot.gen_set["figname"] = "dens_modes_3sim.png"
    o_plot.set_plot_dics = []

    " --- --- --- SIM0 --- --- --- ---  --- --- --- "
    sim = "DD2_M13641364_M0_LK_HR_R04"

    # o_1d_data = COMPUTE_STORE_PAR(sim)

    o_dm_data = LOAD_DENSITY_MODES(sim)
    densmode_m0 = {
        'task': 'line', 'dtype': 'dm', 'ptype': 'cartesian',
        'data': o_dm_data,
        'position': (1, 1),
        'v_n_x': 'times', 'v_n_y': 'int_phi_r abs',
        'xunits': 'ms', #'-t': o_1d_data.get_par("tmerger_gw"),
        'mode': 0, 'norm_to_m': 0,
        'ls': ':', 'color': 'purple', 'lw': 1.,
        'label': 'DD2 LK HR', 'ylabel': r'$C_m/C_0$ Magnitude', 'xlabel': r'time [ms]',
        'xmin': None, 'xmax': None, 'ymin': None, 'ymax': None,
        'xscale': None, 'yscale': 'log', 'legend': True,
        'fancyticks': True, 'minorticks': True
    }
    densmode_m1 = {
        'task': 'line', 'dtype': 'dm', 'ptype': 'cartesian',
        'data': o_dm_data,
        'position': (1, 1),
        'v_n_x': 'times', 'v_n_y': 'int_phi_r abs',
        'xunits': 'ms', #'-t': o_1d_data.get_par("tmerger_gw"),
        'mode': 1, 'norm_to_m': 0,
        'ls': '-', 'color': 'purple', 'lw': 1.,
        'label': None, 'ylabel': r'$C_m/C_0$ Magnitude', 'xlabel': r'time [ms]',
        'xmin': None, 'xmax': None, 'ymin': None, 'ymax': None,
        'xscale': None, 'yscale': 'log', 'legend': False,
        'fancyticks': True, 'minorticks': True
    }
    densmode_m2 = {
        'task': 'line', 'dtype': 'dm', 'ptype': 'cartesian',
        'data': o_dm_data,
        'position': (1, 1),
        'v_n_x': 'times', 'v_n_y': 'int_phi_r abs',
        'xunits': 'ms', #'-t': o_1d_data.get_par("tmerger_gw"),
        'mode': 2, 'norm_to_m': 0,
        'ls': '--', 'color': 'purple', 'lw': 0.6,
        'label': None, 'ylabel': r'$C_m/C_0$ Magnitude', 'xlabel': r'time [ms]',
        'xmin': None, 'xmax': None, 'ymin': None, 'ymax': None,
        'xscale': None, 'yscale': 'log', 'legend': False,
        'fancyticks': True, 'minorticks': True
    }
    o_plot.set_plot_dics.append(densmode_m0)
    o_plot.set_plot_dics.append(densmode_m1)
    o_plot.set_plot_dics.append(densmode_m2)



    " --- --- --- SIM1 --- --- --- "
    sim = "DD2_M13641364_M0_SR"

    o_1d_data = COMPUTE_STORE_PAR(sim)

    o_dm_data = LOAD_DENSITY_MODES(sim)
    densmode_m0 = {
        'task': 'line', 'dtype': 'dm', 'ptype': 'cartesian',
        'data': o_dm_data,
        'position': (1, 1),
        'v_n_x': 'times', 'v_n_y': 'int_phi_r abs',
        'xunits': 'ms', '-t': o_1d_data.get_par("tmerger_gw"),
        'mode': 0, 'norm_to_m': 0,
        'ls': ':', 'color': 'green', 'lw': 1.,
        'label': 'DD2', 'ylabel': r'$C_m/C_0$ Magnitude', 'xlabel': r'time [ms]',
        'xmin': None, 'xmax': None, 'ymin': None, 'ymax': None,
        'xscale': None, 'yscale': 'log', 'legend': True,
        'fancyticks': True, 'minorticks': True
    }
    densmode_m1 = {
        'task': 'line', 'dtype': 'dm', 'ptype': 'cartesian',
        'data': o_dm_data,
        'position': (1, 1),
        'v_n_x': 'times', 'v_n_y': 'int_phi_r abs',
        'xunits': 'ms', '-t': o_1d_data.get_par("tmerger_gw"),
        'mode': 1, 'norm_to_m': 0,
        'ls': '-', 'color': 'green', 'lw': 1.,
        'label': None, 'ylabel': r'$C_m/C_0$ Magnitude', 'xlabel': r'time [ms]',
        'xmin': None, 'xmax': None, 'ymin': None, 'ymax': None,
        'xscale': None, 'yscale': 'log', 'legend': False,
        'fancyticks': True, 'minorticks': True
    }
    densmode_m2 = {
        'task': 'line', 'dtype': 'dm', 'ptype': 'cartesian',
        'data': o_dm_data,
        'position': (1, 1),
        'v_n_x': 'times', 'v_n_y': 'int_phi_r abs',
        'xunits': 'ms', '-t': o_1d_data.get_par("tmerger_gw"),
        'mode': 2, 'norm_to_m': 0,
        'ls': '--', 'color': 'green', 'lw': 0.6,
        'label': None, 'ylabel': r'$C_m/C_0$ Magnitude', 'xlabel': r'time [ms]',
        'xmin': None, 'xmax': None, 'ymin': None, 'ymax': None,
        'xscale': None, 'yscale': 'log', 'legend': False,
        'fancyticks': True, 'minorticks': True
    }
    o_plot.set_plot_dics.append(densmode_m0)
    o_plot.set_plot_dics.append(densmode_m1)
    o_plot.set_plot_dics.append(densmode_m2)



    " --- --- --- SIM1 --- --- --- "
    sim = "LS220_M13641364_M0_SR"
    o_1d_data = COMPUTE_STORE_PAR(sim)
    coll_time_vert = {
        'task': 'vertline', 'dtype': '-', 'ptype': 'cartesian',
        'value': (o_1d_data.get_par("tcoll_gw") - o_1d_data.get_par("tmerger_gw"))*1e3,
        'position': (1, 1),
        'ls': '-.', 'color': 'blue', 'lw': 0.5,
    }
    o_plot.set_plot_dics.append(coll_time_vert)

    o_dm_data = LOAD_DENSITY_MODES(sim)
    densmode_m0 = {
        'task': 'line', 'dtype': 'dm', 'ptype': 'cartesian',
        'data': o_dm_data,
        'position': (1, 1),
        'v_n_x': 'times', 'v_n_y': 'int_phi_r abs',
        'xunits': 'ms', '-t': o_1d_data.get_par("tmerger_gw"),
        'mode': 0, 'norm_to_m': 0,
        'ls': ':', 'color': 'blue', 'lw': 1.,
        'label': 'LS220', 'ylabel': r'$C_m/C_0$ Magnitude', 'xlabel': r'time [ms]',
        'xmin': None, 'xmax': None, 'ymin': None, 'ymax': None,
        'xscale': None, 'yscale': 'log', 'legend': True,
        'fancyticks': True, 'minorticks': True
    }
    densmode_m1 = {
        'task': 'line', 'dtype': 'dm', 'ptype': 'cartesian',
        'data': o_dm_data,
        'position': (1, 1),
        'v_n_x': 'times', 'v_n_y': 'int_phi_r abs',
        'xunits': 'ms', '-t': o_1d_data.get_par("tmerger_gw"),
        'mode': 1, 'norm_to_m': 0,
        'ls': '-', 'color': 'blue', 'lw': 1.,
        'label': None, 'ylabel': r'$C_m/C_0$ Magnitude', 'xlabel': r'time [ms]',
        'xmin': None, 'xmax': None, 'ymin': None, 'ymax': None,
        'xscale': None, 'yscale': 'log', 'legend': False,
        'fancyticks': True, 'minorticks': True
    }
    densmode_m2 = {
        'task': 'line', 'dtype': 'dm', 'ptype': 'cartesian',
        'data': o_dm_data,
        'position': (1, 1),
        'v_n_x': 'times', 'v_n_y': 'int_phi_r abs',
        'xunits': 'ms', '-t': o_1d_data.get_par("tmerger_gw"),
        'mode': 2, 'norm_to_m': 0,
        'ls': '--', 'color': 'blue', 'lw': 0.6,
        'label': None, 'ylabel': r'$C_m/C_0$ Magnitude', 'xlabel': r'time [ms]',
        'xmin': None, 'xmax': None, 'ymin': None, 'ymax': None,
        'xscale': None, 'yscale': 'log', 'legend': False,
        'fancyticks': True, 'minorticks': True
    }
    o_plot.set_plot_dics.append(densmode_m0)
    o_plot.set_plot_dics.append(densmode_m1)
    o_plot.set_plot_dics.append(densmode_m2)

    " --- --- --- SIM3 --- --- --- "
    sim = "SLy4_M13641364_M0_SR"

    o_1d_data = COMPUTE_STORE_PAR(sim)
    coll_time_vert = {
        'task': 'vertline', 'dtype': '-', 'ptype': 'cartesian',
        'value': (o_1d_data.get_par("tcoll_gw") - o_1d_data.get_par("tmerger_gw"))*1e3,
        'position': (1, 1),
        'ls': '-.', 'color': 'red', 'lw': 0.5,
    }
    o_plot.set_plot_dics.append(coll_time_vert)

    o_dm_data = LOAD_DENSITY_MODES(sim)
    densmode_m0 = {
        'task': 'line', 'dtype': 'dm', 'ptype': 'cartesian',
        'data': o_dm_data,
        'position': (1, 1),
        'v_n_x': 'times', 'v_n_y': 'int_phi_r abs',
        'xunits': 'ms', '-t': o_1d_data.get_par("tmerger_gw"),
        'mode': 0, 'norm_to_m': 0,
        'ls': ':', 'color': 'red', 'lw': 1.,
        'label': 'SLy4', 'ylabel': r'$C_m/C_0$ Magnitude', 'xlabel': r'time [ms]',
        'xmin': None, 'xmax': None, 'ymin': None, 'ymax': None,
        'xscale': None, 'yscale': 'log', 'legend': True,
        'fancyticks': True, 'minorticks': True
    }
    densmode_m1 = {
        'task': 'line', 'dtype': 'dm', 'ptype': 'cartesian',
        'data': o_dm_data,
        'position': (1, 1),
        'v_n_x': 'times', 'v_n_y': 'int_phi_r abs',
        'xunits': 'ms', '-t': o_1d_data.get_par("tmerger_gw"),
        'mode': 1, 'norm_to_m': 0,
        'ls': '-', 'color': 'red', 'lw': 1.,
        'label': None, 'ylabel': r'$C_m/C_0$ Magnitude', 'xlabel': r'time [ms]',
        'xmin': None, 'xmax': None, 'ymin': None, 'ymax': None,
        'xscale': None, 'yscale': 'log', 'legend': False,
        'fancyticks': True, 'minorticks': True
    }

    densmode_m2 = {
        'task': 'line', 'dtype': 'dm', 'ptype': 'cartesian',
        'data': o_dm_data,
        'position': (1, 1),
        'v_n_x': 'times', 'v_n_y': 'int_phi_r abs',
        'xunits': 'ms', '-t': o_1d_data.get_par("tmerger_gw"),
        'mode': 2, 'norm_to_m': 0,
        'ls': '--', 'color': 'red', 'lw': 0.6,
        'label': None, 'ylabel': r'$C_m/C_0$ Magnitude', 'xlabel': r'time [ms]',
        'xmin': None, 'xmax': None, 'ymin': None, 'ymax': None,
        'xscale': None, 'yscale': 'log', 'legend': False,
        'fancyticks': True, 'minorticks': True
    }
    o_plot.set_plot_dics.append(densmode_m0)
    o_plot.set_plot_dics.append(densmode_m1)
    o_plot.set_plot_dics.append(densmode_m2)

    o_plot.main()
# plot_m1_in_time_3sims_tcoll_tsub()

def plot_m1_modified_comparison():


    o_plot = PLOT_MANY_TASKS()
    o_plot.gen_set["figdir"] = Paths.plots + '/combine_test/'
    o_plot.gen_set["type"] = "cartesian"
    o_plot.gen_set["figsize"] = (9.0, 9.0)  # <->, |]
    o_plot.gen_set["figname"] = "dens_modes_final.png"
    o_plot.gen_set["sharex"] = True
    o_plot.gen_set["sharey"] = False
    o_plot.set_plot_dics = []

    " --- --- --- --- --- SIM0 --- --- --- ---  --- --- --- "
    sim = "DD2_M13641364_M0_SR"

    o_1d_data = COMPUTE_STORE_PAR(sim)

    o_dm_data = LOAD_DENSITY_MODES(sim)
    o_dm_data.gen_set['fname'] =  Paths.ppr_sims + sim + '/res_3d/' + "density_modes.h5"
    # densmode_m0 = {
    #     'task': 'line', 'dtype': 'dm', 'ptype': 'cartesian',
    #     'data': o_dm_data,
    #     'position': (1, 1),
    #     'v_n_x': 'times', 'v_n_y': 'int_phi_r abs',
    #     'xunits': 'ms', '-t': o_1d_data.get_par("tmerger_gw"),
    #     'mode': 0, 'norm_to_m': 0,
    #     'ls': '-.', 'color': 'red', 'lw': 0.4,
    #     'label': None, 'ylabel': r'$C_m/C_0$ Magnitude', 'xlabel': r'time [ms]',
    #     'xmin': None, 'xmax': None, 'ymin': None, 'ymax': None,
    #     'xscale': None, 'yscale': 'log', 'legend': True,
    #     'fancyticks': True, 'minorticks': True
    # }
    densmode_m1 = {
        'task': 'line', 'dtype': 'dm', 'ptype': 'cartesian',
        'data': o_dm_data,
        'position': (1, 1),
        'v_n_x': 'times', 'v_n_y': 'int_phi_r abs',
        'xunits': 'ms', '-t': o_1d_data.get_par("tmerger_gw"),
        'mode': 1, 'norm_to_m': 0,
        'ls': '-', 'color': 'green', 'lw': 1.,
        'label': 'DD2 (int)', 'ylabel': r'$C_m/C_0$ Magnitude', 'xlabel': r'time [ms]',
        'xmin': None, 'xmax': None, 'ymin': None, 'ymax': None,
        'xscale': None, 'yscale': 'log', 'legend': False,
        'fancyticks': True, 'minorticks': True,
    }
    densmode_m2 = {
        'task': 'line', 'dtype': 'dm', 'ptype': 'cartesian',
        'data': o_dm_data,
        'position': (1, 1),
        'v_n_x': 'times', 'v_n_y': 'int_phi_r abs',
        'xunits': 'ms', '-t': o_1d_data.get_par("tmerger_gw"),
        'mode': 2, 'norm_to_m': 0,
        'ls': ':', 'color': 'green', 'lw': 0.6,
        'label': None, 'ylabel': r'$C_m/C_0$ Magnitude', 'xlabel': r'time [ms]',
        'xmin': None, 'xmax': None, 'ymin': None, 'ymax': None,
        'xscale': None, 'yscale': 'log', 'legend': False,
        'fancyticks': True, 'minorticks': True
    }
    zero_time_vert = {
        'task': 'vertline', 'dtype': '-', 'ptype': 'cartesian',
        'value': 0,
        'position': (1, 1),
        'ls': '-.', 'color': 'gray', 'lw': 0.5, 'label': 'tmerg'
    }
    o_plot.set_plot_dics.append(zero_time_vert)
    # o_plot.set_plot_dics.append(densmode_m0)
    o_plot.set_plot_dics.append(densmode_m1)
    o_plot.set_plot_dics.append(densmode_m2)


    o_dm_data = LOAD_DENSITY_MODES(sim)
    o_dm_data.gen_set['fname'] = Paths.ppr_sims + sim + '/res_3d/' + "density_modes_lap15.h5"
    # densmode_m0 = {
    #     'task': 'line', 'dtype': 'dm', 'ptype': 'cartesian',
    #     'data': o_dm_data,
    #     'position': (1, 1),
    #     'v_n_x': 'times', 'v_n_y': 'int_phi_r abs',
    #     'xunits': 'ms', '-t': o_1d_data.get_par("tmerger_gw"),
    #     'mode': 0, 'norm_to_m': 0,
    #     'ls': '-.', 'color': 'green', 'lw': 0.4,
    #     'label': None, 'ylabel': r'$C_m/C_0$ Magnitude', 'xlabel': r'time [ms]',
    #     'xmin': None, 'xmax': None, 'ymin': None, 'ymax': None,
    #     'xscale': None, 'yscale': 'log', 'legend': True,
    #     'fancyticks': True, 'minorticks': True
    # }
    densmode_m1 = {
        'task': 'line', 'dtype': 'dm', 'ptype': 'cartesian',
        'data': o_dm_data,
        'position': (1, 1),
        'v_n_x': 'times', 'v_n_y': 'int_phi_r abs',
        'xunits': 'ms', '-t': o_1d_data.get_par("tmerger_gw"),
        'mode': 1, 'norm_to_m': 0,
        'ls': '-', 'color': 'lightgreen', 'lw': 1.,
        'label': 'DD2 (CM)', 'ylabel': r'$C_m/C_0$ Magnitude', 'xlabel': r'time [ms]',
        'xmin': None, 'xmax': None, 'ymin': None, 'ymax': None,
        'xscale': None, 'yscale': 'log', 'legend': False,
        'fancyticks': True, 'minorticks': True
    }
    densmode_m2 = {
        'task': 'line', 'dtype': 'dm', 'ptype': 'cartesian',
        'data': o_dm_data,
        'position': (1, 1),
        'v_n_x': 'times', 'v_n_y': 'int_phi_r abs',
        'xunits': 'ms', '-t': o_1d_data.get_par("tmerger_gw"),
        'mode': 2, 'norm_to_m': 0,
        'ls': ':', 'color': 'lightgreen', 'lw': 0.6,
        'label': None, 'ylabel': r'$C_m/C_0$ Magnitude', 'xlabel': r'time [ms]',
        'xmin': None, 'xmax': None, 'ymin': None, 'ymax': None,
        'xscale': None, 'yscale': 'log', 'legend': True,
        'fancyticks': True, 'minorticks': True
    }
    # o_plot.set_plot_dics.append(densmode_m0)
    o_plot.set_plot_dics.append(densmode_m1)
    o_plot.set_plot_dics.append(densmode_m2)

    " --- --- --- --- --- --- ---- --- --- --- --- --- -----"
    # ------------------------------------------------------#
    " --- --- --- --- --- SIM0 --- --- --- ---  --- --- --- "
    sim = "LS220_M13641364_M0_SR"

    o_1d_data = COMPUTE_STORE_PAR(sim)

    o_dm_data = LOAD_DENSITY_MODES(sim)
    o_dm_data.gen_set['fname'] = Paths.ppr_sims + sim + '/res_3d/' + "density_modes.h5"
    # densmode_m0 = {
    #     'task': 'line', 'dtype': 'dm', 'ptype': 'cartesian',
    #     'data': o_dm_data,
    #     'position': (1, 1),
    #     'v_n_x': 'times', 'v_n_y': 'int_phi_r abs',
    #     'xunits': 'ms', '-t': o_1d_data.get_par("tmerger_gw"),
    #     'mode': 0, 'norm_to_m': 0,
    #     'ls': '-.', 'color': 'blue', 'lw': 0.4,
    #     'label': None, 'ylabel': r'$C_m/C_0$ Magnitude', 'xlabel': r'time [ms]',
    #     'xmin': None, 'xmax': None, 'ymin': None, 'ymax': None,
    #     'xscale': None, 'yscale': 'log', 'legend': True,
    #     'fancyticks': True, 'minorticks': True
    # }
    densmode_m1 = {
        'task': 'line', 'dtype': 'dm', 'ptype': 'cartesian',
        'data': o_dm_data,
        'position': (2, 1),
        'v_n_x': 'times', 'v_n_y': 'int_phi_r abs',
        'xunits': 'ms', '-t': o_1d_data.get_par("tmerger_gw"),
        'mode': 1, 'norm_to_m': 0,
        'ls': '-', 'color': 'blue', 'lw': 1.,
        'label': 'LS220 (int)', 'ylabel': r'$C_m/C_0$ Magnitude', 'xlabel': r'time [ms]',
        'xmin': None, 'xmax': None, 'ymin': None, 'ymax': None,
        'xscale': None, 'yscale': 'log', 'legend': False,
        'fancyticks': True, 'minorticks': True
    }
    densmode_m2 = {
        'task': 'line', 'dtype': 'dm', 'ptype': 'cartesian',
        'data': o_dm_data,
        'position': (2, 1),
        'v_n_x': 'times', 'v_n_y': 'int_phi_r abs',
        'xunits': 'ms', '-t': o_1d_data.get_par("tmerger_gw"),
        'mode': 2, 'norm_to_m': 0,
        'ls': ':', 'color': 'blue', 'lw': 0.6,
        'label': None, 'ylabel': r'$C_m/C_0$ Magnitude', 'xlabel': r'time [ms]',
        'xmin': None, 'xmax': None, 'ymin': None, 'ymax': None,
        'xscale': None, 'yscale': 'log', 'legend': False,
        'fancyticks': True, 'minorticks': True
    }
    coll_time_vert = {
        'task': 'vertline', 'dtype': '-', 'ptype': 'cartesian',
        'value': (o_1d_data.get_par("tcoll_gw") - o_1d_data.get_par("tmerger_gw"))*1e3,
        'position': (2, 1),
        'ls': '-.', 'color': 'blue', 'lw': 0.5, 'label': 'tcoll'
    }
    zero_time_vert = {
        'task': 'vertline', 'dtype': '-', 'ptype': 'cartesian',
        'value': 0,
        'position': (2, 1),
        'ls': '-.', 'color': 'gray', 'lw': 0.5, 'label': 'tmerg'
    }
    o_plot.set_plot_dics.append(zero_time_vert)
    o_plot.set_plot_dics.append(coll_time_vert)
    # o_plot.set_plot_dics.append(densmode_m0)
    o_plot.set_plot_dics.append(densmode_m1)
    o_plot.set_plot_dics.append(densmode_m2)


    o_dm_data = LOAD_DENSITY_MODES(sim)
    o_dm_data.gen_set['fname'] = Paths.ppr_sims + sim + '/res_3d/' + "density_modes_lap15.h5"
    # densmode_m0 = {
    #     'task': 'line', 'dtype': 'dm', 'ptype': 'cartesian',
    #     'data': o_dm_data,
    #     'position': (1, 1),
    #     'v_n_x': 'times', 'v_n_y': 'int_phi_r abs',
    #     'xunits': 'ms', '-t': o_1d_data.get_par("tmerger_gw"),
    #     'mode': 0, 'norm_to_m': 0,
    #     'ls': '-.', 'color': 'blue', 'lw': 0.4,
    #     'label': None, 'ylabel': r'$C_m/C_0$ Magnitude', 'xlabel': r'time [ms]',
    #     'xmin': None, 'xmax': None, 'ymin': None, 'ymax': None,
    #     'xscale': None, 'yscale': 'log', 'legend': True,
    #     'fancyticks': True, 'minorticks': True
    # }
    densmode_m1 = {
        'task': 'line', 'dtype': 'dm', 'ptype': 'cartesian',
        'data': o_dm_data,
        'position': (2, 1),
        'v_n_x': 'times', 'v_n_y': 'int_phi_r abs',
        'xunits': 'ms', '-t': o_1d_data.get_par("tmerger_gw"),
        'mode': 1, 'norm_to_m': 0,
        'ls': '-', 'color': 'cyan', 'lw': 1.,
        'label': 'LS220 (CM)', 'ylabel': r'$C_m/C_0$ Magnitude', 'xlabel': r'time [ms]',
        'xmin': None, 'xmax': None, 'ymin': None, 'ymax': None,
        'xscale': None, 'yscale': 'log', 'legend': False,
        'fancyticks': True, 'minorticks': True
    }
    densmode_m2 = {
        'task': 'line', 'dtype': 'dm', 'ptype': 'cartesian',
        'data': o_dm_data,
        'position': (2, 1),
        'v_n_x': 'times', 'v_n_y': 'int_phi_r abs',
        'xunits': 'ms', '-t': o_1d_data.get_par("tmerger_gw"),
        'mode': 2, 'norm_to_m': 0,
        'ls': ':', 'color': 'cyan', 'lw': 0.6,
        'label': None, 'ylabel': r'$C_m/C_0$ Magnitude', 'xlabel': r'time [ms]',
        'xmin': None, 'xmax': None, 'ymin': None, 'ymax': None,
        'xscale': None, 'yscale': 'log', 'legend': True,
        'fancyticks': True, 'minorticks': True
    }
    # o_plot.set_plot_dics.append(densmode_m0)
    o_plot.set_plot_dics.append(densmode_m1)
    o_plot.set_plot_dics.append(densmode_m2)

    " --- --- --- --- --- --- ---- --- --- --- --- --- -----"
    # ------------------------------------------------------#
    " --- --- --- --- --- SIM0 --- --- --- ---  --- --- --- "
    sim = "SLy4_M13641364_M0_SR"

    o_1d_data = COMPUTE_STORE_PAR(sim)

    o_dm_data = LOAD_DENSITY_MODES(sim)
    o_dm_data.gen_set['fname'] = Paths.ppr_sims + sim + '/res_3d/' + "density_modes.h5"
    # densmode_m0 = {
    #     'task': 'line', 'dtype': 'dm', 'ptype': 'cartesian',
    #     'data': o_dm_data,
    #     'position': (1, 1),
    #     'v_n_x': 'times', 'v_n_y': 'int_phi_r abs',
    #     'xunits': 'ms', '-t': o_1d_data.get_par("tmerger_gw"),
    #     'mode': 0, 'norm_to_m': 0,
    #     'ls': '-.', 'color': 'red', 'lw': 0.4,
    #     'label': None, 'ylabel': r'$C_m/C_0$ Magnitude', 'xlabel': r'time [ms]',
    #     'xmin': None, 'xmax': None, 'ymin': None, 'ymax': None,
    #     'xscale': None, 'yscale': 'log', 'legend': True,
    #     'fancyticks': True, 'minorticks': True
    # }
    densmode_m1 = {
        'task': 'line', 'dtype': 'dm', 'ptype': 'cartesian',
        'data': o_dm_data,
        'position': (3, 1),
        'v_n_x': 'times', 'v_n_y': 'int_phi_r abs',
        'xunits': 'ms', '-t': o_1d_data.get_par("tmerger_gw"),
        'mode': 1, 'norm_to_m': 0,
        'ls': '-', 'color': 'red', 'lw': 1.,
        'label': 'Sly4 (int)', 'ylabel': r'$C_m/C_0$ Magnitude', 'xlabel': r'time [ms]',
        'xmin': None, 'xmax': None, 'ymin': None, 'ymax': None,
        'xscale': None, 'yscale': 'log', 'legend': False,
        'fancyticks': True, 'minorticks': True
    }
    densmode_m2 = {
        'task': 'line', 'dtype': 'dm', 'ptype': 'cartesian',
        'data': o_dm_data,
        'position': (3, 1),
        'v_n_x': 'times', 'v_n_y': 'int_phi_r abs',
        'xunits': 'ms', '-t': o_1d_data.get_par("tmerger_gw"),
        'mode': 2, 'norm_to_m': 0,
        'ls': ':', 'color': 'red', 'lw': 0.6,
        'label': None, 'ylabel': r'$C_m/C_0$ Magnitude', 'xlabel': r'time [ms]',
        'xmin': None, 'xmax': None, 'ymin': None, 'ymax': None,
        'xscale': None, 'yscale': 'log', 'legend': False,
        'fancyticks': True, 'minorticks': True
    }
    coll_time_vert = {
        'task': 'vertline', 'dtype': '-', 'ptype': 'cartesian',
        'value': (o_1d_data.get_par("tcoll_gw") - o_1d_data.get_par("tmerger_gw"))*1e3,
        'position': (3, 1),
        'ls': '-.', 'color': 'red', 'lw': 0.5, 'label': 'tcoll'
    }
    zero_time_vert = {
        'task': 'vertline', 'dtype': '-', 'ptype': 'cartesian',
        'value': 0,
        'position': (3, 1),
        'ls': '-.', 'color': 'gray', 'lw': 0.5, 'label': 'tmerg'
    }
    o_plot.set_plot_dics.append(zero_time_vert)
    o_plot.set_plot_dics.append(coll_time_vert)
    # o_plot.set_plot_dics.append(densmode_m0)
    o_plot.set_plot_dics.append(densmode_m1)
    o_plot.set_plot_dics.append(densmode_m2)


    o_dm_data = LOAD_DENSITY_MODES(sim)
    o_dm_data.gen_set['fname'] = Paths.ppr_sims + sim + '/res_3d/' + "density_modes_lap15.h5"
    # densmode_m0 = {
    #     'task': 'line', 'dtype': 'dm', 'ptype': 'cartesian',
    #     'data': o_dm_data,
    #     'position': (1, 1),
    #     'v_n_x': 'times', 'v_n_y': 'int_phi_r abs',
    #     'xunits': 'ms', '-t': o_1d_data.get_par("tmerger_gw"),
    #     'mode': 0, 'norm_to_m': 0,
    #     'ls': '-.', 'color': 'red', 'lw': 0.4,
    #     'label': None, 'ylabel': r'$C_m/C_0$ Magnitude', 'xlabel': r'time [ms]',
    #     'xmin': None, 'xmax': None, 'ymin': None, 'ymax': None,
    #     'xscale': None, 'yscale': 'log', 'legend': True,
    #     'fancyticks': True, 'minorticks': True
    # }
    densmode_m1 = {
        'task': 'line', 'dtype': 'dm', 'ptype': 'cartesian',
        'data': o_dm_data,
        'position': (3, 1),
        'v_n_x': 'times', 'v_n_y': 'int_phi_r abs',
        'xunits': 'ms', '-t': o_1d_data.get_par("tmerger_gw"),
        'mode': 1, 'norm_to_m': 0,
        'ls': '-', 'color': 'orange', 'lw': 1.,
        'label': 'Sly4 (CM)', 'ylabel': r'$C_m/C_0$ Magnitude', 'xlabel': r'time [ms]',
        'xmin': None, 'xmax': None, 'ymin': None, 'ymax': None,
        'xscale': None, 'yscale': 'log', 'legend': True,
        'fancyticks': True, 'minorticks': True
    }
    densmode_m2 = {
        'task': 'line', 'dtype': 'dm', 'ptype': 'cartesian',
        'data': o_dm_data,
        'position': (3, 1),
        'v_n_x': 'times', 'v_n_y': 'int_phi_r abs',
        'xunits': 'ms', '-t': o_1d_data.get_par("tmerger_gw"),
        'mode': 2, 'norm_to_m': 0,
        'ls': ':', 'color': 'orange', 'lw': 0.6,
        'label': None, 'ylabel': r'$C_m/C_0$ Magnitude', 'xlabel': r'time [ms]',
        'xmin': None, 'xmax': None, 'ymin': None, 'ymax': None,
        'xscale': None, 'yscale': 'log', 'legend': False,
        'fancyticks': True, 'minorticks': True
    }
    # o_plot.set_plot_dics.append(densmode_m0)
    o_plot.set_plot_dics.append(densmode_m1)
    o_plot.set_plot_dics.append(densmode_m2)

    o_plot.main()
# plot_m1_modified_comparison()

" --- OUTFLOWED PLOTS --- "

def plot_hist_1sim():

    o_plot = PLOT_MANY_TASKS()
    o_plot.gen_set["figdir"] = Paths.plots + '/combine_test/'
    o_plot.gen_set["type"] = "cartesian"
    o_plot.gen_set["figsize"] = (4.5, 3.5)  # <->, |]
    o_plot.gen_set["figname"] = "hist_theta.png"

    sim = "DD2_M13641364_M0_SR"
    o_data = COMPUTE_STORE_PAR(sim)

    dic_hist = {
        'task': 'hist1d', 'ptype': 'cartesian',
        'position': (1, 1),
        'data': o_data, 'norm': True, 'criterion': '_0',
        'v_n_x': 'hist_theta', 'v_n_y': 'hist_theta_m',
        'color': 'black',
        'xlabel': r"Angle from orbital plane", 'ylabel': 'mass',
        'label': None, 'yscale': 'log',
        'fancyticks': True
    }

    o_plot.set_plot_dics = [dic_hist]
    o_plot.main()
# plot_hist_1sim()

def plot_3_hists_1sim():

    o_plot = PLOT_MANY_TASKS()
    o_plot.gen_set["figdir"] = Paths.plots + '/combine_test/'
    o_plot.gen_set["type"] = "cartesian"
    o_plot.gen_set["figsize"] = (3.5, 7.5)  # <->, |]
    o_plot.gen_set["figname"] = "hist_theta_vel_ye.png"
    o_plot.gen_set["sharex"] = False
    o_plot.gen_set["sharey"] = False
    o_plot.gen_set["subplots_adjust_h"] = 0.3

    sim = "DD2_M13641364_M0_SR"
    o_data = COMPUTE_STORE_PAR(sim)

    dic_hist_theta = {
        'task': 'hist1d', 'ptype': 'cartesian',
        'position': (1, 1),
        'data': o_data, 'norm': True, 'criterion': '_0',
        'v_n_x': 'hist_theta', 'v_n_y': 'hist_theta_m',
        'color': 'black', 'ls':'-', 'lw':1., 'ds':'steps',
        'ymin': 1e-4, 'ymax': 1e0,
        'xlabel': r"Angle from orbital plane", 'ylabel': 'mass',
        'label': None, 'yscale': 'log',
        'fancyticks': True
    }

    dic_hist_vel_inf = {
        'task': 'hist1d', 'ptype': 'cartesian',
        'position': (2, 1),
        'data': o_data, 'norm': True, 'criterion': '_0',
        'v_n_x': 'hist_vel_inf', 'v_n_y': 'hist_vel_inf_m',
        'color': 'black', 'ls':'-', 'lw':1., 'ds':'steps',
        'xmin': 0., 'xmax': 1., 'ymin': 1e-4, 'ymax': 1e0,
        'xlabel': r"$\upsilon_{\infty}$ [c]", 'ylabel': 'mass',
        'label': None, 'yscale': 'log',
        'fancyticks': True
    }

    dic_hist_ye = {
        'task': 'hist1d', 'ptype': 'cartesian',
        'position': (3, 1),
        'data': o_data, 'norm': True, 'criterion': '_0',
        'v_n_x': 'hist_ye', 'v_n_y': 'hist_ye_m',
        'color': 'black', 'ls':'-', 'lw':1., 'ds':'steps',
        'xmin': 0., 'xmax': 0.5, 'ymin': 1e-4, 'ymax': 1e0,
        'xlabel': r"$Y_e$", 'ylabel': 'mass',
        'label': None, 'yscale': 'log',
        'fancyticks': True
    }

    ''' --- BERNOULLI '''

    dic_hist_theta_b = {
        'task': 'hist1d', 'ptype': 'cartesian',
        'position': (1, 1),
        'data': o_data, 'norm': True, 'criterion': '_0_b_w',
        'v_n_x': 'hist_theta', 'v_n_y': 'hist_theta_m',
        'color': 'black', 'ls':'--', 'lw':1., 'ds':'steps',
        'ymin': 1e-4, 'ymax': 1e0,
        'xlabel': r"Angle from orbital plane", 'ylabel': 'mass',
        'label': None, 'yscale': 'log',
        'fancyticks': True
    }

    dic_hist_vel_inf_b = {
        'task': 'hist1d', 'ptype': 'cartesian',
        'position': (2, 1),
        'data': o_data, 'norm': True, 'criterion': '_0_b_w',
        'v_n_x': 'hist_vel_inf', 'v_n_y': 'hist_vel_inf_m',
        'color': 'black', 'ls':'--', 'lw':1., 'ds':'steps',
        'xmin': 0., 'xmax': 1., 'ymin': 1e-4, 'ymax': 1e0,
        'xlabel': r"$\upsilon_{\infty}$ [c]", 'ylabel': 'mass',
        'label': None, 'yscale': 'log',
        'fancyticks': True
    }

    dic_hist_ye_b = {
        'task': 'hist1d', 'ptype': 'cartesian',
        'position': (3, 1),
        'data': o_data, 'norm': True, 'criterion': '_0_b_w',
        'v_n_x': 'hist_ye', 'v_n_y': 'hist_ye_m',
        'color': 'black', 'ls':'--', 'lw':1., 'ds':'steps',
        'xmin': 0., 'xmax': 0.5, 'ymin': 1e-4, 'ymax': 1e0,
        'xlabel': r"$Y_e$", 'ylabel': 'mass',
        'label': None, 'yscale': 'log',
        'fancyticks': True
    }


    o_plot.set_plot_dics = [dic_hist_theta, dic_hist_theta_b,
                            dic_hist_vel_inf, dic_hist_vel_inf_b,
                            dic_hist_ye, dic_hist_ye_b]
    o_plot.main()
# plot_3_hists_1sim()


def plot_3_hists_3sim():

    o_plot = PLOT_MANY_TASKS()
    o_plot.gen_set["figdir"] = Paths.plots + '/combine_test/'
    o_plot.gen_set["type"] = "cartesian"
    o_plot.gen_set["figsize"] = (8.5, 7.5)  # <->, |]
    o_plot.gen_set["figname"] = "hist_theta_vel_ye_3sims.png"
    o_plot.gen_set["sharex"] = False
    o_plot.gen_set["sharey"] = True
    o_plot.gen_set["subplots_adjust_h"] = 0.3

    ' --- --- --- SIM1 --- --- --- '

    sim = "DD2_M13641364_M0_SR"
    o_data = COMPUTE_STORE_PAR(sim)

    dic_hist_theta = {
        'task': 'hist1d', 'ptype': 'cartesian',
        'position': (1, 1),
        'data': o_data, 'norm': True, 'criterion': '_0',
        'v_n_x': 'hist_theta', 'v_n_y': 'hist_theta_m',
        'color': 'black', 'ls':'-', 'lw':1., 'ds':'steps',
        'ymin': 1e-4, 'ymax': 1e0,
        'xlabel': r"Angle from orbital plane", 'ylabel': 'mass',
        'label': None, 'yscale': 'log', 'title': o_data.sim.replace('_', '\_'),
        'fancyticks': True
    }

    dic_hist_vel_inf = {
        'task': 'hist1d', 'ptype': 'cartesian',
        'position': (2, 1),
        'data': o_data, 'norm': True, 'criterion': '_0',
        'v_n_x': 'hist_vel_inf', 'v_n_y': 'hist_vel_inf_m',
        'color': 'black', 'ls':'-', 'lw':1., 'ds':'steps',
        'xmin': 0., 'xmax': 1., 'ymin': 1e-4, 'ymax': 1e0,
        'xlabel': r"$\upsilon_{\infty}$ [c]", 'ylabel': 'mass',
        'label': None, 'yscale': 'log',
        'fancyticks': True
    }

    dic_hist_ye = {
        'task': 'hist1d', 'ptype': 'cartesian',
        'position': (3, 1),
        'data': o_data, 'norm': True, 'criterion': '_0',
        'v_n_x': 'hist_ye', 'v_n_y': 'hist_ye_m',
        'color': 'black', 'ls':'-', 'lw':1., 'ds':'steps',
        'xmin': 0., 'xmax': 0.5, 'ymin': 1e-4, 'ymax': 1e0,
        'xlabel': r"$Y_e$", 'ylabel': 'mass',
        'label': None, 'yscale': 'log',
        'fancyticks': True
    }

    ''' --- BERNOULLI '''

    dic_hist_theta_b = {
        'task': 'hist1d', 'ptype': 'cartesian',
        'position': (1, 1),
        'data': o_data, 'norm': True, 'criterion': '_0_b_w',
        'v_n_x': 'hist_theta', 'v_n_y': 'hist_theta_m',
        'color': 'black', 'ls':'--', 'lw':1., 'ds':'steps',
        'ymin': 1e-4, 'ymax': 1e0,
        'xlabel': r"Angle from orbital plane", 'ylabel': 'mass',
        'label': None, 'yscale': 'log',
        'fancyticks': True
    }

    dic_hist_vel_inf_b = {
        'task': 'hist1d', 'ptype': 'cartesian',
        'position': (2, 1),
        'data': o_data, 'norm': True, 'criterion': '_0_b_w',
        'v_n_x': 'hist_vel_inf', 'v_n_y': 'hist_vel_inf_m',
        'color': 'black', 'ls':'--', 'lw':1., 'ds':'steps',
        'xmin': 0., 'xmax': 1., 'ymin': 1e-4, 'ymax': 1e0,
        'xlabel': r"$\upsilon_{\infty}$ [c]", 'ylabel': 'mass',
        'label': None, 'yscale': 'log',
        'fancyticks': True
    }

    dic_hist_ye_b = {
        'task': 'hist1d', 'ptype': 'cartesian',
        'position': (3, 1),
        'data': o_data, 'norm': True, 'criterion': '_0_b_w',
        'v_n_x': 'hist_ye', 'v_n_y': 'hist_ye_m',
        'color': 'black', 'ls':'--', 'lw':1., 'ds':'steps',
        'xmin': 0., 'xmax': 0.5, 'ymin': 1e-4, 'ymax': 1e0,
        'xlabel': r"$Y_e$", 'ylabel': 'mass',
        'label': None, 'yscale': 'log',
        'fancyticks': True
    }

    ' --- --- --- SIM2 --- --- --- '

    sim = "LS220_M13641364_M0_SR"
    o_data2 = COMPUTE_STORE_PAR(sim)

    dic_hist_theta2 = {
        'task': 'hist1d', 'ptype': 'cartesian',
        'position': (1, 2),
        'data': o_data2, 'norm': True, 'criterion': '_0',
        'v_n_x': 'hist_theta', 'v_n_y': 'hist_theta_m',
        'color': 'blue', 'ls':'-', 'lw':1., 'ds':'steps',
        'ymin': 1e-4, 'ymax': 1e0,
        'xlabel': r"Angle from orbital plane", 'ylabel': 'mass',
        'label': None, 'yscale': 'log', 'title': o_data2.sim.replace('_', '\_'),
        'fancyticks': True
    }

    dic_hist_vel_inf2 = {
        'task': 'hist1d', 'ptype': 'cartesian',
        'position': (2, 2),
        'data': o_data2, 'norm': True, 'criterion': '_0',
        'v_n_x': 'hist_vel_inf', 'v_n_y': 'hist_vel_inf_m',
        'color': 'blue', 'ls':'-', 'lw':1., 'ds':'steps',
        'xmin': 0., 'xmax': 1., 'ymin': 1e-4, 'ymax': 1e0,
        'xlabel': r"$\upsilon_{\infty}$ [c]", 'ylabel': 'mass',
        'label': None, 'yscale': 'log',
        'fancyticks': True
    }

    dic_hist_ye2 = {
        'task': 'hist1d', 'ptype': 'cartesian',
        'position': (3, 2),
        'data': o_data2, 'norm': True, 'criterion': '_0',
        'v_n_x': 'hist_ye', 'v_n_y': 'hist_ye_m',
        'color': 'blue', 'ls':'-', 'lw':1., 'ds':'steps',
        'xmin': 0., 'xmax': 0.5, 'ymin': 1e-4, 'ymax': 1e0,
        'xlabel': r"$Y_e$", 'ylabel': 'mass',
        'label': None, 'yscale': 'log',
        'fancyticks': True
    }

    ''' --- BERNOULLI '''

    dic_hist_theta_b2 = {
        'task': 'hist1d', 'ptype': 'cartesian',
        'position': (1, 2),
        'data': o_data2, 'norm': True, 'criterion': '_0_b_w',
        'v_n_x': 'hist_theta', 'v_n_y': 'hist_theta_m',
        'color': 'blue', 'ls':'--', 'lw':1., 'ds':'steps',
        'ymin': 1e-4, 'ymax': 1e0,
        'xlabel': r"Angle from orbital plane", 'ylabel': 'mass',
        'label': None, 'yscale': 'log',
        'fancyticks': True
    }

    dic_hist_vel_inf_b2 = {
        'task': 'hist1d', 'ptype': 'cartesian',
        'position': (2, 2),
        'data': o_data2, 'norm': True, 'criterion': '_0_b_w',
        'v_n_x': 'hist_vel_inf', 'v_n_y': 'hist_vel_inf_m',
        'color': 'blue', 'ls':'--', 'lw':1., 'ds':'steps',
        'xmin': 0., 'xmax': 1., 'ymin': 1e-4, 'ymax': 1e0,
        'xlabel': r"$\upsilon_{\infty}$ [c]", 'ylabel': 'mass',
        'label': None, 'yscale': 'log',
        'fancyticks': True
    }

    dic_hist_ye_b2 = {
        'task': 'hist1d', 'ptype': 'cartesian',
        'position': (3, 2),
        'data': o_data2, 'norm': True, 'criterion': '_0_b_w',
        'v_n_x': 'hist_ye', 'v_n_y': 'hist_ye_m',
        'color': 'blue', 'ls':'--', 'lw':1., 'ds':'steps',
        'xmin': 0., 'xmax': 0.5, 'ymin': 1e-4, 'ymax': 1e0,
        'xlabel': r"$Y_e$", 'ylabel': 'mass',
        'label': None, 'yscale': 'log',
        'fancyticks': True
    }

    ' --- --- --- SIM3 --- --- --- '

    sim = "SLy4_M13641364_M0_SR"
    o_data3 = COMPUTE_STORE_PAR(sim)

    dic_hist_theta3 = {
        'task': 'hist1d', 'ptype': 'cartesian',
        'position': (1, 3),
        'data': o_data3, 'norm': True, 'criterion': '_0',
        'v_n_x': 'hist_theta', 'v_n_y': 'hist_theta_m',
        'color': 'red', 'ls': '-', 'lw': 1., 'ds': 'steps',
        'ymin': 1e-4, 'ymax': 1e0,
        'xlabel': r"Angle from orbital plane", 'ylabel': 'mass',
        'label': None, 'yscale': 'log', 'title': o_data3.sim.replace('_', '\_'),
        'fancyticks': True
    }

    dic_hist_vel_inf3 = {
        'task': 'hist1d', 'ptype': 'cartesian',
        'position': (2, 3),
        'data': o_data3, 'norm': True, 'criterion': '_0',
        'v_n_x': 'hist_vel_inf', 'v_n_y': 'hist_vel_inf_m',
        'color': 'red', 'ls': '-', 'lw': 1., 'ds': 'steps',
        'xmin': 0., 'xmax': 1., 'ymin': 1e-4, 'ymax': 1e0,
        'xlabel': r"$\upsilon_{\infty}$ [c]", 'ylabel': 'mass',
        'label': None, 'yscale': 'log',
        'fancyticks': True
    }

    dic_hist_ye3 = {
        'task': 'hist1d', 'ptype': 'cartesian',
        'position': (3, 3),
        'data': o_data3, 'norm': True, 'criterion': '_0',
        'v_n_x': 'hist_ye', 'v_n_y': 'hist_ye_m',
        'color': 'red', 'ls': '-', 'lw': 1., 'ds': 'steps',
        'xmin': 0., 'xmax': 0.5, 'ymin': 1e-4, 'ymax': 1e0,
        'xlabel': r"$Y_e$", 'ylabel': 'mass',
        'label': None, 'yscale': 'log',
        'fancyticks': True
    }

    ''' --- BERNOULLI '''

    dic_hist_theta_b3 = {
        'task': 'hist1d', 'ptype': 'cartesian',
        'position': (1, 3),
        'data': o_data3, 'norm': True, 'criterion': '_0_b_w',
        'v_n_x': 'hist_theta', 'v_n_y': 'hist_theta_m',
        'color': 'red', 'ls': '--', 'lw': 1., 'ds': 'steps',
        'ymin': 1e-4, 'ymax': 1e0,
        'xlabel': r"Angle from orbital plane", 'ylabel': 'mass',
        'label': None, 'yscale': 'log',
        'fancyticks': True
    }

    dic_hist_vel_inf_b3 = {
        'task': 'hist1d', 'ptype': 'cartesian',
        'position': (2, 3),
        'data': o_data3, 'norm': True, 'criterion': '_0_b_w',
        'v_n_x': 'hist_vel_inf', 'v_n_y': 'hist_vel_inf_m',
        'color': 'red', 'ls': '--', 'lw': 1., 'ds': 'steps',
        'xmin': 0., 'xmax': 1., 'ymin': 1e-4, 'ymax': 1e0,
        'xlabel': r"$\upsilon_{\infty}$ [c]", 'ylabel': 'mass',
        'label': None, 'yscale': 'log',
        'fancyticks': True
    }

    dic_hist_ye_b3 = {
        'task': 'hist1d', 'ptype': 'cartesian',
        'position': (3, 3),
        'data': o_data3, 'norm': True, 'criterion': '_0_b_w',
        'v_n_x': 'hist_ye', 'v_n_y': 'hist_ye_m',
        'color': 'red', 'ls': '--', 'lw': 1., 'ds': 'steps',
        'xmin': 0., 'xmax': 0.5, 'ymin': 1e-4, 'ymax': 1e0,
        'xlabel': r"$Y_e$", 'ylabel': 'mass',
        'label': None, 'yscale': 'log',
        'fancyticks': True
    }



    o_plot.set_plot_dics = [dic_hist_theta, dic_hist_theta_b,
                            dic_hist_vel_inf, dic_hist_vel_inf_b,
                            dic_hist_ye, dic_hist_ye_b,

                            dic_hist_theta2, dic_hist_theta_b2,
                            dic_hist_vel_inf2, dic_hist_vel_inf_b2,
                            dic_hist_ye2, dic_hist_ye_b2,

                            dic_hist_theta3, dic_hist_theta_b3,
                            dic_hist_vel_inf3, dic_hist_vel_inf_b3,
                            dic_hist_ye3, dic_hist_ye_b3
                            ]
    o_plot.main()
# plot_3_hists_3sim()

def plot_ej_profile_1sim():

    o_plot = PLOT_MANY_TASKS()
    o_plot.gen_set["figdir"] = Paths.plots + '/combine_test/'
    o_plot.gen_set["type"] = "cartesian"
    o_plot.gen_set["figsize"] = (3.5, 3.5)  # <->, |]
    o_plot.gen_set["figname"] = "ej_profile.png"
    o_plot.gen_set["sharex"] = False
    o_plot.gen_set["sharey"] = True
    o_plot.gen_set["subplots_adjust_h"] = 0.3

    sim = "DD2_M13641364_M0_SR"
    o_data = COMPUTE_STORE_PAR(sim)

    dic_ej_prof = {
        'task': 'ejprof', 'ptype': 'cartesian',
        'position': (1, 1),
        'data': o_data, 'criterion': '_0',
        'v_n_x': 't_tot_flux', 'v_n_y': 'mass_tot_flux',
        'color': 'black', 'ls': '-', 'lw': 1.,
        # 'ymin': 1e-4, 'ymax': 1e0,
        'xlabel': r"time [ms]", 'ylabel': r'$M_{ej}$ $[10^{-2}M_{\odot}]$',
        'xunits': 'ms',
        'label': None, 'yscale': None, 'title': o_data.sim.replace('_', '\_'),
        'fancyticks': True, 'minorticks': True
    }
    dic_ej_prof_b = {
        'task': 'ejprof', 'ptype': 'cartesian',
        'position': (1, 1),
        'data': o_data, 'criterion': '_0_b_w',
        'v_n_x': 't_tot_flux', 'v_n_y': 'mass_tot_flux',
        'color': 'black', 'ls': '--', 'lw': 1.,
        # 'ymin': 1e-4, 'ymax': 1e0,
        'xlabel': r"time [ms]", 'ylabel': r'$M_{ej}$ $[10^{-2}M_{\odot}]$',
        'xunits': 'ms',
        'label': None, 'yscale': None, 'title': o_data.sim.replace('_', '\_'),
        'fancyticks': True, 'minorticks': True
    }

    o_plot.set_plot_dics = [dic_ej_prof, dic_ej_prof_b]
    o_plot.main()
# plot_ej_profile_1sim()

def plot_ej_profile_3sim():

    o_plot = PLOT_MANY_TASKS()
    o_plot.gen_set["figdir"] = Paths.plots + '/combine_test/'
    o_plot.gen_set["type"] = "cartesian"
    o_plot.gen_set["figsize"] = (8.5, 3.5)  # <->, |]
    o_plot.gen_set["figname"] = "ej_profile_3sim.png"
    o_plot.gen_set["sharex"] = False
    o_plot.gen_set["sharey"] = True
    o_plot.gen_set["subplots_adjust_h"] = 0.3

    ''' --- --- sim 1 --- --- --- '''

    sim = "DD2_M13641364_M0_SR"
    o_data = COMPUTE_STORE_PAR(sim)

    dic_ej_prof = {
        'task': 'ejprof', 'ptype': 'cartesian',
        'position': (1, 1),
        'data': o_data, 'criterion': '_0',
        'v_n_x': 't_tot_flux', 'v_n_y': 'mass_tot_flux',
        'color': 'black', 'ls': '-', 'lw': 1.,
        # 'ymin': 1e-4, 'ymax': 1e0,
        'xlabel': r"time [ms]", 'ylabel': r'$M_{ej}$ $[10^{-2}M_{\odot}]$',
        'xunits': 'ms',
        'label': None, 'yscale': None, 'title': o_data.sim.replace('_', '\_'),
        'fancyticks': True, 'minorticks': True
    }
    dic_ej_prof_b = {
        'task': 'ejprof', 'ptype': 'cartesian',
        'position': (1, 1),
        'data': o_data, 'criterion': '_0_b_w',
        'v_n_x': 't_tot_flux', 'v_n_y': 'mass_tot_flux',
        'color': 'black', 'ls': '--', 'lw': 1.,
        # 'ymin': 1e-4, 'ymax': 1e0,
        'xlabel': r"time [ms]", 'ylabel': r'$M_{ej}$ $[10^{-2}M_{\odot}]$',
        'xunits': 'ms',
        'label': None, 'yscale': None, 'title': o_data.sim.replace('_', '\_'),
        'fancyticks': True, 'minorticks': True
    }

    ''' --- --- sim 2 --- --- --- '''

    sim = "LS220_M13641364_M0_SR"
    o_data = COMPUTE_STORE_PAR(sim)

    dic_ej_prof2 = {
        'task': 'ejprof', 'ptype': 'cartesian',
        'position': (1, 2),
        'data': o_data, 'criterion': '_0',
        'v_n_x': 't_tot_flux', 'v_n_y': 'mass_tot_flux',
        'color': 'blue', 'ls': '-', 'lw': 1.,
        # 'ymin': 1e-4, 'ymax': 1e0,
        'xlabel': r"time [ms]", 'ylabel': r'$M_{ej}$ $[10^{-2}M_{\odot}]$',
        'xunits': 'ms',
        'label': None, 'yscale': None, 'title': o_data.sim.replace('_', '\_'),
        'fancyticks': True, 'minorticks': True
    }
    dic_ej_prof_b2 = {
        'task': 'ejprof', 'ptype': 'cartesian',
        'position': (1, 2),
        'data': o_data, 'criterion': '_0_b_w',
        'v_n_x': 't_tot_flux', 'v_n_y': 'mass_tot_flux',
        'color': 'blue', 'ls': '--', 'lw': 1.,
        # 'ymin': 1e-4, 'ymax': 1e0,
        'xlabel': r"time [ms]", 'ylabel': r'$M_{ej}$ $[10^{-2}M_{\odot}]$',
        'xunits': 'ms',
        'label': None, 'yscale': None, 'title': o_data.sim.replace('_', '\_'),
        'fancyticks': True, 'minorticks': True
    }

    ''' --- --- sim 3 --- --- --- '''

    sim = "SLy4_M13641364_M0_SR"
    o_data = COMPUTE_STORE_PAR(sim)

    dic_ej_prof3 = {
        'task': 'ejprof', 'ptype': 'cartesian',
        'position': (1, 3),
        'data': o_data, 'criterion': '_0',
        'v_n_x': 't_tot_flux', 'v_n_y': 'mass_tot_flux',
        'color': 'red', 'ls': '-', 'lw': 1.,
        # 'ymin': 1e-4, 'ymax': 1e0,
        'xlabel': r"time [ms]", 'ylabel': r'$M_{ej}$ $[10^{-2}M_{\odot}]$',
        'xunits': 'ms',
        'label': None, 'yscale': None, 'title': o_data.sim.replace('_', '\_'),
        'fancyticks': True, 'minorticks': True
    }
    dic_ej_prof_b3 = {
        'task': 'ejprof', 'ptype': 'cartesian',
        'position': (1, 3),
        'data': o_data, 'criterion': '_0_b_w',
        'v_n_x': 't_tot_flux', 'v_n_y': 'mass_tot_flux',
        'color': 'red', 'ls': '--', 'lw': 1.,
        # 'ymin': 1e-4, 'ymax': 1e0,
        'xlabel': r"time [ms]", 'ylabel': r'$M_{ej}$ $[10^{-2}M_{\odot}]$',
        'xunits': 'ms',
        'label': None, 'yscale': None, 'title': o_data.sim.replace('_', '\_'),
        'fancyticks': True, 'minorticks': True
    }

    o_plot.set_plot_dics = [dic_ej_prof, dic_ej_prof_b,
                            dic_ej_prof2, dic_ej_prof_b2,
                            dic_ej_prof3, dic_ej_prof_b3]
    o_plot.main()
# plot_ej_profile_3sim()


def plot_3hists_and_ej_prof_for_3sims():

    o_plot = PLOT_MANY_TASKS()
    o_plot.gen_set["figdir"] = Paths.plots + '/combine_test/'
    o_plot.gen_set["type"] = "cartesian"
    o_plot.gen_set["figsize"] = (8.5, 9.5)  # <->, |]
    o_plot.gen_set["figname"] = "hists_ej_profiles.png"
    o_plot.gen_set["sharex"] = False
    o_plot.gen_set["sharey"] = True
    o_plot.gen_set["subplots_adjust_h"] = 0.3
    o_plot.gen_set["subplots_adjust_w"] = 0.4

    ''' --- --- sim 1 --- --- --- '''

    sim = "DD2_M13641364_M0_SR"
    o_data = COMPUTE_STORE_PAR(sim)

    dic_ej_prof = {
        'task': 'ejprof', 'ptype': 'cartesian',
        'position': (1, 1),
        'data': o_data, 'criterion': '_0',
        'v_n_x': 't_tot_flux', 'v_n_y': 'mass_tot_flux',
        'color': 'black', 'ls': '-', 'lw': 1.,
        # 'ymin': 1e-4, 'ymax': 1e0,
        'xlabel': r"time [ms]", 'ylabel': r'$M_{ej}$ $[10^{-2}M_{\odot}]$',
        'xunits': 'ms', 'ymin':1e-4, 'ymax':1e-2,
        'label': None, 'yscale': 'log', 'title': o_data.sim.replace('_', '\_'),
        'fancyticks': True, 'minorticks': True
    }
    dic_ej_prof_b = {
        'task': 'ejprof', 'ptype': 'cartesian',
        'position': (1, 1),
        'data': o_data, 'criterion': '_0_b_w',
        'v_n_x': 't_tot_flux', 'v_n_y': 'mass_tot_flux',
        'color': 'black', 'ls': '--', 'lw': 1.,
        # 'ymin': 1e-4, 'ymax': 1e0,
        'xlabel': r"time [ms]", 'ylabel': r'$M_{ej}$ $[10^{-2}M_{\odot}]$',
        'xunits': 'ms','ymin':1e-4, 'ymax':1e-2,
        'label': None, 'yscale': 'log',
        'fancyticks': True, 'minorticks': True
    }


    dic_hist_theta = {
        'task': 'hist1d', 'ptype': 'cartesian',
        'position': (2, 1),
        'data': o_data, 'norm': True, 'criterion': '_0',
        'v_n_x': 'hist_theta', 'v_n_y': 'hist_theta_m',
        'color': 'black', 'ls': '-', 'lw': 1., 'ds': 'steps',
        'ymin': 1e-4, 'ymax': 1e0,
        'xlabel': r"Angle from orbital plane", 'ylabel': 'mass',
        'label': None, 'yscale': 'log',
        'fancyticks': True
    }
    dic_hist_theta_b = {
        'task': 'hist1d', 'ptype': 'cartesian',
        'position': (2, 1),
        'data': o_data, 'norm': True, 'criterion': '_0_b_w',
        'v_n_x': 'hist_theta', 'v_n_y': 'hist_theta_m',
        'color': 'black', 'ls': '--', 'lw': 1., 'ds': 'steps',
        'ymin': 1e-4, 'ymax': 1e0,
        'xlabel': r"Angle from orbital plane", 'ylabel': 'mass',
        'label': None, 'yscale': 'log',
        'fancyticks': True
    }

    dic_hist_vel_inf = {
        'task': 'hist1d', 'ptype': 'cartesian',
        'position': (3, 1),
        'data': o_data, 'norm': True, 'criterion': '_0',
        'v_n_x': 'hist_vel_inf', 'v_n_y': 'hist_vel_inf_m',
        'color': 'black', 'ls': '-', 'lw': 1., 'ds': 'steps',
        'xmin': 0., 'xmax': 1., 'ymin': 1e-4, 'ymax': 1e0,
        'xlabel': r"$\upsilon_{\infty}$ [c]", 'ylabel': 'mass',
        'label': None, 'yscale': 'log',
        'fancyticks': True
    }
    dic_hist_vel_inf_b = {
        'task': 'hist1d', 'ptype': 'cartesian',
        'position': (3, 1),
        'data': o_data, 'norm': True, 'criterion': '_0_b_w',
        'v_n_x': 'hist_vel_inf', 'v_n_y': 'hist_vel_inf_m',
        'color': 'black', 'ls': '--', 'lw': 1., 'ds': 'steps',
        'xmin': 0., 'xmax': 1., 'ymin': 1e-4, 'ymax': 1e0,
        'xlabel': r"$\upsilon_{\infty}$ [c]", 'ylabel': 'mass',
        'label': None, 'yscale': 'log',
        'fancyticks': True
    }

    dic_hist_ye = {
        'task': 'hist1d', 'ptype': 'cartesian',
        'position': (4, 1),
        'data': o_data, 'norm': True, 'criterion': '_0',
        'v_n_x': 'hist_ye', 'v_n_y': 'hist_ye_m',
        'color': 'black', 'ls': '-', 'lw': 1., 'ds': 'steps',
        'xmin': 0., 'xmax': 0.5, 'ymin': 1e-4, 'ymax': 1e0,
        'xlabel': r"$Y_e$", 'ylabel': 'mass',
        'label': None, 'yscale': 'log',
        'fancyticks': True
    }
    dic_hist_ye_b = {
        'task': 'hist1d', 'ptype': 'cartesian',
        'position': (4, 1),
        'data': o_data, 'norm': True, 'criterion': '_0_b_w',
        'v_n_x': 'hist_ye', 'v_n_y': 'hist_ye_m',
        'color': 'black', 'ls': '--', 'lw': 1., 'ds': 'steps',
        'xmin': 0., 'xmax': 0.5, 'ymin': 1e-4, 'ymax': 1e0,
        'xlabel': r"$Y_e$", 'ylabel': 'mass',
        'label': None, 'yscale': 'log',
        'fancyticks': True
    }
    ''' --- --- sim 2 --- --- --- '''

    sim = "LS220_M13641364_M0_SR"
    o_data2 = COMPUTE_STORE_PAR(sim)

    dic_ej_prof2 = {
        'task': 'ejprof', 'ptype': 'cartesian',
        'position': (1, 2),
        'data': o_data2, 'criterion': '_0',
        'v_n_x': 't_tot_flux', 'v_n_y': 'mass_tot_flux',
        'color': 'blue', 'ls': '-', 'lw': 1.,
        # 'ymin': 1e-4, 'ymax': 1e0,
        'xlabel': r"time [ms]", 'ylabel': r'$M_{ej}$ $[10^{-2}M_{\odot}]$',
        'xunits': 'ms','ymin':1e-4, 'ymax':1e-2,
        'label': None, 'yscale': 'log', 'title': o_data2.sim.replace('_', '\_'),
        'fancyticks': True, 'minorticks': True
    }
    dic_ej_prof_b2 = {
        'task': 'ejprof', 'ptype': 'cartesian',
        'position': (1, 2),
        'data': o_data2, 'criterion': '_0_b_w',
        'v_n_x': 't_tot_flux', 'v_n_y': 'mass_tot_flux',
        'color': 'blue', 'ls': '--', 'lw': 1.,
        # 'ymin': 1e-4, 'ymax': 1e0,
        'xlabel': r"time [ms]", 'ylabel': r'$M_{ej}$ $[10^{-2}M_{\odot}]$',
        'xunits': 'ms','ymin':1e-4, 'ymax':1e-2,
        'label': None, 'yscale': 'log',
        'fancyticks': True, 'minorticks': True
    }


    dic_hist_theta2 = {
        'task': 'hist1d', 'ptype': 'cartesian',
        'position': (2, 2),
        'data': o_data2, 'norm': True, 'criterion': '_0',
        'v_n_x': 'hist_theta', 'v_n_y': 'hist_theta_m',
        'color': 'blue', 'ls': '-', 'lw': 1., 'ds': 'steps',
        'ymin': 1e-4, 'ymax': 1e0,
        'xlabel': r"Angle from orbital plane", 'ylabel': 'mass',
        'label': None, 'yscale': 'log',
        'fancyticks': True
    }
    dic_hist_theta_b2 = {
        'task': 'hist1d', 'ptype': 'cartesian',
        'position': (2, 2),
        'data': o_data2, 'norm': True, 'criterion': '_0_b_w',
        'v_n_x': 'hist_theta', 'v_n_y': 'hist_theta_m',
        'color': 'blue', 'ls': '--', 'lw': 1., 'ds': 'steps',
        'ymin': 1e-4, 'ymax': 1e0,
        'xlabel': r"Angle from orbital plane", 'ylabel': 'mass',
        'label': None, 'yscale': 'log',
        'fancyticks': True
    }

    dic_hist_vel_inf2 = {
        'task': 'hist1d', 'ptype': 'cartesian',
        'position': (3, 2),
        'data': o_data2, 'norm': True, 'criterion': '_0',
        'v_n_x': 'hist_vel_inf', 'v_n_y': 'hist_vel_inf_m',
        'color': 'blue', 'ls': '-', 'lw': 1., 'ds': 'steps',
        'xmin': 0., 'xmax': 1., 'ymin': 1e-4, 'ymax': 1e0,
        'xlabel': r"$\upsilon_{\infty}$ [c]", 'ylabel': 'mass',
        'label': None, 'yscale': 'log',
        'fancyticks': True
    }
    dic_hist_vel_inf_b2 = {
        'task': 'hist1d', 'ptype': 'cartesian',
        'position': (3, 2),
        'data': o_data2, 'norm': True, 'criterion': '_0_b_w',
        'v_n_x': 'hist_vel_inf', 'v_n_y': 'hist_vel_inf_m',
        'color': 'blue', 'ls': '--', 'lw': 1., 'ds': 'steps',
        'xmin': 0., 'xmax': 1., 'ymin': 1e-4, 'ymax': 1e0,
        'xlabel': r"$\upsilon_{\infty}$ [c]", 'ylabel': 'mass',
        'label': None, 'yscale': 'log',
        'fancyticks': True
    }

    dic_hist_ye2 = {
        'task': 'hist1d', 'ptype': 'cartesian',
        'position': (4, 2),
        'data': o_data2, 'norm': True, 'criterion': '_0',
        'v_n_x': 'hist_ye', 'v_n_y': 'hist_ye_m',
        'color': 'blue', 'ls': '-', 'lw': 1., 'ds': 'steps',
        'xmin': 0., 'xmax': 0.5, 'ymin': 1e-4, 'ymax': 1e0,
        'xlabel': r"$Y_e$", 'ylabel': 'mass',
        'label': None, 'yscale': 'log',
        'fancyticks': True
    }
    dic_hist_ye_b2 = {
        'task': 'hist1d', 'ptype': 'cartesian',
        'position': (4, 2),
        'data': o_data2, 'norm': True, 'criterion': '_0_b_w',
        'v_n_x': 'hist_ye', 'v_n_y': 'hist_ye_m',
        'color': 'blue', 'ls': '--', 'lw': 1., 'ds': 'steps',
        'xmin': 0., 'xmax': 0.5, 'ymin': 1e-4, 'ymax': 1e0,
        'xlabel': r"$Y_e$", 'ylabel': 'mass',
        'label': None, 'yscale': 'log',
        'fancyticks': True
    }
    ''' --- --- sim 3 --- --- --- '''

    sim = "SLy4_M13641364_M0_SR"
    o_data3 = COMPUTE_STORE_PAR(sim)

    dic_ej_prof3 = {
        'task': 'ejprof', 'ptype': 'cartesian',
        'position': (1, 3),
        'data': o_data3, 'criterion': '_0',
        'v_n_x': 't_tot_flux', 'v_n_y': 'mass_tot_flux',
        'color': 'red', 'ls': '-', 'lw': 1.,
        # 'ymin': 1e-4, 'ymax': 1e0,
        'xlabel': r"time [ms]", 'ylabel': r'$M_{ej}$ $[10^{-2}M_{\odot}]$',
        'xunits': 'ms','ymin':1e-4, 'ymax':1e-2,
        'label': None, 'yscale': 'log', 'title': o_data3.sim.replace('_', '\_'),
        'fancyticks': True, 'minorticks': True
    }
    dic_ej_prof_b3 = {
        'task': 'ejprof', 'ptype': 'cartesian',
        'position': (1, 3),
        'data': o_data3, 'criterion': '_0_b_w',
        'v_n_x': 't_tot_flux', 'v_n_y': 'mass_tot_flux',
        'color': 'red', 'ls': '--', 'lw': 1.,
        # 'ymin': 1e-4, 'ymax': 1e0,
        'xlabel': r"time [ms]", 'ylabel': r'$M_{ej}$ $[10^{-2}M_{\odot}]$',
        'xunits': 'ms','ymin':1e-4, 'ymax':1e-2,
        'label': None, 'yscale': 'log', 'title': o_data3.sim.replace('_', '\_'),
        'fancyticks': True, 'minorticks': True
    }


    dic_hist_theta3 = {
        'task': 'hist1d', 'ptype': 'cartesian',
        'position': (2, 3),
        'data': o_data3, 'norm': True, 'criterion': '_0',
        'v_n_x': 'hist_theta', 'v_n_y': 'hist_theta_m',
        'color': 'red', 'ls': '-', 'lw': 1., 'ds': 'steps',
        'ymin': 1e-4, 'ymax': 1e0,
        'xlabel': r"Angle from orbital plane", 'ylabel': 'mass',
        'label': None, 'yscale': 'log',
        'fancyticks': True
    }
    dic_hist_theta_b3 = {
        'task': 'hist1d', 'ptype': 'cartesian',
        'position': (2, 3),
        'data': o_data3, 'norm': True, 'criterion': '_0_b_w',
        'v_n_x': 'hist_theta', 'v_n_y': 'hist_theta_m',
        'color': 'red', 'ls': '--', 'lw': 1., 'ds': 'steps',
        'ymin': 1e-4, 'ymax': 1e0,
        'xlabel': r"Angle from orbital plane", 'ylabel': 'mass',
        'label': None, 'yscale': 'log',
        'fancyticks': True
    }

    dic_hist_vel_inf3 = {
        'task': 'hist1d', 'ptype': 'cartesian',
        'position': (3, 3),
        'data': o_data3, 'norm': True, 'criterion': '_0',
        'v_n_x': 'hist_vel_inf', 'v_n_y': 'hist_vel_inf_m',
        'color': 'red', 'ls': '-', 'lw': 1., 'ds': 'steps',
        'xmin': 0., 'xmax': 1., 'ymin': 1e-4, 'ymax': 1e0,
        'xlabel': r"$\upsilon_{\infty}$ [c]", 'ylabel': 'mass',
        'label': None, 'yscale': 'log',
        'fancyticks': True
    }
    dic_hist_vel_inf_b3 = {
        'task': 'hist1d', 'ptype': 'cartesian',
        'position': (3, 3),
        'data': o_data3, 'norm': True, 'criterion': '_0_b_w',
        'v_n_x': 'hist_vel_inf', 'v_n_y': 'hist_vel_inf_m',
        'color': 'red', 'ls': '--', 'lw': 1., 'ds': 'steps',
        'xmin': 0., 'xmax': 1., 'ymin': 1e-4, 'ymax': 1e0,
        'xlabel': r"$\upsilon_{\infty}$ [c]", 'ylabel': 'mass',
        'label': None, 'yscale': 'log',
        'fancyticks': True
    }

    dic_hist_ye3 = {
        'task': 'hist1d', 'ptype': 'cartesian',
        'position': (4, 3),
        'data': o_data3, 'norm': True, 'criterion': '_0',
        'v_n_x': 'hist_ye', 'v_n_y': 'hist_ye_m',
        'color': 'red', 'ls': '-', 'lw': 1., 'ds': 'steps',
        'xmin': 0., 'xmax': 0.5, 'ymin': 1e-4, 'ymax': 1e0,
        'xlabel': r"$Y_e$", 'ylabel': 'mass',
        'label': None, 'yscale': 'log',
        'fancyticks': True
    }
    dic_hist_ye_b3 = {
        'task': 'hist1d', 'ptype': 'cartesian',
        'position': (4, 3),
        'data': o_data3, 'norm': True, 'criterion': '_0_b_w',
        'v_n_x': 'hist_ye', 'v_n_y': 'hist_ye_m',
        'color': 'red', 'ls': '--', 'lw': 1., 'ds': 'steps',
        'xmin': 0., 'xmax': 0.5, 'ymin': 1e-4, 'ymax': 1e0,
        'xlabel': r"$Y_e$", 'ylabel': 'mass',
        'label': None, 'yscale': 'log',
        'fancyticks': True
    }

    o_plot.set_plot_dics = [
        dic_ej_prof, dic_ej_prof_b,
        dic_hist_theta, dic_hist_theta_b,
        dic_hist_vel_inf, dic_hist_vel_inf_b,
        dic_hist_ye, dic_hist_ye_b,

        dic_ej_prof2, dic_ej_prof_b2,
        dic_hist_theta2, dic_hist_theta_b2,
        dic_hist_vel_inf2, dic_hist_vel_inf_b2,
        dic_hist_ye2, dic_hist_ye_b2,

        dic_ej_prof3, dic_ej_prof_b3,
        dic_hist_theta3, dic_hist_theta_b3,
        dic_hist_vel_inf3, dic_hist_vel_inf_b3,
        dic_hist_ye3, dic_hist_ye_b3,
    ]
    o_plot.main()
# plot_3hists_and_ej_prof_for_3sims()




def plot_nucleo():

    o_plot = PLOT_MANY_TASKS()
    o_plot.gen_set["figdir"] = Paths.plots + '/combine_test/'
    o_plot.gen_set["type"] = "cartesian"
    o_plot.gen_set["figsize"] = (3.5, 3.5)  # <->, |]
    o_plot.gen_set["figname"] = "nucleo.png"
    o_plot.gen_set["sharex"] = False
    o_plot.gen_set["sharey"] = False
    o_plot.gen_set["subplots_adjust_h"] = 0.3
    o_plot.gen_set["subplots_adjust_w"] = 0.4

    sim = "DD2_M13641364_M0_SR"
    o_data = NORMALIZE_NUCLEO(sim)

    sim_nucleo = {
        'task': 'nucleo', 'ptype': 'cartesian',
        'position': (1, 1),
        'data': o_data, 'criterion': '_0', 'method': 'sum',
        'v_n_x': 'A', 'v_n_y': 'Y_final',
        'color': 'green', 'ls': '--', 'lw': 0.8, 'ds': 'steps', 'alpha': 0.6,
        'ymin': 1e-5, 'ymax': 2e-1, 'xmin': 50, 'xmax': 210,
        'xlabel': r"A", 'ylabel': r'Relative final abundances',
        'label': None, 'yscale': 'log',
        'fancyticks': True, 'minorticks': True
    }

    sim_nucleo_b = {
        'task': 'nucleo', 'ptype': 'cartesian',
        'position': (1, 1),
        'data': o_data, 'criterion': '_0_b_w', 'method': 'sum',
        'v_n_x': 'A', 'v_n_y': 'Y_final',
        'color': 'green', 'ls': '-', 'lw': 1.0, 'ds': 'steps',
        'ymin': 1e-5, 'ymax': 2e-1, 'xmin': 50, 'xmax': 210,
        'xlabel': r"A", 'ylabel': r'Relative final abundances',
        'label': None, 'yscale': 'log',
        'fancyticks': True, 'minorticks': True
    }


    sol_yeilds = {
        'task': 'nucleo', 'ptype': 'cartesian',
        'position': (1, 1),
        'data': o_data, 'criterion': '_0', 'method': 'sum',
        'v_n_x': 'Asun', 'v_n_y': 'Ysun',
        'color': 'gray', 'marker': 'o', 'ms':4,
        'ymin': 1e-5, 'ymax': 2e-1, 'xmin': 50, 'xmax': 210,
        'xlabel': r"A", 'ylabel': r'Relative final abundances',
        'label': None, 'yscale': 'log',
        'fancyticks': True, 'minorticks': True
    }

    o_plot.set_plot_dics = [sim_nucleo, sim_nucleo_b, sol_yeilds]
    o_plot.main()
# plot_nucleo()

def plot_nucleo_3sim():
    o_plot = PLOT_MANY_TASKS()
    o_plot.gen_set["figdir"] = Paths.plots + '/combine_test/'
    o_plot.gen_set["type"] = "cartesian"
    o_plot.gen_set["figsize"] = (9.0, 2.7)  # <->, |]
    o_plot.gen_set["figname"] = "nucleo_3sims.png"
    o_plot.gen_set["sharex"] = False
    o_plot.gen_set["sharey"] = False
    o_plot.gen_set["subplots_adjust_h"] = 0.3
    o_plot.gen_set["subplots_adjust_w"] = 0.0

    ''' --- --- -- SIM 1 --- --- --- '''

    sim = "DD2_M13641364_M0_SR"
    o_data = NORMALIZE_NUCLEO(sim)

    sim_nucleo = {
        'task': 'nucleo', 'ptype': 'cartesian',
        'position': (1, 1),
        'data': o_data, 'criterion': '_0', 'method': 'sum',
        'v_n_x': 'A', 'v_n_y': 'Y_final',
        'color': 'green', 'ls': '--', 'lw': 0.8, 'ds': 'steps', 'alpha': 0.6,
        'ymin': 1e-5, 'ymax': 2e-1, 'xmin': 50, 'xmax': 210,
        'xlabel': r"A", 'ylabel': r'Relative final abundances',
        'label': None, 'yscale': 'log',
        'fancyticks': True, 'minorticks': True
    }

    sim_nucleo_b = {
        'task': 'nucleo', 'ptype': 'cartesian',
        'position': (1, 1),
        'data': o_data, 'criterion': '_0_b_w', 'method': 'sum',
        'v_n_x': 'A', 'v_n_y': 'Y_final',
        'color': 'green', 'ls': '-', 'lw': 1.0, 'ds': 'steps',
        'ymin': 1e-5, 'ymax': 2e-1, 'xmin': 50, 'xmax': 210,
        'xlabel': r"A", 'ylabel': r'Relative final abundances',
        'label': None, 'yscale': 'log',
        'fancyticks': True, 'minorticks': True
    }

    sol_yeilds = {
        'task': 'nucleo', 'ptype': 'cartesian',
        'position': (1, 1),
        'data': o_data, 'criterion': '_0', 'method': 'sum',
        'v_n_x': 'Asun', 'v_n_y': 'Ysun',
        'color': 'gray', 'marker': 'o', 'ms': 4, 'alpha':0.4,
        'ymin': 1e-5, 'ymax': 2e-1, 'xmin': 50, 'xmax': 210,
        'xlabel': r"A", 'ylabel': r'Relative final abundances',
        'label': None, 'yscale': 'log',
        'fancyticks': True, 'minorticks': True
    }

    ''' --- --- -- SIM 2 --- --- --- '''

    sim = "LS220_M13641364_M0_SR"
    o_data = NORMALIZE_NUCLEO(sim)

    sim_nucleo2 = {
        'task': 'nucleo', 'ptype': 'cartesian',
        'position': (1, 2),
        'data': o_data, 'criterion': '_0', 'method': 'sum',
        'v_n_x': 'A', 'v_n_y': 'Y_final',
        'color': 'blue', 'ls': '--', 'lw': 0.8, 'ds': 'steps', 'alpha': 0.6,
        'ymin': 1e-5, 'ymax': 2e-1, 'xmin': 50, 'xmax': 210,
        'xlabel': r"A", 'ylabel': r'Relative final abundances',
        'label': None, 'yscale': 'log',
        'fancyticks': True, 'minorticks': True
    }

    sim_nucleo_b2 = {
        'task': 'nucleo', 'ptype': 'cartesian',
        'position': (1, 2),
        'data': o_data, 'criterion': '_0_b_w', 'method': 'sum',
        'v_n_x': 'A', 'v_n_y': 'Y_final',
        'color': 'blue', 'ls': '-', 'lw': 1.0, 'ds': 'steps',
        'ymin': 1e-5, 'ymax': 2e-1, 'xmin': 50, 'xmax': 210,
        'xlabel': r"A", 'ylabel': r'Relative final abundances',
        'label': None, 'yscale': 'log',
        'fancyticks': True, 'minorticks': True
    }

    sol_yeilds2 = {
        'task': 'nucleo', 'ptype': 'cartesian',
        'position': (1, 2),
        'data': o_data, 'criterion': '_0', 'method': 'sum',
        'v_n_x': 'Asun', 'v_n_y': 'Ysun',
        'color': 'gray', 'marker': 'o', 'ms': 4, 'alpha':0.4,
        'ymin': 1e-5, 'ymax': 2e-1, 'xmin': 50, 'xmax': 210,
        'xlabel': r"A", 'ylabel': r'Relative final abundances',
        'label': None, 'yscale': 'log',
        'fancyticks': True, 'minorticks': True,
        'sharey': True
    }

    ''' --- --- -- SIM 3 --- --- --- '''

    sim = "SLy4_M13641364_M0_SR"
    o_data = NORMALIZE_NUCLEO(sim)

    sim_nucleo3 = {
        'task': 'nucleo', 'ptype': 'cartesian',
        'position': (1, 3),
        'data': o_data, 'criterion': '_0', 'method': 'sum',
        'v_n_x': 'A', 'v_n_y': 'Y_final',
        'color': 'red', 'ls': '--', 'lw': 0.8, 'ds': 'steps', 'alpha': 0.6,
        'ymin': 1e-5, 'ymax': 2e-1, 'xmin': 50, 'xmax': 210,
        'xlabel': r"A", 'ylabel': r'Relative final abundances',
        'label': None, 'yscale': 'log',
        'fancyticks': True, 'minorticks': True
    }

    sim_nucleo_b3 = {
        'task': 'nucleo', 'ptype': 'cartesian',
        'position': (1, 3),
        'data': o_data, 'criterion': '_0_b_w', 'method': 'sum',
        'v_n_x': 'A', 'v_n_y': 'Y_final',
        'color': 'red', 'ls': '-', 'lw': 1.0, 'ds': 'steps',
        'ymin': 1e-5, 'ymax': 2e-1, 'xmin': 50, 'xmax': 210,
        'xlabel': r"A", 'ylabel': r'Relative final abundances',
        'label': None, 'yscale': 'log',
        'fancyticks': True, 'minorticks': True
    }

    sol_yeilds3 = {
        'task': 'nucleo', 'ptype': 'cartesian',
        'position': (1, 3),
        'data': o_data, 'criterion': '_0', 'method': 'sum',
        'v_n_x': 'Asun', 'v_n_y': 'Ysun',
        'color': 'gray', 'marker': 'o', 'ms': 4, 'alpha':0.4,
        'ymin': 1e-5, 'ymax': 2e-1, 'xmin': 50, 'xmax': 210,
        'xlabel': r"A", 'ylabel': r'Relative final abundances',
        'label': None, 'yscale': 'log',
        'fancyticks': True, 'minorticks': True,
        'sharey': True
    }

    o_plot.set_plot_dics = [sim_nucleo, sim_nucleo_b, sol_yeilds,
                            sim_nucleo2, sim_nucleo_b2, sol_yeilds2,
                            sim_nucleo3, sim_nucleo_b3, sol_yeilds3]
    o_plot.main()
# plot_nucleo_3sim()

def plot_nucleo_3sim_one_plot():
    o_plot = PLOT_MANY_TASKS()
    o_plot.gen_set["figdir"] = Paths.plots + '/combine_test/'
    o_plot.gen_set["type"] = "cartesian"
    o_plot.gen_set["figsize"] = (4.6, 3.6)  # <->, |]
    o_plot.gen_set["figname"] = "nucleo_3sims_1plot.png"
    o_plot.gen_set["sharex"] = False
    o_plot.gen_set["sharey"] = False
    o_plot.gen_set["subplots_adjust_h"] = 0.3
    o_plot.gen_set["subplots_adjust_w"] = 0.0
    o_plot.set_plot_dics = []

    ''' --- --- -- SIM 1 --- --- --- '''

    sim = "DD2_M13641364_M0_SR"
    o_data = NORMALIZE_NUCLEO(sim)

    sim_nucleo = {
        'task': 'nucleo', 'ptype': 'cartesian',
        'position': (1, 1),
        'data': o_data, 'criterion': '_0', 'method': 'sum',
        'v_n_x': 'A', 'v_n_y': 'Y_final',
        'color': 'green', 'ls': ':', 'lw': 0.6, 'ds': 'steps', 'alpha': 0.6,
        'ymin': 1e-5, 'ymax': 2e-1, 'xmin': 50, 'xmax': 210,
        'xlabel': r"A", 'ylabel': r'Relative final abundances',
        'label': None, 'yscale': 'log',
        'fancyticks': True, 'minorticks': True
    }

    sim_nucleo_b = {
        'task': 'nucleo', 'ptype': 'cartesian',
        'position': (1, 1),
        'data': o_data, 'criterion': '_0_b_w', 'method': 'sum',
        'v_n_x': 'A', 'v_n_y': 'Y_final',
        'color': 'green', 'ls': '-', 'lw': 1.0, 'ds': 'steps',
        'ymin': 1e-5, 'ymax': 2e-1, 'xmin': 50, 'xmax': 210,
        'xlabel': r"A", 'ylabel': r'Relative final abundances',
        'label': "DD2", 'yscale': 'log',
        'fancyticks': True, 'minorticks': True
    }

    sol_yeilds = {
        'task': 'nucleo', 'ptype': 'cartesian',
        'position': (1, 1),
        'data': o_data, 'criterion': '_0', 'method': 'sum',
        'v_n_x': 'Asun', 'v_n_y': 'Ysun',
        'color': 'gray', 'marker': 'o', 'ms': 4, 'alpha':0.4,
        'ymin': 1e-5, 'ymax': 2e-1, 'xmin': 50, 'xmax': 210,
        'xlabel': r"A", 'ylabel': r'Relative final abundances',
        'label': None, 'yscale': 'log',
        'fancyticks': True, 'minorticks': True
    }

    o_plot.set_plot_dics.append(sim_nucleo)
    o_plot.set_plot_dics.append(sim_nucleo_b)
    o_plot.set_plot_dics.append(sol_yeilds)

    ''' --- --- -- SIM 2 --- --- --- '''

    sim = "LS220_M13641364_M0_SR"
    o_data = NORMALIZE_NUCLEO(sim)

    sim_nucleo = {
        'task': 'nucleo', 'ptype': 'cartesian',
        'position': (1, 1),
        'data': o_data, 'criterion': '_0', 'method': 'sum',
        'v_n_x': 'A', 'v_n_y': 'Y_final',
        'color': 'blue', 'ls': ':', 'lw': 0.6, 'ds': 'steps', 'alpha': 0.6,
        'ymin': 1e-5, 'ymax': 2e-1, 'xmin': 50, 'xmax': 210,
        'xlabel': r"A", 'ylabel': r'Relative final abundances',
        'label': None, 'yscale': 'log',
        'fancyticks': True, 'minorticks': True
    }

    sim_nucleo_b = {
        'task': 'nucleo', 'ptype': 'cartesian',
        'position': (1, 1),
        'data': o_data, 'criterion': '_0_b_w', 'method': 'sum',
        'v_n_x': 'A', 'v_n_y': 'Y_final',
        'color': 'blue', 'ls': '-', 'lw': 1.0, 'ds': 'steps',
        'ymin': 1e-5, 'ymax': 2e-1, 'xmin': 50, 'xmax': 210,
        'xlabel': r"A", 'ylabel': r'Relative final abundances',
        'label': "LS220", 'yscale': 'log',
        'fancyticks': True, 'minorticks': True
    }

    o_plot.set_plot_dics.append(sim_nucleo)
    o_plot.set_plot_dics.append(sim_nucleo_b)

    ''' --- --- -- SIM 3 --- --- --- '''

    sim = "SLy4_M13641364_M0_SR"
    o_data = NORMALIZE_NUCLEO(sim)

    sim_nucleo = {
        'task': 'nucleo', 'ptype': 'cartesian',
        'position': (1, 1),
        'data': o_data, 'criterion': '_0', 'method': 'sum',
        'v_n_x': 'A', 'v_n_y': 'Y_final',
        'color': 'red', 'ls': ':', 'lw': 0.6, 'ds': 'steps', 'alpha': 0.6,
        'ymin': 1e-5, 'ymax': 2e-1, 'xmin': 50, 'xmax': 210,
        'xlabel': r"A", 'ylabel': r'Relative final abundances',
        'label': None, 'yscale': 'log',
        'fancyticks': True, 'minorticks': True
    }

    sim_nucleo_b = {
        'task': 'nucleo', 'ptype': 'cartesian',
        'position': (1, 1),
        'data': o_data, 'criterion': '_0_b_w', 'method': 'sum',
        'v_n_x': 'A', 'v_n_y': 'Y_final',
        'color': 'red', 'ls': '-', 'lw': 1.0, 'ds': 'steps',
        'ymin': 1e-5, 'ymax': 2e-1, 'xmin': 50, 'xmax': 210,
        'xlabel': r"A", 'ylabel': r'Relative final abundances',
        'label': "SLy4", 'yscale': 'log',
        'fancyticks': True, 'minorticks': True,
        'legend': True
    }

    o_plot.set_plot_dics.append(sim_nucleo)
    o_plot.set_plot_dics.append(sim_nucleo_b)

    o_plot.main()
# plot_nucleo_3sim_one_plot()

""" --- --- LIGHTCURVES --- --- --- """

def plot_mkn_3sims_multiple_component_arrangment():


    o_plot = PLOT_MANY_TASKS()
    o_plot.gen_set["figdir"] = Paths.plots + '/combine_test/'
    o_plot.gen_set["type"] = "cartesian"
    o_plot.gen_set["figsize"] = (8.0, 6.0)  # <->, |]
    o_plot.gen_set["figname"] = "mkn_mult_comp_anal.png"
    o_plot.gen_set["sharex"] = False
    o_plot.gen_set["sharey"] = False
    o_plot.gen_set["invert_y"] = False
    o_plot.gen_set["subplots_adjust_h"] = 0.0
    o_plot.gen_set["subplots_adjust_w"] = 0.0
    o_plot.set_plot_dics = []

    """ ===================== SIM 1 ============================"""

    sim = "DD2_M13641364_M0_SR"
    o_lc_data = EXTRACT_LIGHTCURVE(sim)

    def dd2_g(o_lc_data):

        obs_g = {
            'task': 'mkn obs', "ptype": "cartesian",
            'position': (1, 1),
            'data': o_lc_data, 'band': 'g', 'obs': True,
            'v_n_x': 'time', 'v_n_y': 'mag',
            'color': 'gray', 'marker': 'o', 'ms': 5., 'alpha': 0.2,
            'ymin': 25, 'ymax': 15,  'xmin': 3e-1, 'xmax': 3e1,
            'xlabel': r"time [days]", 'ylabel': r"AB magnitude at 40 Mpc",
            'label': None, 'xscale': 'log',
            'fancyticks': True, 'minorticks': True,
            'legend': True, 'title': 'g band',
            'sharex': True
        }
        o_plot.set_plot_dics.append(obs_g)

        dd2_g_1 = {
            'task': 'mkn model', "ptype": "cartesian",
            'position': (1, 1),
            'data': o_lc_data, 'band': 'g', 'fname':'mkn_model1.h5',
            'obs': True,
            'v_n_x': 'time', 'v_n_y': 'mag',
            'color': 'green', 'alpha': 0.1,
            'ymin': 25, 'ymax': 15, 'xmin': 3e-1, 'xmax': 3e1,
            'xlabel': r"time [days]", 'ylabel': r"AB magnitude at 40 Mpc",
            'label': '1 comp', 'xscale': 'log',
            'fancyticks': True, 'minorticks': True,
            'sharex':True
        }
        o_plot.set_plot_dics.append(dd2_g_1)

        dd2_g_2 = copy.deepcopy(dd2_g_1)
        dd2_g_2['fname'] = 'mkn_model2.h5'
        dd2_g_2['alpha'] = 0.3
        dd2_g_2['label'] = '2 comp'
        o_plot.set_plot_dics.append(dd2_g_2)

        dd2_g_3 = copy.deepcopy(dd2_g_1)
        dd2_g_3['fname'] = 'mkn_model3.h5'
        dd2_g_3['alpha'] = 0.6
        dd2_g_3['label'] = '3 comp'
        o_plot.set_plot_dics.append(dd2_g_3)

        dd2_g_4 = copy.deepcopy(dd2_g_1)
        dd2_g_4['fname'] = 'mkn_model4.h5'
        dd2_g_4['alpha'] = 0.9
        dd2_g_4['label'] = '4 comp'
        dd2_g_4['legend'] = True
        o_plot.set_plot_dics.append(dd2_g_4)
    dd2_g(o_lc_data)

    def dd2_z(o_lc_data):
        obs = {
            'task': 'mkn obs', "ptype": "cartesian",
            'position': (1, 2),
            'data': o_lc_data, 'band': 'z', 'obs': True,
            'v_n_x': 'time', 'v_n_y': 'mag',
            'color': 'gray', 'marker': 'o', 'ms': 5., 'alpha': 0.2,
            'ymin': 25, 'ymax': 15, 'xmin': 3e-1, 'xmax': 3e1,
            'xlabel': r"time [days]", 'ylabel': r"AB magnitude at 40 Mpc",
            'label': None, 'xscale': 'log',
            'fancyticks': True, 'minorticks': True,
            'legend': True, 'title': 'z band',
            'sharey': True,
            'sharex': True
        }
        o_plot.set_plot_dics.append(obs)

        dd2_1 = {
            'task': 'mkn model', "ptype": "cartesian",
            'position': (1, 2),
            'data': o_lc_data, 'band': 'z', 'fname': 'mkn_model1.h5',
            'obs': True,
            'v_n_x': 'time', 'v_n_y': 'mag',
            'color': 'green', 'alpha': 0.1,
            'ymin': 25, 'ymax': 15, 'xmin': 3e-1, 'xmax': 3e1,
            'xlabel': r"time [days]", 'ylabel': r"AB magnitude at 40 Mpc",
            'label': '1 comp', 'xscale': 'log',
            'fancyticks': True, 'minorticks': True,
            'sharey': True,
            'sharex': True
        }
        o_plot.set_plot_dics.append(dd2_1)

        dd2_2 = copy.deepcopy(dd2_1)
        dd2_2['fname'] = 'mkn_model2.h5'
        dd2_2['alpha'] = 0.3
        dd2_2['label'] = '2 comp'
        o_plot.set_plot_dics.append(dd2_2)

        dd2_3 = copy.deepcopy(dd2_1)
        dd2_3['fname'] = 'mkn_model3.h5'
        dd2_3['alpha'] = 0.6
        dd2_3['label'] = '3 comp'
        o_plot.set_plot_dics.append(dd2_3)

        dd2_4 = copy.deepcopy(dd2_1)
        dd2_4['fname'] = 'mkn_model4.h5'
        dd2_4['alpha'] = 0.9
        dd2_4['label'] = '4 comp'
        o_plot.set_plot_dics.append(dd2_4)
    dd2_z(o_lc_data)

    def dd2_Ks(o_lc_data):
        obs = {
            'task': 'mkn obs', "ptype": "cartesian",
            'position': (1, 3),
            'data': o_lc_data, 'band': 'Ks', 'obs': True,
            'v_n_x': 'time', 'v_n_y': 'mag',
            'color': 'gray', 'marker': 'o', 'ms': 5., 'alpha': 0.2,
            'ymin': 25, 'ymax': 15, 'xmin': 3e-1, 'xmax': 3e1,
            'xlabel': r"time [days]", 'ylabel': r"AB magnitude at 40 Mpc",
            'label': None, 'xscale': 'log',
            'fancyticks': True, 'minorticks': True,
            'legend': True, 'title': 'Ks band',
            'sharey': True,
            'sharex': True
        }
        o_plot.set_plot_dics.append(obs)

        dd2_1 = {
            'task': 'mkn model', "ptype": "cartesian",
            'position': (1, 3),
            'data': o_lc_data, 'band': 'Ks', 'fname': 'mkn_model1.h5',
            'obs': True,
            'v_n_x': 'time', 'v_n_y': 'mag',
            'color': 'green', 'alpha': 0.1,
            'ymin': 25, 'ymax': 15, 'xmin': 3e-1, 'xmax': 3e1,
            'xlabel': r"time [days]", 'ylabel': r"AB magnitude at 40 Mpc",
            'label': '1 comp', 'xscale': 'log',
            'fancyticks': True, 'minorticks': True,
            'sharey': True,
            'sharex': True
        }
        o_plot.set_plot_dics.append(dd2_1)

        dd2_2 = copy.deepcopy(dd2_1)
        dd2_2['fname'] = 'mkn_model2.h5'
        dd2_2['alpha'] = 0.3
        dd2_2['label'] = '2 comp'
        o_plot.set_plot_dics.append(dd2_2)

        dd2_3 = copy.deepcopy(dd2_1)
        dd2_3['fname'] = 'mkn_model3.h5'
        dd2_3['alpha'] = 0.6
        dd2_3['label'] = '3 comp'
        o_plot.set_plot_dics.append(dd2_3)

        dd2_4 = copy.deepcopy(dd2_1)
        dd2_4['fname'] = 'mkn_model4.h5'
        dd2_4['alpha'] = 0.9
        dd2_4['label'] = '4 comp'
        o_plot.set_plot_dics.append(dd2_4)
    dd2_Ks(o_lc_data)

    """ ===================== SIM 2 ============================"""

    sim = "LS220_M13641364_M0_SR"
    o_lc_data = EXTRACT_LIGHTCURVE(sim)

    def ls220_g(o_lc_data):

        obs_g = {
            'task': 'mkn obs', "ptype": "cartesian",
            'position': (2, 1),
            'data': o_lc_data, 'band': 'g', 'obs': True,
            'v_n_x': 'time', 'v_n_y': 'mag',
            'color': 'gray', 'marker': 'o', 'ms': 5., 'alpha': 0.2,
            'ymin': 25, 'ymax': 15,  'xmin': 3e-1, 'xmax': 3e1,
            'xlabel': r"time [days]", 'ylabel': '',
            'label': None, 'xscale': 'log',
            'fancyticks': True, 'minorticks': True,
            'legend': True, #'title': 'g band',
            # 'sharex': True
        }
        o_plot.set_plot_dics.append(obs_g)

        g_1 = {
            'task': 'mkn model', "ptype": "cartesian",
            'position': (2, 1),
            'data': o_lc_data, 'band': 'g', 'fname':'mkn_model1.h5',
            'obs': True,
            'v_n_x': 'time', 'v_n_y': 'mag',
            'color': 'blue', 'alpha': 0.1,
            'ymin': 25, 'ymax': 15, 'xmin': 3e-1, 'xmax': 3e1,
            'xlabel': r"time [days]", 'ylabel': '',
            'label': '1 comp', 'xscale': 'log',
            'fancyticks': True, 'minorticks': True,
            # 'sharex': True
        }
        o_plot.set_plot_dics.append(g_1)

        g_2 = copy.deepcopy(g_1)
        g_2['fname'] = 'mkn_model2.h5'
        g_2['alpha'] = 0.3
        g_2['label'] = '2 comp'
        o_plot.set_plot_dics.append(g_2)

        g_3 = copy.deepcopy(g_1)
        g_3['fname'] = 'mkn_model3.h5'
        g_3['alpha'] = 0.6
        g_3['label'] = '3 comp'
        o_plot.set_plot_dics.append(g_3)

        g_4 = copy.deepcopy(g_1)
        g_4['fname'] = 'mkn_model4.h5'
        g_4['alpha'] = 0.9
        g_4['label'] = '4 comp'
        g_4['legend'] = True
        o_plot.set_plot_dics.append(g_4)
    ls220_g(o_lc_data)

    def ls220_z(o_lc_data):

        obs_g = {
            'task': 'mkn obs', "ptype": "cartesian",
            'position': (2, 2),
            'data': o_lc_data, 'band': 'z', 'obs': True,
            'v_n_x': 'time', 'v_n_y': 'mag',
            'color': 'gray', 'marker': 'o', 'ms': 5., 'alpha': 0.2,
            'ymin': 25, 'ymax': 15,  'xmin': 3e-1, 'xmax': 3e1,
            'xlabel': r"time [days]", 'ylabel': r"AB magnitude at 40 Mpc",
            'label': None, 'xscale': 'log',
            'fancyticks': True, 'minorticks': True,
            'legend': True, #'title': 'z band',
            'sharey': True,
            # 'sharex': True
        }
        o_plot.set_plot_dics.append(obs_g)

        g_1 = {
            'task': 'mkn model', "ptype": "cartesian",
            'position': (2, 2),
            'data': o_lc_data, 'band': 'z', 'fname':'mkn_model1.h5',
            'obs': True,
            'v_n_x': 'time', 'v_n_y': 'mag',
            'color': 'blue', 'alpha': 0.1,
            'ymin': 25, 'ymax': 15, 'xmin': 3e-1, 'xmax': 3e1,
            'xlabel': r"time [days]", 'ylabel': r"AB magnitude at 40 Mpc",
            'label': '1 comp', 'xscale': 'log',
            'fancyticks': True, 'minorticks': True,
            'sharey': True,
            # 'sharex': True
        }
        o_plot.set_plot_dics.append(g_1)

        g_2 = copy.deepcopy(g_1)
        g_2['fname'] = 'mkn_model2.h5'
        g_2['alpha'] = 0.3
        g_2['label'] = '2 comp'
        o_plot.set_plot_dics.append(g_2)

        g_3 = copy.deepcopy(g_1)
        g_3['fname'] = 'mkn_model3.h5'
        g_3['alpha'] = 0.6
        g_3['label'] = '3 comp'
        o_plot.set_plot_dics.append(g_3)

        g_4 = copy.deepcopy(g_1)
        g_4['fname'] = 'mkn_model4.h5'
        g_4['alpha'] = 0.9
        g_4['label'] = '4 comp'
        # g_4['legend'] = True
        o_plot.set_plot_dics.append(g_4)
    ls220_z(o_lc_data)

    def ls220_Ks(o_lc_data):

        obs_g = {
            'task': 'mkn obs', "ptype": "cartesian",
            'position': (2, 3),
            'data': o_lc_data, 'band': 'Ks', 'obs': True,
            'v_n_x': 'time', 'v_n_y': 'mag',
            'color': 'gray', 'marker': 'o', 'ms': 5., 'alpha': 0.2,
            'ymin': 25, 'ymax': 15,  'xmin': 3e-1, 'xmax': 3e1,
            'xlabel': r"time [days]", 'ylabel': r"AB magnitude at 40 Mpc",
            'label': None, 'xscale': 'log',
            'fancyticks': True, 'minorticks': True,
            'legend': True, #'title': 'Ks band',
            'sharey': True,
            # 'sharex': True
        }
        o_plot.set_plot_dics.append(obs_g)

        g_1 = {
            'task': 'mkn model', "ptype": "cartesian",
            'position': (2, 3),
            'data': o_lc_data, 'band': 'Ks', 'fname':'mkn_model1.h5',
            'obs': True,
            'v_n_x': 'time', 'v_n_y': 'mag',
            'color': 'blue', 'alpha': 0.1,
            'ymin': 25, 'ymax': 15, 'xmin': 3e-1, 'xmax': 3e1,
            'xlabel': r"time [days]", 'ylabel': r"AB magnitude at 40 Mpc",
            'label': '1 comp', 'xscale': 'log',
            'fancyticks': True, 'minorticks': True,
            'sharey': True,
            # 'sharex': True
        }
        o_plot.set_plot_dics.append(g_1)

        g_2 = copy.deepcopy(g_1)
        g_2['fname'] = 'mkn_model2.h5'
        g_2['alpha'] = 0.3
        g_2['label'] = '2 comp'
        o_plot.set_plot_dics.append(g_2)

        g_3 = copy.deepcopy(g_1)
        g_3['fname'] = 'mkn_model3.h5'
        g_3['alpha'] = 0.6
        g_3['label'] = '3 comp'
        o_plot.set_plot_dics.append(g_3)

        g_4 = copy.deepcopy(g_1)
        g_4['fname'] = 'mkn_model4.h5'
        g_4['alpha'] = 0.9
        g_4['label'] = '4 comp'
        # g_4['legend'] = True
        o_plot.set_plot_dics.append(g_4)
    ls220_Ks(o_lc_data)

    """ ===================== SIM 3 ============================"""

    sim = "SLy4_M13641364_M0_SR"
    o_lc_data = EXTRACT_LIGHTCURVE(sim)

    def SLy4_g(o_lc_data):
        obs_g = {
            'task': 'mkn obs', "ptype": "cartesian",
            'position': (3, 1),
            'data': o_lc_data, 'band': 'g', 'obs': True,
            'v_n_x': 'time', 'v_n_y': 'mag',
            'color': 'gray', 'marker': 'o', 'ms': 5., 'alpha': 0.2,
            'ymin': 25, 'ymax': 15, 'xmin': 3e-1, 'xmax': 3e1,
            'xlabel': r"time [days]", 'ylabel': '',
            'label': None, 'xscale': 'log',
            'fancyticks': True, 'minorticks': True,
            'legend': True, #'title': 'g band',
        }
        o_plot.set_plot_dics.append(obs_g)

        g_1 = {
            'task': 'mkn model', "ptype": "cartesian",
            'position': (3, 1),
            'data': o_lc_data, 'band': 'g', 'fname': 'mkn_model1.h5',
            'obs': True,
            'v_n_x': 'time', 'v_n_y': 'mag',
            'color': 'red', 'alpha': 0.1,
            'ymin': 25, 'ymax': 15, 'xmin': 3e-1, 'xmax': 3e1,
            'xlabel': r"time [days]", 'ylabel': '',
            'label': '1 comp', 'xscale': 'log',
            'fancyticks': True, 'minorticks': True,
        }
        o_plot.set_plot_dics.append(g_1)

        g_2 = copy.deepcopy(g_1)
        g_2['fname'] = 'mkn_model2.h5'
        g_2['alpha'] = 0.3
        g_2['label'] = '2 comp'
        o_plot.set_plot_dics.append(g_2)

        g_3 = copy.deepcopy(g_1)
        g_3['fname'] = 'mkn_model3.h5'
        g_3['alpha'] = 0.6
        g_3['label'] = '3 comp'
        o_plot.set_plot_dics.append(g_3)

        g_4 = copy.deepcopy(g_1)
        g_4['fname'] = 'mkn_model4.h5'
        g_4['alpha'] = 0.9
        g_4['label'] = '4 comp'
        g_4['legend'] = True
        o_plot.set_plot_dics.append(g_4)
    SLy4_g(o_lc_data)

    def SLy4_z(o_lc_data):
        obs_g = {
            'task': 'mkn obs', "ptype": "cartesian",
            'position': (3, 2),
            'data': o_lc_data, 'band': 'z', 'obs': True,
            'v_n_x': 'time', 'v_n_y': 'mag',
            'color': 'gray', 'marker': 'o', 'ms': 5., 'alpha': 0.2,
            'ymin': 25, 'ymax': 15, 'xmin': 3e-1, 'xmax': 3e1,
            'xlabel': r"time [days]", 'ylabel': '',
            'label': None, 'xscale': 'log',
            'fancyticks': True, 'minorticks': True,
            'legend': True, #'title': 'z band',
            'sharey': True,
        }
        o_plot.set_plot_dics.append(obs_g)

        g_1 = {
            'task': 'mkn model', "ptype": "cartesian",
            'position': (3, 2),
            'data': o_lc_data, 'band': 'z', 'fname': 'mkn_model1.h5',
            'obs': True,
            'v_n_x': 'time', 'v_n_y': 'mag',
            'color': 'red', 'alpha': 0.1,
            'ymin': 25, 'ymax': 15, 'xmin': 3e-1, 'xmax': 3e1,
            'xlabel': r"time [days]", 'ylabel': '',
            'label': '1 comp', 'xscale': 'log',
            'fancyticks': True, 'minorticks': True,
            'sharey': True,
        }
        o_plot.set_plot_dics.append(g_1)

        g_2 = copy.deepcopy(g_1)
        g_2['fname'] = 'mkn_model2.h5'
        g_2['alpha'] = 0.3
        g_2['label'] = '2 comp'
        o_plot.set_plot_dics.append(g_2)

        g_3 = copy.deepcopy(g_1)
        g_3['fname'] = 'mkn_model3.h5'
        g_3['alpha'] = 0.6
        g_3['label'] = '3 comp'
        o_plot.set_plot_dics.append(g_3)

        g_4 = copy.deepcopy(g_1)
        g_4['fname'] = 'mkn_model4.h5'
        g_4['alpha'] = 0.9
        g_4['label'] = '4 comp'
        # g_4['legend'] = True
        o_plot.set_plot_dics.append(g_4)
    SLy4_z(o_lc_data)

    def SLy4_Ks(o_lc_data):
        obs_g = {
            'task': 'mkn obs', "ptype": "cartesian",
            'position': (3, 3),
            'data': o_lc_data, 'band': 'Ks', 'obs': True,
            'v_n_x': 'time', 'v_n_y': 'mag',
            'color': 'gray', 'marker': 'o', 'ms': 5., 'alpha': 0.2,
            'ymin': 25, 'ymax': 15, 'xmin': 3e-1, 'xmax': 3e1,
            'xlabel': r"time [days]", 'ylabel': '',
            'label': None, 'xscale': 'log',
            'fancyticks': True, 'minorticks': True,
            'legend': True, #'title': 'Ks band',
            'sharey': True,
        }
        o_plot.set_plot_dics.append(obs_g)

        g_1 = {
            'task': 'mkn model', "ptype": "cartesian",
            'position': (3, 3),
            'data': o_lc_data, 'band': 'Ks', 'fname': 'mkn_model1.h5',
            'obs': True,
            'v_n_x': 'time', 'v_n_y': 'mag',
            'color': 'red', 'alpha': 0.1,
            'ymin': 25, 'ymax': 15, 'xmin': 3e-1, 'xmax': 3e1,
            'xlabel': r"time [days]", 'ylabel': '',
            'label': '1 comp', 'xscale': 'log',
            'fancyticks': True, 'minorticks': True,
            'sharey': True,
        }
        o_plot.set_plot_dics.append(g_1)

        g_2 = copy.deepcopy(g_1)
        g_2['fname'] = 'mkn_model2.h5'
        g_2['alpha'] = 0.3
        g_2['label'] = '2 comp'
        o_plot.set_plot_dics.append(g_2)

        g_3 = copy.deepcopy(g_1)
        g_3['fname'] = 'mkn_model3.h5'
        g_3['alpha'] = 0.6
        g_3['label'] = '3 comp'
        o_plot.set_plot_dics.append(g_3)

        g_4 = copy.deepcopy(g_1)
        g_4['fname'] = 'mkn_model4.h5'
        g_4['alpha'] = 0.9
        g_4['label'] = '4 comp'
        # g_4['legend'] = True
        o_plot.set_plot_dics.append(g_4)
    SLy4_Ks(o_lc_data)

    o_plot.main()
    exit(1)





















    """-------------------------------------BAND1--------------------------------"""

    """ --- --- ---- --- --- --- SIM 1 --- -- --- --- --- --- """

    sim = "DD2_M13641364_M0_SR"
    o_lc_data = EXTRACT_LIGHTCURVE(sim)
    sim_lc = {
        'task': 'mkn model', "ptype": "cartesian",
        'position': (1, 1),
        'data': o_lc_data, 'band': 'g', 'obs': True,
        'v_n_x': 'time', 'v_n_y': 'mag',
        'color': 'green', 'alpha': 0.5,
        'ymin': 25, 'ymax': 15, 'xmin': 3e-1, 'xmax': 3e1,
        'xlabel': r"time [days]", 'ylabel': r"AB magnitude at 40 Mpc",
        'label': 'DD2', 'xscale': 'log',
        'fancyticks': True, 'minorticks': True,
    }
    o_plot.set_plot_dics.append(sim_lc)

    """ --- --- ---- --- --- --- SIM 2 --- -- --- --- --- --- """

    sim = "LS220_M13641364_M0_SR"
    o_lc_data = EXTRACT_LIGHTCURVE(sim)
    sim_lc = {
        'task': 'mkn model', "ptype": "cartesian",
        'position': (1, 1),
        'data': o_lc_data, 'band': 'g', 'obs': True,
        'v_n_x': 'time', 'v_n_y': 'mag',
        'color': 'blue', 'alpha': 0.5,
        'ymin': 25, 'ymax': 15,  'xmin': 3e-1, 'xmax': 3e1,
        'xlabel': r"time [days]", 'ylabel': r"AB magnitude at 40 Mpc",
        'label': 'LS220', 'xscale': 'log',
        'fancyticks': True, 'minorticks': True,
    }
    o_plot.set_plot_dics.append(sim_lc)

    """ --- --- ---- --- --- --- SIM 3 --- -- --- --- --- --- """

    sim = "SLy4_M13641364_M0_SR"
    o_lc_data = EXTRACT_LIGHTCURVE(sim)
    sim_lc = {
        'task': 'mkn model', "ptype": "cartesian",
        'position': (1, 1),
        'data': o_lc_data, 'band': 'g', 'obs': True,
        'v_n_x': 'time', 'v_n_y': 'mag',
        'color': 'red', 'alpha': 0.5,
        'ymin': 25, 'ymax': 15,  'xmin': 3e-1, 'xmax': 3e1,
        'xlabel': r"time [days]", 'ylabel': r"AB magnitude at 40 Mpc",
        'label': 'SLy4', 'xscale': 'log',
        'fancyticks': True, 'minorticks': True,
    }
    o_plot.set_plot_dics.append(sim_lc)

    sim_lc_obs = {
        'task': 'mkn obs', "ptype": "cartesian",
        'position': (1, 1),
        'data': o_lc_data, 'band': 'g', 'obs': True,
        'v_n_x': 'time', 'v_n_y': 'mag',
        'color': 'gray', 'marker': 'o', 'ms': 5., 'alpha': 0.6,
        'ymin': 25, 'ymax': 15,  'xmin': 3e-1, 'xmax': 3e1,
        'xlabel': r"time [days]", 'ylabel': r"AB magnitude at 40 Mpc",
        'label': None, 'xscale': 'log',
        'fancyticks': True, 'minorticks': True,
        'legend': True, 'title': 'g band',
    }
    o_plot.set_plot_dics.append(sim_lc_obs)

    """-------------------------------------BAND2--------------------------------"""

    """ --- --- ---- --- --- --- SIM 1 --- -- --- --- --- --- """

    sim = "DD2_M13641364_M0_SR"
    o_lc_data = EXTRACT_LIGHTCURVE(sim)
    sim_lc = {
        'task': 'mkn model', "ptype": "cartesian",
        'position': (1, 2),
        'data': o_lc_data, 'band': 'z', 'obs': True,
        'v_n_x': 'time', 'v_n_y': 'mag',
        'color': 'green', 'alpha': 0.5,
        'ymin': 25, 'ymax': 15,  'xmin': 3e-1, 'xmax': 3e1,
        'xlabel': r"time [days]", 'ylabel': r"AB magnitude at 40 Mpc",
        'label': None, 'xscale': 'log',
        'fancyticks': True, 'minorticks': True,
    }
    o_plot.set_plot_dics.append(sim_lc)

    """ --- --- ---- --- --- --- SIM 2 --- -- --- --- --- --- """

    sim = "LS220_M13641364_M0_SR"
    o_lc_data = EXTRACT_LIGHTCURVE(sim)
    sim_lc = {
        'task': 'mkn model', "ptype": "cartesian",
        'position': (1, 2),
        'data': o_lc_data, 'band': 'z', 'obs': True,
        'v_n_x': 'time', 'v_n_y': 'mag',
        'color': 'blue', 'alpha': 0.5,
        'ymin': 25, 'ymax': 15, 'xmin': 3e-1, 'xmax': 3e1,
        'xlabel': r"time [days]", 'ylabel': r"AB magnitude at 40 Mpc",
        'label': None, 'xscale': 'log',
        'fancyticks': True, 'minorticks': True,
    }
    o_plot.set_plot_dics.append(sim_lc)

    """ --- --- ---- --- --- --- SIM 3 --- -- --- --- --- --- """

    sim = "SLy4_M13641364_M0_SR"
    o_lc_data = EXTRACT_LIGHTCURVE(sim)
    sim_lc = {
        'task': 'mkn model', "ptype": "cartesian",
        'position': (1, 2),
        'data': o_lc_data, 'band': 'z', 'obs': True,
        'v_n_x': 'time', 'v_n_y': 'mag',
        'color': 'red', 'alpha': 0.5,
        'ymin': 25, 'ymax': 15, 'xmin': 3e-1, 'xmax': 3e1,
        'xlabel': r"time [days]", 'ylabel': r"AB magnitude at 40 Mpc",
        'label': None, 'xscale': 'log',
        'fancyticks': True, 'minorticks': True,
        'sharey': True
    }
    o_plot.set_plot_dics.append(sim_lc)

    sim_lc_obs = {
        'task': 'mkn obs', "ptype": "cartesian",
        'position': (1, 2),
        'data': o_lc_data, 'band': 'z', 'obs': True,
        'v_n_x': 'time', 'v_n_y': 'mag',
        'color': 'gray', 'marker': 'o', 'ms': 5., 'alpha': 0.6,
        'ymin': 25, 'ymax': 15, 'xmin': 3e-1, 'xmax': 3e1,
        'xlabel': r"time [days]", 'ylabel': r"AB magnitude at 40 Mpc",
        'label': None, 'xscale': 'log',
        'fancyticks': True, 'minorticks': True,
        'legend': True, 'title': 'z band',
        'sharey': True
    }
    o_plot.set_plot_dics.append(sim_lc_obs)

    """-------------------------------------BAND3--------------------------------"""

    """ --- --- ---- --- --- --- SIM 1 --- -- --- --- --- --- """

    sim = "DD2_M13641364_M0_SR"
    o_lc_data = EXTRACT_LIGHTCURVE(sim)
    sim_lc = {
        'task': 'mkn model', "ptype": "cartesian",
        'position': (1, 3),
        'data': o_lc_data, 'band': 'Ks', 'obs': True,
        'v_n_x': 'time', 'v_n_y': 'mag',
        'color': 'green', 'alpha': 0.5,
        'ymin': 25, 'ymax': 15, 'xmin': 3e-1, 'xmax': 3e1,
        'xlabel': r"time [days]", 'ylabel': r"AB magnitude at 40 Mpc",
        'label': None, 'xscale': 'log',
        'fancyticks': True, 'minorticks': True,
    }
    o_plot.set_plot_dics.append(sim_lc)

    """ --- --- ---- --- --- --- SIM 2 --- -- --- --- --- --- """

    sim = "LS220_M13641364_M0_SR"
    o_lc_data = EXTRACT_LIGHTCURVE(sim)
    sim_lc = {
        'task': 'mkn model', "ptype": "cartesian",
        'position': (1, 3),
        'data': o_lc_data, 'band': 'Ks', 'obs': True,
        'v_n_x': 'time', 'v_n_y': 'mag',
        'color': 'blue', 'alpha': 0.5,
        'ymin': 25, 'ymax': 15, 'xmin': 3e-1, 'xmax': 3e1,
        'xlabel': r"time [days]", 'ylabel': r"AB magnitude at 40 Mpc",
        'label': None, 'xscale': 'log',
        'fancyticks': True, 'minorticks': True,
    }
    o_plot.set_plot_dics.append(sim_lc)

    """ --- --- ---- --- --- --- SIM 3 --- -- --- --- --- --- """

    sim = "SLy4_M13641364_M0_SR"
    o_lc_data = EXTRACT_LIGHTCURVE(sim)
    sim_lc = {
        'task': 'mkn model', "ptype": "cartesian",
        'position': (1, 3),
        'data': o_lc_data, 'band': 'Ks', 'obs': True,
        'v_n_x': 'time', 'v_n_y': 'mag',
        'color': 'red', 'alpha': 0.5,
        'ymin': 25, 'ymax': 15, 'xmin': 3e-1, 'xmax': 3e1,
        'xlabel': r"time [days]", 'ylabel': r"AB magnitude at 40 Mpc",
        'label': None, 'xscale': 'log',
        'fancyticks': True, 'minorticks': True,
        'sharey': True
    }
    o_plot.set_plot_dics.append(sim_lc)

    sim_lc_obs = {
        'task': 'mkn obs', "ptype": "cartesian",
        'position': (1, 3),
        'data': o_lc_data, 'band': 'Ks', 'obs': True,
        'v_n_x': 'time', 'v_n_y': 'mag',
        'color': 'gray', 'marker': 'o', 'ms': 5., 'alpha': 0.6,
        'ymin': 25, 'ymax': 15, 'xmin': 3e-1, 'xmax': 3e1,
        'xlabel': r"time [days]", 'ylabel': r"AB magnitude at 40 Mpc",
        'label': None, 'xscale': 'log',
        'fancyticks': True, 'minorticks': True,
        'legend': True, 'title': 'Ks band',
        'sharey': True
    }
    o_plot.set_plot_dics.append(sim_lc_obs)

    o_plot.main()
    exit(1)
# plot_mkn_3sims_multiple_component_arrangment()