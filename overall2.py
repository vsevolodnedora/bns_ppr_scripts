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
import math
import time
import copy
import h5py
import csv
import os
import functools

from scidata.utils import locate
import scidata.carpet.hdf5 as h5
import scidata.xgraph as xg

from matplotlib.ticker import AutoMinorLocator, FixedLocator, NullFormatter, \
    MultipleLocator
from matplotlib.colors import LogNorm, Normalize
from matplotlib.colors import Normalize, LogNorm
from matplotlib.collections import PatchCollection
from matplotlib.patches import Rectangle
from matplotlib import patches

from general import *
from lists import *
from filework import *
from units import time_constant, volume_constant, energy_constant
from outflowed import SIM_EJ_HIST
from d1analysis import LOAD_ITTIME

from plotting_nexus import PLOT_MANY_TASKS
from preanalysis import LOAD_INIT_DATA
from d1analysis import ADD_METHODS_1D

from d3analysis import LOAD_PROFILE_XYXZ, LOAD_RES_CORR
"""=============================================| LISTS of SIMS |===================================================="""

list_sims_ls220 = [
"LS220_M13641364_M0_SR", # FULL 3D!
"LS220_M14001330_M0_SR", # No 3D
"LS220_M14351298_M0_SR", # No 3D
"LS220_M14691268_M0_SR", # No 3D

"LS220_M13641364_M0_LK_SR", # No 3D
"LS220_M14691268_M0_LK_SR", # 10 profs
]

list_sims_dd2 = [
"DD2_M13641364_M0_SR_R04",      # 3 Profs
"DD2_M13641364_M0_LK_SR_R04",   # 14 profs
"DD2_M13641364_M0_SR",          # 107 profs
"DD2_M14971245_M0_SR",          # No 3D data
"DD2_M15091235_M0_LK_SR"        # 12 profs
]

list_sims_sfho = [
"SFHo_M13641364_M0_SR", # 21ms               # FULL 3D -> extracting
"SFHo_M13641364_M0_LK_SR_2019pizza",  # 55ms # No 3D data
"SFHo_M14521283_M0_LK_SR",# 30ms             # No 3D data
"SFHo_M14521283_M0_LK_SR_2019pizza",# 26ms   # No 3D data
"SFHo_M14521283_M0_SR", # 33ms               # Full 3D (29 profs)
]

list_sims_sly4 = [
"SLy4_M13641364_M0_SR", # 21ms              # Full 3D (31 profs)
# "SLy4_M13641364_M0_LK_SR",  # 55ms #      #
"SLy4_M14521283_M0_SR",# 30ms # Full 3D     # Full 3D ()
]

list_sims_all = list_sims_ls220 + list_sims_dd2 + list_sims_sfho + list_sims_sly4
"""=================================================================================================================="""

def eos_color(eos):
    if eos == 'DD2':
        return 'blue'
    elif eos == 'BHBlp':
        return 'purple'
    elif eos == 'LS220':
        return 'orange'
    elif eos == 'SFHo':
        return 'red'
    elif eos == 'SLy4':
        return 'green'
    else:
        return 'black'

def get_ms(q, qmin=1, qmax = 1.4, msmin = 5., msmax = 10.):

    k = (qmax - qmin) / (msmax - msmin)
    b = qmax - (k * msmax)

    return (q - b) / k

"""=================================================================================================================="""

class TEX_TABLE:

    def __init__(self):

        self.translate_dic = {
            "Mdisk3D": "$M_{\\text{disk}}$",

            "M1": "$M_a$",
            "M2": "$M_b$",

            "tcoll_gw": "$t_{\\text{BH}}$",

            "theta_rms": "$\\theta_{\\text{ej}}$",

            "Mej_tot": "$M_{\\text{ej}}$",
            "Ye_ave": "$\\langle Y_e \\rangle$",
            "vel_inf_ave": "$\\upsilon_{\\text{ej}}$",

            # "Mej_bern": "$M_{\\text{ej}}^b$",
            # "Yeej_bern": "$\\langle Y_e \\rangle^b$",
            # "vej_bern": "$\\upsilon_{\\text{ej}}^b$",
        }

        self.units_dic = {
            "Mdisk3D": "$[M_{\odot}]$",
            "EOS": " ",
            "LK": "  ",
            "theta_rms": " ",
            "tcoll_gw": "[ms]",

            "Mej_tot": "$[10^{-2} M_{\odot}]$",
            "Ye_ave": "  ",
            "vel_inf_ave": "$[c]$",

            # "Mej_bern": "$[M_{\odot}]$",
            # "Yeej_bern": "  ",
            # "vej_bern": "$[c]$",

            "M1": "$[M_{\odot}]$",
            "M2": "$[M_{\odot}]$",
        }

    def get_lbl(self, v_n, criterion="_0"):
        if v_n == "Mdisk3D": return "$M_{\\text{disk}}$"
        elif v_n == "M1": return "$M_a$"
        elif v_n == "M2": return "$M_b$"
        elif v_n == "tcoll_gw": return "$t_{\\text{BH}}$"
        elif v_n == "theta_rms" and criterion=="_0": return "$\\theta_{\\text{ej}}$"
        elif v_n == "Mej_tot" and criterion=="_0": return "$M_{\\text{ej}}$"
        elif v_n == "Ye_ave" and criterion=="_0": return "$\\langle Y_e \\rangle$"
        elif v_n == "vel_inf_ave" and criterion=="_0": return "$\\upsilon_{\\text{ej}}$"
        elif v_n == "theta_rms" and criterion=="_0_b_w": return "$\\theta_{\\text{ej}}^{\\text{w}}$"
        elif v_n == "Mej_tot" and criterion=="_0_b_w": return "$M_{\\text{ej}}^{\\text{w}}$"
        elif v_n == "Ye_ave" and criterion=="_0_b_w": return "$\\langle Y_e ^{\\text{w}}  \\rangle$"
        elif v_n == "vel_inf_ave" and criterion=="_0_b_w": return "$\\upsilon_{\\text{ej}}^{\\text{w}}$"
        elif v_n == "LK": return "LK"
        else: return v_n

    def get_value(self, sim, v_n, criterion="_0"):

        d1class = ADD_METHODS_1D(sim)
        selfclass = LOAD_INIT_DATA(sim)

        if v_n in d1class.list_parameters:
            if v_n == "tcoll_gw":
                tcoll = d1class.get_par("tcoll_gw")
                if not np.isnan(tcoll):
                    tcoll = float(tcoll - d1class.get_par("tmerger_gw"))
                    return str("$ {:.1f}$".format(tcoll * 1e3))  # ms
                else:
                    print("Failed to load 'tcoll' or 'tmerg' for {}".format(sim))
                    tlast = d1class.get_arr("t_unb_mass")[-1] - float(d1class.get_par("tmerger_gw"))
                    return str("$>{:.1f}$".format(tlast * 1e3))  # ms

            elif v_n == "Mej_tot" and criterion != "_0_b_w" :
                return "$ {:.1f}$" % (d1class.get_par(v_n, criterion) * 1e2)
            elif v_n == "Mej_tot" and criterion == "_0_b_w" and np.isnan(d1class.get_par("tcoll_gw")):
                return "$>{:.1f}$" % (d1class.get_par(v_n, criterion) * 1e2)
            elif v_n == "Mej_tot" and criterion == "_0_b_w" and not np.isnan(d1class.get_par("tcoll_gw")):
                return "$\sim{:.1f}$" % (d1class.get_par(v_n, criterion) * 1e2)

            else:
                return d1class.get_par(v_n, criterion=criterion)

        if v_n == "LK":
            if sim.__contains__("LK"):
                return("LK")
            else:
                return("  ")

        if v_n in selfclass.par_dic.keys():
            print(selfclass.par_dic.keys())
            return selfclass.get_par(v_n)

        raise NameError("v_n:{} is not recognized".format(v_n))

    def get_value2(self, sim, v_n, crit, prec):

        d1class = ADD_METHODS_1D(sim)
        selfclass = LOAD_INIT_DATA(sim)

        if v_n in d1class.list_parameters:
            if v_n == "tcoll_gw":
                tcoll = d1class.get_par("tcoll_gw")
                if not np.isnan(tcoll):
                    #tcoll = float(tcoll - d1class.get_par("tmerger_gw"))
                    return str("$ {:.1f}$".format(tcoll * 1e3))
                else:
                    print("Failed to load 'tcoll' or 'tmerg' for {}".format(sim))
                    tlast = d1class.get_arr("t_unb_mass")[-1] - float(d1class.get_par("tmerger_gw"))
                    return str("$>{:.1f}$".format(tlast * 1e3))  # ms
            elif v_n == "Mdisk3D":
                mdisk = d1class.get_par("Mdisk3D")
                if not np.isnan(mdisk):
                    return("%{}f".format(prec) % mdisk)
                else:
                    return("N/A")

            elif v_n == "Mej_tot" and crit == "_0":
                mej = d1class.get_par(v_n, criterion=crit) * 1e2
                return ("%{}f".format(prec) % mej)

            elif v_n == "Mej_tot" and crit != "_0":
                mej = d1class.get_par(v_n, criterion=crit) * 1e2
                if np.isnan(d1class.get_par("tcoll_gw")):
                    return "$>" + "%{}f".format(prec) % mej + "$"
                else:
                    return "$\sim" + "%{}f".format(prec) % mej + "$"

            elif prec != "str":
                return("%{}f".format(prec) % d1class.get_par(v_n, crit))
            elif prec == "str":
                return (str(d1class.get_par(v_n, crit)))

        if v_n in selfclass.par_dic.keys():

            if v_n == "EOS" and selfclass.get_par("EOS") == "SFHo" and sim.__contains__("2019pizza"):
                return("{}$^{}$".format(str(selfclass.get_par(v_n)), "p"))

            if prec == "str":
                return(str(selfclass.get_par(v_n)))
            else:
                return("%{}f".format(prec) % selfclass.get_par(v_n))

            # print(selfclass.par_dic.keys())
            # return selfclass.get_par(v_n)

        if v_n == "LK":
            if sim.__contains__("LK"):
                return("LK")
            else:
                return("  ")

        raise NameError("v_n:{} is not recognized\n{}\n{}"
                        .format(v_n,d1class.list_parameters, selfclass.par_dic.keys()))

    def print_latex_table(self):

        sims = list_sims_all
        v_ns =  ["EOS", "LK", "M1",   "M2", 'tcoll_gw', 'Mdisk3D', 'Mej_tot', 'Ye_ave', 'vel_inf_ave', 'theta_rms',
                 'Mej_tot', 'Ye_ave', 'vel_inf_ave', 'theta_rms']
        crits = ["",     "",   "",    "",    '',         '',        '_0',       '_0',       '_0',
                 '_0', '_0_b_w', '_0_b_w', '_0_b_w', '_0_b_w']
        precs = ["str", "str", "1.2", "1.2", "str",      ".2",     ".2",       ".2",       ".2",
                 ".2", ".2", ".2", ".2", ".0"]


        assert len(v_ns) == len(crits)
        assert len(precs) == len(v_ns)

        rows = []
        for i, sim in enumerate(sims):
            # 1 & 6 & 87837 & 787 \\
            row = ''
            j = 0
            for v_n, prec, crit in zip(v_ns, precs, crits):
                print("\tPrinting {}".format(v_n))
                # if prec != "str":
                #     __val = self.get_value2(sim, v_n, crit=crit, prec=prec)
                #     # if not np.isnan(__val):
                #     #     val = "%{}f".format(prec) % __val
                #     # else:
                #     #     val = "N/A"
                # else:
                val = self.get_value2(sim, v_n, crit=crit, prec=prec)
                row = row + val
                if j != len(v_ns)-1: row = row + ' & '
                j = j + 1
                # if v_n != v_ns[-1]: row = row + ' & '
            row = row + ' \\\\'  # = \\
            rows.append(row)

        # -------------------------------------------
        print("\n")
        size = '{'
        head = ''
        i = 0
        for v_n, crit in zip(v_ns, crits):
            v_n = self.get_lbl(v_n, criterion=crit)
            size = size + 'c'
            head = head + '{}'.format(v_n)
            if v_n != v_ns[-1]: size = size + ' '
            if i != len(v_ns) - 1: head = head + ' & '
            i = i + 1
        size = size + '}'

        unit_bar = ''
        for i, v_n in enumerate(v_ns):
            if v_n in self.units_dic.keys():
                unit = self.units_dic[v_n]
            else:
                unit = v_n
            unit_bar = unit_bar + '{}'.format(unit)
            # if v_ns.index(v_n) != len(v_ns): unit_bar = unit_bar + ' & '
            if i != len(v_ns)-1: unit_bar = unit_bar + ' & '

        head = head + ' \\\\'  # = \\
        unit_bar = unit_bar + ' \\\\ '

        print('\\begin{table*}[t]')
        print('\\begin{center}')
        print('\\begin{tabular}' + '{}'.format(size))
        print('\\hline')
        print(head)
        print(unit_bar)
        print('\\hline\\hline')

        for row in rows:
            print(row)

        # rows = []
        # for i, sim in enumerate(sims):
        #
        #     # 1 & 6 & 87837 & 787 \\
        #     row = ''
        #     for v_n, prec, crit in zip(v_ns, precs, crits):
        #         print("\tPrinting {}".format(v_n))
        #         if prec != "str":
        #             val = "%{}f".format(prec) % self.get_value(sim, v_n, crit)
        #         else:
        #             val = self.get_value(sim, v_n)
        #         row = row + val
        #         if v_n != v_ns[-1]: row = row + ' & '
        #     row = row + ' \\\\'  # = \\
        #     rows.append(row)
            # print(row)


        print('\\hline')
        print('\\end{tabular}')
        print('\\end{center}')
        print('\\caption{I am your table! }')
        print('\\label{tbl:1}')
        print('\\end{table*}')

        exit(0)


# plot ejecta properties

# def get_color(list_of_sims):





tex = TEX_TABLE()
# tex.print_latex_table()

def plot_ejecta_properites():

    o_plot = PLOT_MANY_TASKS()
    o_plot.gen_set["figdir"] = Paths.plots + "all/"
    o_plot.gen_set["type"] = "cartesian"
    o_plot.gen_set["figsize"] = (4.2, 9.2)  # <->, |]
    o_plot.gen_set["figname"] = "dyn_ej_ave.png"
    o_plot.gen_set["sharex"] = True
    o_plot.gen_set["sharey"] = False
    o_plot.gen_set["subplots_adjust_h"] = 0.0
    o_plot.gen_set["subplots_adjust_w"] = 0.0
    o_plot.set_plot_dics = []

    # label1 = "Dyn"
    # label2 = "Wind"
    label1v = "Dyn."
    label2v = "Wind"
    fontsize = 12
    labelsize = 12

    lams, mejs, eoss = [], [], []
    for sim in list_sims_all:
        print(sim),
        self_prop = LOAD_INIT_DATA(sim)
        lam = self_prop.get_par("Lambda")
        lams.append(lam)
        print("lam:{}".format(lam)),

        d1class = ADD_METHODS_1D(sim)
        mej = d1class.get_par("Mej_tot", criterion="_0")
        mejs.append(mej)
        print("mej:{}".format(mej)),

        print('')

        if sim.__contains__("LK"):
            marker = "d"
        else:
            marker = "x"

        if sim.split('_')[0] in eoss:
            lbl = None
        else:
            eoss.append(sim.split('_')[0])
            lbl = sim.split('_')[0]

        ms = get_ms(float(self_prop.get_par("q")))

        # ================ Mej
        dic_ej_prof = {
            'task': 'marker', 'ptype': 'cartesian',
            'position': (1, 1),
            'x': np.array(lam, dtype=float), 'y': np.array(mej * 1e2, dtype=float),
            'v_n_x': 'Lambda', 'v_n_y': None,
            'color': eos_color(sim.split('_')[0]), 'marker': marker, 'ms': ms, 'alpha': 1.,
            # 'xmin': 400, 'xmax': 850,
            'ymin': 0.07, 'ymax': 0.7,
            'xlabel': r"$\tilde{\Lambda}$", 'ylabel': r'$M_{\rm ej}$ $[10^{-2} M_{\odot}]$',
            'label': lbl, 'yscale': None, 'title': {},
            'fancyticks': True, 'minorticks': True,
            'fontsize': fontsize,
            'labelsize': labelsize,
            'legend': {'loc': 'upper right', 'ncol': 2, 'fontsize': 12}
        }
        o_plot.set_plot_dics.append(dic_ej_prof)

        # ================================================= Ye
        ye = d1class.get_par("Ye_ave", criterion="_0")
        dic_ye_prof = {
            'task': 'marker', 'ptype': 'cartesian',
            'position': (2, 1),
            'x': np.array(lam, dtype=float), 'y': np.array(ye, dtype=float),
            'v_n_x': 'Lambda', 'v_n_y': None,
            'color': eos_color(sim.split('_')[0]), 'marker': marker, 'ms': ms, 'alpha': 1.,
            # 'xmin': 400, 'xmax': 850,
            'ymin': 0.10, 'ymax': 0.34,
            'xlabel': r"$\tilde{\Lambda}$", 'ylabel': r'$Y_e$',
            'label': lbl, 'yscale': None, 'title': {},
            'fancyticks': True, 'minorticks': True,
            'fontsize': fontsize,
            'labelsize': labelsize,
            # 'legend': {'loc': 'upper right', 'ncol': 2, 'fontsize': 12}
        }
        o_plot.set_plot_dics.append(dic_ye_prof)

        # =========================================== vinf
        vinf = d1class.get_par("vel_inf_ave", criterion="_0")
        dic_ye_prof = {
            'task': 'marker', 'ptype': 'cartesian',
            'position': (3, 1),
            'x': np.array(lam, dtype=float), 'y': np.array(vinf, dtype=float),
            'v_n_x': 'Lambda', 'v_n_y': None,
            'color': eos_color(sim.split('_')[0]), 'marker': marker, 'ms': ms, 'alpha': 1.,
            # 'xmin': 400, 'xmax': 850,
            'ymin': 0.10, 'ymax': 0.27,
            'xlabel': r"$\tilde{\Lambda}$", 'ylabel': r'$\upsilon_{\infty}$ $[c]$',
            'label': lbl, 'yscale': None, 'title': {},
            'fancyticks': True, 'minorticks': True,
            'fontsize': fontsize,
            'labelsize': labelsize,
            # 'legend': {'loc': 'upper right', 'ncol': 2, 'fontsize': 12}
        }
        o_plot.set_plot_dics.append(dic_ye_prof)

        # =========================================== DiskMass3D
        # mdisk = d1class.get_par("Mdisk3D")
        # if not np.isnan(mdisk):
        #     dic_ye_prof = {
        #         'task': 'marker', 'ptype': 'cartesian',
        #         'position': (4, 1),
        #         'x': np.array(lam, dtype=float), 'y': np.array(mdisk, dtype=float),
        #         'v_n_x': 'Lambda', 'v_n_y': None,
        #         'color': eos_color(sim.split('_')[0]), 'marker': marker, 'ms': ms, 'alpha': 1.,
        #         # 'xmin': 400, 'xmax': 850,
        #         'ymin': 0., 'ymax': 0.32,
        #         'xlabel': r"$\tilde{\Lambda}$", 'ylabel': r'$M_{\rm{disk}}$ $[M_{\odot}]$',
        #         'label': lbl, 'yscale': None, 'title': {},
        #         'fancyticks': True, 'minorticks': True,
        #         'fontsize': fontsize,
        #         'labelsize': labelsize,
        #         # 'legend': {'loc': 'upper right', 'ncol': 2, 'fontsize': 12}
        #     }
        #     o_plot.set_plot_dics.append(dic_ye_prof)


        # print('')
        # break

    o_plot.main()
    exit(1)
# plot_ejecta_properites()

def plot_disk_mass_evol():

    o_plot = PLOT_MANY_TASKS()
    o_plot.gen_set["figdir"] = Paths.plots + "all/"
    o_plot.gen_set["type"] = "cartesian"
    o_plot.gen_set["figsize"] = (4.2, 3.6)  # <->, |]
    o_plot.gen_set["figname"] = "disk_mass_ev.png"
    o_plot.gen_set["sharex"] = True
    o_plot.gen_set["sharey"] = False
    o_plot.gen_set["subplots_adjust_h"] = 0.0
    o_plot.gen_set["subplots_adjust_w"] = 0.0
    o_plot.set_plot_dics = []

    # label1 = "Dyn"
    # label2 = "Wind"
    label1v = "Dyn."
    label2v = "Wind"
    fontsize = 12
    labelsize = 12

    lams, mejs, eoss = [], [], []
    list_sims_all = list_sims_ls220 + list_sims_sfho + list_sims_sly4
    for sim in list_sims_all:
        print(sim),
        self_prop = LOAD_INIT_DATA(sim)
        lam = self_prop.get_par("Lambda")
        lams.append(lam)
        print("lam:{}".format(lam)),

        d1class = ADD_METHODS_1D(sim)
        mej = d1class.get_par("Mej_tot", criterion="_0")
        mejs.append(mej)
        print("mej:{}".format(mej)),

        print('')

        if sim.__contains__("LK"):
            marker = "d"
        else:
            marker = "x"

        if sim.split('_')[0] in eoss:
            lbl = None
        else:
            eoss.append(sim.split('_')[0])
            lbl = sim.split('_')[0]

        lw = get_ms(float(self_prop.get_par("q")), msmin=0.5, msmax=2.)

        # =========================================== DiskMass3D
        mdisk = d1class.get_arr("disk_mass")
        tdisk = (d1class.get_arr("t_disk_mass") - d1class.get_par("tmerger_gw")) * 1e3
        tcoll = d1class.get_par("tcoll_gw") * 1e3

        if len(mdisk) > 0:
            dic_ye_prof = {
                'task': 'marker', 'ptype': 'cartesian',
                'position': (1, 1),
                'x': np.array(tdisk, dtype=float), 'y': np.array(mdisk, dtype=float),
                'v_n_x': 'Lambda', 'v_n_y': None,
                'color': eos_color(sim.split('_')[0]), 'ls': '-', 'lw': lw, 'alpha': 1., 'ds': 'default',
                # 'xmin': 400, 'xmax': 850,
                'ymin': 0., 'ymax': 0.32,
                'xlabel': r"$t-t_{\rm{merg}}$ [ms]", 'ylabel': r'$M_{\rm{disk}}$ $[M_{\odot}]$',
                'label': lbl, 'yscale': None, 'title': {},
                'fancyticks': True, 'minorticks': True,
                'fontsize': fontsize,
                'labelsize': labelsize,
                'legend': {'loc': 'upper right', 'ncol': 1, 'fontsize': 12}
            }
            o_plot.set_plot_dics.append(dic_ye_prof)

            # tcoll = d1class.get_par("tcoll_gw")
            if not np.isnan(tcoll):
                # print( tdisk[np.where(tdisk==tcoll)])
                dic_ye_prof = {
                    'task': 'marker', 'ptype': 'cartesian',
                    'position': (1, 1),
                    'x': tdisk[find_nearest_index(tdisk, tcoll)], 'y': mdisk[find_nearest_index(tdisk, tcoll)],
                    'v_n_x': 'Lambda', 'v_n_y': None,
                    'color': eos_color(sim.split('_')[0]), 'marker': 'o', 'ms': 2, 'alpha': 1.,
                    'xmin': -5, 'xmax': 30,
                    'ymin': 0., 'ymax': 0.2,
                    'xlabel': r"$t-t_{\rm{merg}}$ [ms]", 'ylabel': r'$M_{\rm{disk}}$ $[M_{\odot}]$',
                    'label': None, 'yscale': None, 'title': {},
                    'fancyticks': True, 'minorticks': True,
                    'fontsize': fontsize,
                    'labelsize': labelsize,
                    # 'legend': {'loc': 'upper right', 'ncol': 2, 'fontsize': 12}
                }
                o_plot.set_plot_dics.append(dic_ye_prof)


                # coll_time_vert = {
                #     'task': 'vertline', 'dtype': '-', 'ptype': 'cartesian',
                #     'value': (tcoll) * 1e3,
                #     'position': (1, 1),
                #     'ls': '-.', 'color': eos_color(sim.split('_')[0]), 'lw': 0.8,
                # }
                # o_plot.set_plot_dics.append(coll_time_vert)




        # print('')
        # break

    o_plot.main()
    exit(1)
# plot_disk_mass_evol()

def plot_ejecta_properites_bern():

    o_plot = PLOT_MANY_TASKS()
    o_plot.gen_set["figdir"] = Paths.plots + "all/"
    o_plot.gen_set["type"] = "cartesian"
    o_plot.gen_set["figsize"] = (4.2, 9.2)  # <->, |]
    o_plot.gen_set["figname"] = "bern_ej_ave.png"
    o_plot.gen_set["sharex"] = True
    o_plot.gen_set["sharey"] = False
    o_plot.gen_set["subplots_adjust_h"] = 0.0
    o_plot.gen_set["subplots_adjust_w"] = 0.0
    o_plot.set_plot_dics = []

    # label1 = "Dyn"
    # label2 = "Wind"
    label1v = "Dyn."
    label2v = "Wind"
    fontsize = 12
    labelsize = 12

    lams, mejs, eoss = [], [], []
    for sim in list_sims_all:
        print(sim),
        self_prop = LOAD_INIT_DATA(sim)
        lam = self_prop.get_par("Lambda")
        lams.append(lam)
        print("lam:{}".format(lam)),

        d1class = ADD_METHODS_1D(sim)
        mej = d1class.get_par("Mej_tot", criterion="_0_b_w")
        mejs.append(mej)
        print("mej:{}".format(mej)),

        print('')

        if sim.__contains__("LK"):
            marker = "d"
        else:
            marker = "x"

        if sim.split('_')[0] in eoss:
            lbl = None
        else:
            eoss.append(sim.split('_')[0])
            lbl = sim.split('_')[0]

        ms = get_ms(float(self_prop.get_par("q")))

        tcoll = d1class.get_par("tcoll_gw")
        if not np.isnan(tcoll):
            arrow = None
        else:
            arrow = "up"
        print("arrow: {}".format(arrow))
        # ================ Mej
        dic_ej_prof = {
            'task': 'marker', 'ptype': 'cartesian',
            'position': (1, 1),
            'x': np.array(lam, dtype=float), 'y': np.array(mej * 1e2, dtype=float),
            'arrow': arrow,
            'v_n_x': 'Lambda', 'v_n_y': None,
            'color': eos_color(sim.split('_')[0]), 'marker': marker, 'ms': ms, 'alpha': 1.,
            # 'xmin': 400, 'xmax': 850,
            'ymin': 0.07, 'ymax': 2.9,
            'xlabel': r"$\tilde{\Lambda}$", 'ylabel': r'$M_{\rm ej}$ $[10^{-2} M_{\odot}]$',
            'label': lbl, 'yscale': None, 'title': {},
            'fancyticks': True, 'minorticks': True,
            'fontsize': fontsize,
            'labelsize': labelsize,
            'legend': {'loc': 'upper center', 'ncol': 2, 'fontsize': 12}
        }
        o_plot.set_plot_dics.append(dic_ej_prof)

        # ================================================= Ye
        ye = d1class.get_par("Ye_ave", criterion="_0_b_w")
        dic_ye_prof = {
            'task': 'marker', 'ptype': 'cartesian',
            'position': (2, 1),
            'x': np.array(lam, dtype=float), 'y': np.array(ye, dtype=float),
            'v_n_x': 'Lambda', 'v_n_y': None,
            'color': eos_color(sim.split('_')[0]), 'marker': marker, 'ms': ms, 'alpha': 1.,
            # 'xmin': 400, 'xmax': 850,
            'ymin': 0.10, 'ymax': 0.34,
            'xlabel': r"$\tilde{\Lambda}$", 'ylabel': r'$Y_e$',
            'label': lbl, 'yscale': None, 'title': {},
            'fancyticks': True, 'minorticks': True,
            'fontsize': fontsize,
            'labelsize': labelsize,
            # 'legend': {'loc': 'upper right', 'ncol': 2, 'fontsize': 12}
        }
        o_plot.set_plot_dics.append(dic_ye_prof)

        # =========================================== vinf
        vinf = d1class.get_par("vel_inf_ave", criterion="_0_b_w")
        dic_ye_prof = {
            'task': 'marker', 'ptype': 'cartesian',
            'position': (3, 1),
            'x': np.array(lam, dtype=float), 'y': np.array(vinf, dtype=float),
            'v_n_x': 'Lambda', 'v_n_y': None,
            'color': eos_color(sim.split('_')[0]), 'marker': marker, 'ms': ms, 'alpha': 1.,
            # 'xmin': 400, 'xmax': 850,
            'ymin': 0.10, 'ymax': 0.27,
            'xlabel': r"$\tilde{\Lambda}$", 'ylabel': r'$\upsilon_{\infty}$ $[c]$',
            'label': lbl, 'yscale': None, 'title': {},
            'fancyticks': True, 'minorticks': True,
            'fontsize': fontsize,
            'labelsize': labelsize,
            # 'legend': {'loc': 'upper right', 'ncol': 2, 'fontsize': 12}
        }
        o_plot.set_plot_dics.append(dic_ye_prof)

        # =========================================== DiskMass3D
        # mdisk = d1class.get_par("Mdisk3D")
        # if not np.isnan(mdisk):
        #     dic_ye_prof = {
        #         'task': 'marker', 'ptype': 'cartesian',
        #         'position': (4, 1),
        #         'x': np.array(lam, dtype=float), 'y': np.array(mdisk, dtype=float),
        #         'v_n_x': 'Lambda', 'v_n_y': None,
        #         'color': eos_color(sim.split('_')[0]), 'marker': marker, 'ms': ms, 'alpha': 1.,
        #         # 'xmin': 400, 'xmax': 850,
        #         'ymin': 0., 'ymax': 0.32,
        #         'xlabel': r"$\tilde{\Lambda}$", 'ylabel': r'$M_{\rm{disk}}$ $[M_{\odot}]$',
        #         'label': lbl, 'yscale': None, 'title': {},
        #         'fancyticks': True, 'minorticks': True,
        #         'fontsize': fontsize,
        #         'labelsize': labelsize,
        #         # 'legend': {'loc': 'upper right', 'ncol': 2, 'fontsize': 12}
        #     }
        #     o_plot.set_plot_dics.append(dic_ye_prof)


        # print('')
        # break

    o_plot.main()
    exit(1)
# plot_ejecta_properites_bern()

# plot disk masses

def plot_disk_masses():

    o_plot = PLOT_MANY_TASKS()
    o_plot.gen_set["figdir"] = Paths.plots + "all/"
    o_plot.gen_set["type"] = "cartesian"
    o_plot.gen_set["figsize"] = (4.2, 3.6)  # <->, |]
    o_plot.gen_set["figname"] = "disk_masses.png"
    o_plot.gen_set["sharex"] = True
    o_plot.gen_set["sharey"] = False
    o_plot.gen_set["subplots_adjust_h"] = 0.0
    o_plot.gen_set["subplots_adjust_w"] = 0.0
    o_plot.set_plot_dics = []

    # label1 = "Dyn"
    # label2 = "Wind"
    label1v = "Dyn."
    label2v = "Wind"
    fontsize = 12
    labelsize = 12

    lams, mejs, eoss = [], [], []
    for sim in list_sims_all:
        print(sim),
        self_prop = LOAD_INIT_DATA(sim)
        lam = self_prop.get_par("Lambda")
        lams.append(lam)
        print("lam:{}".format(lam)),

        d1class = ADD_METHODS_1D(sim)
        mej = d1class.get_par("Mej_tot", criterion="_0_b_w")
        mejs.append(mej)
        print("mej:{}".format(mej)),

        print('')

        if sim.__contains__("LK"):
            marker = "d"
        else:
            marker = "x"

        if sim.split('_')[0] in eoss:
            lbl = None
        else:
            eoss.append(sim.split('_')[0])
            lbl = sim.split('_')[0]

        ms = get_ms(float(self_prop.get_par("q")))


        # =========================================== DiskMass3D
        mdisk = d1class.get_par("Mdisk3D")
        if not np.isnan(mdisk):
            dic_ye_prof = {
                'task': 'marker', 'ptype': 'cartesian',
                'position': (4, 1),
                'x': np.array(lam, dtype=float), 'y': np.array(mdisk, dtype=float),
                'v_n_x': 'Lambda', 'v_n_y': None,
                'color': eos_color(sim.split('_')[0]), 'marker': marker, 'ms': ms, 'alpha': 1.,
                # 'xmin': 400, 'xmax': 850,
                'ymin': 0., 'ymax': 0.32,
                'xlabel': r"$\tilde{\Lambda}$", 'ylabel': r'$M_{\rm{disk}}$ $[M_{\odot}]$',
                'label': lbl, 'yscale': None, 'title': {},
                'fancyticks': True, 'minorticks': True,
                'fontsize': fontsize,
                'labelsize': labelsize,
                # 'legend': {'loc': 'upper right', 'ncol': 2, 'fontsize': 12}
            }
            o_plot.set_plot_dics.append(dic_ye_prof)


        # print('')
        # break

    o_plot.main()
    exit(1)
# plot_ejecta_properites_bern()

"""=================================================| 2D slices |===================================================="""

def get_xmin_xmax_ymin_ymax_zmin_zmax(rl):
    if rl == 6:
        xmin, xmax = -14, 14
        ymin, ymax = -14, 14
        zmin, zmax = 0, 14
    elif rl == 5:
        xmin, xmax = -28, 28
        ymin, ymax = -28, 28
        zmin, zmax = 0, 28
    elif rl == 4:
        xmin, xmax = -48, 48
        ymin, ymax = -48, +48
        zmin, zmax = 0, 48
    elif rl == 3:
        xmin, xmax = -88, 88
        ymin, ymax = -88, 88
        zmin, zmax = 0, 88
    elif rl == 2:
        xmin, xmax = -178, 178
        ymin, ymax = -178, +178
        zmin, zmax = 0, 178
    elif rl == 1:
        xmin, xmax = -354, 354
        ymin, ymax = -354, +354
        zmin, zmax = 0, 354
    elif rl == 0:
        xmin, xmax = -1044, 1044
        ymin, ymax = -1044, 1044
        zmin, zmax = 0, 1044
    else:
        # pass
        raise IOError("Set limits for rl:{}".format(rl))

    return xmin, xmax, ymin, ymax, zmin, zmax

def plot_den_unb__vel_z():

    # tmp = d3class.get_data(688128, 3, "xy", "ang_mom_flux")
    # print(tmp.min(), tmp.max())
    # print(tmp)
    # exit(1) # dens_unb_geo

    """ --- --- --- """


    o_plot = PLOT_MANY_TASKS()
    o_plot.gen_set["figdir"] = Paths.plots + 'all/'
    o_plot.gen_set["type"] = "cartesian"
    o_plot.gen_set["figsize"] = (16.8, 8.0)  # <->, |] # to match hists with (8.5, 2.7)
    o_plot.gen_set["figname"] = "DD2_ls220_slices.png"#"DD2_1512_slices.png" # LS_1412_slices
    o_plot.gen_set["sharex"] = False
    o_plot.gen_set["sharey"] = True
    o_plot.gen_set["dpi"] = 128
    o_plot.gen_set["subplots_adjust_h"] = -0.35
    o_plot.gen_set["subplots_adjust_w"] = 0.1
    o_plot.set_plot_dics = []

    ''' -----------------------------------------| DD2_M15091235_M0_LK_SR |----------------------------------------- '''

    # sim = "LS220_M14691268_M0_LK_SR" # "DD2_M15091235_M0_LK_SR"
    sim = "DD2_M15091235_M0_LK_SR"
    d3class = LOAD_PROFILE_XYXZ(sim)
    d1class = ADD_METHODS_1D(sim)
    rl = 3
    it = 1490944

    t = d3class.get_time_for_it(it, d1d2d3prof="prof")
    tmerg = d1class.get_par("tmerger_gw")
    xmin, xmax, ymin, ymax, zmin, zmax = get_xmin_xmax_ymin_ymax_zmin_zmax(rl)


    # --------------------------------------------
    mask = "x>0"

    v_n = "rho"

    rho_dic_xz = {'task': 'slice', 'dtype': '3d rl', 'ptype': 'cartesian', 'aspect': 1.,
                  'data': d3class, 'it': it, 'plane': 'xz', 'rl': rl,
                  'position': (1, 1),  # 'title': '[{:.1f} ms]'.format(time_),
                  'cbar': {},
                  'v_n_x': 'x', 'v_n_y': 'z', 'v_n': v_n,
                  'xmin': xmin, 'xmax': xmax, 'ymin': zmin, 'ymax': zmax, 'vmin': 1e-8, 'vmax': 1e-4,
                  'fill_vmin': False,  # fills the x < vmin with vmin
                  'xscale': None, 'yscale': None,
                  'mask': mask, 'cmap': 'viridis', 'norm': "log",
                  'fancyticks': True,
                  'title': {"text": r'$t-t_{merg}:$' + r'${:.1f}$ [ms]'.format((t - tmerg)*1e3), 'fontsize': 14},
                  'sharex': True,  # removes angular citkscitks
                  'fontsize': 14,
                  'labelsize': 14
                  }
    rho_dic_xy = {'task': 'slice', 'dtype': '3d rl', 'ptype': 'cartesian', 'aspect': 1.,
                  'data': d3class, 'it': it, 'plane': 'xy', 'rl': rl,
                  'position': (2, 1),  # 'title': '[{:.1f} ms]'.format(time_),
                  'cbar': {},
                  'v_n_x': 'x', 'v_n_y': 'y', 'v_n': v_n,
                  'xmin': xmin, 'xmax': xmax, 'ymin': ymin, 'ymax': ymax, 'vmin': 1e-8, 'vmax': 1e-4,
                  'fill_vmin': False,  # fills the x < vmin with vmin
                  'xscale': None, 'yscale': None,
                  'mask': mask, 'cmap': 'viridis', 'norm': "log",
                  'fancyticks': True,
                  'title': {},
                  'sharex': False,  # removes angular citkscitks
                  'fontsize': 14,
                  'labelsize': 14
                  }

    o_plot.set_plot_dics.append(rho_dic_xz)
    o_plot.set_plot_dics.append(rho_dic_xy)

    v_n = "dens_unb_bern"

    dunb_dic_xz = {'task': 'slice', 'dtype': '3d rl', 'ptype': 'cartesian', 'aspect': 1.,
                  'data': d3class, 'it': it, 'plane': 'xz', 'rl': rl,
                  'position': (1, 1),  # 'title': '[{:.1f} ms]'.format(time_),
                  'cbar': {},
                  'v_n_x': 'x', 'v_n_y': 'z', 'v_n': v_n,
                  'xmin': xmin, 'xmax': xmax, 'ymin': zmin, 'ymax': zmax, 'vmin': 1e-10, 'vmax': 1e-7,
                  'fill_vmin': False,  # fills the x < vmin with vmin
                  'xscale': None, 'yscale': None,
                  'mask': mask, 'cmap': 'Blues', 'norm': "log",
                  'fancyticks': True,
                  'title': {"text": r'$t-t_{merg}:$' + r'${:.1f}$ [ms]'.format((t - tmerg) * 1e3), 'fontsize': 14},
                  'sharex': True,  # removes angular citkscitks
                  'fontsize': 14,
                  'labelsize': 14
                  }
    dunb_dic_xy = {'task': 'slice', 'dtype': '3d rl', 'ptype': 'cartesian', 'aspect': 1.,
                  'data': d3class, 'it': it, 'plane': 'xy', 'rl': rl,
                  'position': (2, 1),  # 'title': '[{:.1f} ms]'.format(time_),
                  'cbar': {},
                  'v_n_x': 'x', 'v_n_y': 'y', 'v_n': v_n,
                  'xmin': xmin, 'xmax': xmax, 'ymin': ymin, 'ymax': ymax, 'vmin': 1e-10, 'vmax': 1e-7,
                  'fill_vmin': False,  # fills the x < vmin with vmin
                  'xscale': None, 'yscale': None,
                  'mask': mask, 'cmap': 'Blues', 'norm': "log",
                  'fancyticks': True,
                  'title': {},
                  'sharex': False,  # removes angular citkscitks
                  'fontsize': 14,
                  'labelsize': 14
                  }

    o_plot.set_plot_dics.append(dunb_dic_xz)
    o_plot.set_plot_dics.append(dunb_dic_xy)

    # ---------------------------------------------
    mask = "x<0"

    v_n = "Ye"

    ye_dic_xz = {'task': 'slice', 'dtype': '3d rl', 'ptype': 'cartesian', 'aspect': 1.,
                  'data': d3class, 'it': it, 'plane': 'xz', 'rl': rl,
                  'position': (1, 1),  # 'title': '[{:.1f} ms]'.format(time_),
                  'cbar': {},
                  'v_n_x': 'x', 'v_n_y': 'z', 'v_n': v_n,
                  'xmin': xmin, 'xmax': xmax, 'ymin': zmin, 'ymax': zmax, 'vmin': 0.05, 'vmax': 0.40,
                  'fill_vmin': False,  # fills the x < vmin with vmin
                  'xscale': None, 'yscale': None,
                  'mask': mask, 'cmap': 'RdBu', 'norm': None,
                  'fancyticks': True,
                  'title': {"text": r'$t-t_{merg}:$' + r'${:.1f}$ [ms]'.format((t - tmerg)*1e3), 'fontsize': 14},
                  'sharex': True,  # removes angular citkscitks
                  'fontsize': 14,
                  'labelsize': 14
                  }
    ye_dic_xy = {'task': 'slice', 'dtype': '3d rl', 'ptype': 'cartesian', 'aspect': 1.,
                  'data': d3class, 'it': it, 'plane': 'xy', 'rl': rl,
                  'position': (2, 1),  # 'title': '[{:.1f} ms]'.format(time_),
                  'cbar':{},
                  # 'cbar': {'location': 'right .06 .05', 'label': r'$Y_e$',  # 'fmt': '%.1e',
                  #          'labelsize': 14,
                  #          'fontsize': 14},
                  'v_n_x': 'x', 'v_n_y': 'y', 'v_n': v_n,
                  'xmin': xmin, 'xmax': xmax, 'ymin': ymin, 'ymax': ymax, 'vmin': 0.05, 'vmax': 0.40,
                  'fill_vmin': False,  # fills the x < vmin with vmin
                  'xscale': None, 'yscale': None,
                  'mask': mask, 'cmap': 'RdBu', 'norm': None,
                  'fancyticks': True,
                  'title': {},
                  'sharex': False,  # removes angular citkscitks
                  'fontsize': 14,
                  'labelsize': 14
                  }

    o_plot.set_plot_dics.append(ye_dic_xz)
    o_plot.set_plot_dics.append(ye_dic_xy)



    ''' -----------------------------------------| LS220_M14691268_M0_LK_SR |----------------------------------------- '''

    # sim = "LS220_M14691268_M0_LK_SR" # "DD2_M15091235_M0_LK_SR"
    sim = "LS220_M14691268_M0_LK_SR"
    d3class = LOAD_PROFILE_XYXZ(sim)
    d1class = ADD_METHODS_1D(sim)
    rl = 3
    it = 1490944

    t = d3class.get_time_for_it(it, d1d2d3prof="prof")
    tmerg = d1class.get_par("tmerger_gw")
    xmin, xmax, ymin, ymax, zmin, zmax = get_xmin_xmax_ymin_ymax_zmin_zmax(rl)

    # --------------------------------------------
    mask = "x>0"

    v_n = "rho"

    rho_dic_xz = {'task': 'slice', 'dtype': '3d rl', 'ptype': 'cartesian', 'aspect': 1.,
                  'data': d3class, 'it': it, 'plane': 'xz', 'rl': rl,
                  'position': (1, 2),  # 'title': '[{:.1f} ms]'.format(time_),
                  'cbar': {},
                  'v_n_x': 'x', 'v_n_y': 'z', 'v_n': v_n,
                  'xmin': xmin, 'xmax': xmax, 'ymin': zmin, 'ymax': zmax, 'vmin': 1e-8, 'vmax': 1e-4,
                  'fill_vmin': False,  # fills the x < vmin with vmin
                  'xscale': None, 'yscale': None,
                  'mask': mask, 'cmap': 'viridis', 'norm': "log",
                  'fancyticks': True,
                  'title': {"text": r'$t-t_{merg}:$' + r'${:.1f}$ [ms]'.format((t - tmerg) * 1e3), 'fontsize': 14},
                  'sharex': True,  # removes angular citkscitks
                  'fontsize': 14,
                  'labelsize': 14
                  }
    rho_dic_xy = {'task': 'slice', 'dtype': '3d rl', 'ptype': 'cartesian', 'aspect': 1.,
                  'data': d3class, 'it': it, 'plane': 'xy', 'rl': rl,
                  'position': (2, 2),  # 'title': '[{:.1f} ms]'.format(time_),
                  'cbar': {},
                  'v_n_x': 'x', 'v_n_y': 'y', 'v_n': v_n,
                  'xmin': xmin, 'xmax': xmax, 'ymin': ymin, 'ymax': ymax, 'vmin': 1e-8, 'vmax': 1e-4,
                  'fill_vmin': False,  # fills the x < vmin with vmin
                  'xscale': None, 'yscale': None,
                  'mask': mask, 'cmap': 'viridis', 'norm': "log",
                  'fancyticks': True,
                  'title': {},
                  'sharex': False,  # removes angular citkscitks
                  'fontsize': 14,
                  'labelsize': 14,
                  }

    o_plot.set_plot_dics.append(rho_dic_xz)
    o_plot.set_plot_dics.append(rho_dic_xy)

    v_n = "dens_unb_bern"

    dunb_dic_xz = {'task': 'slice', 'dtype': '3d rl', 'ptype': 'cartesian', 'aspect': 1.,
                   'data': d3class, 'it': it, 'plane': 'xz', 'rl': rl,
                   'position': (1, 2),  # 'title': '[{:.1f} ms]'.format(time_),
                   'cbar': {},
                   'v_n_x': 'x', 'v_n_y': 'z', 'v_n': v_n,
                   'xmin': xmin, 'xmax': xmax, 'ymin': zmin, 'ymax': zmax, 'vmin': 1e-10, 'vmax': 1e-7,
                   'fill_vmin': False,  # fills the x < vmin with vmin
                   'xscale': None, 'yscale': None,
                   'mask': mask, 'cmap': 'Blues', 'norm': "log",
                   'fancyticks': True,
                   'title': {"text": r'$t-t_{merg}:$' + r'${:.1f}$ [ms]'.format((t - tmerg) * 1e3), 'fontsize': 14},
                   'sharex': True,  # removes angular citkscitks
                   'fontsize': 14,
                   'labelsize': 14
                   }
    dunb_dic_xy = {'task': 'slice', 'dtype': '3d rl', 'ptype': 'cartesian', 'aspect': 1.,
                   'data': d3class, 'it': it, 'plane': 'xy', 'rl': rl,
                   'position': (2, 2),  # 'title': '[{:.1f} ms]'.format(time_),
                   'cbar': {},
                   'v_n_x': 'x', 'v_n_y': 'y', 'v_n': v_n,
                   'xmin': xmin, 'xmax': xmax, 'ymin': ymin, 'ymax': ymax, 'vmin': 1e-10, 'vmax': 1e-7,
                   'fill_vmin': False,  # fills the x < vmin with vmin
                   'xscale': None, 'yscale': None,
                   'mask': mask, 'cmap': 'Blues', 'norm': "log",
                   'fancyticks': True,
                   'title': {},
                   'sharex': False,  # removes angular citkscitks
                   'fontsize': 14,
                   'labelsize': 14
                   }

    o_plot.set_plot_dics.append(dunb_dic_xz)
    o_plot.set_plot_dics.append(dunb_dic_xy)

    # ---------------------------------------------
    mask = "x<0"

    v_n = "Ye"

    ye_dic_xz = {'task': 'slice', 'dtype': '3d rl', 'ptype': 'cartesian', 'aspect': 1.,
                 'data': d3class, 'it': it, 'plane': 'xz', 'rl': rl,
                 'position': (1, 2),  # 'title': '[{:.1f} ms]'.format(time_),
                 'cbar': {},
                 'v_n_x': 'x', 'v_n_y': 'z', 'v_n': v_n,
                 'xmin': xmin, 'xmax': xmax, 'ymin': zmin, 'ymax': zmax, 'vmin': 0.05, 'vmax': 0.40,
                 'fill_vmin': False,  # fills the x < vmin with vmin
                 'xscale': None, 'yscale': None,
                 'mask': mask, 'cmap': 'RdBu', 'norm': None,
                 'fancyticks': True,
                 'title': {"text": r'$t-t_{merg}:$' + r'${:.1f}$ [ms]'.format((t - tmerg) * 1e3), 'fontsize': 14},
                 'sharex': True,  # removes angular citkscitks
                 'sharey': True,
                 'fontsize': 14,
                 'labelsize': 14
                 }
    ye_dic_xy = {'task': 'slice', 'dtype': '3d rl', 'ptype': 'cartesian', 'aspect': 1.,
                 'data': d3class, 'it': it, 'plane': 'xy', 'rl': rl,
                 'position': (2, 2),  # 'title': '[{:.1f} ms]'.format(time_),
                 'cbar': {},
                 # 'cbar': {'location': 'right .06 .05', 'label': r'$Y_e$',  # 'fmt': '%.1e',
                 #          'labelsize': 14,
                 #          'fontsize': 14},
                 'v_n_x': 'x', 'v_n_y': 'y', 'v_n': v_n,
                 'xmin': xmin, 'xmax': xmax, 'ymin': ymin, 'ymax': ymax, 'vmin': 0.05, 'vmax': 0.40,
                 'fill_vmin': False,  # fills the x < vmin with vmin
                 'xscale': None, 'yscale': None,
                 'mask': mask, 'cmap': 'RdBu', 'norm': None,
                 'fancyticks': True,
                 'title': {},
                 'sharex': False,  # removes angular citkscitks
                 'sharey': True,
                 'fontsize': 14,
                 'labelsize': 14
                 }

    o_plot.set_plot_dics.append(ye_dic_xz)
    o_plot.set_plot_dics.append(ye_dic_xy)


    ''' -----------------------------------------| LS220_M13641364_M0_LK_SR |----------------------------------------- '''

    # sim = "LS220_M14691268_M0_LK_SR" # "DD2_M15091235_M0_LK_SR"
    sim = "LS220_M13641364_M0_LK_SR_restart"
    d3class = LOAD_PROFILE_XYXZ(sim)
    d1class = ADD_METHODS_1D(sim)
    rl = 3
    it = 696320

    t = d3class.get_time_for_it(it, d1d2d3prof="prof")
    tmerg = d1class.get_par("tmerger_gw")
    xmin, xmax, ymin, ymax, zmin, zmax = get_xmin_xmax_ymin_ymax_zmin_zmax(rl)

    # --------------------------------------------
    mask = "x>0"

    v_n = "rho"

    rho_dic_xz = {'task': 'slice', 'dtype': '3d rl', 'ptype': 'cartesian', 'aspect': 1.,
                  'data': d3class, 'it': it, 'plane': 'xz', 'rl': rl,
                  'position': (1, 3),  # 'title': '[{:.1f} ms]'.format(time_),
                  'cbar': {},
                  'v_n_x': 'x', 'v_n_y': 'z', 'v_n': v_n,
                  'xmin': xmin, 'xmax': xmax, 'ymin': zmin, 'ymax': zmax, 'vmin': 1e-8, 'vmax': 1e-4,
                  'fill_vmin': False,  # fills the x < vmin with vmin
                  'xscale': None, 'yscale': None,
                  'mask': mask, 'cmap': 'viridis', 'norm': "log",
                  'fancyticks': True,
                  'title': {"text": r'$t-t_{merg}:$' + r'${:.1f}$ [ms]'.format((t - tmerg) * 1e3), 'fontsize': 14},
                  'sharex': True,  # removes angular citkscitks
                  'fontsize': 14,
                  'labelsize': 14
                  }
    rho_dic_xy = {'task': 'slice', 'dtype': '3d rl', 'ptype': 'cartesian', 'aspect': 1.,
                  'data': d3class, 'it': it, 'plane': 'xy', 'rl': rl,
                  'position': (2, 3),  # 'title': '[{:.1f} ms]'.format(time_),
                  'cbar': {},
                  'v_n_x': 'x', 'v_n_y': 'y', 'v_n': v_n,
                  'xmin': xmin, 'xmax': xmax, 'ymin': ymin, 'ymax': ymax, 'vmin': 1e-8, 'vmax': 1e-4,
                  'fill_vmin': False,  # fills the x < vmin with vmin
                  'xscale': None, 'yscale': None,
                  'mask': mask, 'cmap': 'viridis', 'norm': "log",
                  'fancyticks': True,
                  'title': {},
                  'sharex': False,  # removes angular citkscitks
                  'fontsize': 14,
                  'labelsize': 14
                  }

    o_plot.set_plot_dics.append(rho_dic_xz)
    o_plot.set_plot_dics.append(rho_dic_xy)

    v_n = "dens_unb_bern"

    dunb_dic_xz = {'task': 'slice', 'dtype': '3d rl', 'ptype': 'cartesian', 'aspect': 1.,
                   'data': d3class, 'it': it, 'plane': 'xz', 'rl': rl,
                   'position': (1, 3),  # 'title': '[{:.1f} ms]'.format(time_),
                   'cbar': {},
                   'v_n_x': 'x', 'v_n_y': 'z', 'v_n': v_n,
                   'xmin': xmin, 'xmax': xmax, 'ymin': zmin, 'ymax': zmax, 'vmin': 1e-10, 'vmax': 1e-7,
                   'fill_vmin': False,  # fills the x < vmin with vmin
                   'xscale': None, 'yscale': None,
                   'mask': mask, 'cmap': 'Blues', 'norm': "log",
                   'fancyticks': True,
                   'title': {"text": r'$t-t_{merg}:$' + r'${:.1f}$ [ms]'.format((t - tmerg) * 1e3), 'fontsize': 14},
                   'sharex': True,  # removes angular citkscitks
                   'fontsize': 14,
                   'labelsize': 14
                   }
    dunb_dic_xy = {'task': 'slice', 'dtype': '3d rl', 'ptype': 'cartesian', 'aspect': 1.,
                   'data': d3class, 'it': it, 'plane': 'xy', 'rl': rl,
                   'position': (2, 3),  # 'title': '[{:.1f} ms]'.format(time_),
                   'cbar': {'location': 'bottom -.05 .00', 'label': r'$D_{\rm{unb}}$',  # 'fmt': '%.1e',
                            'labelsize': 14,
                            'fontsize': 14},
                   'v_n_x': 'x', 'v_n_y': 'y', 'v_n': v_n,
                   'xmin': xmin, 'xmax': xmax, 'ymin': ymin, 'ymax': ymax, 'vmin': 1e-10, 'vmax': 1e-7,
                   'fill_vmin': False,  # fills the x < vmin with vmin
                   'xscale': None, 'yscale': None,
                   'mask': mask, 'cmap': 'Blues', 'norm': "log",
                   'fancyticks': True,
                   'title': {},
                   'sharex': False,  # removes angular citkscitks
                   'fontsize': 14,
                   'labelsize': 14
                   }

    o_plot.set_plot_dics.append(dunb_dic_xz)
    o_plot.set_plot_dics.append(dunb_dic_xy)

    # ---------------------------------------------
    mask = "x<0"

    v_n = "Ye"

    ye_dic_xz = {'task': 'slice', 'dtype': '3d rl', 'ptype': 'cartesian', 'aspect': 1.,
                 'data': d3class, 'it': it, 'plane': 'xz', 'rl': rl,
                 'position': (1, 3),  # 'title': '[{:.1f} ms]'.format(time_),
                 'cbar': {},
                 'v_n_x': 'x', 'v_n_y': 'z', 'v_n': v_n,
                 'xmin': xmin, 'xmax': xmax, 'ymin': zmin, 'ymax': zmax, 'vmin': 0.05, 'vmax': 0.40,
                 'fill_vmin': False,  # fills the x < vmin with vmin
                 'xscale': None, 'yscale': None,
                 'mask': mask, 'cmap': 'RdBu', 'norm': None,
                 'fancyticks': True,
                 'title': {"text": r'$t-t_{merg}:$' + r'${:.1f}$ [ms]'.format((t - tmerg) * 1e3), 'fontsize': 14},
                 'sharex': True,  # removes angular citkscitks
                 'sharey': True,
                 'fontsize': 14,
                 'labelsize': 14
                 }
    ye_dic_xy = {'task': 'slice', 'dtype': '3d rl', 'ptype': 'cartesian', 'aspect': 1.,
                 'data': d3class, 'it': it, 'plane': 'xy', 'rl': rl,
                 'position': (2, 3),  # 'title': '[{:.1f} ms]'.format(time_),
                 'cbar': {},
                 # 'cbar': {'location': 'right .06 .05', 'label': r'$Y_e$',  # 'fmt': '%.1e',
                 #          'labelsize': 14,
                 #          'fontsize': 14},
                 'v_n_x': 'x', 'v_n_y': 'y', 'v_n': v_n,
                 'xmin': xmin, 'xmax': xmax, 'ymin': ymin, 'ymax': ymax, 'vmin': 0.05, 'vmax': 0.40,
                 'fill_vmin': False,  # fills the x < vmin with vmin
                 'xscale': None, 'yscale': None,
                 'mask': mask, 'cmap': 'RdBu', 'norm': None,
                 'fancyticks': True,
                 'title': {},
                 'sharex': False,  # removes angular citkscitk
                 'sharey': True,
                 'fontsize': 14,
                 'labelsize': 14
                 }

    o_plot.set_plot_dics.append(ye_dic_xz)
    o_plot.set_plot_dics.append(ye_dic_xy)

    ''' -----------------------------------------| DD2_M13641364_M0_LK_SR |----------------------------------------- '''

    # sim = "LS220_M14691268_M0_LK_SR" # "DD2_M15091235_M0_LK_SR"
    sim = "DD2_M13641364_M0_LK_SR_R04"
    d3class = LOAD_PROFILE_XYXZ(sim)
    d1class = ADD_METHODS_1D(sim)
    rl = 3
    it = 1581030

    t = d3class.get_time_for_it(it, d1d2d3prof="prof")
    tmerg = d1class.get_par("tmerger_gw")
    xmin, xmax, ymin, ymax, zmin, zmax = get_xmin_xmax_ymin_ymax_zmin_zmax(rl)

    # --------------------------------------------
    mask = "x>0"

    v_n = "rho"

    rho_dic_xz = {'task': 'slice', 'dtype': '3d rl', 'ptype': 'cartesian', 'aspect': 1.,
                  'data': d3class, 'it': it, 'plane': 'xz', 'rl': rl,
                  'position': (1, 4),  # 'title': '[{:.1f} ms]'.format(time_),
                  'cbar': {},
                  'v_n_x': 'x', 'v_n_y': 'z', 'v_n': v_n,
                  'xmin': xmin, 'xmax': xmax, 'ymin': zmin, 'ymax': zmax, 'vmin': 1e-8, 'vmax': 1e-4,
                  'fill_vmin': False,  # fills the x < vmin with vmin
                  'xscale': None, 'yscale': None,
                  'mask': mask, 'cmap': 'viridis', 'norm': "log",
                  'fancyticks': True,
                  'title': {"text": r'$t-t_{merg}:$' + r'${:.1f}$ [ms]'.format((t - tmerg) * 1e3), 'fontsize': 14},
                  'sharex': True,  # removes angular citkscitks
                  'fontsize': 14,
                  'labelsize': 14
                  }
    rho_dic_xy = {'task': 'slice', 'dtype': '3d rl', 'ptype': 'cartesian', 'aspect': 1.,
                  'data': d3class, 'it': it, 'plane': 'xy', 'rl': rl,
                  'position': (2, 4),  # 'title': '[{:.1f} ms]'.format(time_),
                  'cbar': {},
                  'v_n_x': 'x', 'v_n_y': 'y', 'v_n': v_n,
                  'xmin': xmin, 'xmax': xmax, 'ymin': ymin, 'ymax': ymax, 'vmin': 1e-8, 'vmax': 1e-4,
                  'fill_vmin': False,  # fills the x < vmin with vmin
                  'xscale': None, 'yscale': None,
                  'mask': mask, 'cmap': 'viridis', 'norm': "log",
                  'fancyticks': True,
                  'title': {},
                  'sharex': False,  # removes angular citkscitks
                  'fontsize': 14,
                  'labelsize': 14
                  }

    o_plot.set_plot_dics.append(rho_dic_xz)
    o_plot.set_plot_dics.append(rho_dic_xy)

    v_n = "dens_unb_bern"

    dunb_dic_xz = {'task': 'slice', 'dtype': '3d rl', 'ptype': 'cartesian', 'aspect': 1.,
                   'data': d3class, 'it': it, 'plane': 'xz', 'rl': rl,
                   'position': (1, 4),  # 'title': '[{:.1f} ms]'.format(time_),
                   'cbar': {},
                   'v_n_x': 'x', 'v_n_y': 'z', 'v_n': v_n,
                   'xmin': xmin, 'xmax': xmax, 'ymin': zmin, 'ymax': zmax, 'vmin': 1e-10, 'vmax': 1e-7,
                   'fill_vmin': False,  # fills the x < vmin with vmin
                   'xscale': None, 'yscale': None,
                   'mask': mask, 'cmap': 'Blues', 'norm': "log",
                   'fancyticks': True,
                   'title': {"text": r'$t-t_{merg}:$' + r'${:.1f}$ [ms]'.format((t - tmerg) * 1e3), 'fontsize': 14},
                   'sharex': True,  # removes angular citkscitks
                   'fontsize': 14,
                   'labelsize': 14
                   }
    dunb_dic_xy = {'task': 'slice', 'dtype': '3d rl', 'ptype': 'cartesian', 'aspect': 1.,
                   'data': d3class, 'it': it, 'plane': 'xy', 'rl': rl,
                   'position': (2, 4),  # 'title': '[{:.1f} ms]'.format(time_),
                   'cbar': {},
                   'v_n_x': 'x', 'v_n_y': 'y', 'v_n': v_n,
                   'xmin': xmin, 'xmax': xmax, 'ymin': ymin, 'ymax': ymax, 'vmin': 1e-10, 'vmax': 1e-7,
                   'fill_vmin': False,  # fills the x < vmin with vmin
                   'xscale': None, 'yscale': None,
                   'mask': mask, 'cmap': 'Blues', 'norm': "log",
                   'fancyticks': True,
                   'title': {},
                   'sharex': False,  # removes angular citkscitks
                   'fontsize': 14,
                   'labelsize': 14
                   }

    o_plot.set_plot_dics.append(dunb_dic_xz)
    o_plot.set_plot_dics.append(dunb_dic_xy)

    # ---------------------------------------------
    mask = "x<0"

    v_n = "Ye"

    ye_dic_xz = {'task': 'slice', 'dtype': '3d rl', 'ptype': 'cartesian', 'aspect': 1.,
                 'data': d3class, 'it': it, 'plane': 'xz', 'rl': rl,
                 'position': (1, 4),  # 'title': '[{:.1f} ms]'.format(time_),
                 'cbar': {},
                 'v_n_x': 'x', 'v_n_y': 'z', 'v_n': v_n,
                 'xmin': xmin, 'xmax': xmax, 'ymin': zmin, 'ymax': zmax, 'vmin': 0.05, 'vmax': 0.40,
                 'fill_vmin': False,  # fills the x < vmin with vmin
                 'xscale': None, 'yscale': None,
                 'mask': mask, 'cmap': 'RdBu', 'norm': None,
                 'fancyticks': True,
                 'title': {"text": r'$t-t_{merg}:$' + r'${:.1f}$ [ms]'.format((t - tmerg) * 1e3), 'fontsize': 14},
                 'sharex': True,  # removes angular citkscitks
                 'sharey': True,
                 'fontsize': 14,
                 'labelsize': 14
                 }
    ye_dic_xy = {'task': 'slice', 'dtype': '3d rl', 'ptype': 'cartesian', 'aspect': 1.,
                 'data': d3class, 'it': it, 'plane': 'xy', 'rl': rl,
                 'position': (2, 4),  # 'title': '[{:.1f} ms]'.format(time_),
                 # 'cbar': {},
                 'cbar': {'location': 'bottom -.05 .00', 'label': r'$Y_e$',  # 'fmt': '%.1e',
                          'labelsize': 14,
                          'fontsize': 14},
                 'v_n_x': 'x', 'v_n_y': 'y', 'v_n': v_n,
                 'xmin': xmin, 'xmax': xmax, 'ymin': ymin, 'ymax': ymax, 'vmin': 0.05, 'vmax': 0.40,
                 'fill_vmin': False,  # fills the x < vmin with vmin
                 'xscale': None, 'yscale': None,
                 'mask': mask, 'cmap': 'RdBu', 'norm': None,
                 'fancyticks': True,
                 'title': {},
                 'sharex': False,  # removes angular citkscitks
                 'sharey': True,
                 'fontsize': 14,
                 'labelsize': 14
                 }

    o_plot.set_plot_dics.append(ye_dic_xz)
    o_plot.set_plot_dics.append(ye_dic_xy)


    # o_plot.set_plot_dics.append(def_dic_xz)
    # o_plot.set_plot_dics.append(def_dic_xy)

    o_plot.main()

    exit(0)
# plot_den_unb__vel_z()


def plot_den_unb__vel_z2():

    # tmp = d3class.get_data(688128, 3, "xy", "ang_mom_flux")
    # print(tmp.min(), tmp.max())
    # print(tmp)
    # exit(1) # dens_unb_geo

    """ --- --- --- """


    o_plot = PLOT_MANY_TASKS()
    o_plot.gen_set["figdir"] = Paths.plots + 'all/'
    o_plot.gen_set["type"] = "cartesian"
    o_plot.gen_set["figsize"] = (16.8, 8.0)  # <->, |] # to match hists with (8.5, 2.7)
    o_plot.gen_set["figname"] = "DD2_ls220_slices.png"#"DD2_1512_slices.png" # LS_1412_slices
    o_plot.gen_set["sharex"] = False
    o_plot.gen_set["sharey"] = True
    o_plot.gen_set["dpi"] = 128
    o_plot.gen_set["subplots_adjust_h"] = -0.35
    o_plot.gen_set["subplots_adjust_w"] = 0.05
    o_plot.set_plot_dics = []

    ''' -----------------------------------------------| SETTINGS |--------------------------------------------------'''

    sims = ["DD2_M15091235_M0_LK_SR", "DD2_M13641364_M0_LK_SR_R04", "LS220_M14691268_M0_LK_SR", "LS220_M13641364_M0_LK_SR_restart"]
    iterations = [1490944, 1581030, 1490944, 696320]
    rls = [3, 3, 3, 3]



    ''' ------------------------------------------------------------------------------------------------------------ '''
    n = 1
    for sim, it, rl in zip(sims, iterations, rls):
        d3class = LOAD_PROFILE_XYXZ(sim)
        d1class = ADD_METHODS_1D(sim)
        t = d3class.get_time_for_it(it, d1d2d3prof="prof")
        tmerg = d1class.get_par("tmerger_gw")
        xmin, xmax, ymin, ymax, zmin, zmax = get_xmin_xmax_ymin_ymax_zmin_zmax(rl)

        # ------------------ left side -----------------
        mask = "x>0"

        v_n = "rho"

        rho_dic_xz = {'task': 'slice', 'dtype': '3d rl', 'ptype': 'cartesian', 'aspect': 1.,
                      'data': d3class, 'it': it, 'plane': 'xz', 'rl': rl,
                      'position': (1, n),  # 'title': '[{:.1f} ms]'.format(time_),
                      'cbar': {},
                      'v_n_x': 'x', 'v_n_y': 'z', 'v_n': v_n,
                      'xmin': xmin, 'xmax': xmax, 'ymin': zmin, 'ymax': zmax, 'vmin': 1e-8, 'vmax': 1e-4,
                      'fill_vmin': False,  # fills the x < vmin with vmin
                      'xscale': None, 'yscale': None,
                      'mask': mask, 'cmap': 'viridis', 'norm': "log",
                      'fancyticks': True,
                      'title': {"text": r'$t-t_{merg}:$' + r'${:.1f}$ [ms]'.format((t - tmerg) * 1e3), 'fontsize': 14},
                      'sharex': True,  # removes angular citkscitks
                      'sharey': False,
                      'fontsize': 14,
                      'labelsize': 14
                      }
        rho_dic_xy = {'task': 'slice', 'dtype': '3d rl', 'ptype': 'cartesian', 'aspect': 1.,
                      'data': d3class, 'it': it, 'plane': 'xy', 'rl': rl,
                      'position': (2, n),  # 'title': '[{:.1f} ms]'.format(time_),
                      'cbar': {},
                      'v_n_x': 'x', 'v_n_y': 'y', 'v_n': v_n,
                      'xmin': xmin, 'xmax': xmax, 'ymin': ymin, 'ymax': ymax, 'vmin': 1e-8, 'vmax': 1e-4,
                      'fill_vmin': False,  # fills the x < vmin with vmin
                      'xscale': None, 'yscale': None,
                      'mask': mask, 'cmap': 'viridis', 'norm': "log",
                      'fancyticks': True,
                      'title': {},
                      'sharex': False,  # removes angular citkscitks
                      'sharey': False,
                      'fontsize': 14,
                      'labelsize': 14
                      }
        o_plot.set_plot_dics.append(rho_dic_xz)
        o_plot.set_plot_dics.append(rho_dic_xy)

        v_n = "dens_unb_bern"

        dunb_dic_xz = {'task': 'slice', 'dtype': '3d rl', 'ptype': 'cartesian', 'aspect': 1.,
                       'data': d3class, 'it': it, 'plane': 'xz', 'rl': rl,
                       'position': (1, n),  # 'title': '[{:.1f} ms]'.format(time_),
                       'cbar': {},
                       'v_n_x': 'x', 'v_n_y': 'z', 'v_n': v_n,
                       'xmin': xmin, 'xmax': xmax, 'ymin': zmin, 'ymax': zmax, 'vmin': 1e-10, 'vmax': 1e-7,
                       'fill_vmin': False,  # fills the x < vmin with vmin
                       'xscale': None, 'yscale': None,
                       'mask': mask, 'cmap': 'Blues', 'norm': "log",
                       'fancyticks': True,
                       'title': {"text": r'$t-t_{merg}:$' + r'${:.1f}$ [ms]'.format((t - tmerg) * 1e3), 'fontsize': 14},
                       'sharex': True,  # removes angular citkscitks
                       'sharey': False,
                       'fontsize': 14,
                       'labelsize': 14
                       }
        dunb_dic_xy = {'task': 'slice', 'dtype': '3d rl', 'ptype': 'cartesian', 'aspect': 1.,
                       'data': d3class, 'it': it, 'plane': 'xy', 'rl': rl,
                       'position': (2, n),  # 'title': '[{:.1f} ms]'.format(time_),
                       'cbar': {},
                       'v_n_x': 'x', 'v_n_y': 'y', 'v_n': v_n,
                       'xmin': xmin, 'xmax': xmax, 'ymin': ymin, 'ymax': ymax, 'vmin': 1e-10, 'vmax': 1e-7,
                       'fill_vmin': False,  # fills the x < vmin with vmin
                       'xscale': None, 'yscale': None,
                       'mask': mask, 'cmap': 'Blues', 'norm': "log",
                       'fancyticks': True,
                       'title': {},
                       'sharex': False,  # removes angular citkscitks
                       'sharey': False,
                       'fontsize': 14,
                       'labelsize': 14
                       }

        o_plot.set_plot_dics.append(dunb_dic_xz)
        o_plot.set_plot_dics.append(dunb_dic_xy)

        # ------------------- right side -------------
        mask = "x<0"

        v_n = "Ye"

        ye_dic_xz = {'task': 'slice', 'dtype': '3d rl', 'ptype': 'cartesian', 'aspect': 1.,
                     'data': d3class, 'it': it, 'plane': 'xz', 'rl': rl,
                     'position': (1, n),  # 'title': '[{:.1f} ms]'.format(time_),
                     'cbar': {},
                     'v_n_x': 'x', 'v_n_y': 'z', 'v_n': v_n,
                     'xmin': xmin, 'xmax': xmax, 'ymin': zmin, 'ymax': zmax, 'vmin': 0.1, 'vmax': 0.5,
                     'fill_vmin': False,  # fills the x < vmin with vmin
                     'xscale': None, 'yscale': None,
                     'mask': mask, 'cmap': 'Reds_r', 'norm': None,
                     'fancyticks': True,
                     'title': {"text": r'$t-t_{merg}:$' + r'${:.1f}$ [ms]'.format((t - tmerg) * 1e3), 'fontsize': 14},
                     'sharex': True,  # removes angular citkscitks
                     'sharey': False,
                     'fontsize': 14,
                     'labelsize': 14
                     }
        ye_dic_xy = {'task': 'slice', 'dtype': '3d rl', 'ptype': 'cartesian', 'aspect': 1.,
                     'data': d3class, 'it': it, 'plane': 'xy', 'rl': rl,
                     'position': (2, n),  # 'title': '[{:.1f} ms]'.format(time_),
                     'cbar': {},
                     # 'cbar': {'location': 'right .06 .05', 'label': r'$Y_e$',  # 'fmt': '%.1e',
                     #          'labelsize': 14,
                     #          'fontsize': 14},
                     'v_n_x': 'x', 'v_n_y': 'y', 'v_n': v_n,
                     'xmin': xmin, 'xmax': xmax, 'ymin': ymin, 'ymax': ymax, 'vmin': 0.1, 'vmax': 0.5,
                     'fill_vmin': False,  # fills the x < vmin with vmin
                     'xscale': None, 'yscale': None,
                     'mask': mask, 'cmap': 'Reds_r', 'norm': None,
                     'fancyticks': True,
                     'title': {},
                     'sharex': False,  # removes angular citkscitks
                     'sharey': False,
                     'fontsize': 14,
                     'labelsize': 14
                     }

        o_plot.set_plot_dics.append(ye_dic_xz)
        o_plot.set_plot_dics.append(ye_dic_xy)

        # --------------------------------------------------
        if n > 1:
            rho_dic_xz['sharey'] = True
            rho_dic_xy['sharey'] = True
            dunb_dic_xz['sharey'] = True
            dunb_dic_xy['sharey'] = True
            ye_dic_xz['sharey'] = True
            ye_dic_xy['sharey'] = True

        if n == 1:
            rho_dic_xy['cbar'] = \
                {'location': 'bottom -.05 .00', 'label': r'$\rho$ [geo]',  # 'fmt': '%.1e',
                 'labelsize': 14, 'fontsize': 14}

        if n ==2:
            dunb_dic_xy['cbar'] = \
                {'location': 'bottom -.05 .00', 'label': r'$D_{\rm{unb}}$ [geo]',  # 'fmt': '%.1e',
                 'labelsize': 14, 'fontsize': 14}

        if n ==3:
            ye_dic_xy['cbar'] = \
                {'location': 'bottom -.05 .00', 'label': r'$Y_e$ [geo]',  'fmt': '%.1f',
                 'labelsize': 14, 'fontsize': 14}

        n = n + 1

    ''' ------------------------------------------------------------------------------------------------------------ '''
    o_plot.main()

    exit(0)
# plot_den_unb__vel_z2()

def plot_den_unb__vel_z__ls220_LK_MO():

    # tmp = d3class.get_data(688128, 3, "xy", "ang_mom_flux")
    # print(tmp.min(), tmp.max())
    # print(tmp)
    # exit(1) # dens_unb_geo

    """ --- --- --- """


    o_plot = PLOT_MANY_TASKS()
    o_plot.gen_set["figdir"] = Paths.plots + 'all/'
    o_plot.gen_set["type"] = "cartesian"
    o_plot.gen_set["figsize"] = (8.4, 8.0)  # <->, |] # to match hists with (8.5, 2.7)
    o_plot.gen_set["figname"] = "ls220_lk_m0.png"#"DD2_1512_slices.png" # LS_1412_slices
    o_plot.gen_set["sharex"] = False
    o_plot.gen_set["sharey"] = True
    o_plot.gen_set["dpi"] = 128
    o_plot.gen_set["subplots_adjust_h"] = -0.35
    o_plot.gen_set["subplots_adjust_w"] = 0.05
    o_plot.set_plot_dics = []

    ''' -----------------------------------------------| SETTINGS |--------------------------------------------------'''


    sims = ["LS220_M140120_LK", "LS220_M140120_M0"]
    iterations = [692398, 720294]
    rls = [2, 2]



    ''' ------------------------------------------------------------------------------------------------------------ '''
    n = 1
    for sim, it, rl in zip(sims, iterations, rls):
        d3class = LOAD_PROFILE_XYXZ(sim)
        d1class = ADD_METHODS_1D(sim)
        t = d3class.get_time_for_it(it, d1d2d3prof="prof")
        tmerg = 0#d1class.get_par("tmerger_gw")
        xmin, xmax, ymin, ymax, zmin, zmax = get_xmin_xmax_ymin_ymax_zmin_zmax(rl)

        # ------------------ left side -----------------
        mask = "x>0"

        v_n = "rho"

        rho_dic_xz = {'task': 'slice', 'dtype': '3d rl', 'ptype': 'cartesian', 'aspect': 1.,
                      'data': d3class, 'it': it, 'plane': 'xz', 'rl': rl,
                      'position': (1, n),  # 'title': '[{:.1f} ms]'.format(time_),
                      'cbar': {},
                      'v_n_x': 'x', 'v_n_y': 'z', 'v_n': v_n,
                      'xmin': xmin, 'xmax': xmax, 'ymin': zmin, 'ymax': zmax, 'vmin': 1e-8, 'vmax': 1e-4,
                      'fill_vmin': False,  # fills the x < vmin with vmin
                      'xscale': None, 'yscale': None,
                      'mask': mask, 'cmap': 'viridis', 'norm': "log",
                      'fancyticks': True,
                      'title': {"text": r'$t-t_{merg}:$' + r'${:.1f}$ [ms]'.format((t - tmerg) * 1e3), 'fontsize': 14},
                      'sharex': True,  # removes angular citkscitks
                      'sharey': False,
                      'fontsize': 14,
                      'labelsize': 14
                      }
        rho_dic_xy = {'task': 'slice', 'dtype': '3d rl', 'ptype': 'cartesian', 'aspect': 1.,
                      'data': d3class, 'it': it, 'plane': 'xy', 'rl': rl,
                      'position': (2, n),  # 'title': '[{:.1f} ms]'.format(time_),
                      'cbar': {},
                      'v_n_x': 'x', 'v_n_y': 'y', 'v_n': v_n,
                      'xmin': xmin, 'xmax': xmax, 'ymin': ymin, 'ymax': ymax, 'vmin': 1e-8, 'vmax': 1e-4,
                      'fill_vmin': False,  # fills the x < vmin with vmin
                      'xscale': None, 'yscale': None,
                      'mask': mask, 'cmap': 'viridis', 'norm': "log",
                      'fancyticks': True,
                      'title': {},
                      'sharex': False,  # removes angular citkscitks
                      'sharey': False,
                      'fontsize': 14,
                      'labelsize': 14
                      }
        # o_plot.set_plot_dics.append(rho_dic_xz)
        # o_plot.set_plot_dics.append(rho_dic_xy)

        v_n = "dens_unb_bern"

        dunb_dic_xz = {'task': 'slice', 'dtype': '3d rl', 'ptype': 'cartesian', 'aspect': 1.,
                       'data': d3class, 'it': it, 'plane': 'xz', 'rl': rl,
                       'position': (1, n),  # 'title': '[{:.1f} ms]'.format(time_),
                       'cbar': {},
                       'v_n_x': 'x', 'v_n_y': 'z', 'v_n': v_n,
                       'xmin': xmin, 'xmax': xmax, 'ymin': zmin, 'ymax': zmax, 'vmin': 1e-10, 'vmax': 1e-7,
                       'fill_vmin': False,  # fills the x < vmin with vmin
                       'xscale': None, 'yscale': None,
                       'mask': mask, 'cmap': 'Blues', 'norm': "log",
                       'fancyticks': True,
                       'title': {"text": r'$t-t_{merg}:$' + r'${:.1f}$ [ms]'.format((t - tmerg) * 1e3), 'fontsize': 14},
                       'sharex': True,  # removes angular citkscitks
                       'sharey': False,
                       'fontsize': 14,
                       'labelsize': 14
                       }
        dunb_dic_xy = {'task': 'slice', 'dtype': '3d rl', 'ptype': 'cartesian', 'aspect': 1.,
                       'data': d3class, 'it': it, 'plane': 'xy', 'rl': rl,
                       'position': (2, n),  # 'title': '[{:.1f} ms]'.format(time_),
                       'cbar': {},
                       'v_n_x': 'x', 'v_n_y': 'y', 'v_n': v_n,
                       'xmin': xmin, 'xmax': xmax, 'ymin': ymin, 'ymax': ymax, 'vmin': 1e-10, 'vmax': 1e-7,
                       'fill_vmin': False,  # fills the x < vmin with vmin
                       'xscale': None, 'yscale': None,
                       'mask': mask, 'cmap': 'Blues', 'norm': "log",
                       'fancyticks': True,
                       'title': {},
                       'sharex': False,  # removes angular citkscitks
                       'sharey': False,
                       'fontsize': 14,
                       'labelsize': 14
                       }

        o_plot.set_plot_dics.append(dunb_dic_xz)
        o_plot.set_plot_dics.append(dunb_dic_xy)

        # ------------------- right side -------------
        mask = "x<0"

        v_n = "Ye"

        ye_dic_xz = {'task': 'slice', 'dtype': '3d rl', 'ptype': 'cartesian', 'aspect': 1.,
                     'data': d3class, 'it': it, 'plane': 'xz', 'rl': rl,
                     'position': (1, n),  # 'title': '[{:.1f} ms]'.format(time_),
                     'cbar': {},
                     'v_n_x': 'x', 'v_n_y': 'z', 'v_n': v_n,
                     'xmin': xmin, 'xmax': xmax, 'ymin': zmin, 'ymax': zmax, 'vmin': 0.1, 'vmax': 0.5,
                     'fill_vmin': False,  # fills the x < vmin with vmin
                     'xscale': None, 'yscale': None,
                     'mask': mask, 'cmap': 'Reds_r', 'norm': None,
                     'fancyticks': True,
                     'title': {"text": r'$t-t_{merg}:$' + r'${:.1f}$ [ms]'.format((t - tmerg) * 1e3), 'fontsize': 14},
                     'sharex': True,  # removes angular citkscitks
                     'sharey': False,
                     'fontsize': 14,
                     'labelsize': 14
                     }
        ye_dic_xy = {'task': 'slice', 'dtype': '3d rl', 'ptype': 'cartesian', 'aspect': 1.,
                     'data': d3class, 'it': it, 'plane': 'xy', 'rl': rl,
                     'position': (2, n),  # 'title': '[{:.1f} ms]'.format(time_),
                     'cbar': {},
                     # 'cbar': {'location': 'right .06 .05', 'label': r'$Y_e$',  # 'fmt': '%.1e',
                     #          'labelsize': 14,
                     #          'fontsize': 14},
                     'v_n_x': 'x', 'v_n_y': 'y', 'v_n': v_n,
                     'xmin': xmin, 'xmax': xmax, 'ymin': ymin, 'ymax': ymax, 'vmin': 0.1, 'vmax': 0.5,
                     'fill_vmin': False,  # fills the x < vmin with vmin
                     'xscale': None, 'yscale': None,
                     'mask': mask, 'cmap': 'Reds_r', 'norm': None,
                     'fancyticks': True,
                     'title': {},
                     'sharex': False,  # removes angular citkscitks
                     'sharey': False,
                     'fontsize': 14,
                     'labelsize': 14
                     }

        o_plot.set_plot_dics.append(ye_dic_xz)
        o_plot.set_plot_dics.append(ye_dic_xy)

        # --------------------------------------------------
        if n > 1:
            rho_dic_xz['sharey'] = True
            rho_dic_xy['sharey'] = True
            dunb_dic_xz['sharey'] = True
            dunb_dic_xy['sharey'] = True
            ye_dic_xz['sharey'] = True
            ye_dic_xy['sharey'] = True

        # if n == 1:
        #     rho_dic_xy['cbar'] = \
        #         {'location': 'bottom -.05 .00', 'label': r'$\rho$ [geo]',  # 'fmt': '%.1e',
        #          'labelsize': 14, 'fontsize': 14}

        if n == 2:
            dunb_dic_xy['cbar'] = \
                {'location': 'bottom -.05 .00', 'label': r'$D_{\rm{unb}}$ [geo]',  # 'fmt': '%.1e',
                 'labelsize': 14, 'fontsize': 14}

        if n == 2:
            ye_dic_xy['cbar'] = \
                {'location': 'right .03 .00', 'label': r'$Y_e$ [geo]',  'fmt': '%.1f',
                 'labelsize': 14, 'fontsize': 14}

        n = n + 1

    ''' ------------------------------------------------------------------------------------------------------------ '''
    o_plot.main()

    exit(0)
# plot_den_unb__vel_z__ls220_LK_MO()

def plot_corr_ls220_LK_MO():

    o_plot = PLOT_MANY_TASKS()
    o_plot.gen_set["figdir"] = Paths.plots + '/all/'
    o_plot.gen_set["type"] = "cartesian"
    o_plot.gen_set["figsize"] = (8.0, 12.0)  # <->, |]
    o_plot.gen_set["figname"] = "corr_ls22_LK_MO_test2.png"
    o_plot.gen_set["sharex"] = False
    o_plot.gen_set["sharey"] = False
    o_plot.gen_set["dpi"] = 128
    o_plot.gen_set["subplots_adjust_h"] = 0.3
    o_plot.gen_set["subplots_adjust_w"] = 0.4
    o_plot.set_plot_dics = []
    ''' -----------------------------------------------| SETTINGS |--------------------------------------------------'''

    sims = ["LS220_M140120_LK", "LS220_M140120_M0"]
    iterations = [692398, 720294]

    ''' ------------------------------------------------------------------------------------------------------------ '''

    p = 0
    for sim, it in zip(sims, iterations):
        o_corr_data = LOAD_RES_CORR(sim)
        p = p + 1

        corr_dic = {  # relies on the "get_res_corr(self, it, v_n): " method of data object
            'task': '2d colormesh', 'dtype': 'corr', 'ptype': 'cartesian',
            'data': o_corr_data, 'it': it,
            'position': (1, p),
            'v_n_x': 'velz', 'v_n_y': 'Ye', 'v_n': 'mass',
            'xlabel': r"$\upsilon_{z}$", 'ylabel': r"$Y_e$",
            'cbar' : {},
            'xmin': -0.25, 'xmax': 0.25, 'ymin': 0.05, 'ymax': 0.5, 'vmin': 1e-9, 'vmax': 1e-6,
            'xscale': "linear", 'yscale': "linear",
            # 'mask_below': None, 'mask_above': None,
            'cmap': 'viridis_r',
            'norm': 'log', #'todo': None,
            # 'title': o_corr_data.sim.replace('_', '\_'),
            'title': {"text": o_corr_data.sim.replace('_', '\_'), 'fontsize': 14},
            'fancyticks': True,
            'minorticks': True,
            # 'sharex': False,  # removes angular citkscitks
            # 'sharey': False,
            'fontsize': 14,
            'labelsize': 14
        }
        # o_plot.set_plot_dics.append(corr_dic)

        corr_dic = {  # relies on the "get_res_corr(self, it, v_n): " method of data object
            'task': '2d colormesh', 'dtype': 'corr', 'ptype': 'cartesian',
            'data': o_corr_data, 'it': it,
            'position': (2, p),
            'v_n_x': 'velz', 'v_n_y': 'theta', 'v_n': 'mass',
            'xlabel': r"$\upsilon_{z}$", 'ylabel': r"$\theta$",
            'cbar' : {},
            'xmin': -0.25, 'xmax': 0.25, 'ymin': 0.0, 'ymax': 90, 'vmin': 1e-9, 'vmax': 1e-6,
            'xscale': "linear", 'yscale': "linear",
            # 'mask_below': None, 'mask_above': None,
            'cmap': 'viridis_r',
            'norm': 'log', #'todo': None,
            # 'title': o_corr_data.sim.replace('_', '\_'),
            'title': {"text": o_corr_data.sim.replace('_', '\_'), 'fontsize': 14},
            'fancyticks': True,
            'minorticks': True,
            # 'sharex': False,  # removes angular citkscitks
            # 'sharey': False,
            'fontsize': 14,
            'labelsize': 14
        }
        # o_plot.set_plot_dics.append(corr_dic)

        corr_dic2 = {  # relies on the "get_res_corr(self, it, v_n): " method of data object
            'task': '2d colormesh', 'dtype': 'corr', 'ptype': 'cartesian',
            'data': o_corr_data, 'it': it,
            'position': (1, p),
            'v_n_x': 'dens_unb_bern', 'v_n_y': 'Ye', 'v_n': 'mass',
            'cbar': {},
            'xlabel': r"$D_{\rm{unb}}$", 'ylabel': r"$Y_e$",
            'xmin': 1e-8, 'xmax': 1e-6, 'ymin': 0.05, 'ymax': 0.5, 'vmin': 1e-9, 'vmax': 1e-6,
            'xscale': "log", 'yscale': "linear",
            # 'mask_below': None, 'mask_above': None,
            'cmap': 'viridis_r',
            'norm': 'log', #'todo': None,
            # 'title': o_corr_data.sim.replace('_', '\_'),
            'title': {"text": o_corr_data.sim.replace('_', '\_'), 'fontsize': 14}, # "text": o_corr_data.sim.replace('_', '\_'), 'fontsize': 14
            'fancyticks': True,
            'minorticks': True,
            # 'sharex': False,  # removes angular citkscitks
            # 'sharey': False,
            'fontsize': 14,
            'labelsize': 14
        }
        o_plot.set_plot_dics.append(corr_dic2)

        corr_dic3 = {  # relies on the "get_res_corr(self, it, v_n): " method of data object
            'task': '2d colormesh', 'dtype': 'corr', 'ptype': 'cartesian',
            'data': o_corr_data, 'it': it,
            'position': (2, p),
            'v_n_x': 'dens_unb_bern', 'v_n_y': 'velz', 'v_n': 'mass',
            'cbar': {},
            'xlabel': r"$D_{\rm{unb}}$", 'ylabel': r"$\upsilon_{z}$",
            'xmin': 1e-8, 'xmax': 1e-6, 'ymin': -0.25, 'ymax': 0.25, 'vmin': 1e-9, 'vmax': 1e-6,
            'xscale': "log", 'yscale': "linear",
            # 'mask_below': None, 'mask_above': None,
            'cmap': 'viridis_r',
            'norm': 'log', #'todo': None,
            # 'title': o_corr_data.sim.replace('_', '\_'),
            'title': {}, # "text": o_corr_data.sim.replace('_', '\_'), 'fontsize': 14
            'fancyticks': True,
            'minorticks': True,
            # 'sharex': False,  # removes angular citkscitks
            # 'sharey': False,
            'fontsize': 14,
            'labelsize': 14
        }
        o_plot.set_plot_dics.append(corr_dic3)

        corr_dic3 = {  # relies on the "get_res_corr(self, it, v_n): " method of data object
            'task': '2d colormesh', 'dtype': 'corr', 'ptype': 'cartesian',
            'data': o_corr_data, 'it': it,
            'position': (3, p),
            'v_n_x': 'dens_unb_bern', 'v_n_y': 'theta', 'v_n': 'mass',
            'cbar': {},
            'xlabel': r"$D_{\rm{unb}}$", 'ylabel': r"$\theta$",
            'xmin': 1e-8, 'xmax': 1e-6, 'ymin': 0, 'ymax': 90, 'vmin': 1e-9, 'vmax': 1e-6,
            'xscale': "log", 'yscale': "linear",
            # 'mask_below': None, 'mask_above': None,
            'cmap': 'viridis_r',
            'norm': 'log', #'todo': None,
            # 'title': o_corr_data.sim.replace('_', '\_'),
            'title': {}, # "text": o_corr_data.sim.replace('_', '\_'), 'fontsize': 14
            'fancyticks': True,
            'minorticks': True,
            # 'sharex': False,  # removes angular citkscitks
            # 'sharey': False,
            'fontsize': 14,
            'labelsize': 14
        }
        o_plot.set_plot_dics.append(corr_dic3)

    o_plot.main()
plot_corr_ls220_LK_MO()
