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

class MODELS:

    def __init__(self):
        self.simulations = pd.read_csv(Paths.output + Files.models)
        self.simulations = self.simulations.set_index("name")
        # self.simulations["res"] = self.simulations["name"]
        self.add_res()
        self.add_q()
        self.add_viscosity()
        # self.add_dyn_phase_status()

    def add_res(self):
        """
        because of the way how pandas load dictionary, one has to loop over
        all the first row entries to select the needed resolution
        then if the resolution SR, LR or HR is found in the simulation
        name it is added, otherwise -- it complains and adds SR.
        :return: Nothing
        """
        Printcolor.blue("...adding resolutions from sim. names")
        resolutions = []
        for sim_name in [sim[0] for sim in self.simulations.iterrows()]:
            appended = False
            for res in Lists.res:
                # print(sim_name)
                if str(sim_name).__contains__(str(res)):
                    resolutions.append(res)
                    appended = True
            if not appended:
                Printcolor.yellow("Warning: No 'res' found in {} name. Using 'SR' instead".format(sim_name))
                resolutions.append('SR')
                # raise NameError("for sim:{} resolution not found".format(sim_name))
        self.simulations["res"] = resolutions

    def add_viscosity(self):
        """
        because of the way how pandas load dictionary, one has to loop over
        all the first row entries to select the needed resolution
        then if the resolution 'LK' is found in the simulation
        name it is added, otherwise - it complains and adds '-' (no viscosity).
        :return: Nothing
        """
        Printcolor.blue("...adding viscosity from sim. names")
        viscosities = []
        for sim_name in [sim[0] for sim in self.simulations.iterrows()]:
            appended = False
            for vis in Lists.visc:
                # print(sim_name)
                if str(sim_name).__contains__(str(vis)):
                    viscosities.append(vis)
                    appended = True
            if not appended:
                print("Note: No 'visc' found in {} name. Using '-' instead".format(sim_name))
                viscosities.append('-')
                # raise NameError("for sim:{} resolution not found".format(sim_name))
        self.simulations["visc"] = viscosities

    def add_q(self):
        Printcolor.blue("...adding q = M1/M2")
        self.simulations["q"] = self.simulations["M1"] / self.simulations["M2"]

    # def add_dyn_phase_status(self):
    #     Printcolor.blue("...adding dynamical phase info from static list")
    #     passed_not_passed = []
    #
    #     for name_sim, sim in self.simulations.iterrows():
    #         if name_sim in Lists.dyn_not_pas:
    #             passed_not_passed.append("passed")
    #         else:
    #             passed_not_passed.append("not passed")
    #     self.simulations["dyn_phase"] = passed_not_passed

    def get_all(self):
        return self.simulations

    def get_selected_models(self, lim_dic):
        """
        To use, provide a dictioany like :
        {'EOS':['DD2'], 'res':['SR'], 'q':[1.], 'visc':['-'], 'maximum':'tend'}
        Where the desired preperties are in the lists [], as they can be
        selected in a several.
        the last entry - 'maximum' will result in returning 1
        simulation with the maximum (minimum) of this value

        :param lim_dic:
        :return:
        """

        cropped_sims = copy.deepcopy(self.simulations)
        for v_n in lim_dic.keys():

            if not v_n in cropped_sims.keys() and v_n != 'names':
                raise NameError("key: {} not in table.keys()\n{}"
                                .format(v_n, cropped_sims.keys()))
            elif v_n == 'names':
                # .loc[[sim1, sim2 ...]] returns a dataframe with these simulations
                cropped_sims = cropped_sims.loc[lim_dic[v_n]]

            # if v_n == 'name':
            #     for value in lim_dic[v_n]:
            #         cropped_sims = cropped_sims[cropped_sims[value] == value]

            elif v_n == 'EOS':
                for value in lim_dic[v_n]:
                    cropped_sims = cropped_sims[cropped_sims[v_n] == value]
            elif v_n == 'res':
                for value in lim_dic[v_n]:
                    cropped_sims = cropped_sims[cropped_sims[v_n] == value]
            elif v_n == 'q':
                for value in lim_dic[v_n]:
                    cropped_sims = cropped_sims[cropped_sims[v_n] == value]
            elif v_n == 'visc':
                for value in lim_dic[v_n]:
                    cropped_sims = cropped_sims[cropped_sims[v_n] == value]
            #
            # else:
            #     for value in lim_dic[v_n]:
            #         cropped_sims = cropped_sims[cropped_sims[v_n] == value]

            # if v_n == 'name':
            #     for value in lim_dic[v_n]:
            #         cropped_sims = cropped_sims[cropped_sims[v_n] == value]
            # else:
            #     raise NameError("limit dic entry: {} is not reognized".format(v_n))

        if 'maximum' in lim_dic.keys():
            return cropped_sims[cropped_sims[lim_dic['maximum']] ==
                                cropped_sims[lim_dic['maximum']].max()]

        if 'minimum' in lim_dic.keys():
            return cropped_sims[cropped_sims[lim_dic['minimum']] ==
                                cropped_sims[lim_dic['minimum']].min()]

        return cropped_sims

    def save_new_table(self, new_simulations, fname="../output/summary.csv"):

        header = []
        table = []

        for sim_name, sim in new_simulations.iterrows():
            for v_n in sim.keys():
                header.append(v_n)

        table = []
        for sim_name, sim in new_simulations.iterrows():
            run = {}
            for v_n in sim.keys():
                value = sim[v_n]
                if value == None:
                    run[str(v_n)] = ''
                else:
                    run[str(v_n)] = value
            table.append(run)

                # print(sim_name)
                # print(sim[v_n])
                # exit(1)

        with open(fname, "w") as csvfile:
            writer = csv.DictWriter(csvfile, header)
            writer.writeheader()
            for run in table:
                writer.writerow(run)


    def get_value(self, sim, v_n):
        return self.simulations.get_value(sim, v_n)


class TEX_TABLE(MODELS):

    def __init__(self):
        MODELS.__init__(self)

        self.translate_dic = {
            "Mdisk3D": "$M_{\\text{disk}}$",

            "Mej" : "$M_{ej}$",
            "Yeej": "$\\langle Y_e \\rangle$",
            "vej": "$\\upsilon_{\\text{ej}}$",

            "Mej_bern" : "$M_{ej}^b$",
            "Yeej_bern": "$\\langle Y_e \\rangle^b$",
            "vej_bern": "$\\upsilon_{\\text{ej}}^b$",
        }

        self.units_dic = {
            "Mdisk3D": "$M_{\odot}$",
            "EOS": " ",

            "Mej" : "$[M_{\odot}]$",
            "Yeej": "  ",
            "vej": "$[c]$",

            "Mej_bern": "$[M_{\odot}]$",
            "Yeej_bern": "  ",
            "vej_bern": "$[c]$",

            "M1": "$[M_{\odot}]$",
            "M2": "$[M_{\odot}]$",
        }

    def print_latex_table(self):


        sims = ['DD2_M13641364_M0_SR', 'LS220_M13641364_M0_SR', 'SLy4_M13641364_M0_SR']
        v_ns = ["EOS", "M1", "M2", 'Mdisk3D', 'Mej', 'Yeej', 'vej', 'Mej_bern', 'Yeej_bern', 'vej_bern']
        precs= ["str", "1.2","1.2", ".4", ".4", ".4", ".4", ".4", ".4", ".4"]

        size = '{'
        head = ''
        for i, v_n in enumerate(v_ns):
            if v_n in self.translate_dic.keys(): v_n = self.translate_dic[v_n]
            size = size + 'c'
            head = head + '{}'.format(v_n)
            if v_n != v_ns[-1]: size = size + ' '
            if i != len(v_ns) -1 : head = head + ' & '
        size = size + '}'

        unit_bar = ''
        for v_n in v_ns:
            if v_n in self.units_dic.keys():
                unit = self.units_dic[v_n]
            else:
                unit = v_n
            unit_bar = unit_bar + '{}'.format(unit)
            if v_n != v_ns[-1]: unit_bar = unit_bar + ' & '


        head = head + ' \\\\'  # = \\
        unit_bar = unit_bar + ' \\\\ '




        print('\\begin{table*}[t]')
        print('\\begin{center}')
        print('\\begin{tabular}' + '{}'.format(size))
        print('\\hline')
        print(head)
        print(unit_bar)
        print('\\hline\\hline')

        for i, sim in enumerate(sims):

            # 1 & 6 & 87837 & 787 \\
            row = ''
            for v_n, prec in zip(v_ns, precs):

                if prec != "str":
                    val = "%{}f".format(prec) % self.get_value(sim, v_n)
                else:
                    val = self.get_value(sim, v_n)
                row = row + val
                if v_n != v_ns[-1]: row = row + ' & '
            row = row + ' \\\\'  # = \\
            print(row)

        print('\\hline')
        print('\\end{tabular}')
        print('\\end{center}')
        print('\\caption{I am your table! }')
        print('\\label{tbl:1}')
        print('\\end{table*}')


class PLOT_TASK:

    def __init__(self, o_table_class):

        self.table = o_table_class


    def plot_limits_labels_stuff(self, ax, dic):

        if dic["xmin"] != None and dic["xmax"] != None:
            ax.set_xlim(dic["xmin"], dic["xmax"])
        if dic["ymin"] != None and dic["ymax"] != None:
            ax.set_ylim(dic["ymin"], dic["ymax"])
        if dic["xscale"] == 'log':
            ax.set_xscale("log")
        if dic["yscale"] == 'log':
            ax.set_yscale("log")
        ax.set_ylabel(dic["v_n_y"].replace('_', '\_'))
        ax.set_xlabel(dic["v_n_x"].replace('_', '\_'))

        if dic["legend"]:
            ax.legend(fancybox=False, loc='best', shadow=False, fontsize=8)

    def plot_points(self, ax, dic):

        # dic = {
        #     'name': 'points', 'position': (1, 1), 'title': None, 'cbar': None,
        #     'it': None, 'v_n_x': 'Lambda', 'v_n_y': 'Mej', 'v_n': None,
        #     'xmin': None, 'xmax': None, 'ymin': None, 'ymax': None,
        #     'xscale': None, 'yscale': "log",
        #     'color': 'EOS', 'marker': 'x', 'annotate': 'name'
        # }

        for isim, sim in self.table.iterrows():

            if not dic["v_n_x"] in sim.keys():
                raise NameError("v_n_x:{} not in sim.keys() \n{}".format(dic["v_n_x"], sim.keys()))

            if not dic["v_n_y"] in sim.keys():
                raise NameError("v_n_x:{} not in sim.keys() \n{}".format(dic["v_n_y"], sim.keys()))




            if dic["color"] == "EOS":
                color = eos_color(sim["EOS"])
            else:
                color = dic["color"]
            ax.plot(sim[dic["v_n_x"]], sim[dic["v_n_y"]], dic["marker"], color=color)

            if dic["annotate"] != None and dic["annotate"] != '':
                if dic["annotate"] == 'name':
                    annotation = isim
                else:
                    annotation = sim[dic["annotation"]]
                ax.annotate(r'${}$'.format(annotation.replace('_', '\_')),
                            xy=(sim[dic["v_n_x"]], sim[dic["v_n_y"]]),
                            textcoords='data')

            if dic["arrows_up"] != None and dic["arrows_up"] != '':
                if dic["arrows_up"] == "bern_phase":
                    if sim["bern_phase"] == "not passed":
                        print("sim: {} arrow".format(isim))
                        ax.arrow(sim[dic["v_n_x"]],  sim[dic["v_n_y"]], 0, 0.05 * sim[dic["v_n_y"]],
                                 )


        self.plot_limits_labels_stuff(ax, dic)

    def plot_line(self, ax, dic):

        x_coords = []
        y_coords = []

        for isim, sim in  self.table.iterrows():

            if not dic["v_n_x"] in sim.keys():
                raise NameError("v_n_x:{} not in sim.keys() \n{}".format(dic["v_n_x"], sim.keys()))

            if not dic["v_n_y"] in sim.keys():
                raise NameError("v_n_x:{} not in sim.keys() \n{}".format(dic["v_n_y"], sim.keys()))

            x_coords.append(sim[dic["v_n_x"]])
            y_coords.append(sim[dic["v_n_y"]])

        x_coords = np.array(x_coords)
        y_coords = np.array(y_coords)

        try:
            x_coords, y_coords = x_y_z_sort(x_coords, y_coords)
        except ValueError:
            Printcolor.yellow("Warning: sorting x, y failed. (x:{} y:{})"
                              .format(dic["v_n_x"], dic["v_n_y"]))

        ax.plot(x_coords, y_coords, dic["ls"], color=dic["color"])

        self.plot_limits_labels_stuff(ax, dic)

    def secial_plot(self, ax, dic):

        m_bern_ej = []
        m_dyn_ej = []
        m_disk = []
        lam = []

        for isim, sim in self.table.iterrows():

            m_bern_ej.append(sim["Mej_bern"])
            m_dyn_ej.append(sim["Mej"])
            m_disk.append(sim["Mdisk3D"])
            lam.append(sim["Lambda"])

        m_bern_ej = np.array(m_bern_ej)
        m_dyn_ej = np.array(m_dyn_ej)
        m_disk = np.array(m_disk)
        lam = np.array(lam)



        # try:
        #     lam, m_bern_ej = x_y_z_sort(lam, m_bern_ej)
        #     lam, m_dyn_ej = x_y_z_sort(lam, m_dyn_ej)
        #     lam, m_disk = x_y_z_sort(lam, m_disk)
        # except ValueError:
        #     Printcolor.yellow("Warning: sorting failed")


        ax.plot(lam, m_dyn_ej, '-', color="black")
        ax.plot(lam, m_bern_ej, '-', color="red")
        ax.plot(lam, 0.4 * m_disk, '--', color="gray")

        ax.plot(lam, m_dyn_ej, 'x', color="black", label="Dyn")
        ax.plot(lam, m_bern_ej, 'o', color="red", label="Bern")
        ax.plot(lam, 0.4 * m_disk, '.', color="gray", label="Disk")


        ylim = np.zeros(len(m_disk))
        ylim.fill(ax.get_ylim()[-1])
        ax.fill_between(lam, 0.4 * m_disk, ylim, color='gray')

        self.plot_limits_labels_stuff(ax, dic)

        for isim, sim in self.table.iterrows():
            print("\n--- --- KILONOVA PROPERTIES --- --- ")
            print("sim: {}".format(isim))
            print("\tglob_vars:")
            print("\t\tm_disk:       {} ".format(sim["Mdisk3D"]))

            print("\tDynamical: ")
            print("\t\tm_ej:              {} # average mass".format(sim["Mej"]))
            print("\t\tcentral_vel:       {} # average vel_inf".format(sim["vej"]))

            print("\tPseudodynamical: ")
            print("\t\tm_ej:              {} # average mass".format(sim["Mej_bern"]))
            print("\t\tcentral_vel:       {} # average vel_inf".format(sim["vej_bern"]))
            print("--- --- --- -- DONE -- --- --- --- \n")

class PLOT_MANY_TASKS(PLOT_TASK):

    def __init__(self, sims):

        PLOT_TASK.__init__(self, sims)

        self.data = sims

        self.gen_set = {
            "figdir": Paths.plots + 'overall/',
            "figname": "inv_ang_mom_flux.png",
            # "figsize": (13.5, 3.5), # <->, |
            "figsize": (3.8, 3.5),  # <->, |
            "type": "cartesian",
            "subplots_adjust_h": 0.2,
            "subplots_adjust_w": 0.3,
            "fancy_ticks": True,
            "minorticks_on": True
        }

        self.set_plot_dics = []


    def set_ncols_nrows(self):

        tmp_rows = []
        tmp_cols = []

        for dic in self.set_plot_dics:
            tmp_cols.append(dic['position'][1])
            tmp_rows.append(dic['position'][0])

        max_row = max(tmp_rows)
        max_col = max(tmp_cols)

        for row in range(1, max_row):
            if not row in tmp_rows:
                raise NameError("Please set vertical plot position in a subsequent order: 1,2,3... not 1,3...")

        for col in range(1, max_col):
            if not col in tmp_cols:
                raise NameError("Please set horizontal plot position in a subsequent order: 1,2,3... not 1,3...")

        print("Set {} rows {} columns (total {}) of plots".format(max_row, max_col, len(self.set_plot_dics)))

        return int(max_row), int(max_col)

    def set_plot_dics_matrix(self):

        plot_dic_matrix = [[0
                             for x in range(self.n_rows)]
                             for y in range(self.n_cols)]

        # get a matrix of dictionaries describing plots (for ease of representation)
        for dic in self.set_plot_dics:
            col, row = int(dic['position'][1]-1), int(dic['position'][0]-1) # -1 as position starts with 1
            # print(col, row)
            for n_row in range(self.n_rows):
                for n_col in range(self.n_cols):
                    if int(col) == int(n_col) and int(row) == int(n_row):
                        plot_dic_matrix[n_col][n_row] = dic
                        # print('adding {} {}'.format(col, row))

            if isinstance(plot_dic_matrix[col][row], int):
                raise ValueError("Dictionary to found for n_row {} n_col {} in "
                                 "creating matrix of dictionaries".format(col, row))

        return plot_dic_matrix

    def set_plot_matrix(self):

        fig = plt.figure(figsize=self.gen_set['figsize'])  # (<->; v)


        if self.gen_set['type'] == 'cartesian':
            # initializing the matrix with dummy axis objects
            sbplot_matrix = [[fig.add_subplot(self.n_rows, self.n_cols, 1)
                                  for x in range(self.n_rows)]
                                  for y in range(self.n_cols)]

            i = 1
            for n_row in range(self.n_rows):
                for n_col in range(self.n_cols):

                    if n_col == 0 and n_row == 0:
                        sbplot_matrix[n_col][n_row] = fig.add_subplot(self.n_rows, self.n_cols, i)
                    elif n_col == 0 and n_row > 0:
                        if self.gen_set['sharex']:
                            sbplot_matrix[n_col][n_row] = fig.add_subplot(self.n_rows, self.n_cols, i,
                                                                          harex=self.sbplot_matrix[n_col][0])
                        else:
                           sbplot_matrix[n_col][n_row] = fig.add_subplot(self.n_rows, self.n_cols, i)
                    elif n_col > 0 and n_row == 0:
                        if self.gen_set['sharey']:
                            sbplot_matrix[n_col][n_row] = fig.add_subplot(self.n_rows, self.n_cols, i,
                                                                          sharey=self.sbplot_matrix[0][n_row])
                        else:
                            sbplot_matrix[n_col][n_row] = fig.add_subplot(self.n_rows, self.n_cols, i)
                    else:
                        if self.gen_set['sharex'] and not self.gen_set['sharey']:
                            sbplot_matrix[n_col][n_row] = fig.add_subplot(self.n_rows, self.n_cols, i,
                                                                          sharex=self.sbplot_matrix[n_col][0])
                        elif not self.gen_set['sharex'] and self.gen_set['sharey']:
                            sbplot_matrix[n_col][n_row] = fig.add_subplot(self.n_rows, self.n_cols, i,
                                                                          sharey=self.sbplot_matrix[0][n_row])
                        else:
                           sbplot_matrix[n_col][n_row] = fig.add_subplot(self.n_rows, self.n_cols, i,
                                                                         sharex=self.sbplot_matrix[n_col][0],
                                                                         sharey=self.sbplot_matrix[0][n_row])

                        # sbplot_matrix[n_col][n_row].axes.get_yaxis().set_visible(False)
                    # sbplot_matrix[n_col][n_row] = fig.add_subplot(n_rows, n_cols, i)
                    i += 1

        elif self.gen_set['type'] == 'polar':
            # initializing the matrix with dummy axis objects
            sbplot_matrix = [[fig.add_subplot(self.n_rows, self.n_cols, 1, projection='polar')
                                  for x in range(self.n_rows)]
                                  for y in range(self.n_cols)]

            i = 1
            for n_row in range(self.n_rows):
                for n_col in range(self.n_cols):

                    if n_col == 0 and n_row == 0:
                        sbplot_matrix[n_col][n_row] = fig.add_subplot(self.n_rows, self.n_cols, i, projection='polar')
                    elif n_col == 0 and n_row > 0:
                        sbplot_matrix[n_col][n_row] = fig.add_subplot(self.n_rows, self.n_cols, i, projection='polar')
                                                                      # sharex=self.sbplot_matrix[n_col][0])
                    elif n_col > 0 and n_row == 0:
                        sbplot_matrix[n_col][n_row] = fig.add_subplot(self.n_rows, self.n_cols, i, projection='polar')
                                                                      # sharey=self.sbplot_matrix[0][n_row])
                    else:
                        sbplot_matrix[n_col][n_row] = fig.add_subplot(self.n_rows, self.n_cols, i, projection='polar')
                                                                      # sharex=self.sbplot_matrix[n_col][0],
                                                                      # sharey=self.sbplot_matrix[0][n_row])

                        # sbplot_matrix[n_col][n_row].axes.get_yaxis().set_visible(False)
                    # sbplot_matrix[n_col][n_row] = fig.add_subplot(n_rows, n_cols, i)
                    i += 1
        else:
            raise NameError("type of the plot is not recognized. Use 'polar' or 'cartesian' ")

        return fig, sbplot_matrix

    def plot_one_task(self, ax, dic):


        if dic["name"] == "points":
            self.plot_points(ax, dic)

        elif dic["name"] == "special":
            self.secial_plot(ax, dic)

        im = 0
        return 0


        # if dic["name"] == "corr" or dic["name"] == "int":
        #     im = self.plot_x_y_z2d_colormesh(ax, dic)
        # elif dic["name"] == "count":
        #     self.plot_countours(ax, dic)
        #     im = 0
        # elif dic["name"] == 'densmode':
        #     self.plot_density_mode(ax, dic)
        #     im = 0
        # else:
        #     raise NameError("name:{} is not recognized"
        #                     .format(dic["name"]))
        #
        # # self.time
        # return im

    def set_plot_title(self, ax, plot_dic):
        if plot_dic["title"] != '' and plot_dic["title"] != None:

            if plot_dic["title"] == 'it':
                title = plot_dic["it"]
            elif plot_dic["title"] == 'time [s]' or \
                plot_dic["title"] == 'time':
                title = "%.3f" % self.data.get_time(plot_dic["it"]) + " [s]"
            elif plot_dic["title"] == 'time [ms]':
                title = "%.1f" % (self.data.get_time(plot_dic["it"]) * 1000) + " [ms]"
            else:
                title = plot_dic["title"]
            ax.title.set_text(r'{}'.format(title))

    def plot_images(self):

        # initializing the matrix of images for colorbars (if needed)
        image_matrix = [[0
                        for x in range(self.n_rows)]
                        for y in range(self.n_cols)]


        for n_row in range(self.n_rows):
            for n_col in range(self.n_cols):
                for dic in self.set_plot_dics:
                    if n_col + 1 == int(dic['position'][1]) and n_row + 1 == int(dic['position'][0]):
                        print("Plotting n_row:{} n_col:{}".format(n_row, n_col))
                        ax = self.sbplot_matrix[n_col][n_row]
                        # dic = self.plot_dic_matrix[n_col][n_row]
                        if isinstance(dic, int):
                            Printcolor.yellow("Dictionary for row:{} col:{} not set".format(n_row, n_col))
                            self.fig.delaxes(ax)  # delets the axis for empty plot
                        else:
                            dic = dict(dic)
                            im = self.plot_one_task(ax, dic)
                            self.set_plot_title(ax, dic)
                            if not isinstance(im, int):
                                image_matrix[n_col][n_row] = im

                            self.add_fancy_to_ax(ax)

        return image_matrix

    def add_fancy_to_ax(self, ax):

        if self.gen_set["fancy_ticks"]:
            ax.tick_params(axis='both', which='both', labelleft=True,
                           labelright=False, tick1On=True, tick2On=True,
                           labelsize=12, direction='in')
        if self.gen_set["minorticks_on"]:
            ax.minorticks_on()

    def plot_one_cbar(self, im, dic, n_row, n_col):

        if dic["cbar"] != None and dic["cbar"] != '':

            location = dic["cbar"].split(' ')[0]
            shift_h = float(dic["cbar"].split(' ')[1])
            shift_w = float(dic["cbar"].split(' ')[2])
            cbar_width = 0.02


            if location == 'right':
                ax_to_use = self.sbplot_matrix[-1][n_row]
                pos1 = ax_to_use.get_position()
                pos2 = [pos1.x0 + pos1.width + shift_h,
                        pos1.y0,
                        cbar_width,
                        pos1.height]
            elif location == 'left':
                ax_to_use = self.sbplot_matrix[-1][n_row]
                pos1 = ax_to_use.get_position()
                pos2 = [pos1.x0 - pos1.width - shift_h,
                        pos1.y0,
                        cbar_width,
                        pos1.height]
            elif location == 'bottom':
                ax_to_use = self.sbplot_matrix[n_col][-1]
                pos1 = ax_to_use.get_position()
                pos2 = [pos1.x0,
                        pos1.y0 - pos1.height + shift_w,
                        cbar_width,
                        pos1.height]
            else:
                raise NameError("cbar location {} not recognized. Use 'right' or 'bottom' "
                                .format(location))

            cax1 = self.fig.add_axes(pos2)
            if location == 'right':
                cbar = plt.colorbar(im, cax=cax1, extend='both')#, format='%.1e')
            elif location == 'left':
                cbar = plt.colorbar(im, cax=cax1, extend='both')#, format='%.1e')
                cax1.yaxis.set_ticks_position('left')
                cax1.yaxis.set_label_position('left')
            else:
                raise NameError("cbar location {} not recognized. Use 'right' or 'bottom' "
                                .format(location))
            cbar.ax.set_title(r"{}".format(str(dic["v_n"]).replace('_', '\_')))

    def plot_colobars(self):

        for n_row in range(self.n_rows):
            for n_col in range(self.n_cols):
                for dic in self.set_plot_dics:
                    if n_col + 1 == int(dic['position'][1]) and n_row + 1 == int(dic['position'][0]):
                        print("Colobar for n_row:{} n_col:{}".format(n_row, n_col))
                        # ax  = self.sbplot_matrix[n_col][n_row]
                        # dic = self.plot_dic_matrix[n_col][n_row]
                        im  = self.image_matrix[n_col][n_row]
                        if isinstance(dic, int):
                            Printcolor.yellow("Dictionary for row:{} col:{} not set".format(n_row, n_col))
                        else:
                            self.plot_one_cbar(im, dic, n_row, n_col)


        # for n_row in range(self.n_rows):
        #     for n_col in range(self.n_cols):
        #         print("Colobar for n_row:{} n_col:{}".format(n_row, n_col))
        #         # ax  = self.sbplot_matrix[n_col][n_row]
        #         dic = self.plot_dic_matrix[n_col][n_row]
        #         im  = self.image_matrix[n_col][n_row]
        #         if isinstance(dic, int):
        #             Printcolor.yellow("Dictionary for row:{} col:{} not set".format(n_row, n_col))
        #         else:
        #             self.plot_one_cbar(im, dic, n_row, n_col)

    def save_plot(self):

        plt.subplots_adjust(hspace=self.gen_set["subplots_adjust_h"])
        plt.subplots_adjust(wspace=self.gen_set["subplots_adjust_w"])
        # plt.tight_layout()
        plt.savefig('{}{}'.format(self.gen_set["figdir"], self.gen_set["figname"]),
                    bbox_inches='tight', dpi=128)
        plt.close()

    def main(self):

        # initializing the n_cols, n_rows
        self.n_rows, self.n_cols = self.set_ncols_nrows()
        # initializing the matrix of dictionaries of the
        self.plot_dic_matrix = self.set_plot_dics_matrix()
        # initializing the axis matrix (for all subplots) and image matrix fo colorbars
        self.fig, self.sbplot_matrix = self.set_plot_matrix()
        # plotting
        self.image_matrix = self.plot_images()
        # adding colobars
        self.plot_colobars()


        # saving the result
        self.save_plot()


""" --- old class to be putged --- --- """
class PLOT_MODELS:

    def __init__(self):
        pass

    @staticmethod
    def plot_vn1_vn2(sims, v_n1="Lambda", v_n2="Yeej"):

        # sims = sims[sims["Mej"] > 0.0015]
        # sims = sims[sims["tend"] > 0.0020]
        print(len(sims))


        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.minorticks_on()
        for isim, sim in sims.iterrows():
            ax.plot(sim[v_n1], sim[v_n2], 'x',
                    color=eos_color(sim["EOS"]))

        unique_eos = get_uniqe_eoss(sims["EOS"])

        for eos in unique_eos:
            for isim, sim in sims.iterrows():
                if sim["EOS"] == eos:
                    ax.plot(sim[v_n1], sim[v_n2], 'x',
                            color=eos_color(sim["EOS"]),
                            label='{}'.format(sim["EOS"]))
                    break

        ax.tick_params(axis='both', which='both', labelleft=True,
                       labelright=False, tick1On=True, tick2On=True,
                       labelsize=12, direction='in')
        ax.legend(loc='lower center', bbox_to_anchor=(0.5, 0.80),
                  fancybox=False, shadow=False, ncol=2, fontsize=10)
        ax.set_xlabel(v_n1)
        # ax.set_xlim(390, 450)
        # ax.set_yscale("log")
        ax.set_ylabel(v_n2)
        plt.savefig('fig1', bbox_inches='tight', dpi=128)
        plt.legend()
        plt.close()

    @staticmethod
    def plot_vn1_vn2_flux_evolution(sims, figname):

        # sims = sims[sims["Mej"] > 0.0015]
        # sims = sims[sims["tend"] > 0.0020]
        print(len(sims))
        v_n1 = 'tend'
        v_n2 = 'Mej'

        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.minorticks_on()
        for isim, sim in sims.iterrows():
            ax.plot(sim[v_n1], sim[v_n2], 'x',
                    color=eos_color(sim["EOS"]))

        for isim, sim in sims.iterrows():

            ej = SIM_EJ_HIST(isim)

            ax.plot(ej.time_total_flux, ej.mass_total_flux, '-', color=eos_color(sim["EOS"]))
            ax.plot(sim[v_n1], sim[v_n2], 'x', color=eos_color(sim["EOS"]))
            ax.annotate(r'${}$'.format(str(isim).replace('_', '\_')), xy=(sim[v_n1], sim[v_n2]),
                        textcoords='data')


        unique_eos = get_uniqe_eoss(sims["EOS"])

        for eos in unique_eos:
            for isim, sim in sims.iterrows():
                if sim["EOS"] == eos:
                    ax.plot(sim[v_n1], sim[v_n2], 'x',
                            color=eos_color(sim["EOS"]),
                            label='{}'.format(sim["EOS"]))
                    break

        ax.tick_params(axis='both', which='both', labelleft=True,
                       labelright=False, tick1On=True, tick2On=True,
                       labelsize=12, direction='in')
        ax.legend(loc='lower center', bbox_to_anchor=(0.5, 0.80),
                  fancybox=False, shadow=False, ncol=2, fontsize=10)
        ax.set_xlabel(v_n1)
        # ax.set_xlim(390, 450)
        # ax.set_yscale("log")
        ax.set_ylabel(v_n2)
        plt.savefig(Paths.plots + figname, bbox_inches='tight', dpi=128)
        plt.legend()
        plt.close()




    @staticmethod
    def plot_one_ejecta_prof(ax, sims, v_n1='tend', v_n2='Mej', ext="_0_b"):

        # for isim, sim in sims.iterrows():
        #     ax.plot(sim[v_n1], sim[v_n2], 'x',
        #             color=eos_color(sim["EOS"]))
        # if ext in ["_0_b", "_1_b", "_0_b_w", "_1_b_w", "_0_b_d", "_1_b_d"]:
        #     v_n2 += "_bern"
        idx = 1
        for isim, sim in sims.iterrows():
            from scipy import interpolate
            data_status = LOAD_ITTIME(isim)
            print(isim)
            color = get_color_for_q(sim["M1"] / sim["M2"])

            # geodesic
            ej = SIM_EJ_HIST(isim, extension="_0", norm=True)

            ax.plot(ej.time_total_flux, ej.mass_total_flux, '-', color=color, lw=0.9,
                    label=r'${}:{}$'.format(idx, str(isim).replace('_', '\_')))


            isd3data, d3it, d3times = data_status.get_ittime(output="overall", d1d2d3prof="d3")
            if isd3data:
                ax.plot(d3times, interpolate.interp1d(ej.time_total_flux, ej.mass_total_flux, bounds_error=False)(d3times),
                        ':', color=color, lw=1.8)

            are_profs, profit, proft = data_status.get_ittime(output="profiles", d1d2d3prof="prof")
            if are_profs:
                ax.plot(proft, interpolate.interp1d(ej.time_total_flux, ej.mass_total_flux, bounds_error=False)(proft),
                        marker='|', color=color, markersize=3.3)


            ax.annotate(r'${}$'.format(idx), xy=(ej.time_total_flux[-1], ej.mass_total_flux[-1]), textcoords='data')
            # bernoulli
            ej2 = SIM_EJ_HIST(isim, extension=ext, norm=True)
            try:
                ax.plot(ej2.time_total_flux, ej2.mass_total_flux, '--', color=color, lw=0.5)
                        # label=r'${}:{}$'.format(idx, str(isim).replace('_', '\_')))
                # print(ej2.time_total_flux[-1], ej2.mass_total_flux[-1])
                ax.annotate(r'${}$'.format(idx), xy=(ej2.time_total_flux[-1], ej2.mass_total_flux[-1]),
                            textcoords='data', fontsize=8, color="gray")
            except:
                Printcolor.yellow("bernoulli ejecta failed for {}".format(isim))

            if not np.isinf(sim["tcoll"]):
                ax.scatter(ej.time_total_flux[-1], ej.mass_total_flux[-1], s=3.0,color=color)#marker='o', color=color, markersize=0.5)

            idx += 1

        ax.legend(fancybox=False, loc='best',shadow=False, fontsize=8)

        # unique_eos = get_uniqe_eoss(sims["EOS"])

        # for eos in unique_eos:
        #     for isim, sim in sims.iterrows():
        #         if sim["EOS"] == eos:
        #             ax.text(0.9, 0.8, '{}'.format(sim["EOS"]),
        #                     verticalalignment='bottom', horizontalalignment='right',
        #                     transform=ax.transAxes, color='black', fontsize=12)
                    # ax.text(0.01, 0.5, '{}'.format(sim["EOS"]), fontsize=14)
                    # # ax.plot(sim[v_n1], sim[v_n2], 'x',
                    # #         color=eos_color(sim["EOS"]),
                    # #         label='{}'.format(sim["EOS"]))
                    # break

        ax.tick_params(axis='both', which='both', labelleft=True,
                       labelright=False, tick1On=True, tick2On=True,
                       labelsize=12, direction='in')
        ax.minorticks_on()

        # ax.set_xlabel(v_n1)
        # ax.set_xlim(390, 450)
        # ax.set_yscale("log")
        # ax.set_ylabel(v_n2)
        # plt.savefig(Paths.plots + figname, bbox_inches='tight', dpi=128)
        # plt.legend()
        # plt.close()


    @staticmethod
    def plot_vn1_vn2_flux_evolution2(extension="_0"):

        list_EOSs = ['DD2', 'LS220', 'SLy4', 'SFHo']
        list_of_res=['LR', 'SR', 'HR']

        sim_cl = MODELS()

        fig = plt.figure(figsize=(4.5*2.5, 4.5*3.))  # (<->; v)

        n_rows = len(list_EOSs)
        n_cols = len(list_of_res)

        # initializing the matrix with dummy axis objects
        sbplot_matrix = [[fig.add_subplot(n_rows, n_cols, 1)
                          for x in range(n_rows)] for y in range(n_cols)]

        i = 1
        for n_row in range(n_rows):
            for n_col in range(n_cols):

                # if n_col == 0:
                #     sbplot_matrix[n_col][n_row] = fig.add_subplot(n_rows, n_cols, i)
                # else:
                #     sbplot_matrix[n_col][n_row] = fig.add_subplot(n_rows, n_cols, i, sharey=sbplot_matrix[0][n_row])
                #     sbplot_matrix[n_col][n_row].axes.get_yaxis().set_visible(False)

                if n_col == 0 and n_row == 0:
                    sbplot_matrix[n_col][n_row] = fig.add_subplot(n_rows, n_cols, i)

                elif n_col == 0 and n_row > 0:
                    sbplot_matrix[n_col][n_row] = fig.add_subplot(n_rows, n_cols, i,
                                                                  sharex=sbplot_matrix[n_col][0])

                elif n_cols > 0 and n_row == 0:
                    sbplot_matrix[n_col][n_row] = fig.add_subplot(n_rows, n_cols, i,
                                                                  sharey=sbplot_matrix[0][n_row])

                else:
                    sbplot_matrix[n_col][n_row] = fig.add_subplot(n_rows, n_cols, i,
                                                                  sharex=sbplot_matrix[n_col][0],
                                                                  sharey=sbplot_matrix[0][n_row])

                    # sbplot_matrix[n_col][n_row].axes.get_yaxis().set_visible(False)
                # sbplot_matrix[n_col][n_row] = fig.add_subplot(n_rows, n_cols, i)
                i+=1

        for n_row in range(n_rows):
            for n_col in range(n_cols):
                lim_dic_in = {'EOS': [list_EOSs[n_row]],
                              'res': [list_of_res[n_col]],
                              'q': [], 'visc': []}
                sims = sim_cl.get_selected_models(lim_dic_in)
                ax = sbplot_matrix[n_col][n_row]
                PLOT_MODELS.plot_one_ejecta_prof(ax, sims, v_n1='tend', v_n2='Mej', ext=extension)
                # ax.set_xlim(left=1e-2)
                # ax.set_xscale("log")
                # ax.set_xlim(left=1e-2,right=5e-1)

                if n_col == 0:
                    ax.text(0.9, 0.2, '{}'.format(list_EOSs[n_row]),
                            verticalalignment='bottom', horizontalalignment='right',
                            transform=ax.transAxes, color='black', fontsize=12)

                if n_col == 0:
                    ax.set_ylabel(r'Mass Flux $[M_{\odot}]$')
                if n_row == n_rows-1:
                    ax.set_xlabel(r'Time [s]')

                if n_col > 0 and n_row < n_rows-1:
                    sbplot_matrix[n_col][n_row].tick_params(labelbottom=False)
                if n_col > 0 and n_row < n_rows:
                    sbplot_matrix[n_col][n_row].tick_params(labelleft=False)

                # ax.set_yscale('log')
                # if n_col == 0:
                #     ax.legend(loc='lower center', bbox_to_anchor=(0.5, 0.80),
                #               fancybox=False, shadow=False, ncol=2, fontsize=10)

                # if n_col == n_cols-1:
                    # text = ""
                    # ax.text(0.02, 0.5, text, fontsize=14, transform=plt.gcf().transFigure)

        plt.subplots_adjust(hspace=0.0)
        plt.subplots_adjust(wspace=0.0)

        plt.savefig(Paths.plots + 'flux_evolution{}.pdf'.format(extension), bbox_inches='tight', dpi=128)
        plt.legend()
        plt.close()


    @staticmethod
    def plot_one_ave_final(ax, sims, v_nx, v_ny):

        # for isim, sim in sims.iterrows():
        #     # if sim["resolution"]
        #     ax.plot(sim[v_nx], sim[v_ny], 'x',
        #             color=eos_color(sim["EOS"]))


        for isim, sim in sims.iterrows():
            ax.plot(sim[v_nx], sim[v_ny], 'x',
                    color=eos_color(sim["EOS"]))



        for isim, sim in sims.iterrows():
            if v_ny + "_bern" in sim.keys():
                ax.plot(sim[v_nx], sim[v_ny+"_bern"], 'o',
                        color=eos_color(sim["EOS"]))

        # if v_ny == "Mdisk" or v_ny == "Mdisk3D":
        #     ax.annotate(r'${}$'.format(idx), xy=(ej.time_total_flux[-1], ej.mass_total_flux[-1]), textcoords='data')

            # except:
            #     Printcolor.yellow("Warning: {} {} not found for bernouli".format(isim, v_ny))
        # for isim, sim in sims.iterrows():
        #
        #     ej = SIM_EJ_HIST(isim)
        #
        #     ax.plot(ej.time_total_flux, ej.mass_total_flux, '-', label=r'${}$'.format(str(isim).replace('_', '\_')))
        #     ax.plot(sim[v_nx], sim[v_ny], 'x')
            # ax.annotate(r'${}$'.format(str(isim).replace('_', '\_')), xy=(sim[v_n1], sim[v_n2]),
            #             textcoords='data')



        unique_eos = get_uniqe_eoss(sims["EOS"])

        for eos in unique_eos:
            for isim, sim in sims.iterrows():
                if sim["EOS"] == eos:
                    # ax.text(0.9, 0.8, '{}'.format(sim["EOS"]),
                    #         verticalalignment='bottom', horizontalalignment='right',
                    #         transform=ax.transAxes, color='black', fontsize=12)
                    # ax.text(0.01, 0.5, '{}'.format(sim["EOS"]), fontsize=14)
                    ax.plot(sim[v_nx], sim[v_ny], 'x',
                            color=eos_color(sim["EOS"]),
                            label='{}'.format(sim["EOS"]))
                    break

        ax.tick_params(axis='both', which='both', labelleft=True,
                       labelright=False, tick1On=True, tick2On=True,
                       labelsize=12, direction='in')
        ax.minorticks_on()
        ax.legend(fancybox=False, loc='best', shadow=False, fontsize=8)

        # ax.set_xlabel(v_n1)
        # ax.set_xlim(390, 450)
        # ax.set_yscale("log")
        # ax.set_ylabel(v_n2)
        # plt.savefig(Paths.plots + figname, bbox_inches='tight', dpi=128)
        # plt.legend()
        # plt.close()



    @staticmethod
    def plot_final_average(v_nx, v_ny_arr):

        sim_cl = MODELS()

        fig = plt.figure(figsize=(3.5, 8.6))  # (<->; v) 3.6

        n_rows = len(v_ny_arr)
        n_cols = len(v_nx)

        # initializing the matrix with dummy axis objects
        sbplot_matrix = [[fig.add_subplot(n_rows, n_cols, 1)
                          for x in range(n_rows)] for y in range(n_cols)]

        i = 1
        for n_row in range(n_rows):
            for n_col in range(n_cols):

                if n_col == 0 and n_row == 0:
                    sbplot_matrix[n_col][n_row] = fig.add_subplot(n_rows, n_cols, i)
                elif n_col == 0 and n_row > 0:
                    sbplot_matrix[n_col][n_row] = fig.add_subplot(n_rows, n_cols, i,
                                                                  sharex=sbplot_matrix[n_col][0])
                elif n_cols > 0 and n_row == 0:
                    sbplot_matrix[n_col][n_row] = fig.add_subplot(n_rows, n_cols, i,
                                                                  sharey=sbplot_matrix[0][n_row])
                else:
                    sbplot_matrix[n_col][n_row] = fig.add_subplot(n_rows, n_cols, i,
                                                                  sharex=sbplot_matrix[n_col][0],
                                                                  sharey=sbplot_matrix[0][n_row])
                i += 1



        lim_dic_in = {'dyn_phase': 'passed'}
        sims = sim_cl.get_selected_models(lim_dic_in)
        for n_row in range(n_rows):
            for n_col in range(n_cols):

                ax = sbplot_matrix[n_col][n_row]
                PLOT_MODELS.plot_one_ave_final(ax, sims, v_nx[n_col], v_ny_arr[n_row])

                # if n_col == 0:
                #     ax.text(0.9, 0.8, '{}'.format(list_EOSs[n_row]),
                #             verticalalignment='bottom', horizontalalignment='right',
                #             transform=ax.transAxes, color='black', fontsize=12)

                if v_ny_arr[n_row] == 'Yeej':
                    ax.axhline(y=0.25, linestyle='dashed', color='black', linewidth=.5)
                    ax.set_ylim(0.05, 0.4)
                if v_ny_arr[n_row] == 'Mej':
                    ax.set_yscale('log')
                if v_nx[n_col] == 'Mej':
                    ax.set_xscale('log')

                if v_nx[n_col] == 'Mej' and v_ny_arr[n_row] == 'vej':
                    vej_blue_err = np.array([0.2, 0.3])
                    mej_blue_err = np.array([1*10**-2, 2*10**-2])

                    vej_red_err = np.array([0.07, 0.14])
                    mej_red_err = np.array([4*10**-2, 6*10**-2])

                    vej_blue = 0.27
                    mej_blue = 0.016

                    vej_red = 0.1
                    mej_red = 0.05

                    ax.plot(mej_blue, vej_blue, 'o', color='blue')
                    ax.plot([mej_blue, mej_blue], [vej_blue_err[0], vej_blue_err[1]], '-', color='blue')
                    ax.plot([mej_blue_err[0], mej_blue_err[1]], [vej_blue, vej_blue], '-', color='blue')

                    ax.plot(mej_red, vej_red, 'o', color='red')
                    ax.plot([mej_red, mej_red], [vej_red_err[0], vej_red_err[1]], '-', color='red')
                    ax.plot([mej_red_err[0], mej_red_err[1]], [vej_red, vej_red], '-', color='red')

                    mej_blue_err_s = [mej_blue_err[0], mej_blue_err[1], mej_blue_err[1], mej_blue_err[0]]
                    vej_blue_err_s = [vej_blue_err[0], vej_blue_err[0], vej_blue_err[1], vej_blue_err[1]]
                    ax.add_patch(patches.Polygon(xy=list(zip(mej_blue_err_s, vej_blue_err_s)), fill=True, alpha=.4,
                                                 color='blue',label='Blue (Siegal 19)'))

                    mej_red_err_s = [mej_red_err[0], mej_red_err[1], mej_red_err[1], mej_red_err[0]]
                    vej_red_err_s = [vej_red_err[0], vej_red_err[0], vej_red_err[1], vej_red_err[1]]
                    ax.add_patch(patches.Polygon(xy=list(zip(mej_red_err_s, vej_red_err_s)), fill=True, alpha=.4,
                                                 color='red', label='Red (Siegal 19)'))

                    ax.set_xlim(4*10**-4, 10**-1)
                    ax.set_ylim(0.05, 0.32)

                    # dd2_mej_tot = dd2.mass_total_flux[-1] + dd2_w.mass_total_flux[-1]
                    # print(dd2.get_par('vel_inf'))
                    # print(dd2_w.get_par('vel_inf'))
                    # vel_ave = (dd2.mass_total_flux[-1] * dd2.get_par('vel_inf') +
                    #            dd2_w.mass_total_flux[-1] * dd2_w.get_par('vel_inf_bern')) / \
                    #            dd2_mej_tot

                    # viller
                    kappa_b = 0.5
                    mej_blue_vill = 0.02
                    vej_blue_vill = 0.27

                    kappa_purp = 3
                    mej_purp_vill = 0.047
                    vej_purp_vill = 0.15

                    kappa_red = 10
                    mej_red_vill = 0.011
                    vej_red_vill = 0.14

                    # ax.plot(mej_blue_vill, vej_blue_vill, 's', color='gray', markersize=12, label='Villar+17')
                    # ax.plot(mej_blue_vill, vej_blue_vill, 's', color='blue',markersize=12)
                    # ax.plot(mej_purp_vill, vej_purp_vill, 's', color='purple',markersize=12)
                    # ax.plot(mej_red_vill, vej_red_vill, 's', color='red',markersize=12)

                    # ax.plot(dd2_mej_tot, vel_ave, '*', color='black', label=r'DD2\_M13641364\_M0')

                    # ax.legend()
                    ax.legend(bbox_to_anchor=(1, 1), loc='lower right', ncol=2, fontsize=11)

                if n_col == 0:
                    if v_nx[n_col] == 'Mej':
                        ax.set_ylabel(r'$M_{ej}$', fontsize=12)
                    else:
                        ax.set_ylabel(r'{}'.format(v_ny_arr[n_row]))

                if n_row == n_rows - 1:
                    if v_ny_arr[n_row] == 'vej':
                        ax.set_xlabel(r'$v_{ej}$', fontsize=12)
                    else:
                        ax.set_xlabel(r'{}'.format(v_nx[n_col]))

                if n_col > 0 and n_row < n_rows - 1:
                    sbplot_matrix[n_col][n_row].tick_params(labelbottom=False)
                    sbplot_matrix[n_col][n_row].tick_params(labelleft=False)


        plt.subplots_adjust(hspace=0.0)
        plt.subplots_adjust(wspace=0.0)

        plt.savefig(Paths.plots + 'overall.png', bbox_inches='tight', dpi=128)
        # plt.legend()
        plt.close()



""" --- --- TASK SPECIFIC --- --- """

def plot_table_quantities():

    tbl = TEX_TABLE()
    tbl.print_latex_table()
    exit(1)


    o_models = MODELS()
    lim_dic_in = {'dyn_phase': 'passed', 'names': ['DD2_M13641364_M0_SR', 'LS220_M13641364_M0_SR', 'SLy4_M13641364_M0_SR']}


    sims = o_models.get_selected_models(lim_dic_in)
    o_plots = PLOT_MANY_TASKS(sims)
    o_plots.gen_set["figdir"] = Paths.plots + 'overall/'
    o_plots.gen_set["figname"]  = "Mej.png"
    o_plots.gen_set["figsize"] =  (3.4, 6.8)  # <->, |]

    points_Mej_dic = {
        'name': 'points', 'position': (1, 1), 'title': None, 'cbar': None,
        'it': None, 'v_n_x': 'Lambda', 'v_n_y': 'Mej', 'v_n': None,
        'xmin': None, 'xmax': None, 'ymin': None, 'ymax': None,
        'xscale': None, 'yscale': "log",
        'color': 'EOS', 'marker': 'x', 'annotate': 'name', "arrows_up": None
    }
    points_Mej_bern_dic = {
        'name': 'points', 'position': (1, 1), 'title': None, 'cbar': None,
        'it': None, 'v_n_x': 'Lambda', 'v_n_y': 'Mej_bern', 'v_n': None,
        'xmin': None, 'xmax': None, 'ymin': 1e-3, 'ymax': 5e-2,
        'xscale': None, 'yscale': "log",
        'color': 'EOS', 'marker': 'o', 'annotate': 'name', 'arrows_up': "bern_phase"
    }

    points_Mdisk_dic = {
        'name': 'points', 'position': (2, 1), 'title': None, 'cbar': None,
        'it': None, 'v_n_x': 'Lambda', 'v_n_y': 'Mdisk3D', 'v_n': None,
        'xmin': None, 'xmax': None, 'ymin': None, 'ymax': None,
        'xscale': None, 'yscale': None,
        'color': 'EOS', 'marker': 'x', 'annotate': 'name', "arrows_up": None
    }


    special_MejMdisk = {
        'name': 'special', 'position': (1, 1), 'title': "Masses", 'cbar': None,
        'v_n_x': 'Lambda', 'v_n_y': r'Masses [$M_{\odot}$]',
        'xmin': None, 'xmax': None, 'ymin': None, 'ymax': None,
        'xscale': None, 'yscale': "log", 'legend': True,
    }


    o_plots.set_plot_dics = [special_MejMdisk]

    o_plots.main()

    exit(1)


if __name__ == '__main__':

    # plot_table_quantities()

    # sim_cl = MODELS()

    # lim_dic_in = {'EOS':['SFHo'], 'res':['HR'], 'q':[], 'visc':[]} # 'maximum':'tend'    '-'
    #
    # sims = sim_cl.get_selected_models(lim_dic_in)

    PLOT_MODELS.plot_vn1_vn2_flux_evolution2(extension="_0_b")
    # PLOT_MODELS.plot_final_average(['Lambda'], ["q"]) # 'Mej', 'vej', 'Yeej', 'Mdisk3D',  'vej', 'Yeej'



    # PLOT_MODELS.plot_vn1_vn2_flux_evolution(sims, "SFHo_HR")


    # for isim, sim in sims.iterrows():
    #     print(sim[isim])


# cropped_sims = get_selected_models(lim_dic_in)
#
# for sim in cropped_sims.iterrows():
#     print sim