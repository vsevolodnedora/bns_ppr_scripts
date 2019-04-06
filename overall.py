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
            if v_n == 'EOS':
                for value in lim_dic[v_n]:
                    cropped_sims = cropped_sims[cropped_sims[v_n] == value]
            if v_n == 'res':
                for value in lim_dic[v_n]:
                    cropped_sims = cropped_sims[cropped_sims[v_n] == value]
            if v_n == 'q':
                for value in lim_dic[v_n]:
                    cropped_sims = cropped_sims[cropped_sims[v_n] == value]
            if v_n == 'visc':
                for value in lim_dic[v_n]:
                    cropped_sims = cropped_sims[cropped_sims[v_n] == value]

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


    # print(cropped_sims)

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
    def plot_one_ejecta_prof(ax, sims, v_n1='tend', v_n2='Mej'):

        # for isim, sim in sims.iterrows():
        #     ax.plot(sim[v_n1], sim[v_n2], 'x',
        #             color=eos_color(sim["EOS"]))

        idx = 1
        for isim, sim in sims.iterrows():

            ej = SIM_EJ_HIST(isim)

            color=get_color_for_q(sim["M1"]/sim["M2"])

            ax.plot(ej.time_total_flux, ej.mass_total_flux, '-', color=color,
                    label=r'${}:{}$'.format(idx, str(isim).replace('_', '\_')))
            if not np.isinf(sim["tcoll"]):

                ax.plot(sim[v_n1], sim[v_n2], 'o', color=color)

            # ax.plot(sim[v_n1], sim[v_n2], 'x', color=color)
            ax.annotate(r'${}$'.format(idx), xy=(sim[v_n1], sim[v_n2]), textcoords='data')
            idx += 1

        ax.legend(fancybox=False, loc='best',shadow=False, fontsize=8)

        unique_eos = get_uniqe_eoss(sims["EOS"])

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
    def plot_vn1_vn2_flux_evolution2():

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
                PLOT_MODELS.plot_one_ejecta_prof(ax, sims)
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

        plt.savefig(Paths.plots + 'flux_evolution.pdf', bbox_inches='tight', dpi=128)
        plt.legend()
        plt.close()


    @staticmethod
    def plot_one_ave_final(ax, sims, v_nx, v_ny):

        for isim, sim in sims.iterrows():
            ax.plot(sim[v_nx], sim[v_ny], 'x',
                    color=eos_color(sim["EOS"]))

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

        fig = plt.figure(figsize=(4.5, 3.6))  # (<->; v)

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

        dd2 = SIM_EJ_HIST('DD2_M13641364_M0_SR', extension='_0')
        dd2_w = SIM_EJ_HIST('DD2_M13641364_M0_SR', extension='_0_b_w')

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

                    dd2_mej_tot = dd2.mass_total_flux[-1] + dd2_w.mass_total_flux[-1]
                    # print(dd2.get_par('vel_inf'))
                    # print(dd2_w.get_par('vel_inf'))
                    vel_ave = (dd2.mass_total_flux[-1] * dd2.get_par('vel_inf') +
                               dd2_w.mass_total_flux[-1] * dd2_w.get_par('vel_inf_bern')) / \
                               dd2_mej_tot

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

                    ax.plot(dd2_mej_tot, vel_ave, '*', color='black', label=r'DD2\_M13641364\_M0')

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

        plt.savefig(Paths.plots + 'siegal_like_plot', bbox_inches='tight', dpi=128)
        # plt.legend()
        plt.close()


if __name__ == '__main__':

    # sim_cl = MODELS()

    # lim_dic_in = {'EOS':['SFHo'], 'res':['HR'], 'q':[], 'visc':[]} # 'maximum':'tend'    '-'
    #
    # sims = sim_cl.get_selected_models(lim_dic_in)

    PLOT_MODELS.plot_vn1_vn2_flux_evolution2()
    # PLOT_MODELS.plot_final_average(['Mej'], ['vej']) # 'vej', 'Yeej'



    # PLOT_MODELS.plot_vn1_vn2_flux_evolution(sims, "SFHo_HR")


    # for isim, sim in sims.iterrows():
    #     print(sim[isim])


# cropped_sims = get_selected_models(lim_dic_in)
#
# for sim in cropped_sims.iterrows():
#     print sim