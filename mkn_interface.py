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

path.append(Paths.mkn)
from mkn import MKN


class COMPUTE_LIGHTCURVE():

    def __init__(self, sim):


        # if criteria == '' or criteria == 'geo':
        #     self.path_to_outflow_dir = LISTS.loc_of_sims + sim + '/outflow_{}/'.format(det)
        # elif criteria == 'bern' or criteria == ' bern':
        #     self.path_to_outflow_dir = LISTS.loc_of_sims + sim + '/outflow_{}_b/'.format(det)
        # elif criteria == 'bern dyn' or criteria == ' bern dyn':
        #     self.path_to_outflow_dir = LISTS.loc_of_sims + sim + '/outflow_{}_b_d/'.format(det)
        # elif criteria == 'bern wind' or criteria == ' bern wind':
        #     self.path_to_outflow_dir = LISTS.loc_of_sims + sim + '/outflow_{}_b_w/'.format(det)
        # else:
        #     raise NameError("Criteria '{}' is not recongnized".format(criteria))
        self.path_to_outflow_dir = MakePath.outflow(sim, '_0')
        self.path_to_outflow_dir_psdyn = MakePath.outflow(sim, '_0_b_w')
        dyn_ejecta_profile_fpath = self.path_to_outflow_dir + Files.ejecta_profile
        psdyn_ejecta_profile_fpath=self.path_to_outflow_dir_psdyn + Files.ejecta_profile


        set_dyn_iso_aniso       = "aniso"
        set_psdyn_iso_aniso     = "aniso"
        set_wind_iso_aniso      = ""
        set_secular_iso_aniso   = "aniso"

        self.glob_params    = {}
        self.glob_vars      = {}
        self.ejecta_params  = {}
        self.ejecta_vars    = {}
        self.source_name    = {}

        self.set_glob_par_var_source(True, dyn_ejecta_profile_fpath,
                                     True, psdyn_ejecta_profile_fpath)
        self.set_dyn_par_var(set_dyn_iso_aniso)
        self.set_psdyn_par_var(set_psdyn_iso_aniso)
        self.set_wind_par_war(set_wind_iso_aniso)
        self.set_secular_par_war(set_secular_iso_aniso)

        # self.compute_save_lightcurve(write_output=True)

    '''------------------------------------------'''

    def set_glob_par_var_source(self, NR_data=True, dyn_ej_prof_fpath="",NR2_data=True, psdyn_ej_prof_fpath=""):

        self.glob_params = {'lc model'   : 'grossman',  # model for the lightcurve (grossman or villar)
                       #              'mkn model': 'aniso1comp',  # possible choices: iso1comp, iso2comp, iso3comp, aniso1comp, aniso2comp, aniso3comp
                       'omega frac':0.5,      #
                       'rad shell': False,    #
                       'v_min':     1.e-7,    # minimal velocity for the Grossman model
                       'n_v':       400,      # number of points for the Grossman model
                       'vscale':    'linear', # scale for the velocity in the Grossman model
                       'sigma0':    0.11,     # parameter for the nuclear heating rate
                       'alpha':     1.3,      # parameter for the nuclear heating rate
                       't0eps':     1.3,      # parameter for the nuclear heating rate
                       'cnst_eff':  0.3333,   # parameter for the constant heating efficiency
                       'n slices':  24,       # number for the number of slices along the polar angle [12,18,24,30]
                       'dist slices': 'cos_uniform',  # discretization law for the polar angle [uniform or cos_uniform]
                       'time min':  3600.,    # minimum time [s]
                       'time max':  2000000., # maximum time [s]
                       'n time':    200,      # integer number of bins in time
                       'scale for t': 'log',    # kind of spacing in time [log - linear - measures]
                       'NR_data':   NR_data,     # use (True) or not use (False) NR profiles
                       'NR2_data':  NR2_data,
                       'NR_filename': dyn_ej_prof_fpath,
                       'NR2_filename': psdyn_ej_prof_fpath
                       # path of the NR profiles, necessary if NR_data is True
                       }

        self.source_name = 'AT2017gfo'

        self.glob_vars = {'m_disk':     0.20,
                         'eps0':        1.5e19,
                         'T_floor_LA':  1000.,
                         'T_floor_Ni':  3500.,
                         'a_eps_nuc':   0.5,
                         'b_eps_nuc':   2.5,
                         't_eps_nuc':   1.0}

        # return glob_params, glob_vars, source_name

    def set_dyn_par_var(self, iso_or_aniso):

        if iso_or_aniso == 'iso':
            self.ejecta_params['dynamics'] = {'mass_dist': 'uniform', 'vel_dist': 'uniform', 'op_dist': 'uniform',
                                             'therm_model': 'BKWM', 'eps_ye_dep': True, 'v_law': 'uniform'}
            self.ejecta_vars['dynamics'] = {'xi_disk': None,
                                           'm_ej': 0.003,
                                           'step_angle_mass': None,
                                           'high_lat_flag': None,
                                           'central_vel': 0.24,
                                           'high_lat_vel': None,
                                           'low_lat_vel': None,
                                           'step_angle_vel': None,
                                           'central_op': 30.,
                                           'high_lat_op': None,
                                           'low_lat_op': None,
                                           'step_angle_op': None}
        elif iso_or_aniso == 'aniso':
            self.ejecta_params['dynamics'] = {'mass_dist': 'sin2', 'vel_dist': 'uniform', 'op_dist': 'step',
                                              'therm_model': 'BKWM', 'eps_ye_dep': True, 'v_law': 'power'}
            self.ejecta_vars['dynamics'] = {'xi_disk':           None,
                                           'm_ej':              0.003,
                                           'step_angle_mass':   None,
                                           'high_lat_flag':     None,
                                           'central_vel':       0.33,
                                           'high_lat_vel':      None,
                                           'low_lat_vel':       None,
                                           'step_angle_vel':    None,
                                           'central_op':        None,
                                           'high_lat_op':       1.,
                                           'low_lat_op':        10., # does not work for NR
                                           'step_angle_op': math.radians(45.)}
        elif iso_or_aniso == "":
            pass
        else:
            raise NameError('only iso or aniso')

        # return dyn_ej_pars, dyn_ej_vars

    def set_psdyn_par_var(self, iso_or_aniso):

        if iso_or_aniso == 'iso':
            self.ejecta_params['psdynamics'] = {'mass_dist': 'uniform', 'vel_dist': 'uniform', 'op_dist': 'uniform',
                                             'therm_model': 'BKWM', 'eps_ye_dep': True, 'v_law': 'uniform'}
            self.ejecta_vars['psdynamics'] = {'xi_disk':        None,
                                           'm_ej':              0.01,
                                           'step_angle_mass':   None,
                                           'high_lat_flag':     None,
                                           'central_vel':       0.24,
                                           'high_lat_vel':      None,
                                           'low_lat_vel':       None,
                                           'step_angle_vel':    None,
                                           'central_op':        30.,
                                           'high_lat_op':       None,
                                           'low_lat_op':        None,
                                           'step_angle_op':     None}
        elif iso_or_aniso == 'aniso':
            self.ejecta_params['psdynamics'] = {'mass_dist': 'sin2', 'vel_dist': 'uniform', 'op_dist': 'step',
                                              'therm_model': 'BKWM', 'eps_ye_dep': True, 'v_law': 'power'}
            self.ejecta_vars['psdynamics'] = {'xi_disk':           None,
                                           'm_ej':              0.008,
                                           'step_angle_mass':   None,
                                           'high_lat_flag':     None,
                                           'central_vel':       0.33,
                                           'high_lat_vel':      None,
                                           'low_lat_vel':       None,
                                           'step_angle_vel':    None,
                                           'central_op':        None,
                                           'high_lat_op':       1.,
                                           'low_lat_op':        10.,
                                           'step_angle_op': math.radians(45.)}
        elif iso_or_aniso == "":
            pass
        else:
            raise NameError('only iso or aniso')

        # return dyn_ej_pars, dyn_ej_vars

    def set_wind_par_war(self, iso_or_aniso):

        if iso_or_aniso == 'iso':
            self.ejecta_params['wind'] = {'mass_dist':'uniform','vel_dist':'uniform','op_dist':'uniform',
                         'therm_model':'BKWM','eps_ye_dep':True,'v_law':'power'}
            self.ejecta_vars['wind'] = {'xi_disk' :None,
                         'm_ej': 0.02,
                         'step_angle_mass': None,
                         'high_lat_flag': True,
                         'central_vel': 0.08,
                         'high_lat_vel': None,
                         'low_lat_vel': None,
                         'step_angle_vel': None,
                         'central_op': 1.0,
                         'high_lat_op': None,
                         'low_lat_op': None,
                         'step_angle_op': None}
        elif iso_or_aniso == 'aniso':
            self.ejecta_params['wind'] = {'mass_dist':'step', 'vel_dist':'uniform', 'op_dist':'step',
                         'therm_model':'BKWM','eps_ye_dep':True,'v_law':'power'}
            self.ejecta_vars['wind'] = {'xi_disk':         0.05,
                         'm_ej':            None,
                         'step_angle_mass': math.radians(60.),
                         'high_lat_flag':   True,
                         'central_vel':     0.08,
                         'high_lat_vel':    None,
                         'low_lat_vel':     None,
                         'step_angle_vel':  None,
                         'central_op':      None,
                         'high_lat_op':     0.5,
                         'low_lat_op':      5.0,
                         'step_angle_op':   math.radians(30.)}
        elif iso_or_aniso == "":
            pass
        else:
            raise NameError("iso_or_aniso: {} is not recognized".format(iso_or_aniso))

    def set_secular_par_war(self, iso_or_aniso):

        if iso_or_aniso == 'iso':
            self.ejecta_params['secular'] =  {'mass_dist':'uniform','vel_dist':'uniform','op_dist':'uniform',
                             'therm_model':'BKWM','eps_ye_dep':True,'v_law':'power'}
            self.ejecta_vars['secular'] = {'xi_disk'        :0.4,
                            'm_ej': None,
                            'step_angle_mass': None,
                            'high_lat_flag': None,
                            'central_vel': 0.06,
                            'high_lat_vel': None,
                            'low_lat_vel': None,
                            'step_angle_vel': None,
                            'central_op': 5.0,
                            'low_lat_op': None,
                            'high_lat_op': None,
                            'step_angle_op': None}
        elif iso_or_aniso == 'aniso':
            self.ejecta_params['secular'] = {'mass_dist':'sin2', 'vel_dist':'uniform', 'op_dist':'uniform',
                            'therm_model':'BKWM','eps_ye_dep':True,'v_law':'power'}
            self.ejecta_vars['secular'] = {'xi_disk':          0.2,
                            'm_ej':             None,
                            'step_angle_mass':  None,
                            'high_lat_flag':    None,
                            'central_vel':      0.06,
                            'high_lat_vel':     None,
                            'low_lat_vel':      None,
                            'step_angle_vel':   None,
                            'central_op':       5.0,
                            'low_lat_op':       None,
                            'high_lat_op':      None,
                            'step_angle_op':    None}
        elif iso_or_aniso == "":
            pass
        else:
            raise NameError("iso_or_aniso: {} is not recognized".format(iso_or_aniso))

    def modify_input(self, place, v_n, value):

        ''' Replaces the default value with the given '''

        if place == 'glob_params':
            if not v_n in self.glob_params.keys():
                raise NameError('v_n:{} is not in glob_params:{}'
                                .format(v_n, self.glob_params.keys()))
            self.glob_params[v_n] = value

        if place == 'glob_vars':
            if not v_n in self.glob_vars.keys():
                raise NameError('v_n:{} is not in glob_vars:{}'
                                .format(v_n, self.glob_vars.keys()))
            self.glob_vars[v_n] = value

        # ejecta_params[]
        if place == 'ejecta_params[dynamics]':
            if not v_n in self.ejecta_params['dynamics'].keys():
                raise NameError(
                    'v_n:{} is not in ejecta_params[dynamics]:{}'
                        .format(v_n, self.ejecta_params['dynamics'].keys()))
            self. ejecta_params['dynamics'][v_n] = value

        if place == 'ejecta_params[wind]':
            if not v_n in self.ejecta_params['wind'].keys():
                raise NameError('v_n:{} is not in ejecta_params[wind]:{}'
                                .format(v_n, self.ejecta_params['wind'].keys()))
            self.ejecta_params['wind'][v_n] = value

        if place == 'ejecta_params[secular]':
            if not v_n in self.ejecta_params['secular'].keys():
                raise NameError(
                    'v_n:{} is not in ejecta_params[secular]:{}'
                        .format(v_n, self.ejecta_params['secular'].keys()))
            self.ejecta_params['secular'][v_n] = value

        # shell_vars[]
        if place == 'shell_vars[dynamics]':
            if not v_n in self.ejecta_vars['dynamics'].keys():
                raise NameError('v_n:{} is not in shell_vars[dynamics]:{}'
                                .format(v_n, self.ejecta_vars['dynamics'].keys()))
            self.ejecta_vars['dynamics'][v_n] = value

        if place == 'shell_vars[wind]':
            if not v_n in self.ejecta_vars['wind'].keys():
                raise NameError('v_n:{} is not in shell_vars[wind]:{}'
                                .format(v_n, self.ejecta_vars['wind'].keys()))
            self.ejecta_vars['wind'][v_n] = value

        if place == 'shell_vars[secular]':
            if not v_n in self.ejecta_vars['secular'].keys():
                raise NameError('v_n:{} is not in shell_vars[wind]:{}'
                                .format(v_n, self.ejecta_vars['secular'].keys()))
            self.ejecta_vars['secular'][v_n] = value

    def compute_save_lightcurve(self, write_output = True):
        # glob_params, glob_vars, ejecta_params, shell_vars, source_name_d

        print('I am initializing the model')
        # glob_params, glob_vars, ejecta_params, shell_vars, source_name = self.mkn_parameters()

        # go into the fold with all classes of mkn
        os.chdir(Paths.mkn)
        # from mkn import MKN
        model = MKN(self.glob_params, self.glob_vars, self.ejecta_params, self.ejecta_vars, self.source_name)

        print('I am computing the light curves')
        #    r_ph,L_bol,T_eff = model.lightcurve(ejecta_vars,glob_params['NR_data'],glob_params['NR_filename'])
        r_ph, L_bol, T_eff = model.E.lightcurve(model.angular_distribution,
                                                model.omega_distribution,
                                                model.time,
                                                model.ejecta_vars,
                                                model.ejecta_params,
                                                model.glob_vars,
                                                model.glob_params)

        print('I am computing the likelihood')
        logL = model.log_likelihood(r_ph, T_eff)

        if (write_output):
            print('I am printing out the output')
            model.write_output_h5(r_ph, T_eff, L_bol)
            model.write_filters_h5()

        # copy the result into sim folder and go back into the main script folder

        if (write_output):
            # from shutil import move
            from shutil import copyfile
            # move('./mkn_model.txt', self.path_to_outflow_dir + 'mkn_model.txt')
            copyfile('./mkn_model.h5', self.path_to_outflow_dir + 'mkn_model.h5')

        os.chdir(Paths.home)
        return logL

    # table methods
    def print_latex_table_of_glob_pars(self):

        '''
        \begin{table}
            \footnotesize
                \begin{tabular}[t]{|p{3.2cm}|c|}
                    \hline
                    bla   & 1\\ \hline
                    blubb & 2 \\ \hline
                    bla   & 1\\ \hline
                    blubb & 2 \\ \hline
                    bla   & 1\\ \hline
                    blubb & 2 \\ \hline
                    xxx   & x \\ \hline
                \end{tabular}
                % \hfill
                \begin{tabular}[t]{|p{3.2cm}|c|}
                    \hline
                    bla&1\\ \hline
                    blubb&2 \\ \hline
                    bla&1\\ \hline
                    blubb&2 \\ \hline
                    bla&1\\ \hline
                    blubb&2 \\ \hline
                \end{tabular}
            \hfill
            \caption{99 most frequent hashtags in the data set.}
        \end{table}
        :return:
        '''

        # glob_params, glob_vars, source_name = self.mkn_parameters_glob()

        print('\n')
        print('\\begin{table}[!ht]')
        print('\\footnotesize')

        # table of glob. parameters
        print('\\begin{tabular}[t]{ p{3.2cm} c }')
        print('\\hline')

        # printing rows
        for v_n, value in zip(self.glob_params.keys(), self.glob_params.values()):
            print(' {}  &  {} \\\\'.format(v_n.replace('_', '\\_'), value))
        print('\\hline')

        # table of glob. vars
        print('\\end{tabular}')
        print('\\begin{tabular}[t]{ p{3.2cm} c }')
        print('\\hline')

        for v_n, value in zip(self.glob_vars.keys(), self.glob_vars.values()):
            print(' {}  &  {} \\\\'.format(v_n.replace('_', '\\_'), value))

        print('\\hline')
        print('\\end{tabular}')
        print('\\\caption{Global parameters (left) and global variables (right)}')
        print('\\end{table}')

    def print_latex_table_of_ejecta_pars(self, dynamics=True, wind=True, secular=True):

        print('\n')
        print('\\begin{table}[!ht]')
        print('\\footnotesize')

        if dynamics:

            # dyn_ej_pars, dyn_ej_vars = self.mkn_parameters_dynamics()

            print('\\begin{tabular}[t]{ p{3.cm} c }')
            print('Dynamic & \\\\')
            print('\\hline')

            for v_n, value in zip(self.ejecta_params["dynamical"].keys(), self.ejecta_params["dynamical"].values()):
                print(' {}  &  {} \\\\'.format(v_n.replace('_', '\\_'), value))
            print('\\hline')

            print('\\hline')

            for v_n, value in zip(self.ejecta_vars["dynamical"].keys(), self.ejecta_vars["dynamical"].values()):
                print(' {}  &  {} \\\\'.format(v_n.replace('_', '\\_'), value))
            print('\\hline')

            print('\\end{tabular}')

        if wind:

            # wind_pars, wind_vars = self.mkn_parameters_wind()

            print('\\begin{tabular}[t]{ p{3.cm} c }')
            print('Wind & \\\\')
            print('\\hline')

            for v_n, value in zip(self.ejecta_params["wind"].keys(), self.ejecta_params["wind"].values()):
                print(' {}  &  {} \\\\'.format(v_n.replace('_', '\\_'), value))
            print('\\hline')

            print('\\hline')

            for v_n, value in zip(self.ejecta_vars["wind"].keys(), self.ejecta_vars["wind"].values()):
                print(' {}  &  {} \\\\'.format(v_n.replace('_', '\\_'), value))
            print('\\hline')

            print('\\end{tabular}')

        if secular:

            # secular_pars, secular_vars = self.mkn_parameters_secular()

            print('\\begin{tabular}[t]{ p{3.cm} c }')
            print('Secualr & \\\\')
            print('\\hline')

            for v_n, value in zip(self.ejecta_params["secular"].keys(), self.ejecta_params["secular"].values()):
                print(' {}  &  {} \\\\'.format(v_n.replace('_', '\\_'), value))
            print('\\hline')

            print('\\hline')

            for v_n, value in zip(self.ejecta_vars["secular"].keys(), self.ejecta_vars["secular"].values()):
                print(' {}  &  {} \\\\'.format(v_n.replace('_', '\\_'), value))
            print('\\hline')

            print('\\end{tabular}')

        print('\\caption{Ejecta parameters}')
        print('\\end{table}')


class PLOT_LIGHTCURVE():

    def __init__(self, sim, extension):

        self.model_fpath = MakePath.outflow(sim, extension) + Files.mkn_model

        self.filter_fpath = Paths.mkn + Files.filt_at2017gfo


    def load_mkn_model(self):

        dict_model = {}

        model = h5py.File(self.model_fpath, "r")
        filters_model = []
        for it in model:
            filters_model.append(it)
            dict_model[str(it)] = np.array(model[it])

        # print('\t Following filters are available in mkn_model.h5: \n\t  {}'.format(filters_model))

        return dict_model

    def load_obs_filters(self):

        dict_obs_filters = {}

        obs_filters = h5py.File(self.filter_fpath, "r")

        filters_model = []
        for it in obs_filters:
            filters_model.append(it)
            arr = np.array(obs_filters[it])
            # print(arr.shape)
            dict_obs_filters[str(it)] = np.array(obs_filters[it])

        # print('\t Following filters are available in AT2017gfo.h5: \n\t  {}'.format(filters_model))

        return dict_obs_filters

    @staticmethod
    def get_min_max_for_all_band(dict_model, band):

        arr = np.zeros(len(dict_model['time'])) # should always be there

        for filter in dict_model.keys():
            if filter.split('_')[0] == band:
                arr = np.vstack((arr, dict_model[filter]))
        #
        # if len(arr) <= len(dict_model['time']):
        #     raise ValueError("No min/max found for a band: {}\n Use on of the {}".format(band, dict_model.keys()))

        arr = np.delete(arr, 0, 0)

        maxs = []
        for i in range(len(arr[0,:])):
            maxs = np.append(maxs, arr[:,i].max())

        mins = []
        for i in range(len(arr[0,:])):
            mins = np.append(mins, arr[:,i].min())

        return mins, maxs

    def plot_for_all_band(self, band, plot_as_band = True, fname = 'tst10'):

        dict_obs_filters = self.load_obs_filters()

        dict_model = self.load_mkn_model()

        # plot models as a continouse band of gray color
        if plot_as_band == True:
            mins, maxs = self.get_min_max_for_all_band(dict_model, band)
            plt.fill_between(dict_model['time'], maxs, mins, alpha=.5, color='gray')

            for filter in dict_model.keys():
                if filter.split('_')[0] == band and filter in dict_obs_filters.keys():
                    plt.errorbar(dict_obs_filters[filter][:, 0], dict_obs_filters[filter][:, 1],
                                 yerr=dict_obs_filters[filter][:, 2],
                                 fmt='o', color='black')
        # plot models as separate lines with different colors and labels
        else:
            for filter in dict_model.keys():
                if filter.split('_')[0] == band and filter in dict_obs_filters.keys():
                    plt.errorbar(dict_obs_filters[filter][:, 0], dict_obs_filters[filter][:, 1],
                                 yerr=dict_obs_filters[filter][:, 2],
                                 fmt='o', color='black', label='filter')

                    plt.plot(dict_model['time'], dict_model[filter], '-', label=filter)

        plt.minorticks_on()
        plt.xticks(fontsize=12)
        plt.yticks(fontsize=12)
        plt.tight_layout()
        plt.tick_params(axis='both', which='both', labelleft=True,
                        labelright=False, tick1On=True, tick2On=True,
                        labelsize=12, direction='in')  # labeltop
        plt.gca().invert_yaxis()
        plt.xscale("log")
        # plt.ylim(ymin=1e-5, ymax=2e-1)
        # plt.xlim(xmin=50, xmax=210)
        plt.ylabel("AB magnitude at 40 Mpc", fontsize=12)
        plt.xlabel("time [s]", fontsize=12)
        plt.legend(loc='best', numpoints=1)

        plt.savefig('{}.png'.format(fname), bbox_inches='tight', dpi=128)
        plt.close()

    def plot_for_several_bands(self, bands, fname):

        dict_obs_filters = self.load_obs_filters()

        dict_model = self.load_mkn_model()

        def plot_for_one_band(ax, band, color):

            mins, maxs = self.get_min_max_for_all_band(dict_model, band)
            plt.fill_between(dict_model['time'], maxs, mins, alpha=.5, color=color, label=band)

            for filter in dict_model.keys():
                # print("filter:{} ")
                if filter.split('_')[0] == band and filter in dict_obs_filters.keys():
                    ax.errorbar(dict_obs_filters[filter][:, 0], dict_obs_filters[filter][:, 1],
                                 yerr=dict_obs_filters[filter][:, 2],
                                 fmt='o', color=color)

        plt.figure()
        ax = plt.subplot(111)

        for band, color in zip(bands, clrs(bands)):
            plot_for_one_band(ax, band, color)

        plt.minorticks_on()
        plt.xticks(fontsize=12)
        plt.yticks(fontsize=12)
        plt.tight_layout()
        plt.tick_params(axis='both', which='both', labelleft=True, labelright=False,
                        tick1On=True, tick2On=True, labelsize=12, direction='in')  # labeltop
        plt.gca().invert_yaxis()
        plt.xscale("log")
        plt.ylim(ymin=25)#, ymax=2e-1)
        # plt.xlim(xmin=50, xmax=210)
        plt.ylabel(r"AB magnitude at 40 Mpc", fontsize=12)
        plt.xlabel(r"time [days]", fontsize=12)
        plt.legend(loc='best', numpoints=1)

        plt.savefig('{}{}.png'.format(Paths.plots, fname), bbox_inches='tight', dpi=128)
        plt.close()


class SIM_KILONOVA():

    def __init__(self, sim, criteria, det):

        self.sim = sim
        self.det = det

        if criteria == '' or criteria == 'geo':
            self.path_to_outflow_dir = LISTS.loc_of_sims + sim + '/outflow_{}/'.format(det)
        elif criteria == 'bern' or criteria == ' bern':
            self.path_to_outflow_dir = LISTS.loc_of_sims + sim + '/outflow_{}_b/'.format(det)
        elif criteria == 'bern dyn' or criteria == ' bern dyn':
            self.path_to_outflow_dir = LISTS.loc_of_sims + sim + '/outflow_{}_b_d/'.format(det)
        elif criteria == 'bern wind' or criteria == ' bern wind':
            self.path_to_outflow_dir = LISTS.loc_of_sims + sim + '/outflow_{}_b_w/'.format(det)
        else:
            raise NameError("Criteria '{}' is not recongnized".format(criteria))

        self.set_dynamics  = 'iso' # only central parameters to be changed [everything is UNIFORM]
        self.set_wind      = 'iso'
        self.set_seculiar  = 'iso'

        self.set_use_NR_dyn_ej_distrib = False
        self.set_use_NR_dyn_ej_tot_m   = False
        self.set_use_NR_disk_mass      = False



        SIM.__init__(self, sim)


        self.outflow_dir = "{}{}/outflow_{}/".format(LISTS.loc_of_sims, self.sim, int(self.det))

        self.ejecta_profile = "ejecta_profile.dat"
        self.ejecta_profile_fpath = self.outflow_dir + self.ejecta_profile

        self.kilonova_fname = 'mkn_model.h5'
        self.kilonovae_path = self.outflow_dir + self.kilonova_fname

        self.mkn_dir = '/data01/numrel/vsevolod.nedora/sim_anal_tst/macrokilonova_bayes/mk_source/source/'
        self.scr_dir = '/data01/numrel/vsevolod.nedora/sim_anal_tst/Postprocessing/'

        self.filter_fname = 'AT2017gfo.h5'
        self.obs_filter_path = self.mkn_dir + self.filter_fname


    def load_ej_profile_for_mkn(self):
        th, mass, vel, ye = np.loadtxt(self.ejecta_profile_fpath, unpack=True, usecols=(0, 1, 2, 3))
        return th, mass, vel, ye

    # tech func to check how the smoothing actually done
    def smooth_profile(self, mass):

        lmass = np.log10(mass)

        mass_smooth = []
        for i in range(len(mass)):
            if (i == 0):
                mass_smooth.append(10. ** ((lmass[0])))
            elif (i == 1):
                mass_smooth.append(10. ** ((lmass[i - 1] + lmass[i] + lmass[i + 1]) / 3.))
            elif (i == 2):
                mass_smooth.append(10. ** ((lmass[i - 2] + lmass[i - 1] + lmass[i] + lmass[i + 1] + lmass[i + 2]) / 5.))
            elif (i == len(mass) - 3):
                mass_smooth.append(10. ** ((lmass[i - 2] + lmass[i - 1] + lmass[i] + lmass[i + 1] + lmass[i + 2]) / 5.))
            elif (i == len(mass) - 2):
                mass_smooth.append(10. ** ((lmass[i - 1] + lmass[i] + lmass[i + 1]) / 3.))
            elif (i == len(mass) - 1):
                mass_smooth.append(10. ** ((lmass[i])))
            else:
                mass_smooth.append(10. ** ((lmass[i - 3] + lmass[i - 2] + lmass[i - 1] + lmass[i] + lmass[i + 1] +
                                            lmass[i + 2] + lmass[i + 3]) / 7.))
        mass_smooth = np.asarray(mass_smooth)
        # tmp1 = np.sum(mass)
        # tmp2 = np.sum(mass_smooth)
        # mass_smooth = tmp1 / tmp2 * mass_smooth

        return mass_smooth

    # tech func to show how the smoothing actually done
    def plot_test_smooth_profile(self):

        # loading original profiles
        th, mass, vel, ye = self.load_ej_profile_for_mkn()

        th *= 180 / np.pi
        # th -= 90

        # HAVE NO IDEA WHI THIS EXISTS
        for i in range(len(th)):
            if (mass[i] < 1.e-9):
                mass[i] = 1.e-9
                vel[i] = 0.1
                ye[i] = 0.1


        # smoothing data
        mass_smooth = self.smooth_profile(mass)
        tmp1 = np.sum(mass)
        tmp2 = np.sum(mass_smooth)
        mass_smooth = tmp1 / tmp2 * mass_smooth

        ye_smooth = self.smooth_profile(ye)
        tmp1 = np.sum(ye * mass)
        tmp2 = np.sum(ye_smooth * mass)
        print(ye_smooth)
        ye_smooth = tmp1 / tmp2 * ye_smooth

        vel_smooth = self.smooth_profile(vel)
        tmp1 = np.sum(vel * mass)
        tmp2 = np.sum(vel_smooth)
        vel_smooth = (tmp1 / tmp2 * vel_smooth) / mass_smooth

        # central angle (plane of orbit)
        # th_central = []
        # for a in ang_dist:
        #     th_central.append(0.5 * (a[1] + a[0]))

        # cropping everything for just above orbital plane (0-> 90)
        idx = Math.find_nearest_index(th, 90)
        th = th[idx:] - 90 # offsetting to orbital plane
        mass = mass[idx:]
        mass_smooth = mass_smooth[idx:]

        ye = ye[idx:]
        ye_smooth = ye_smooth[idx:]

        vel = vel[idx:]
        vel_smooth = vel_smooth[idx:]

        # Plot Results of the smoothing for comparison
        f, (ax1, ax2, ax3) = plt.subplots(3, sharex=True) # sharey=True

        ax1.plot(th, mass, '-.', color='gray')
        ax1.plot(th, mass_smooth, '-', color='black', label='smoothed')
        ax1.set_ylabel(LABELS.labels('mass'))
        ax1.tick_params(labelsize=12)
        ax1.legend(loc='best', numpoints=1)

        ax2.plot(th, ye, '-.', color='gray')
        ax2.plot(th, ye_smooth, '-', color='black', label='smoothed')
        ax2.set_ylabel(LABELS.labels("ye"))
        ax2.tick_params(labelsize=12)
        ax2.legend(loc='best', numpoints=1)

        ax3.plot(th, vel, '-.', color='gray')
        ax3.plot(th, vel_smooth, '-', color='black', label='smoothed')
        ax3.set_ylabel(LABELS.labels('vel'))
        ax3.tick_params(labelsize=12)
        ax3.legend(loc='best', numpoints=1)


        f.subplots_adjust(hspace=0)

        plt.setp([a.get_xticklabels() for a in f.axes[:-1]], visible=False)
        plt.xlabel(LABELS.labels("theta"))
        plt.tick_params(axis='both', which='both', labelleft=True, labelright=False, tick1On=True, tick2On=True,
                        labelsize=12)  # labeltop

        plt.minorticks_on()
        plt.savefig('./tst3.png', bbox_inches='tight', dpi=128)
        plt.close()

    def modify_input(self, place, v_n, value, glob_params, glob_vars, ejecta_params, ejecta_vars):

        # load default parameters
        # glob_params, glob_vars, ejecta_params, shell_vars, source_name = self.mkn_parameters()

        if place == 'glob_params':
            if not v_n in glob_params.keys():
                raise NameError('v_n:{} is not in glob_params:{}'.format(v_n, glob_params.keys()))
            glob_params[v_n] = value

        if place == 'glob_vars':
            if not v_n in glob_vars.keys():
                raise NameError('v_n:{} is not in glob_vars:{}'.format(v_n, glob_vars.keys()))
            glob_vars[v_n] = value

        # ejecta_params[]
        if place == 'ejecta_params[dynamics]':
            if not v_n in ejecta_params['dynamics'].keys():
                raise NameError(
                    'v_n:{} is not in ejecta_params[dynamics]:{}'.format(v_n, ejecta_params['dynamics'].keys()))
            ejecta_params['dynamics'][v_n] = value

        if place == 'ejecta_params[wind]':
            if not v_n in ejecta_params['wind'].keys():
                raise NameError('v_n:{} is not in ejecta_params[wind]:{}'.format(v_n, ejecta_params['wind'].keys()))
            ejecta_params['wind'][v_n] = value

        if place == 'ejecta_params[secular]':
            if not v_n in ejecta_params['secular'].keys():
                raise NameError(
                    'v_n:{} is not in ejecta_params[secular]:{}'.format(v_n, ejecta_params['secular'].keys()))
            ejecta_params['secular'][v_n] = value

        # shell_vars[]
        if place == 'shell_vars[dynamics]':
            if not v_n in ejecta_vars['dynamics'].keys():
                raise NameError('v_n:{} is not in shell_vars[dynamics]:{}'.format(v_n, ejecta_vars['dynamics'].keys()))
            ejecta_vars['dynamics'][v_n] = value

        if place == 'shell_vars[wind]':
            if not v_n in ejecta_vars['wind'].keys():
                raise NameError('v_n:{} is not in shell_vars[wind]:{}'.format(v_n, ejecta_vars['wind'].keys()))
            ejecta_vars['wind'][v_n] = value

        if place == 'shell_vars[secular]':
            if not v_n in ejecta_vars['secular'].keys():
                raise NameError('v_n:{} is not in shell_vars[wind]:{}'.format(v_n, ejecta_vars['secular'].keys()))
            ejecta_vars['secular'][v_n] = value

        return glob_params, glob_vars, ejecta_params, ejecta_vars
    def modify_mkn_parameters(self, glob_params, glob_vars, ejecta_params, ejecta_vars):
        '''
            Changes MKN input parameters (one by one)
        '''

        print ("Error! Modification of the mkn parameters is not finished... FINISH")
        return glob_params, glob_vars, ejecta_params, ejecta_vars

        if self.set_use_NR_dyn_ej_distrib or self.set_use_NR_disk_mass or self.set_use_NR_dyn_ej_tot_m:
            ej = SIM_EJECTA_PARS(self.sim, self.det)
        else:
            ej = None

        # set to use the Num Rel Dyn Ej destribution
        if self.set_use_NR_dyn_ej_distrib:
            # print(glob_params['NR_data'])
            glob_params, glob_vars, ejecta_params, ejecta_vars = \
                self.modify_input('glob_params', 'NR_data', True,
                                  glob_params, glob_vars, ejecta_params, ejecta_vars)
            # print(glob_params['NR_data'])
            # exit(1)

        # set to use Num Rel mass of the disk
        if self.set_use_NR_disk_mass:
            glob_params, glob_vars, ejecta_params, ejecta_vars = \
                self.modify_input('glob_vars', 'm_disk', np.float(ej.dic['Mdisk']) / 1e2,
                                  glob_params, glob_vars, ejecta_params, ejecta_vars)

        # set to use Num Rel total mass of Dyn. Ej.
        if self.set_use_NR_dyn_ej_tot_m:
            glob_params, glob_vars, ejecta_params, ejecta_vars = \
                self.modify_input('ejecta_vars[dynamics]', 'm_ej', np.float(ej.dic['Mej_tot']) / 1e2,
                                  glob_params, glob_vars, ejecta_params, ejecta_vars)

        return glob_params, glob_vars, ejecta_params, ejecta_vars

    # basic set of parameters for mkn class --- --- ---
    def mkn_parameters_glob(self):

        glob_params = {'lc model'   : 'grossman',  # model for the lightcurve (grossman or villar)
                       #              'mkn model': 'aniso1comp',  # possible choices: iso1comp, iso2comp, iso3comp, aniso1comp, aniso2comp, aniso3comp
                       'omega frac':0.5,      #
                       'rad shell': False,    #
                       'v_min':     1.e-7,    # minimal velocity for the Grossman model
                       'n_v':       400,      # number of points for the Grossman model
                       'vscale':    'linear', # scale for the velocity in the Grossman model
                       'sigma0':    0.11,     # parameter for the nuclear heating rate
                       'alpha':     1.3,      # parameter for the nuclear heating rate
                       't0eps':     1.3,      # parameter for the nuclear heating rate
                       'cnst_eff':  0.3333,   # parameter for the constant heating efficiency
                       'n slices':  25,       # number for the number of slices along the polar angle [12,18,24,30]
                       'dist slices': 'cos_uniform',  # discretization law for the polar angle [uniform or cos_uniform]
                       'time min':  3600.,    # minimum time [s]
                       'time max':  1000000., # maximum time [s]
                       'n time':    200,      # integer number of bins in time
                       'scale for t': 'log',    # kind of spacing in time [log - linear - measures]
                       'NR_data':   True,     # use (True) or not use (False) NR profiles
                       'NR_filename': self.ejecta_profile_fpath
                       # path of the NR profiles, necessary if NR_data is True
                       }

        source_name = 'AT2017gfo'

        glob_vars = {'m_disk':      0.12,
                     'eps0':        1.5e19,
                     'T_floor_LA':  1000.,
                     'T_floor_Ni':  3500.,
                     'a_eps_nuc':   0.5,
                     'b_eps_nuc':   2.5,
                     't_eps_nuc':   1.0}

        return glob_params, glob_vars, source_name
    def mkn_parameters_dynamics(self):

        if self.set_dynamics == 'iso':
            dyn_ej_pars = {'mass_dist': 'uniform', 'vel_dist': 'uniform', 'op_dist': 'uniform',
                           'therm_model': 'BKWM', 'eps_ye_dep': True, 'v_law': 'uniform'}
            dyn_ej_vars = {'xi_disk': None,
                           'm_ej': 0.04,
                           'step_angle_mass': None,
                           'high_lat_flag': None,
                           'central_vel': 0.24,
                           'high_lat_vel': None,
                           'low_lat_vel': None,
                           'step_angle_vel': None,
                           'central_op': 30.,
                           'high_lat_op': None,
                           'low_lat_op': None,
                           'step_angle_op': None}
        elif self.set_dynamics == 'aniso':
            dyn_ej_pars = {'mass_dist': 'sin2', 'vel_dist': 'uniform', 'op_dist': 'step',
                           'therm_model': 'BKWM', 'eps_ye_dep': True, 'v_law': 'power'}
            dyn_ej_vars = {'xi_disk':           None,
                           'm_ej':              0.3*1e-2,
                           'step_angle_mass':   None,
                           'high_lat_flag':     None,
                           'central_vel':       0.33,
                           'high_lat_vel':      None,
                           'low_lat_vel':       None,
                           'step_angle_vel':    None,
                           'central_op':        None,
                           'high_lat_op':       1.,
                           'low_lat_op':        30.,
                           'step_angle_op': math.radians(45.)}
        else:
            raise NameError('only iso or aniso')

        return dyn_ej_pars, dyn_ej_vars
    def mkn_parameters_wind(self):

        if self.set_wind == 'iso':
            wind_pars = {'mass_dist':'uniform','vel_dist':'uniform','op_dist':'uniform',
                         'therm_model':'BKWM','eps_ye_dep':True,'v_law':'power'}
            wind_vars = {'xi_disk' :None,
                         'm_ej': 0.02,
                         'step_angle_mass': None,
                         'high_lat_flag': True,
                         'central_vel': 0.08,
                         'high_lat_vel': None,
                         'low_lat_vel': None,
                         'step_angle_vel': None,
                         'central_op': 1.0,
                         'high_lat_op': None,
                         'low_lat_op': None,
                         'step_angle_op': None}
        else:
            wind_pars = {'mass_dist':'step', 'vel_dist':'uniform', 'op_dist':'step',
                         'therm_model':'BKWM','eps_ye_dep':True,'v_law':'power'}
            wind_vars = {'xi_disk':         0.05,
                         'm_ej':            None,
                         'step_angle_mass': math.radians(60.),
                         'high_lat_flag':   True,
                         'central_vel':     0.08,
                         'high_lat_vel':    None,
                         'low_lat_vel':     None,
                         'step_angle_vel':  None,
                         'central_op':      None,
                         'high_lat_op':     0.5,
                         'low_lat_op':      5.0,
                         'step_angle_op':   math.radians(30.)}

        return wind_pars, wind_vars
    def mkn_parameters_secular(self):

        if self.set_wind == 'iso':
            secular_pars =  {'mass_dist':'uniform','vel_dist':'uniform','op_dist':'uniform',
                             'therm_model':'BKWM','eps_ye_dep':True,'v_law':'power'}
            secular_vars = {'xi_disk'        :0.4,
                            'm_ej': None,
                            'step_angle_mass': None,
                            'high_lat_flag': None,
                            'central_vel': 0.06,
                            'high_lat_vel': None,
                            'low_lat_vel': None,
                            'step_angle_vel': None,
                            'central_op': 5.0,
                            'low_lat_op': None,
                            'high_lat_op': None,
                            'step_angle_op': None}
        else:
            secular_pars = {'mass_dist':'sin2', 'vel_dist':'uniform', 'op_dist':'uniform',
                            'therm_model':'BKWM','eps_ye_dep':True,'v_law':'power'}
            secular_vars = {'xi_disk':          0.2,
                            'm_ej':             None,
                            'step_angle_mass':  None,
                            'high_lat_flag':    None,
                            'central_vel':      0.06,
                            'high_lat_vel':     None,
                            'low_lat_vel':      None,
                            'step_angle_vel':   None,
                            'central_op':       5.0,
                            'low_lat_op':       None,
                            'high_lat_op':      None,
                            'step_angle_op':    None}

        return secular_pars, secular_vars
    def mkn_parameters(self):

        glob_params, glob_vars, source_name = self.mkn_parameters_glob()

        ejecta_params = {}
        ejecta_vars = {}

        if self.set_dynamics != None:
            dyn_ej_pars, dyn_ej_vars = self.mkn_parameters_dynamics()
            ejecta_params['dynamics'] = dyn_ej_pars
            ejecta_vars['dynamics'] = dyn_ej_vars

        if self.set_wind != None:
            wind_pars, wind_vars = self.mkn_parameters_wind()
            ejecta_params['wind'] = wind_pars
            ejecta_vars['wind'] = wind_vars

        if self.set_seculiar != None:
            secular_pars, secular_vars = self.mkn_parameters_secular()
            ejecta_params['secular'] = secular_pars
            ejecta_vars['secular'] = secular_vars

        glob_params, glob_vars, ejecta_params, ejecta_vars = \
            self.modify_mkn_parameters(glob_params, glob_vars, ejecta_params, ejecta_vars)

        return glob_params, glob_vars, ejecta_params, ejecta_vars, source_name

    # parses all parameters into the MKN --- --- ---
    def mkn_interface(self, glob_params, glob_vars, ejecta_params, shell_vars, source_name, write_output = True):
        # glob_params, glob_vars, ejecta_params, shell_vars, source_name_d

        print('I am initializing the model')
        # glob_params, glob_vars, ejecta_params, shell_vars, source_name = self.mkn_parameters()

        # go into the fold with all classes of mkn
        os.chdir(self.mkn_dir)
        model = MKN(glob_params, glob_vars, ejecta_params, shell_vars, source_name)


        print('I am computing the light curves')
        #    r_ph,L_bol,T_eff = model.lightcurve(ejecta_vars,glob_params['NR_data'],glob_params['NR_filename'])
        r_ph, L_bol, T_eff = model.E.lightcurve(model.angular_distribution,
                                                model.omega_distribution,
                                                model.time,
                                                model.ejecta_vars,
                                                model.ejecta_params,
                                                model.glob_vars,
                                                model.glob_params)

        print('I am computing the likelihood')
        logL = model.log_likelihood(r_ph, T_eff)


        if (write_output):
            print('I am printing out the output')
            model.write_output_h5(r_ph, T_eff, L_bol)
            model.write_filters_h5()

        # copy the result into sim folder and go back into the main script folder

        if (write_output):
            from shutil import copyfile
            copyfile('./mkn_model.txt', self.outflow_dir + 'mkn_model.txt')
            copyfile('./mkn_model.h5', self.outflow_dir + 'mkn_model.h5')

        os.chdir(self.scr_dir)
        return logL

    def mkn_parameters_old(self):

        '''
        It is a COPY of a method of calling MKN() class from @aperego code "macrokilonova_bayes"
        Modifications:
            - addtd __init__.py into the MAIN FOLDER so I can use it as a package
            - copy results into simulation folder
            - save data from all filters into h5 file to be used for quick plotting
            - save model in h5 file (as normal output has more filters than columns, (do not know why)
        :return:
        '''

        ej = SIM_EJECTA(self.sim, self.det)

        glob_params = {'lc model': 'grossman',  # model for the lightcurve (grossman or villar)
        #              'mkn model': 'aniso1comp',  # possible choices: iso1comp, iso2comp, iso3comp, aniso1comp, aniso2comp, aniso3comp
                       'omega frac': 0.5,       #
                       'rad shell': False,      #
                       'v_min': 1.e-7,          # minimal velocity for the Grossman model
                       'n_v': 400,              # number of points for the Grossman model
                       'vscale': 'linear',      # scale for the velocity in the Grossman model
                       'sigma0': 0.11,          # parameter for the nuclear heating rate
                       'alpha': 1.3,            # parameter for the nuclear heating rate
                       't0eps': 1.3,            # parameter for the nuclear heating rate
                       'cnst_eff': 0.3333,      # parameter for the constant heating efficiency
                       'n slices': 12,          # number for the number of slices along the polar angle [12,18,24,30]
                       'dist slices': 'cos_uniform',  # discretization law for the polar angle [uniform or cos_uniform]
                       'time min': 3600.,       # minimum time [s]
                       'time max': 2000000.,    # maximum time [s]
                       'n time': 200,           # integer number of bins in time
                       'scale for t': 'log',    # kind of spacing in time [log - linear - measures]
                       'NR_data': True,         # use (True) or not use (False) NR profiles
                       'NR_filename': '../../../'+self.outflow_dir+'ejecta_profile.dat' # path of the NR profiles, necessary if NR_data is True
                       }
        source_name = 'AT2017gfo'

        glob_vars = {'m_disk': np.float(ej.dic['Mdisk'])/1e2, # 0.08
                     'eps0': 1.5e19,
                     'T_floor_LA': 1000.,
                     'T_floor_Ni': 3500.,
                     'a_eps_nuc': 0.5,
                     'b_eps_nuc': 2.5,
                     't_eps_nuc': 1.0}

        ejecta_params = {}
        ejecta_params['dynamics'] = {'mass_dist': 'sin2', 'vel_dist': 'uniform', 'op_dist': 'step', 'therm_model': 'BKWM', 'eps_ye_dep': True, 'v_law': 'power'}  # fixed
        ejecta_params['wind']     = {'mass_dist': 'step', 'vel_dist': 'uniform', 'op_dist': 'step', 'therm_model': 'BKWM', 'eps_ye_dep': True, 'v_law': 'power'}  # fixed
        ejecta_params['secular']  = {'mass_dist': 'sin2', 'vel_dist': 'uniform', 'op_dist': 'uniform', 'therm_model': 'BKWM', 'eps_ye_dep': True, 'v_law': 'power'}  # fixed

        shell_vars={}
        shell_vars['dynamics'] = {'xi_disk': None,
                                  'm_ej': np.float(ej.dic['Mej_tot'])/1e2, #0.1,         # dynamical ejecta mass [Msun], to be changed!
                                  'step_angle_mass': None,
                                  'high_lat_flag': None,
                                  'central_vel': 0.2,  # dynamical ejecta mean velocity [c], to be changed!
                                  'low_lat_vel': None,
                                  'high_lat_vel': None,
                                  'step_angle_vel': None,
                                  'central_op': None,
                                  'low_lat_op': 30.,    # equatorial opacity, usually for Lanthanide rich matter, > 10. [cm^2/g]
                                  'high_lat_op': 1.0,   # polar opacity, usually Lanthanide-free, ~ 1. [cm^2/g]
                                  'step_angle_op': math.radians(45.)}  # fixed
        shell_vars['wind'] = {'xi_disk': 0.001,          # nu-driven wind simulations says 0.01-0.05, but B field can increase it (<~ 0.1)
                              'm_ej': None,
                              'step_angle_mass': math.radians(60.),  # fixed
                              'high_lat_flag': True,    # fixed
                              'central_vel': 0.067,      # 0.04 c < v < 0.1
                              'low_lat_vel': None,
                              'high_lat_vel': None,
                              'step_angle_vel': None,
                              'central_op': None,
                              'low_lat_op': 5.0,        # 1-5
                              'high_lat_op': 0.5,       # 0.1-1
                              'step_angle_op': math.radians(30.)}  # fixed
        shell_vars['secular'] = {'xi_disk': 0.4,        # 0.1-0.4
                                 'm_ej': None,
                                 'step_angle_mass': None,
                                 'high_lat_flag': None,
                                 'central_vel': 0.06,   # 0.04 < v < 0.1
                                 'low_lat_vel': None,
                                 'high_lat_vel': None,
                                 'step_angle_vel': None,
                                 'central_op': 1.0,     # 1-10
                                 'low_lat_op': None,
                                 'high_lat_op': None,
                                 'step_angle_op': None}

        return glob_params, glob_vars, ejecta_params, shell_vars, source_name

    # def mkn_interface_old(self):
    #     '''
    #     It is a COPY of a method of calling MKN() class from @aperego code "macrokilonova_bayes"
    #     Modifications:
    #         - addtd __init__.py into the MAIN FOLDER so I can use it as a package
    #         - copy results into simulation folder
    #         - save data from all filters into h5 file to be used for quick plotting
    #         - save model in h5 file (as normal output has more filters than columns, (do not know why)
    #     :return:
    #     '''
    #
    #
    #
    #     # dictionary with the global parameters of the model
    #     glob_params = {'lc model': 'grossman',  # model for the lightcurve (grossman or villar)
    #                    'mkn model': 'aniso1comp',
    #                    # possible choices: iso1comp, iso2comp, iso3comp, aniso1comp, aniso2comp, aniso3comp
    #                    'omega frac': 0.5,  #
    #                    'rad shell': False,  #
    #                    'v_min': 1.e-7,  # minimal velocity for the Grossman model
    #                    'n_v': 400,  # number of points for the Grossman model
    #                    'vscale': 'linear',  # scale for the velocity in the Grossman model
    #                    'sigma0': 0.11,  # parameter for the nuclear heating rate
    #                    'alpha': 1.3,  # parameter for the nuclear heating rate
    #                    't0eps': 1.3,  # parameter for the nuclear heating rate
    #                    'cnst_eff': 0.3333,  # parameter for the constant heating efficiency
    #                    'n slices': 30,  # number for the number of slices along the polar angle [12,18,24,30]
    #                    'dist slices': 'cos_uniform',  # discretization law for the polar angle [uniform or cos_uniform]
    #                    'time min': 3600.,  # minimum time [s]
    #                    'time max': 2000000.,  # maximum time [s]
    #                    'n time': 200,  # integer number of bins in time
    #                    'scale for t': 'log',  # kind of spacing in time [log - linear - measures]
    #                    'NR_data': True,  # use (True) or not use (False) NR profiles
    #                    'NR_filename': '../../../'+self.outflow_dir+'ejecta_profile.dat'#  '../example_NR_data/DD2_M125125_LK/outflow_1/ejecta_profile.dat'
    #                    # path of the NR profiles, necessary if NR_data is True
    #                    }
    #
    #     source_name = 'AT2017gfo'   # name of the source or "default"
    #     # source_name = 'default'  # name of the source or "default"
    #     save_source_h5 = True
    #
    #     # dictionary for the global variables
    #     glob_vars = {'m_disk': 0.08,
    #                  'eps0': 1.5e19,
    #                  'T_floor_LA': 1000.,
    #                  'T_floor_Ni': 3500.,
    #                  'a_eps_nuc': 0.5,
    #                  'b_eps_nuc': 2.5,
    #                  't_eps_nuc': 1.0}
    #
    #
    #     ###############################
    #     # Template for isotropic case #
    #     ###############################
    #
    #     # hardcoded ejecta geometric and thermal parameters for the spherical case
    #     ejecta_params_iso = {}
    #     ejecta_params_iso['dynamics'] = {'mass_dist': 'uniform', 'vel_dist': 'uniform', 'op_dist': 'uniform', 'therm_model': 'BKWM', 'eps_ye_dep': True, 'v_law': 'uniform'}
    #     ejecta_params_iso['wind'] = {'mass_dist': 'uniform', 'vel_dist': 'uniform', 'op_dist': 'uniform', 'therm_model': 'BKWM', 'eps_ye_dep': True, 'v_law': 'power'}
    #     ejecta_params_iso['secular'] = {'mass_dist': 'uniform', 'vel_dist': 'uniform', 'op_dist': 'uniform', 'therm_model': 'BKWM', 'eps_ye_dep': True, 'v_law': 'power'}
    #
    #     # set of shell parameters to be sampled on
    #     ejecta_vars_iso = {}
    #
    #     ejecta_vars_iso['dynamics'] = {'xi_disk': None,
    #                                    'm_ej': 0.04,
    #                                    'step_angle_mass': None,
    #                                    'high_lat_flag': None,
    #                                    'central_vel': 0.24,
    #                                    'high_lat_vel': None,
    #                                    'low_lat_vel': None,
    #                                    'step_angle_vel': None,
    #                                    'central_op': 30.,
    #                                    'high_lat_op': None,
    #                                    'low_lat_op': None,
    #                                    'step_angle_op': None}
    #
    #     ejecta_vars_iso['secular'] = {'xi_disk': None,
    #                                   'm_ej': 0.05,
    #                                   'step_angle_mass': None,
    #                                   'high_lat_flag': None,
    #                                   'central_vel': 0.06,
    #                                   'high_lat_vel': None,
    #                                   'low_lat_vel': None,
    #                                   'step_angle_vel': None,
    #                                   'central_op': 5.0,
    #                                   'low_lat_op': None,
    #                                   'high_lat_op': None,
    #                                   'step_angle_op': None}
    #
    #     ejecta_vars_iso['wind'] = {'xi_disk': None,
    #                                'm_ej': 0.02,
    #                                'step_angle_mass': None,
    #                                'high_lat_flag': True,
    #                                'central_vel': 0.08,
    #                                'high_lat_vel': None,
    #                                'low_lat_vel': None,
    #                                'step_angle_vel': None,
    #                                'central_op': 1.0,
    #                                'high_lat_op': None,
    #                                'low_lat_op': None,
    #                                'step_angle_op': None}
    #
    #     #################################
    #     # Template for anisotropic case #
    #     #################################
    #
    #     # hardcoded ejecta geometric and thermal parameters for the aspherical case
    #     ejecta_params_aniso = {}
    #     ejecta_params_aniso['dynamics'] = {'mass_dist': 'sin2', 'vel_dist': 'uniform', 'op_dist': 'step', 'therm_model': 'BKWM', 'eps_ye_dep': True, 'v_law': 'power'}
    #     ejecta_params_aniso['wind'] = {'mass_dist': 'step', 'vel_dist': 'uniform', 'op_dist': 'step', 'therm_model': 'BKWM', 'eps_ye_dep': True, 'v_law': 'power'}
    #     ejecta_params_aniso['secular'] = {'mass_dist': 'sin2', 'vel_dist': 'uniform', 'op_dist': 'uniform', 'therm_model': 'BKWM', 'eps_ye_dep': True, 'v_law': 'power'}
    #
    #     # set of shell parameters to be sampled on
    #     ejecta_vars_aniso = {}
    #
    #     shell_vars = {}
    #
    #     ejecta_vars_aniso['dynamics'] = {'xi_disk': None,
    #                                      'm_ej': 0.01,
    #                                      'step_angle_mass': None,
    #                                      'high_lat_flag': None,
    #                                      'central_vel': 0.33,
    #                                      'high_lat_vel': None,
    #                                      'low_lat_vel': None,
    #                                      'step_angle_vel': None,
    #                                      'central_op': None,
    #                                      'high_lat_op': 10.,
    #                                      'low_lat_op': 30.,
    #                                      'step_angle_op': math.radians(45.)}
    #
    #     ejecta_vars_aniso['secular'] = {'xi_disk': 0.4,
    #                                     'm_ej': None,
    #                                     'step_angle_mass': None,
    #                                     'high_lat_flag': None,
    #                                     'central_vel': 0.06,
    #                                     'high_lat_vel': None,
    #                                     'low_lat_vel': None,
    #                                     'step_angle_vel': None,
    #                                     'central_op': 5.0,
    #                                     'low_lat_op': None,
    #                                     'high_lat_op': None,
    #                                     'step_angle_op': None}
    #
    #
    #     ejecta_vars_aniso['wind'] = {'xi_disk': 0.001,
    #                                  'm_ej': None,
    #                                  'step_angle_mass': math.radians(60.),
    #                                  'high_lat_flag': True,
    #                                  'central_vel': 0.08,
    #                                  'high_lat_vel': None,
    #                                  'low_lat_vel': None,
    #                                  'step_angle_vel': None,
    #                                  'central_op': None,
    #                                  'high_lat_op': 0.1,
    #                                  'low_lat_op': 5.0,
    #                                  'step_angle_op': math.radians(30.)}
    #
    #
    #     ##########################################################
    #     # choose the appropriate set of parameters and variables #
    #     ##########################################################
    #
    #     if (glob_params['mkn model'] == 'iso1comp'):
    #         ejecta_params = {}
    #         ejecta_vars = {}
    #         ejecta_params['dynamics'] = ejecta_params_iso['dynamics']
    #         ejecta_vars['dynamics'] = ejecta_vars_iso['dynamics']
    #     elif (glob_params['mkn model'] == 'iso2comp'):
    #         ejecta_params = {}
    #         ejecta_vars = {}
    #         ejecta_params['dynamics'] = ejecta_params_iso['dynamics']
    #         ejecta_vars['dynamics'] = ejecta_vars_iso['dynamics']
    #         ejecta_params['secular'] = ejecta_params_iso['secular']
    #         ejecta_vars['secular'] = ejecta_vars_iso['secular']
    #     elif (glob_params['mkn model'] == 'iso3comp'):
    #         ejecta_params = ejecta_params_iso
    #         ejecta_vars = ejecta_vars_iso
    #     elif (glob_params['mkn model'] == 'aniso1comp'):
    #         ejecta_params = {}
    #         ejecta_vars = {}
    #         ejecta_params['dynamics'] = ejecta_params_aniso['dynamics']
    #         ejecta_vars['dynamics'] = ejecta_vars_aniso['dynamics']
    #     elif (glob_params['mkn model'] == 'aniso2comp'):
    #         ejecta_params = {}
    #         ejecta_vars = {}
    #         ejecta_params['dynamics'] = ejecta_params_aniso['dynamics']
    #         ejecta_vars['dynamics'] = ejecta_vars_aniso['dynamics']
    #         ejecta_params['secular'] = ejecta_params_aniso['secular']
    #         ejecta_vars['secular'] = ejecta_vars_aniso['secular']
    #     elif (glob_params['mkn model'] == 'aniso3comp'):
    #         ejecta_params = ejecta_params_aniso
    #         ejecta_vars = ejecta_vars_aniso
    #
    #     print('I am initializing the model')
    #
    #     # go into the fold with all classes of mkn
    #     os.chdir(self.mkn_dir)
    #     model = mkn.MKN(glob_params, glob_vars, ejecta_params, ejecta_vars, source_name)
    #
    #     print('I am computing the light curves')
    #     #    r_ph,L_bol,T_eff = model.lightcurve(ejecta_vars,glob_params['NR_data'],glob_params['NR_filename'])
    #     r_ph, L_bol, T_eff = model.E.lightcurve(model.angular_distribution,
    #                                             model.omega_distribution,
    #                                             model.time,
    #                                             model.ejecta_vars,
    #                                             model.ejecta_params,
    #                                             model.glob_vars,
    #                                             model.glob_params)
    #
    #     print('I am computing the likelihood')
    #     logL = model.log_likelihood(r_ph, T_eff)
    #
    #     write_output = True
    #     if (write_output):
    #         print('I am printing out the output')
    #         # model.write_output(r_ph, T_eff, L_bol)
    #         model.write_output_h5(r_ph, T_eff, L_bol)
    #
    #     plot_separately = False  # Choose to plot all lightcurves in different bands on the same plot
    #     plot_together = False  # or to plot lightcurve and data in each band on different plots
    #
    #     if (plot_separately):
    #         model.plot_output_sep()
    #     elif (plot_together):
    #         model.plot_output_tog()
    #
    #     # added ability to save all filters as h5 file with rescaled magnitudes
    #     if save_source_h5:
    #         model.write_filters_h5()
    #
    #
    #     # copy the result into sim folder and go back into the main script folder
    #     os.chdir(self.scr_dir)
    #
    #     from shutil import copyfile
    #     copyfile('./macrokilonova_bayes/mk_source/source/mkn_model.txt', self.outflow_dir + 'mkn_model.txt')
    #     copyfile('./macrokilonova_bayes/mk_source/source/mkn_model.h5', self.outflow_dir + 'mkn_model.h5')
    #     # os.chdir(self.scr_dir)

    def load_obs_filters(self):

        dict_obs_filters = {}

        obs_filters = h5py.File(self.obs_filter_path, "r")

        filters_model = []
        for it in obs_filters:
            filters_model.append(it)
            arr = np.array(obs_filters[it])
            # print(arr.shape)
            dict_obs_filters[str(it)] = np.array(obs_filters[it])

        # print('\t Following filters are available in AT2017gfo.h5: \n\t  {}'.format(filters_model))

        return dict_obs_filters

    def save_one_model(self):

        glob_params, glob_vars, ejecta_params, shell_vars, source_name = self.mkn_parameters()
        self.mkn_interface(glob_params, glob_vars, ejecta_params, shell_vars, source_name, True)

    def load_mkn_model(self):

        dict_model = {}

        model = h5py.File(self.kilonovae_path, "r")
        filters_model = []
        for it in model:
            filters_model.append(it)
            dict_model[str(it)] = np.array(model[it])

        # print('\t Following filters are available in mkn_model.h5: \n\t  {}'.format(filters_model))

        return dict_model

    def plot_one_model(self, filter):

        dict_obs_filters = self.load_obs_filters()

        dict_model = self.load_mkn_model()

        plt.errorbar(dict_obs_filters[filter][:,0], dict_obs_filters[filter][:,1], yerr=dict_obs_filters[filter][:,2],
                     color = 'black', fmt='o')

        plt.plot(dict_model['time'], dict_model[filter], '-', color='gray')

        plt.minorticks_on()
        plt.xticks(fontsize=12)
        plt.yticks(fontsize=12)
        plt.tight_layout()
        plt.tick_params(axis='both', which='both', labelleft=True, labelright=False, tick1On=True, tick2On=True, labelsize=12, direction='in')  # labeltop
        plt.gca().invert_yaxis()
        plt.xscale("log")
        # plt.ylim(ymin=1e-5, ymax=2e-1)
        # plt.xlim(xmin=50, xmax=210)
        plt.ylabel("AB magnitude at 40 Mpc", fontsize=12)
        plt.xlabel("time [s]", fontsize=12)
        plt.legend(loc='best', numpoints=1)

        plt.savefig('./tst4.png', bbox_inches='tight', dpi=128)
        plt.close()

    @staticmethod
    def get_min_max_for_all_band(dict_model, band):

        arr = np.zeros(len(dict_model['time'])) # should always be there

        for filter in dict_model.keys():
            if filter.split('_')[0] == band:
                arr = np.vstack((arr, dict_model[filter]))
        arr = np.delete(arr, 0, 0)

        maxs = []
        for i in range(len(arr[0,:])):
            maxs = np.append(maxs, arr[:,i].max())

        mins = []
        for i in range(len(arr[0,:])):
            mins = np.append(mins, arr[:,i].min())

        return mins, maxs

    def plot_for_all_band(self, band, plot_as_band = True, fname = 'tst10'):

        dict_obs_filters = self.load_obs_filters()

        dict_model = self.load_mkn_model()

        # plot models as a continouse band of gray color
        if plot_as_band == True:
            mins, maxs = self.get_min_max_for_all_band(dict_model, band)
            plt.fill_between(dict_model['time'], maxs, mins, alpha=.5, color='gray')

            for filter in dict_model.keys():
                if filter.split('_')[0] == band and filter in dict_obs_filters.keys():
                    plt.errorbar(dict_obs_filters[filter][:, 0], dict_obs_filters[filter][:, 1],
                                 yerr=dict_obs_filters[filter][:, 2],
                                 fmt='o', color='black')
        # plot models as separate lines with different colors and labels
        else:
            for filter in dict_model.keys():
                if filter.split('_')[0] == band and filter in dict_obs_filters.keys():
                    plt.errorbar(dict_obs_filters[filter][:, 0], dict_obs_filters[filter][:, 1],
                                 yerr=dict_obs_filters[filter][:, 2],
                                 fmt='o', color='black', label='filter')

                    plt.plot(dict_model['time'], dict_model[filter], '-', label=filter)

        plt.minorticks_on()
        plt.xticks(fontsize=12)
        plt.yticks(fontsize=12)
        plt.tight_layout()
        plt.tick_params(axis='both', which='both', labelleft=True, labelright=False, tick1On=True, tick2On=True, labelsize=12, direction='in')  # labeltop
        plt.gca().invert_yaxis()
        plt.xscale("log")
        # plt.ylim(ymin=1e-5, ymax=2e-1)
        # plt.xlim(xmin=50, xmax=210)
        plt.ylabel("AB magnitude at 40 Mpc", fontsize=12)
        plt.xlabel("time [s]", fontsize=12)
        plt.legend(loc='best', numpoints=1)

        plt.savefig('{}.png'.format(fname), bbox_inches='tight', dpi=128)
        plt.close()

    def plot_for_several_bands(self, bands, fname):

        dict_obs_filters = self.load_obs_filters()

        dict_model = self.load_mkn_model()

        def plot_for_one_band(ax, band, color):

            mins, maxs = self.get_min_max_for_all_band(dict_model, band)
            plt.fill_between(dict_model['time'], maxs, mins, alpha=.5, color=color, label=band)

            for filter in dict_model.keys():
                if filter.split('_')[0] == band and filter in dict_obs_filters.keys():
                    ax.errorbar(dict_obs_filters[filter][:, 0], dict_obs_filters[filter][:, 1],
                                 yerr=dict_obs_filters[filter][:, 2],
                                 fmt='o', color=color)

        plt.figure()
        ax = plt.subplot(111)

        for band, color in zip(bands, Math.clrs(bands)):
            plot_for_one_band(ax, band, color)

        plt.minorticks_on()
        plt.xticks(fontsize=12)
        plt.yticks(fontsize=12)
        plt.tight_layout()
        plt.tick_params(axis='both', which='both', labelleft=True, labelright=False, tick1On=True, tick2On=True, labelsize=12, direction='in')  # labeltop
        plt.gca().invert_yaxis()
        plt.xscale("log")
        plt.ylim(ymin=25)#, ymax=2e-1)
        # plt.xlim(xmin=50, xmax=210)
        plt.ylabel(r"AB magnitude at 40 Mpc", fontsize=12)
        plt.xlabel(r"time [days]", fontsize=12)
        plt.legend(loc='best', numpoints=1)

        plt.savefig('./{}.png'.format(fname), bbox_inches='tight', dpi=128)
        plt.close()


    def print_latex_table_of_glob_pars(self):

        '''
        \begin{table}
            \footnotesize
                \begin{tabular}[t]{|p{3.2cm}|c|}
                    \hline
                    bla   & 1\\ \hline
                    blubb & 2 \\ \hline
                    bla   & 1\\ \hline
                    blubb & 2 \\ \hline
                    bla   & 1\\ \hline
                    blubb & 2 \\ \hline
                    xxx   & x \\ \hline
                \end{tabular}
                % \hfill
                \begin{tabular}[t]{|p{3.2cm}|c|}
                    \hline
                    bla&1\\ \hline
                    blubb&2 \\ \hline
                    bla&1\\ \hline
                    blubb&2 \\ \hline
                    bla&1\\ \hline
                    blubb&2 \\ \hline
                \end{tabular}
            \hfill
            \caption{99 most frequent hashtags in the data set.}
        \end{table}
        :return:
        '''

        glob_params, glob_vars, source_name = self.mkn_parameters_glob()

        print('\n')
        print('\\begin{table}[!ht]')
        print('\\footnotesize')

        # table of glob. parameters
        print('\\begin{tabular}[t]{ p{3.2cm} c }')
        print('\\hline')

        # printing rows
        for v_n, value in zip(glob_params.keys(), glob_params.values()):
            print(' {}  &  {} \\\\'.format(v_n.replace('_', '\\_'), value))
        print('\\hline')

        # table of glob. vars
        print('\\end{tabular}')
        print('\\begin{tabular}[t]{ p{3.2cm} c }')
        print('\\hline')

        for v_n, value in zip(glob_vars.keys(), glob_vars.values()):
            print(' {}  &  {} \\\\'.format(v_n.replace('_', '\\_'), value))

        print('\\hline')
        print('\\end{tabular}')
        print('\\\caption{Global parameters (left) and global variables (right)}')
        print('\\end{table}')

    def print_latex_table_of_ejecta_pars(self, dynamics=True, wind=True, secular=True):

        print('\n')
        print('\\begin{table}[!ht]')
        print('\\footnotesize')

        if dynamics:

            dyn_ej_pars, dyn_ej_vars = self.mkn_parameters_dynamics()

            print('\\begin{tabular}[t]{ p{3.cm} c }')
            print('Dynamic & \\\\')
            print('\\hline')

            for v_n, value in zip(dyn_ej_pars.keys(), dyn_ej_pars.values()):
                print(' {}  &  {} \\\\'.format(v_n.replace('_', '\\_'), value))
            print('\\hline')

            print('\\hline')

            for v_n, value in zip(dyn_ej_vars.keys(), dyn_ej_vars.values()):
                print(' {}  &  {} \\\\'.format(v_n.replace('_', '\\_'), value))
            print('\\hline')

            print('\\end{tabular}')

        if wind:

            wind_pars, wind_vars = self.mkn_parameters_wind()

            print('\\begin{tabular}[t]{ p{3.cm} c }')
            print('Wind & \\\\')
            print('\\hline')

            for v_n, value in zip(wind_pars.keys(), wind_pars.values()):
                print(' {}  &  {} \\\\'.format(v_n.replace('_', '\\_'), value))
            print('\\hline')

            print('\\hline')

            for v_n, value in zip(wind_vars.keys(), wind_vars.values()):
                print(' {}  &  {} \\\\'.format(v_n.replace('_', '\\_'), value))
            print('\\hline')

            print('\\end{tabular}')

        if secular:

            secular_pars, secular_vars = self.mkn_parameters_secular()

            print('\\begin{tabular}[t]{ p{3.cm} c }')
            print('Secualr & \\\\')
            print('\\hline')

            for v_n, value in zip(secular_pars.keys(), secular_pars.values()):
                print(' {}  &  {} \\\\'.format(v_n.replace('_', '\\_'), value))
            print('\\hline')

            print('\\hline')

            for v_n, value in zip(secular_vars.keys(), secular_vars.values()):
                print(' {}  &  {} \\\\'.format(v_n.replace('_', '\\_'), value))
            print('\\hline')

            print('\\end{tabular}')

        print('\\caption{Ejecta parameters}')
        print('\\end{table}')


class BEST_SIM_KILONOVA(SIM_KILONOVA):

    def __init__(self, sim, det):

        SIM_KILONOVA.__init__(self, sim, det)

    @staticmethod
    def convert_input(raw_input_arr):

        '''
        Requires in input as \/   and converts it into \/
        input = [                            result = [
            ['a', [1, 2, 3]],           ->      'a = 1', 'a = 2', 'a = 3',
            ['b', [0.001, 0.1]],        ->      'b = 0.001', 'b = 0.1',
            ['c', [2.1, 3.1, 4.1]]      ->      'c = 2.1', 'c = 3.1', 'c = 4.1'
        ]                                    ]
        '''

        res = []
        for i in range(len(raw_input_arr[:])):
            if len(raw_input_arr[i]) > 2: raise ValueError(
                'Incorrect entry in input array: {}'.format(raw_input_arr[i]))
            v_n = raw_input_arr[i][0]
            values = raw_input_arr[i][1]
            for val in values:
                res.append('{} = {}'.format(v_n, val))

        return res

    @staticmethod
    def get_permutations(input_arr, limit_order = None):

        '''
        This is a clear example what a shitty progrrrammer I am. However.
        It takes input in a following form:

        tst_conv = [
            'a = 1', 'a = 2', 'a = 3',
            'b = 0.001', 'b = 0.1',
            'c = 2.1', 'c = 3.1', 'c = 4.1'
        ]

        and computes all permutations among (a, b, c) removing repetitions and duplicates.

        '''

        res = []

        if limit_order != None:
            n_of_v_ns = limit_order
        else:
            n_of_v_ns = len(input_arr[:])

        if n_of_v_ns < 1: raise ValueError('n_of_v_ns must be bigger than 0. Given: {}'.format(n_of_v_ns))
        for interation in range(1, n_of_v_ns):

            arr = list(itertools.permutations(input_arr, interation))

            print('\t permutations:{} for it:{}'.format(len(arr), interation))

            # getting rid of subarrays with repeated variable
            arr2 = []
            for i in range(len(arr[:])):
                tmp = list(x.split(' = ')[0] for x in arr[i][:])
                if len(tmp) == len(set(tmp)):
                    arr2.append(arr[i])

            # getting rid of equivalent permutations (a=1,b=2) (b=2,a=1)
            arr3 = []
            arr_ex = []
            for i in range(len(arr2)):
                if i == 0:
                    arr3.append(arr2[i])
                else:
                    flag = True
                    for j in range(len(arr3)):
                        if sorted(arr2[i]) == sorted(list(arr3[j])):
                            flag = False
                    if flag:
                        arr3.append(arr2[i])
                    else:
                        arr_ex.append(arr2[i])
            # print(len(arr3))

            res += arr3

        print('total permutations: {}'.format(len(res)))

        return res

    def vary_one_paramter(self, v_n_place, v_n, values):

        # set default parameters
        glob_params_d, glob_vars_d, ejecta_params_d, shell_vars_d, source_name_d = self.mkn_parameters()

        logL = []
        for value in values:
            glob_params, glob_vars, ejecta_params, shell_vars = \
                self.modify_input(v_n_place, v_n, value,
                             glob_params_d, glob_vars_d, ejecta_params_d, shell_vars_d)

            logL_ = self.mkn_interface(glob_params, glob_vars, ejecta_params, shell_vars, source_name_d, False)
            logL = np.append(logL, logL_)

        print('\n')
        print(logL)

        return logL

    def vary_many_paramters(self, input, glob_params, glob_vars, ejecta_params, shell_vars):

        '''
        :param input: (a = 1, b = 2) or (a = 1, b = 2, c = 3)
        :return:
        '''

        glob_params_, glob_vars_, ejecta_params_, shell_vars_ = \
            copy.deepcopy(glob_params), copy.deepcopy(glob_vars), copy.deepcopy(ejecta_params), copy.deepcopy(shell_vars)

        for j in range(len(input)):
            if len(input[j].split(' = ')) != 2:
                raise ValueError('Entry of input_arr must be [place v_n = value] Given:{}'.format(input[j]))

            v_n_place = input[j].split(' = ')[0]
            if len(v_n_place.split()) != 2:
                raise ValueError(
                    'Entry of input_arr[i].split( = )(0) must be (place v_n) Given:{}'.format(v_n_place))

            place = v_n_place.split()[0]
            v_n = v_n_place.split()[1]
            value = np.float(input[j].split(' = ')[1]) # WARNING ! Only floats are supported so far.

            # print('... j:{} setting {} {} to {}'.format(j, place, v_n, value))

            # iteratively, by one, changing input parameters
            glob_params_, glob_vars_, ejecta_params_, shell_vars_ = self.modify_input(place, v_n, value,
                                                                                  glob_params_, glob_vars_,
                                                                                  ejecta_params_, shell_vars_)

        return glob_params_, glob_vars_, ejecta_params_, shell_vars_

    def vary_multiple_parameters(self, input_arr):
        '''
        Use following form:
        input_arr = [
            ('a = 1',), ('a = 2',), ... (a = 1, b = 2) ... (a = 1, b = 2, c = 3) ...
        ]
        where a is 'place v_n' (with ' ' between them)
        SOURCE_NAME is not changed here
        '''

        # get default parameters, -- same at the beginning of every iteration
        glob_params_d, glob_vars_d, ejecta_params_d, shell_vars_d, source_name_d = self.mkn_parameters()

        start_t = time.time()
        print("\n\t Starting Iterations \n"),

        logLs = []
        for i in range(len(input_arr)):
            print('\t Iteration: {}'.format(i))
            print('\t Computing for {}'.format(input_arr[i]))
            # get new input parameters where changes were made to all v_ns listed in input_arr[i]
            glob_params_, glob_vars_, ejecta_params_, shell_vars_ = \
                self.vary_many_paramters(input_arr[i],
                                         copy.deepcopy(glob_params_d),
                                         copy.deepcopy(glob_vars_d),
                                         copy.deepcopy(ejecta_params_d.copy()),
                                         copy.deepcopy(shell_vars_d.copy()))

            # run MKN model for final set of parameters (no saves are done)

            logL_ = self.mkn_interface(glob_params_, glob_vars_, ejecta_params_, shell_vars_, source_name_d, False)

            # record the modified parameters only and logL
            logLs = np.append(logLs, (-1.) * logL_)

        print("\t Iterations are done! (%.2f sec)" % (time.time() - start_t))

        i_best = np.argmin(logLs)
        logL = logLs[int(i_best)]
        best = input_arr[int(i_best)]

        print("\t Best model has logL: {}".format(logL))
        print('\t Best set of parameters: {}'.format(best))

        # computing the best model once again (from scratch) and saving the output
        glob_params, glob_vars, ejecta_params, shell_vars = \
            self.vary_many_paramters(input_arr[int(i_best)], glob_params_d, glob_vars_d, ejecta_params_d, shell_vars_d)
        self.mkn_interface(glob_params, glob_vars, ejecta_params, shell_vars, source_name_d, True)


        return best, logL



                # string = input_arr[i][j]
                # print(string)


            # if len(input_arr[i]) > 2:
            #     raise ValueError('Entry of input_arr must be [str(place v_n), [values]] Given:{}'.format(input_arr[i]))






        # if len(v_n_places) != len(v_ns):
        #     raise IOError('len(v_n_places) != len(v_ns) : {}!={}'.format(len(v_n_places), len(v_ns)))
        #
        # if len(v_n_places) != len(values[:][0]):
        #     raise IOError('len(v_n_places) != len(values[:][0]) : {}!={}'.format(len(v_n_places), len(values[:][0])))
        #
        # glob_params, glob_vars, ejecta_params, shell_vars, source_name = self.mkn_parameters()
        #
        # for v_n_place, v_n, value_arr in zip(v_n_places, v_ns, values):
        #     for value in value_arr:
        #         glob_params, glob_vars, ejecta_params, shell_vars = \
        #             self.modify_input(v_n_place, v_n, value,
        #                               glob_params, glob_vars, ejecta_params, shell_vars)
        #
        #         logL_ = self.mkn_interface(glob_params, glob_vars, ejecta_params, shell_vars, source_name,
        #                                    False)
        #         logL = np.append(logL, logL_)

    def master(self, raw_input, fname):
        '''
        Provide what parameters and where to iterate as:
            raw_input = [
                ['shell_vars[dynamics] high_lat_op', [1., 2., 3.]],
                ['shell_vars[wind]     xi_disk',     [0.001, 0.1]],
                ['shell_vars[secular]  central_op',  [2.1, 3.1, 4.1]]
            ]

        '''

        # modify input for easier use for permutations
        input = self.convert_input(raw_input)

        # compute permutations for all input parameters (no repetitions or equivalent inputs)
        permutations = self.get_permutations(input)

        # print(permutations)
        # exit(1)

        # do
        self.vary_multiple_parameters(permutations)

        # plot
        self.plot_for_several_bands(['g', 'z', 'Ks'], fname)


if __name__ == '__main__':

    make_model = COMPUTE_LIGHTCURVE("DD2_M13641364_M0_SR")
    make_model.compute_save_lightcurve(write_output=True)

    plot_model = PLOT_LIGHTCURVE("DD2_M13641364_M0_SR", '_0')
    plot_model.plot_for_several_bands(['g', 'z'], "tst_mkn")

# if __name__ == '__main__':
#
#
#     mkn_int = SIM_KILONOVA("DD2_M13641364_M0_SR", 0)
#     mkn_int.set_dynamics    = 'aniso'  # only central parameters to be changed [everything is UNIFORM]
#     mkn_int.set_wind        = 'aniso'  # 'iso'
#     mkn_int.set_seculiar    = 'aniso'  # 'iso'
#     mkn_int.set_use_NR_dyn_ej_distrib   = False
#     mkn_int.set_use_NR_dyn_ej_tot_m     = False
#     mkn_int.set_use_NR_disk_mass        = False
#     # mkn_int.print_latex_table_of_glob_pars(); exit(1)
#     # mkn_int.print_latex_table_of_ejecta_pars(); exit(1)
#
#     # mkn.plot_test_smooth_profile()
#     # mkn_int.mkn_interface()
#     # mkn_int.load_obs_filters()
#     mkn_int.save_one_model()
#     # mkn_int.load_mkn_model()
#     # mkn_int.plot_one_model('Ks_Gemini')
#     # mkn_int.plot_for_all_band('z')
#     mkn_int.plot_for_several_bands(['g', 'z', 'Ks'], 'tst_mkn'); exit(1)
#
#     mkn_set = BEST_SIM_KILONOVA("DD2_M13641364_M0_SR", 0)
#     mkn_set.set_dynamics    = 'aniso'  # 'iso'
#     mkn_set.set_wind        = 'aniso'  # 'iso'
#     mkn_set.set_seculiar    = 'aniso'  # 'iso'
#     mkn_set.set_use_NR_dyn_ej_distrib   = True
#     mkn_set.set_use_NR_dyn_ej_tot_m     = True
#     mkn_set.set_use_NR_disk_mass        = False
#
#     raw_input = [
#         ['glob_vars     m_disk',         [0.08, 0.12, 0.2]],
#         ['shell_vars[wind]     xi_disk', [0.01, 0.05, 0.1]],
#         ['shell_vars[wind]     central_vel', [0.04, 0.08]],
#         # ['shell_vars[wind]     low_lat_op',   [1.0, 5.0]],
#         ['shell_vars[wind]     high_lat_op', [.1, 1.]]
#     ]
#
#     # raw_input = [
#     #         ['glob_vars     m_disk', [0.08, 0.12, 0.2]],
#     #         ['shell_vars[secular]     xi_disk',      [0.1, 0.4, 0.6]],
#     #         ['shell_vars[secular]     central_vel',  [0.04, 0.08]],
#     #         ['shell_vars[secular]     central_op',   [5., 10.]],
#     #         # ['glob_vars     eps0',  [4e18, 8e18, 1e19, 2e19, 4e19, 8e19]]
#     #     ]
#
#     # raw_input = [
#     #         ['glob_vars     m_disk', [0.08, 0.12]],
#     #         ['shell_vars[wind]     xi_disk',       [0.1]],
#     #         ['shell_vars[wind]     central_vel',   [0.04]],
#     #         # ['shell_vars[wind]     low_lat_op',   [2.5, 6.5]],
#     #         ['shell_vars[wind]     high_lat_op',   [1.]],
#     #         ['shell_vars[secular]     xi_disk',    [0.6]],
#     #         ['shell_vars[secular]     central_vel',[0.04]],
#     #         ['shell_vars[secular]     central_op',   [10.]],
#     #     ]
#     # mkn_set.master(raw_input, 'tst12')