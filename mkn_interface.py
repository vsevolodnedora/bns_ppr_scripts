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
from glob import glob

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

from overall import MODELS
from plotting_nexus import PLOT_MANY_TASKS

path.append(Paths.mkn)
from mkn import MKN


class COMPUTE_LIGHTCURVE():

    def __init__(self, sim):

        self.table = MODELS()

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
        self.sim = sim
        self.output_fname = 'mkn_model.h5'
        self.path_to_outflow_dir = MakePath.outflow(sim, '_0')
        self.path_to_outflow_dir_psdyn = MakePath.outflow(sim, '_0_b_w')
        self.dyn_ejecta_profile_fpath = self.path_to_outflow_dir + Files.ejecta_profile
        self.psdyn_ejecta_profile_fpath=self.path_to_outflow_dir_psdyn + Files.ejecta_profile_bern


        self.set_use_dyn_NR = True
        self.set_use_bern_NR = True
        self.set_dyn_iso_aniso       = "aniso"
        self.set_psdyn_iso_aniso     = "aniso"
        self.set_wind_iso_aniso      = "aniso"
        self.set_secular_iso_aniso   = "aniso"

        self.glob_params    = {}
        self.glob_vars      = {}
        self.ejecta_params  = {}
        self.ejecta_vars    = {}
        self.source_name    = {}

        # self.set_glob_par_var_source(True, dyn_ejecta_profile_fpath,
        #                              True, psdyn_ejecta_profile_fpath)
        # self.set_dyn_par_var(self.set_dyn_iso_aniso)
        # self.set_psdyn_par_var(self.set_psdyn_iso_aniso)
        # self.set_wind_par_war(self.set_wind_iso_aniso)
        # self.set_secular_par_war(self.set_secular_iso_aniso)

        # self.compute_save_lightcurve(write_output=True)

    '''------------------------------------------'''

    def set_par_war(self):


        self.set_glob_par_var_source(self.set_use_dyn_NR, self.dyn_ejecta_profile_fpath,
                                     self.set_use_bern_NR, self.psdyn_ejecta_profile_fpath)
        self.set_dyn_par_var(self.set_dyn_iso_aniso)
        self.set_psdyn_par_var(self.set_psdyn_iso_aniso)
        self.set_wind_par_war(self.set_wind_iso_aniso)
        self.set_secular_par_war(self.set_secular_iso_aniso)

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
        # self.source_name = 'AT2017gfo view_angle=180/12.' # change the source properties

        self.glob_vars = {'m_disk':     self.table.get_value(self.sim, "Mdisk3D"), # 0.070, # LS220 | 0.20 - DD2
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
            self.ejecta_vars['dynamics'] = {'xi_disk':          None,
                                           'm_ej':              self.table.get_value(self.sim, "Mej"), # 0.00169045, # - LS220 | 0.00263355 - DD2
                                           'step_angle_mass':   None,
                                           'high_lat_flag':     None,
                                           'central_vel':       0.20, # changed from 0.33
                                           'high_lat_vel':      None,
                                           'low_lat_vel':       None,
                                           'step_angle_vel':    None,
                                           'central_op':        None,
                                           'high_lat_op':       1.,  # F:1
                                           'low_lat_op':        10., # F:30    # does not work for NR
                                           'step_angle_op': math.radians(45.)} # F:30
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
            self.ejecta_params['psdynamics'] = {'mass_dist': 'sin2', 'vel_dist': 'uniform', 'op_dist': 'uniform',
                                              'therm_model': 'BKWM', 'eps_ye_dep': True, 'v_law': 'power'}
            self.ejecta_vars['psdynamics'] = {'xi_disk':           None,
                                           'm_ej':              self.table.get_value(self.sim, "Mej_bern"), # WARNING! If set and NR it will overwrite
                                           'step_angle_mass':   None,
                                           'high_lat_flag':     None,
                                           'central_vel':       0.27,
                                           'high_lat_vel':      None,
                                           'low_lat_vel':       None,
                                           'step_angle_vel':    None,
                                           'central_op':        1.,
                                           'high_lat_op':       None,
                                           'low_lat_op':        None,
                                           'step_angle_op':     math.radians(45.),

                                           'override_m_ej':     False, # for manual import
                                           'kappa_low':         1.,    # for Import_NR data
                                           'kappa_high':        30.

                                        }
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
            self.ejecta_vars['wind'] = {
                         'xi_disk':         0.02, # F: 0.2    | 0.001 to 0.2
                         'm_ej':            None,
                         'step_angle_mass': math.radians(60.),
                         'high_lat_flag':   True,
                         'central_vel':     0.068, # F: 0.068 V:0.08
                         'high_lat_vel':    None,
                         'low_lat_vel':     None,
                         'step_angle_vel':  None,
                         'central_op':      None,
                         'high_lat_op':     0.5, # F
                         'low_lat_op':      5.0, # F
                         'step_angle_op':   math.radians(45.)} # F: 45 | might need # N:30
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
            self.ejecta_vars['secular'] = {
                            'xi_disk':          0.1, # F: 0.4
                            'm_ej':             None,
                            'step_angle_mass':  None,
                            'high_lat_flag':    None,
                            'central_vel':      0.04, # F: 0.04 V:0.06
                            'high_lat_vel':     None,
                            'low_lat_vel':      None,
                            'step_angle_vel':   None,
                            'central_op':       5.0, #
                            'low_lat_op':       None,
                            'high_lat_op':      None,
                            'step_angle_op':    None
            }
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

        if len(self.glob_params.keys()) == 0:
            raise ValueError("parameters are not set. Use 'set_par_war()' for that")

        if not os.path.isdir(Paths.ppr_sims+self.sim+'/res_mkn/'):
            print("making directory {}".format(Paths.ppr_sims+self.sim+'/res_mkn/'))
            os.mkdir(Paths.ppr_sims+self.sim+'/res_mkn/')

        print('I am initializing the model')
        # glob_params, glob_vars, ejecta_params, shell_vars, source_name = self.mkn_parameters()

        # go into the fold with all classes of mkn
        os.chdir(Paths.mkn)
        # from mkn import MKN

        # print(self.ejecta_vars['psdynamics']['m_ej'])
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
            copyfile('./mkn_model.h5', Paths.ppr_sims+self.sim+'/res_mkn/' + self.output_fname)

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


    def load_ej_profile_for_mkn(self, fpath):
        th, mass, vel, ye = np.loadtxt(fpath,
                                       unpack=True, usecols=(0, 1, 2, 3))
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
    def plot_test_smooth_profile(self, extension):

        fpath = MakePath.outflow(self.sim, extension) + Files.ejecta_profile

        # loading original profiles
        th, mass, vel, ye = self.load_ej_profile_for_mkn(fpath)

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
        # print(ye_smooth)
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
        idx = find_nearest_index(th, 90)
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
        ax1.set_ylabel('mass')
        ax1.tick_params(labelsize=12)
        ax1.legend(loc='best', numpoints=1)

        ax2.plot(th, ye, '-.', color='gray')
        ax2.plot(th, ye_smooth, '-', color='black', label='smoothed')
        ax2.set_ylabel('ye')
        ax2.tick_params(labelsize=12)
        ax2.legend(loc='best', numpoints=1)

        ax3.plot(th, vel, '-.', color='gray')
        ax3.plot(th, vel_smooth, '-', color='black', label='smoothed')
        ax3.set_ylabel('vel')
        ax3.tick_params(labelsize=12)
        ax3.legend(loc='best', numpoints=1)

        f.subplots_adjust(hspace=0)

        plt.setp([a.get_xticklabels() for a in f.axes[:-1]], visible=False)
        plt.xlabel('theta')
        plt.tick_params(axis='both', which='both', labelleft=True, labelright=False, tick1On=True, tick2On=True,
                        labelsize=12)  # labeltop

        plt.minorticks_on()
        plt.savefig(Paths.plots+'smoothed_profs.png', bbox_inches='tight', dpi=128)
        plt.close()


class LOAD_LIGHTCURVE():

    def __init__(self, sim):

        self.sim = sim
        self.default_fname = "mkn_model.h5"

        self.filter_fpath = Paths.mkn + Files.filt_at2017gfo

        # self.list_model_fnames = ["mkn_model.h5", "mkn_model1.h5", "mkn_model2.h5", "mkn_model3.h5", "mkn_model4.h5",
        #                           'mkn_model2_t50.h5', 'mkn_model2_t100.h5', 'mkn_model2_t150.h5', 'mkn_model2_t200.h5']

        fpaths = glob(Paths.ppr_sims + sim + "/res_mkn/mkn_model*")
        flist = []
        if len(fpaths) == 0: raise IOError("No mkn files found {}"
                                           .format(Paths.ppr_sims + sim + "/res_mkn/mkn_model*"))
        for file_ in fpaths:
            flist.append(file_.split('/')[-1])
        self.list_model_fnames = flist


        self.list_obs_filt_fnames = ["AT2017gfo.h5"]
        self.list_fnames = self.list_model_fnames + self.list_obs_filt_fnames

        self.list_attrs = ["psdynamics", "dynamics", "wind", "secular"]
        self.attrs_matrix = [[{}
                              for z in range(len(self.list_attrs))]
                              for y in range(len(self.list_fnames))]

        self.data_matrix = [{}
                             for y in range(len(self.list_fnames))]

        self.filters = {}

    def check_fname(self, fname=''):
        if not fname in self.list_fnames:
            raise NameError("fname: {} not in list_fnames:\n{}"
                            .format(fname, self.list_fnames))

    def check_attr(self, attr):
        if not attr in self.list_attrs:
            raise NameError("attr:{} not in list of attrs:{}"
                            .format(attr, self.list_attrs))

    def i_attr(self, attr):
        return int(self.list_attrs.index(attr))

    def get_attr(self, attr, fname=''):

        self.check_fname(fname)
        self.check_attr(attr)
        self.is_mkn_file_loaded(fname)

        return self.attrs_matrix[self.i_fname(fname)][self.i_attr(attr)]




    def i_fname(self, fname=''):
        return int(self.list_fnames.index(fname))

    def load_mkn_model(self, fname=''):

        if fname == '': fname = self.default_fname
        model_fpath = Paths.ppr_sims + self.sim + "/res_mkn/" + fname

        dict_model = {}

        model = h5py.File(model_fpath, "r")
        filters_model = []
        for it in model:
            if it in self.list_attrs:
                dic = {}
                for v_n in model[it].attrs:
                    dic[v_n] = model[it].attrs[v_n]
                self.attrs_matrix[self.i_fname(fname)][self.i_attr(it)] = dic
            else:
                filters_model.append(it)
                dict_model[str(it)] = np.array(model[it])

        # print('\t Following filters are available in mkn_model.h5: \n\t  {}'.format(filters_model))

        self.data_matrix[self.i_fname(fname)] = dict_model

    def load_obs_filters(self, fname=''):

        dict_obs_filters = {}

        obs_filters = h5py.File(self.filter_fpath, "r")

        filters_model = []
        for it in obs_filters:
            filters_model.append(it)
            arr = np.array(obs_filters[it])
            # print(arr.shape)
            dict_obs_filters[str(it)] = np.array(obs_filters[it])

        # print('\t Following filters are available in AT2017gfo.h5: \n\t  {}'.format(filters_model))

        self.filters = dict_obs_filters

    def is_filters_loaded(self, fname):

        if not bool(self.filters):
            self.load_obs_filters(fname)

    def is_mkn_file_loaded(self, fname=''):

        if not bool(self.data_matrix[self.i_fname(fname)]):
            self.load_mkn_model(fname)

    def get_mkn_model(self, fname=''):

        self.check_fname(fname)
        self.is_mkn_file_loaded(fname)

        return self.data_matrix[self.i_fname(fname)]

    def get_filters(self, fname):
        self.is_filters_loaded(fname)
        return self.filters


class EXTRACT_LIGHTCURVE(LOAD_LIGHTCURVE):

    def __init__(self, sim):
        LOAD_LIGHTCURVE.__init__(self, sim)
        self.list_bands = ["g", "z", "Ks"]
        self.do_extract_parameters = True

        self.model_params = [[{"psdynamics":{}, "dynamics":{}, "wind":{}, "secular":{}}
                                  for y in range(len(self.list_bands))]
                                 for z in range(len(self.list_fnames))]

        self.model_mag_matrix = [[ []
                                     for y in range(len(self.list_bands))]
                                     for z in range(len(self.list_fnames))]
        self.obs_mag_matrix = [[ []
                                     for y in range(len(self.list_bands))]
                                     for z in range(len(self.list_fnames))]


    def check_band(self, band):
        if not band in self.list_bands:
            raise NameError("band:{} not in tha band list:{}"
                            .format(band, self.list_bands))


    def i_band(self, band):
        return int(self.list_bands.index(band))

    def extract_lightcurve(self, band, fname=''):

        dict_model = self.get_mkn_model(fname)
        # arr = np.zeros(len(dict_model['time']))
        time_ = np.array(dict_model['time'])
        # if not band in dict_model.keys():
        #     raise NameError("band:{} is not in the loaded model:\n{}"
        #                     .format(band, dict_model.keys()))

        res = []
        for filter in dict_model.keys():
            if filter.split('_')[0] == band:
                # arr = np.vstack((arr, dict_model[filter]))
                res.append(np.vstack((time_, np.array(dict_model[filter]))).T)
        # times = arr[:, 0]
        # arr = np.delete(arr, 0, 0)

        if len(res) == 0:
            raise NameError("band:{} is not found in the loaded model:\n{}"
                                .format(band, dict_model.keys()))

        self.model_mag_matrix[self.i_fname(fname)][self.i_band(band)] = res

        # ''' extract parameters '''
        # if self.do_extract_parameters:
        #     if "psdynamics" in


    def extract_obs_data(self, band, fname):

        dict_obs_filters = self.get_filters(fname)
        # dict_model = self.get_mkn_model(fname)

        sub_bands = []
        for filter in dict_obs_filters.keys():
            if filter.split('_')[0] == band:# and filter in dict_obs_filters.keys():
                sub_bands.append(dict_obs_filters[filter])


        if len(sub_bands) == 0:
            raise NameError("band:{} is not found in the loaded obs filters:\n{}"
                            .format(band, dict_obs_filters.keys()))

        self.obs_mag_matrix[self.i_fname(fname)][self.i_band(band)] = sub_bands

    def is_extracted(self, band, fname=''):

        data = self.model_mag_matrix[self.i_fname(fname)][self.i_band(band)]

        if len(data)  == 0 and fname in self.list_model_fnames:
            self.extract_lightcurve(band, fname)

        if len(data) == 0 and fname in self.list_obs_filt_fnames:
            self.extract_obs_data(band, fname)



    def get_obs_data(self, band, fname="AT2017gfo.h5"):
        """

        :param band:
        :param fname:
        :return:     list of [:times:, :magnitudes:, :errors:] 3D array for every subband in band
        """
        self.check_fname(fname)
        self.check_band(band)

        self.is_extracted(band, fname)


        return self.obs_mag_matrix[self.i_fname(fname)][self.i_band(band)]


    def get_model(self, band, fname="mkn_model.h5"):
        self.check_band(band)
        self.check_fname(fname)

        self.is_extracted(band, fname)

        return self.model_mag_matrix[self.i_fname(fname)][self.i_band(band)]

    def get_model_min_max(self, band, fname="mkn_model.h5"):

        band_list = self.get_model(band, fname)

        maxs = []
        mins = []
        times = []
        mags = []
        for i_band, band in enumerate(band_list):
            times = band[:, 0]
            mags = np.append(mags, band[:, 1])

        mags = np.reshape(mags, (len(band_list), len(times)))

        for i in range(len(times)):
            maxs.append(mags[:,i].max())
            mins.append(mags[:,i].min())

        return times, maxs, mins
        #
        #
        #
        #
        # time_ = arr[0, :]
        # # arr = np.delete(arr, 0, 0)
        #
        # print(arr.shape)
        # print(arr)
        #
        # maxs = []
        # for i in range(len(arr[0, :])):
        #     maxs = np.append(maxs, arr[1:,i].max())
        #
        # mins = []
        # for i in range(len(arr[0, :])):
        #     mins = np.append(mins, arr[1:,i].min())
        #
        # if len(time_) != len(mins):
        #     raise ValueError("len(time_) {} != {} len(mins)"
        #                      .format(len(time_) ,len(mins)))
        #
        #
        #
        # return time_, mins, maxs

    def get_model_median(self, band, fname="mkn_model.h5"):

        m_times, m_maxs, m_mins = self.get_model_min_max(band, fname)

        m_times = np.array(m_times)
        m_maxs = np.array(m_maxs)
        m_mins = np.array(m_mins)

        return m_times, m_mins + ((m_maxs - m_mins) / 2)

    def get_mismatch(self, band, fname="mkn_model.h5"):

        from scipy import interpolate

        m_times, m_maxs, m_mins = self.get_model_min_max(band, fname)
        obs_data = self.get_obs_data(band)


        all_obs_times = []
        all_obs_maxs = []
        all_obs_mins = []
        for sumbband in obs_data:

            obs_time = sumbband[:, 0]
            obs_maxs = sumbband[:, 1] + sumbband[:, 2] # data + error bar
            obs_mins = sumbband[:, 1] - sumbband[:, 2]  # data - error bar

            all_obs_times = np.append(all_obs_times, obs_time)
            all_obs_maxs = np.append(all_obs_maxs, obs_maxs)
            all_obs_mins = np.append(all_obs_mins, obs_mins)

        all_obs_times, all_obs_maxs, all_obs_mins = \
            x_y_z_sort(all_obs_times, all_obs_maxs, all_obs_mins)

        # interpolate for observationa times

        int_m_times = all_obs_times
        if all_obs_times.max() > m_times.max():
            int_m_times = all_obs_times[all_obs_times < m_times.max()]
        int_m_maxs = interpolate.interp1d(m_times, m_maxs, kind='linear')(int_m_times)
        int_m_mins = interpolate.interp1d(m_times, m_mins, kind='linear')(int_m_times)

        min_mismatch = []
        max_mismatch = []

        for i in range(len(int_m_times)):
            m_max = int_m_maxs[i]
            m_min = int_m_mins[i]
            o_min = all_obs_mins[i]
            o_max = all_obs_maxs[i]

            if o_max > m_max and o_min < m_min:
                min_mismatch = np.append(min_mismatch, 0)
            elif o_min <= m_max and o_min >= m_min:
                min_mismatch = np.append(min_mismatch, 0)
            elif o_max <= m_max and o_max >= m_min:
                min_mismatch = np.append(min_mismatch, 0)
            elif (o_min > m_max):
                min_mismatch = np.append(min_mismatch, o_min - m_max)
            elif (o_max < m_min):
                min_mismatch = np.append(min_mismatch, o_max - m_min)
            else:
                raise ValueError("mismatched failed m_max:{} m_min:{} o_max:{} o_min:{}"
                                 .format(m_max, m_min, o_max, o_min))
            #
            #
            # min_mismatch = np.append(min_mismatch, min([o_min - m_min, o_min - m_max,
            #                                             m_max - m_min, o_max - m_max]))
            # max_mismatch = np.append(max_mismatch, max([o_min - m_min, o_min - m_max,
            #                                             m_max - m_min, o_max - m_max]))

        # print(min_mismatch)

        return int_m_times, min_mismatch, max_mismatch


        # print(obs_data)

    def get_model_peak(self, band, fname="mkn_model.h5"):
        t, mag = self.get_model_median(band, fname)
        idx = find_nearest_index(mag, mag.min())
        return t[idx], mag[idx]

    def get_obs_peak(self, band, fname = "AT2017gfo.h5"):

        from scipy import interpolate

        obs_data = self.get_obs_data(band, fname)
        obs_times = []
        obs_mags = []

        for sumbband in obs_data:
            obs_times = np.append(obs_times, sumbband[:, 0])
            obs_mags = np.append(obs_mags, sumbband[:, 1])

        obs_times, obs_mags = x_y_z_sort(obs_times, obs_mags)

        int_obs_times = np.mgrid[obs_times[0]:obs_times[-2]:100j]

        assert len(int_obs_times) == 100

        assert obs_times.min() <= int_obs_times.min()
        assert obs_times.max() >= int_obs_times.max()

        int_obs_mags = interpolate.interp1d(obs_times, obs_mags, kind='linear')(int_obs_times)
        print(int_obs_mags)
        idx = find_nearest_index(int_obs_mags, int_obs_mags.min())
        return int_obs_times[idx], int_obs_mags[idx]


        # obs_data = self.get_obs_data(band)

        # all_obs_times = []
        # all_obs_maxs = []
        # all_obs_mins = []
        # for sumbband in obs_data:
        #     obs_time = sumbband[:, 0]
        #     obs_maxs = sumbband[:, 1] + sumbband[:, 2]  # data + error bar
        #     obs_mins = sumbband[:, 1] - sumbband[:, 2]  # data - error bar
        #
        #     all_obs_times = np.append(all_obs_times, obs_time)
        #     all_obs_maxs = np.append(all_obs_maxs, obs_maxs)
        #     all_obs_mins = np.append(all_obs_mins, obs_mins)
        #
        # all_obs_times, all_obs_maxs, all_obs_mins = \
        #     x_y_z_sort(all_obs_times, all_obs_maxs, all_obs_mins)
        #
        #
        # #
        # # print(m_times)
        # # print(all_obs_times)
        # #
        # # mask1 = (m_times < all_obs_times.max())
        # # mask2 = (m_times > all_obs_times.min())
        # # print(mask1)
        # # print(mask2)
        # # int_obs_times = m_times[mask1 & mask2]
        # int_obs_times = np.mgrid[all_obs_times.min():all_obs_times.max():100j]
        # print(np.log10(all_obs_times))
        # int_obs_maxs = interpolate.interp1d(all_obs_times, all_obs_maxs, kind='linear')(int_obs_times)
        # int_obs_mins = interpolate.interp1d(all_obs_times, all_obs_mins, kind='linear')(int_obs_times)
        #
        # idx = find_nearest_index(int_obs_maxs, int_obs_maxs.min())
        #
        # return int_obs_times[idx], int_obs_maxs[idx], int_obs_mins[idx]

        #
        #
        #
        #
        # # interpolate for observationa times
        #
        # int_m_times = all_obs_times
        # if all_obs_times.max() > m_times.max():
        #     int_m_times = all_obs_times[all_obs_times < m_times.max()]
        # int_m_maxs = interpolate.interp1d(all_obs_times, all_obs_maxs, kind='cubic')(int_m_times)
        # int_m_mins = interpolate.interp1d(all_obs_times, all_obs_mins, kind='cubic')(int_m_times)
        #
        # min_mismatch = []
        # max_mismatch = []
        #
        # for i in range(len(int_m_times)):
        #     m_max = int_m_maxs[i]
        #     m_min = int_m_mins[i]
        #     o_min = all_obs_mins[i]
        #     o_max = all_obs_maxs[i]
        #
        #     if o_max > m_max and o_min < m_min:
        #         min_mismatch = np.append(min_mismatch, 0)
        #     elif o_min <= m_max and o_min >= m_min:
        #         min_mismatch = np.append(min_mismatch, 0)
        #     elif o_max <= m_max and o_max >= m_min:
        #         min_mismatch = np.append(min_mismatch, 0)
        #     elif (o_min > m_max):
        #         min_mismatch = np.append(min_mismatch, o_min - m_max)
        #     elif (o_max < m_min):
        #         min_mismatch = np.append(min_mismatch, o_max - m_min)
        #     else:
        #         raise ValueError("mismatched failed m_max:{} m_min:{} o_max:{} o_min:{}"
        #                          .format(m_max, m_min, o_max, o_min))
        #     #
        #     #
        #     # min_mismatch = np.append(min_mismatch, min([o_min - m_min, o_min - m_max,
        #     #                                             m_max - m_min, o_max - m_max]))
        #     # max_mismatch = np.append(max_mismatch, max([o_min - m_min, o_min - m_max,
        #     #                                             m_max - m_min, o_max - m_max]))
        #
        # # print(min_mismatch)
        #
        # return int_m_times, min_mismatch, max_mismatch

    def get_obs_peak_duration(self, band, limit=1.,  fname = "AT2017gfo.h5"):

        from scipy import interpolate

        obs_data = self.get_obs_data(band,  fname)
        obs_times = []
        obs_mags = []

        for sumbband in obs_data:
            obs_times = np.append(obs_times, sumbband[:, 0])
            obs_mags = np.append(obs_mags, sumbband[:, 1])

        obs_times, obs_mags = x_y_z_sort(obs_times, obs_mags)

        int_obs_times = np.mgrid[obs_times[0]:obs_times[-2]:100j]

        assert len(int_obs_times) == 100

        assert obs_times.min() <= int_obs_times.min()
        assert obs_times.max() >= int_obs_times.max()

        int_obs_mags = interpolate.interp1d(obs_times, obs_mags, kind='linear')(int_obs_times)
        print(int_obs_mags)
        idx = find_nearest_index(int_obs_mags, int_obs_mags.min())

        peaktime = int_obs_times[idx]
        peakmag = int_obs_mags[idx]

        mask = (obs_times >= peaktime) & (obs_mags < peakmag + limit)
        assert len(mask) > 1
        post_peak_times = obs_times[mask]
        post_peak_mags = obs_mags[mask]

        assert len(post_peak_times) > 1


        return post_peak_times[-1] - peaktime, post_peak_mags[-1]

    def get_model_peak_duration(self, band, fname="mkn_model.h5", limit = 1.):

        t, mag = self.get_model_median(band, fname)
        idx = find_nearest_index(mag, mag.min())
        tpeak = t[idx]
        magpeak = mag[idx]

        mask = (t >= tpeak) & (mag < magpeak + limit)
        assert len(mask) > 1
        post_peak_times = t[mask]
        post_peak_mags = mag[mask]

        assert len(post_peak_times) > 1

        return post_peak_times[-1] - tpeak, post_peak_mags[-1]

class COMBINE_LIGHTCURVES(EXTRACT_LIGHTCURVE):

    def __init__(self, sim):

        EXTRACT_LIGHTCURVE.__init__(self, sim)


    def get_model_peaks(self, band, files_name_gen=r"mkn_model2_m*.h5"):

        files = glob(Paths.ppr_sims+self.sim+'/res_mkn/'+files_name_gen)
        # print(files)

        tpeaks = []
        mpeaks = []
        attrs = []

        for file_ in files:
            tpeak, mpeak = self.get_model_peak(band, file_.split('/')[-1])
            attr = self.get_attr("psdynamics", file_.split('/')[-1])["m_ej"]

            tpeaks = np.append(tpeaks, tpeak)
            mpeaks = np.append(mpeaks, mpeak)
            attrs = np.append(attrs, attr)

        attrs, tpeaks, mpeaks = x_y_z_sort(attrs, tpeaks, mpeaks)

        return attrs, tpeaks, mpeaks

    def get_model_peak_durations(self, band, files_name_gen=r"mkn_model2_m*.h5"):

        files = glob(Paths.ppr_sims + self.sim + '/res_mkn/' + files_name_gen)
        # print(files)

        tdurs = []
        attrs = []

        for file_ in files:
            tdur, _ = self.get_model_peak_duration(band, file_.split('/')[-1], limit=1.)
            attr = self.get_attr("psdynamics", file_.split('/')[-1])["m_ej"]

            tdurs = np.append(tdurs, tdur)
            attrs = np.append(attrs, attr)

        attrs, tdurs = x_y_z_sort(attrs, tdurs)

        return attrs, tdurs

    def get_table(self, band='g', files_name_gen=r"mkn_model2_m*.h5"):

        files = glob(Paths.ppr_sims+self.sim+'/res_mkn/'+files_name_gen)
        # print(files)

        t_arr = []
        mag_arr = []
        attr_arr = []


        def get_atr(file_):
            return self.get_attr("psdynamics", file_.split('/')[-1])["m_ej"]

        files = sorted(files, key=get_atr)


        for file_ in files:

            m_time, m_mag = self.get_model_median(band, file_.split('/')[-1])
            attr = self.get_attr("psdynamics", file_.split('/')[-1])["m_ej"]

            print('\t processing {} atr: {}'.format(file_.split('/')[-1], attr))

            t_arr = m_time
            mag_arr = np.append(mag_arr, m_mag)
            attr_arr.append(attr)

        mag_table = np.reshape(mag_arr, (len(attr_arr), len(t_arr)))

        t_grid, attr_grid = np.meshgrid(t_arr, attr_arr)

        return  t_grid, attr_grid, mag_table

        #
        # dfile = h5py.File(files[0], "r")
        #
        #
        #
        #
        #
        # ejecta_type = "psdynamics"
        # print(dfile[ejecta_type])
        #
        # # dfile[ejecta_type].attrs[""]
        #
        # v_ns = []
        # values = []
        # for v_n in dfile[ejecta_type].attrs:
        #     v_ns.append(v_n)
        #     values.append(dfile[ejecta_type].attrs[v_n])
        #
        # print(v_ns, values)
        #
        # pass




''' ---------------------------------------------- '''


class LOAD_LIGHTCURVES():

    def __init__(self, sim):

        self.sim = sim
        self.default_fname = "mkn_model.h5"

        self.filter_fpath = Paths.mkn + Files.filt_at2017gfo

        self.list_fnames = ['', "mkn_model.h5"]
        self.list_of_sims = os.listdir(Paths.ppr_sims)

        self.data_matrix = [[{}
                             for x in range(len(self.list_fnames))]
                             for y in range(len(self.list_of_sims))]

        self.filters = {}

    def check_sim(self, sim):
        if not sim in self.list_of_sims:
            raise NameError("sim: {} not in list_of_sims:\n{}"
                            .format(sim, self.list_of_sims))

    def check_fname(self, fname):
        if not fname in self.list_fnames:
            raise NameError("fname: {} not in list_fnames:\n{}"
                            .format(fname, self.list_fnames))

    def i_sim(self, sim):
        return int(self.list_of_sims.index(sim))

    def i_fname(self, fname):
        return int(self.list_fnames.index(fname))


    def load_mkn_model(self, sim, fname):

        if fname == '': fname = self.default_fname
        model_fpath = Paths.ppr_sims + sim + "/res_mkn/" + fname

        dict_model = {}

        model = h5py.File(model_fpath, "r")
        filters_model = []
        for it in model:
            filters_model.append(it)
            dict_model[str(it)] = np.array(model[it])

        # print('\t Following filters are available in mkn_model.h5: \n\t  {}'.format(filters_model))

        self.data_matrix[self.i_sim(sim)][self.i_fname(fname)] = dict_model

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

        self.filters = dict_obs_filters

    def is_filters_loaded(self):

        if not bool(self.filters):
            self.load_obs_filters()

    def is_mkn_file_loaded(self, sim, fname):

        if not bool(self.data_matrix[self.i_sim(sim)][self.i_fname(fname)]):
            self.load_mkn_model(sim, fname)

    def get_mkn_model(self, sim, fname):

        self.check_sim(sim)
        self.check_fname(fname)
        self.is_mkn_file_loaded(sim, fname)

        return self.data_matrix[self.i_sim(sim)][self.i_fname(fname)]

    def get_filters(self):
        self.is_filters_loaded()
        return self.filters

class PLOT_LIGHTCURVE():

    def __init__(self):
        pass
        # self.model_fpath = Paths.ppr_sims+self.sim+'/res_mkn/' + 'mkn_model.h5'
        #
        # self.filter_fpath = Paths.mkn + Files.filt_at2017gfo


    # def load_mkn_model(self):
    #find SFHo_13641364_M0_LK_HR -exec touch '{}' \;


    #     dict_model = {}
    #
    #     model = h5py.File(self.model_fpath, "r")
    #     filters_model = []
    #     for it in model:
    #         filters_model.append(it)
    #         dict_model[str(it)] = np.array(model[it])
    #
    #     # print('\t Following filters are available in mkn_model.h5: \n\t  {}'.format(filters_model))
    #
    #     return dict_model
    #
    # def load_obs_filters(self):
    #
    #     dict_obs_filters = {}
    #
    #     obs_filters = h5py.File(self.filter_fpath, "r")
    #
    #     filters_model = []
    #     for it in obs_filters:
    #         filters_model.append(it)
    #         arr = np.array(obs_filters[it])
    #         # print(arr.shape)
    #         dict_obs_filters[str(it)] = np.array(obs_filters[it])
    #
    #     # print('\t Following filters are available in AT2017gfo.h5: \n\t  {}'.format(filters_model))
    #
    #     return dict_obs_filters
    #
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
    #
    # def plot_for_all_band(self, band, plot_as_band = True, fname = 'tst10'):
    #
    #     dict_obs_filters = self.load_obs_filters()
    #
    #     dict_model = self.load_mkn_model()
    #
    #     # plot models as a continouse band of gray color
    #     if plot_as_band == True:
    #         mins, maxs = self.get_min_max_for_all_band(dict_model, band)
    #         plt.fill_between(dict_model['time'], maxs, mins, alpha=.5, color='gray')
    #
    #         for filter in dict_model.keys():
    #             if filter.split('_')[0] == band and filter in dict_obs_filters.keys():
    #                 plt.errorbar(dict_obs_filters[filter][:, 0], dict_obs_filters[filter][:, 1],
    #                              yerr=dict_obs_filters[filter][:, 2],
    #                              fmt='o', color='black')
    #     # plot models as separate lines with different colors and labels
    #     else:
    #         for filter in dict_model.keys():
    #             if filter.split('_')[0] == band and filter in dict_obs_filters.keys():
    #                 plt.errorbar(dict_obs_filters[filter][:, 0], dict_obs_filters[filter][:, 1],
    #                              yerr=dict_obs_filters[filter][:, 2],
    #                              fmt='o', color='black', label='filter')
    #
    #                 plt.plot(dict_model['time'], dict_model[filter], '-', label=filter)
    #
    #     plt.minorticks_on()
    #     plt.xticks(fontsize=12)
    #     plt.yticks(fontsize=12)
    #     plt.tight_layout()
    #     plt.tick_params(axis='both', which='both', labelleft=True,
    #                     labelright=False, tick1On=True, tick2On=True,
    #                     labelsize=12, direction='in')  # labeltop
    #     plt.gca().invert_yaxis()
    #     plt.xscale("log")
    #     # plt.ylim(ymin=1e-5, ymax=2e-1)
    #     # plt.xlim(xmin=50, xmax=210)
    #     plt.ylabel("AB magnitude at 40 Mpc", fontsize=12)
    #     plt.xlabel("time [s]", fontsize=12)
    #     plt.legend(loc='best', numpoints=1)
    #
    #     plt.savefig('{}.png'.format(fname), bbox_inches='tight', dpi=128)
    #     plt.close()
    #
    # def plot_for_several_bands(self, bands, fname):
    #
    #     dict_obs_filters = self.load_obs_filters()
    #
    #     dict_model = self.load_mkn_model()
    #
    #     def plot_for_one_band(ax, band, color):
    #
    #         mins, maxs = self.get_min_max_for_all_band(dict_model, band)
    #         plt.fill_between(dict_model['time'], maxs, mins, alpha=.5, color=color)
    #
    #         labeled_bands = []
    #         for filter in dict_model.keys():
    #             # print("filter:{} ")
    #             if filter.split('_')[0] == band and filter in dict_obs_filters.keys():
    #                 if not band in labeled_bands:
    #                     ax.errorbar(dict_obs_filters[filter][:, 0], dict_obs_filters[filter][:, 1],
    #                                 yerr=dict_obs_filters[filter][:, 2],
    #                                 fmt='o', color=color, label=band)
    #                 else:
    #                     ax.errorbar(dict_obs_filters[filter][:, 0], dict_obs_filters[filter][:, 1],
    #                                 yerr=dict_obs_filters[filter][:, 2],
    #                                 fmt='o', color=color)
    #                 labeled_bands.append(band)
    #
    #
    #     plt.figure(figsize=(4.5, 3.6))
    #     ax = plt.subplot(111)
    #
    #     for band in bands:
    #         # print('band')
    #         color = color_for_mkn_band(band)
    #         plot_for_one_band(ax, band, color)
    #
    #     plt.minorticks_on()
    #     plt.xticks(fontsize=12)
    #     plt.yticks(fontsize=12)
    #     plt.tight_layout()
    #     plt.tick_params(axis='both', which='both', labelleft=True, labelright=False,
    #                     tick1On=True, tick2On=True, labelsize=12, direction='in')  # labeltop
    #     plt.gca().invert_yaxis()
    #     plt.xscale("log")
    #     plt.xlim(xmin=0.2)
    #     plt.ylim(ymin=25)#, ymax=2e-1)
    #     # plt.xlim(xmin=50, xmax=210)
    #     plt.ylabel(r"AB magnitude at 40 Mpc", fontsize=12)
    #     plt.xlabel(r"time [days]", fontsize=12)
    #     plt.legend(loc='best', numpoints=1, fontsize=12)
    #
    #
    #     plt.savefig(Paths.ppr_sims+self.sim+'/res_mkn/' + self.gen_set['figname'], bbox_inches='tight', dpi=128)
    #     plt.close()


    def plot_obs_data(self, ax, dic, data):

        dict_obs_filters = data.get_filters()
        dict_model = data.get_mkn_model(dic['sim'], dic['fname'])

        if dic['obscolor'] == '' or dic['obscolor'] == None:
            color = color_for_mkn_band(dic["band"])
        else:
            color = dic['obscolor']

        labeled_bands = []
        for filter in dict_model.keys():
            # print("filter:{} ")dic["obscolor"]
            if filter.split('_')[0] == dic["band"] and filter in dict_obs_filters.keys():
                if not dic["band"] in labeled_bands:
                    ax.errorbar(dict_obs_filters[filter][:, 0], dict_obs_filters[filter][:, 1],
                                yerr=dict_obs_filters[filter][:, 2],
                                fmt=dic['marker'], color=color, label=dic["band"])
                else:
                    ax.errorbar(dict_obs_filters[filter][:, 0], dict_obs_filters[filter][:, 1],
                                yerr=dict_obs_filters[filter][:, 2],
                                fmt=dic['marker'], color=color)
                labeled_bands.append(dic["band"])

    def plot_simple_lightcurve(self, ax, dic, data):

        # dic = {
        #     'sim': "DD2_M13641364_M0_LK_SR", 'fname': 'mkn_model.h5',
        #     'obs': True,
        #     'band': 'g', 'fillcolor': 'gray', 'obscolor': 'black'
        # }

        if dic['fillcolor'] == '' or dic['fillcolor'] == None:
            color = color_for_mkn_band(dic["band"])
        else:
            color = dic['fillcolor']

        dict_model = data.get_mkn_model(dic['sim'], dic['fname'])
        mins, maxs = self.get_min_max_for_all_band(dict_model, dic["band"])

        ax.fill_between(dict_model['time'], maxs, mins, alpha=.5, color=color)

    def plot_main(self, ax, dic, data):
        if dic["name"] == '-':
            self.plot_simple_lightcurve(ax, dic, data)
            if dic['obs']:
                self.plot_obs_data(ax, dic, data)

# class PLOT_MANY_TASKS:
#
#     def __init__(self, o_data, o_base_plots):
#
#
#         self.data = o_data
#         self.base_plots = o_base_plots
#
#         self.gen_set = {
#             "figdir": Paths.plots + 'overall/',
#             "figname": "inv_ang_mom_flux.png",
#             # "figsize": (13.5, 3.5), # <->, |
#             "figsize": (3.8, 3.5),  # <->, |
#             "type": "cartesian",
#             "subplots_adjust_h": 0.2,
#             "subplots_adjust_w": 0.3,
#             "fancy_ticks": True,
#             "minorticks_on": True
#         }
#
#         self.set_plot_dics = []
#
#
#     def set_ncols_nrows(self):
#
#         tmp_rows = []
#         tmp_cols = []
#
#         for dic in self.set_plot_dics:
#             tmp_cols.append(dic['position'][1])
#             tmp_rows.append(dic['position'][0])
#
#         max_row = max(tmp_rows)
#         max_col = max(tmp_cols)
#
#         for row in range(1, max_row):
#             if not row in tmp_rows:
#                 raise NameError("Please set vertical plot position in a subsequent order: 1,2,3... not 1,3...")
#
#         for col in range(1, max_col):
#             if not col in tmp_cols:
#                 raise NameError("Please set horizontal plot position in a subsequent order: 1,2,3... not 1,3...")
#
#         print("Set {} rows {} columns (total {}) of plots".format(max_row, max_col, len(self.set_plot_dics)))
#
#         return int(max_row), int(max_col)
#
#     def set_plot_dics_matrix(self):
#
#         plot_dic_matrix = [[0
#                              for x in range(self.n_rows)]
#                              for y in range(self.n_cols)]
#
#         # get a matrix of dictionaries describing plots (for ease of representation)
#         for dic in self.set_plot_dics:
#             col, row = int(dic['position'][1]-1), int(dic['position'][0]-1) # -1 as position starts with 1
#             # print(col, row)
#             for n_row in range(self.n_rows):
#                 for n_col in range(self.n_cols):
#                     if int(col) == int(n_col) and int(row) == int(n_row):
#                         plot_dic_matrix[n_col][n_row] = dic
#                         # print('adding {} {}'.format(col, row))
#
#             if isinstance(plot_dic_matrix[col][row], int):
#                 raise ValueError("Dictionary to found for n_row {} n_col {} in "
#                                  "creating matrix of dictionaries".format(col, row))
#
#         return plot_dic_matrix
#
#     def set_plot_matrix(self):
#
#         fig = plt.figure(figsize=self.gen_set['figsize'])  # (<->; v)
#
#         if self.gen_set['type'] == 'cartesian':
#             # initializing the matrix with dummy axis objects
#             sbplot_matrix = [[fig.add_subplot(self.n_rows, self.n_cols, 1)
#                               for x in range(self.n_rows)]
#                              for y in range(self.n_cols)]
#
#             i = 1
#             for n_row in range(self.n_rows):
#                 for n_col in range(self.n_cols):
#
#                     if n_col == 0 and n_row == 0:
#                         sbplot_matrix[n_col][n_row] = fig.add_subplot(self.n_rows, self.n_cols, i)
#                     elif n_col == 0 and n_row > 0:
#                         if self.gen_set['sharex']:
#                             sbplot_matrix[n_col][n_row] = fig.add_subplot(self.n_rows, self.n_cols, i,
#                                                                           sharex=sbplot_matrix[n_col][0])
#                         else:
#                             sbplot_matrix[n_col][n_row] = fig.add_subplot(self.n_rows, self.n_cols, i)
#                     elif n_col > 0 and n_row == 0:
#                         if self.gen_set['sharey']:
#                             sbplot_matrix[n_col][n_row] = fig.add_subplot(self.n_rows, self.n_cols, i,
#                                                                           sharey=sbplot_matrix[0][n_row])
#                         else:
#                             sbplot_matrix[n_col][n_row] = fig.add_subplot(self.n_rows, self.n_cols, i)
#                     else:
#                         if self.gen_set['sharex'] and not self.gen_set['sharey']:
#                             sbplot_matrix[n_col][n_row] = fig.add_subplot(self.n_rows, self.n_cols, i,
#                                                                           sharex=sbplot_matrix[n_col][0])
#                         elif not self.gen_set['sharex'] and self.gen_set['sharey']:
#                             sbplot_matrix[n_col][n_row] = fig.add_subplot(self.n_rows, self.n_cols, i,
#                                                                           sharey=sbplot_matrix[0][n_row])
#                         else:
#                             sbplot_matrix[n_col][n_row] = fig.add_subplot(self.n_rows, self.n_cols, i,
#                                                                           sharex=sbplot_matrix[n_col][0],
#                                                                           sharey=sbplot_matrix[0][n_row])
#
#                         # sbplot_matrix[n_col][n_row].axes.get_yaxis().set_visible(False)
#                     # sbplot_matrix[n_col][n_row] = fig.add_subplot(n_rows, n_cols, i)
#                     i += 1
#
#         elif self.gen_set['type'] == 'polar':
#             # initializing the matrix with dummy axis objects
#             sbplot_matrix = [[fig.add_subplot(self.n_rows, self.n_cols, 1, projection='polar')
#                                   for x in range(self.n_rows)]
#                                   for y in range(self.n_cols)]
#
#             i = 1
#             for n_row in range(self.n_rows):
#                 for n_col in range(self.n_cols):
#
#                     if n_col == 0 and n_row == 0:
#                         sbplot_matrix[n_col][n_row] = fig.add_subplot(self.n_rows, self.n_cols, i, projection='polar')
#                     elif n_col == 0 and n_row > 0:
#                         sbplot_matrix[n_col][n_row] = fig.add_subplot(self.n_rows, self.n_cols, i, projection='polar')
#                                                                       # sharex=self.sbplot_matrix[n_col][0])
#                     elif n_col > 0 and n_row == 0:
#                         sbplot_matrix[n_col][n_row] = fig.add_subplot(self.n_rows, self.n_cols, i, projection='polar')
#                                                                       # sharey=self.sbplot_matrix[0][n_row])
#                     else:
#                         sbplot_matrix[n_col][n_row] = fig.add_subplot(self.n_rows, self.n_cols, i, projection='polar')
#                                                                       # sharex=self.sbplot_matrix[n_col][0],
#                                                                       # sharey=self.sbplot_matrix[0][n_row])
#
#                         # sbplot_matrix[n_col][n_row].axes.get_yaxis().set_visible(False)
#                     # sbplot_matrix[n_col][n_row] = fig.add_subplot(n_rows, n_cols, i)
#                     i += 1
#         else:
#             raise NameError("type of the plot is not recognized. Use 'polar' or 'cartesian' ")
#
#         return fig, sbplot_matrix
#
#     def plot_one_task(self, ax, dic):
#
#
#         if dic["class"] == "lightcurve":
#             self.base_plots.plot_main(ax, dic, self.data)
#         else:
#             raise NameError("dic['class']: '{}' is not recognized"
#                             .format(dic["class"]))
#
#
#
#         im = 0
#         return 0
#
#
#         # if dic["name"] == "corr" or dic["name"] == "int":
#         #     im = self.plot_x_y_z2d_colormesh(ax, dic)
#         # elif dic["name"] == "count":
#         #     self.plot_countours(ax, dic)
#         #     im = 0
#         # elif dic["name"] == 'densmode':
#         #     self.plot_density_mode(ax, dic)
#         #     im = 0
#         # else:
#         #     raise NameError("name:{} is not recognized"
#         #                     .format(dic["name"]))
#         #
#         # # self.time
#         # return im
#
#     def set_plot_title(self, ax, plot_dic):
#         if plot_dic["title"] != '' and plot_dic["title"] != None:
#
#             if plot_dic["title"] == 'it':
#                 title = plot_dic["it"]
#             elif plot_dic["title"] == 'time [s]' or \
#                 plot_dic["title"] == 'time':
#                 title = "%.3f" % self.data.get_time(plot_dic["it"]) + " [s]"
#             elif plot_dic["title"] == 'time [ms]':
#                 title = "%.1f" % (self.data.get_time(plot_dic["it"]) * 1000) + " [ms]"
#             else:
#                 title = plot_dic["title"]
#             ax.title.set_text(r'{}'.format(title))
#
#     def plot_images(self):
#
#         # initializing the matrix of images for colorbars (if needed)
#         image_matrix = [[0
#                         for x in range(self.n_rows)]
#                         for y in range(self.n_cols)]
#
#
#         for n_row in range(self.n_rows):
#             for n_col in range(self.n_cols):
#                 for dic in self.set_plot_dics:
#                     if n_col + 1 == int(dic['position'][1]) and n_row + 1 == int(dic['position'][0]):
#                         print("Plotting n_row:{} n_col:{}".format(n_row, n_col))
#                         ax = self.sbplot_matrix[n_col][n_row]
#                         # dic = self.plot_dic_matrix[n_col][n_row]
#                         if isinstance(dic, int):
#                             Printcolor.yellow("Dictionary for row:{} col:{} not set".format(n_row, n_col))
#                             self.fig.delaxes(ax)  # delets the axis for empty plot
#                         else:
#                             dic = dict(dic)
#                             im = self.plot_one_task(ax, dic)
#                             self.set_plot_title(ax, dic)
#                             if not isinstance(im, int):
#                                 image_matrix[n_col][n_row] = im
#
#                             self.add_fancy_to_ax(ax, dic)
#                             self.account_for_shared(ax, n_row, n_col)
#
#         return image_matrix
#
#     def account_for_shared(self, ax, n_row, n_col):
#
#         ax.axes.xaxis.set_ticklabels([])
#         ax.axes.yaxis.set_ticklabels([])
#
#         if n_col > 0 and n_row < self.n_rows:
#             # ax.tick_params(labelbottom=False)
#             # ax.tick_params(labelleft=False)
#             #
#             # ax.set_yticklabels([])
#             # ax.set_xticklabels([])
#             #
#             # ax.get_yaxis().set_ticks([])
#             #
#             # ax.set_yticks([])
#             # ax.set_yticklabels(labels=[])
#             # # ax.set_yticklabels([]).remove()
#
#             # ax.axes.get_xaxis().set_visible(False)
#             # ax.axes.get_yaxis().set_visible(False)
#
#             ax.axes.xaxis.set_ticklabels([])
#             ax.axes.yaxis.set_ticklabels([])
#
#             # ax.tick_params(labelbottom=False)
#             # ax.tick_params(labelleft=False)
#             # ax.tick_params(labelright=False)
#
#             # ax.get_yaxis().set_visible(False)
#
#         if n_col > 0:
#             ax.set_ylabel('')
#
#         if n_row != self.n_rows-1:
#             ax.set_xlabel('')
#
#
#     def add_fancy_to_ax(self, ax, dic):
#
#         if self.gen_set["fancy_ticks"]:
#             ax.tick_params(axis='both', which='both', labelleft=True,
#                            labelright=False, tick1On=True, tick2On=True,
#                            labelsize=12, direction='in')
#
#         if "rmylbls" in dic.keys():
#             if dic["rmylbls"]:
#                 ax.set_yticklabels([])
#
#         if "rmxlbls" in dic.keys():
#             if dic["rmxlbls"]:
#                 ax.set_xticklabels([])
#
#         if dic["inverty"]:
#             ax.invert_yaxis()
#
#         if dic["xmin"] != None and dic["xmax"] != None:
#             ax.set_xlim(dic["xmin"], dic["xmax"])
#         if dic["ymin"] != None and dic["ymax"] != None:
#             ax.set_ylim(dic["ymin"], dic["ymax"])
#         if dic["xscale"] == 'log':
#             ax.set_xscale("log")
#         if dic["yscale"] == 'log':
#             ax.set_yscale("log")
#
#         ax.set_ylabel(dic["v_n_y"].replace('_', '\_'))
#         ax.set_xlabel(dic["v_n_x"].replace('_', '\_'))
#
#         if dic["legend"]:
#             ax.legend(fancybox=False, loc='best', shadow=False, fontsize=8)
#
#         if self.gen_set["minorticks_on"]:
#             ax.minorticks_on()
#
#     def plot_one_cbar(self, im, dic, n_row, n_col):
#
#         if dic["cbar"] != None and dic["cbar"] != '':
#
#             location = dic["cbar"].split(' ')[0]
#             shift_h = float(dic["cbar"].split(' ')[1])
#             shift_w = float(dic["cbar"].split(' ')[2])
#             cbar_width = 0.02
#
#
#             if location == 'right':
#                 ax_to_use = self.sbplot_matrix[-1][n_row]
#                 pos1 = ax_to_use.get_position()
#                 pos2 = [pos1.x0 + pos1.width + shift_h,
#                         pos1.y0,
#                         cbar_width,
#                         pos1.height]
#             elif location == 'left':
#                 ax_to_use = self.sbplot_matrix[-1][n_row]
#                 pos1 = ax_to_use.get_position()
#                 pos2 = [pos1.x0 - pos1.width - shift_h,
#                         pos1.y0,
#                         cbar_width,
#                         pos1.height]
#             elif location == 'bottom':
#                 ax_to_use = self.sbplot_matrix[n_col][-1]
#                 pos1 = ax_to_use.get_position()
#                 pos2 = [pos1.x0,
#                         pos1.y0 - pos1.height + shift_w,
#                         cbar_width,
#                         pos1.height]
#             else:
#                 raise NameError("cbar location {} not recognized. Use 'right' or 'bottom' "
#                                 .format(location))
#
#             cax1 = self.fig.add_axes(pos2)
#             if location == 'right':
#                 cbar = plt.colorbar(im, cax=cax1, extend='both')#, format='%.1e')
#             elif location == 'left':
#                 cbar = plt.colorbar(im, cax=cax1, extend='both')#, format='%.1e')
#                 cax1.yaxis.set_ticks_position('left')
#                 cax1.yaxis.set_label_position('left')
#             else:
#                 raise NameError("cbar location {} not recognized. Use 'right' or 'bottom' "
#                                 .format(location))
#             cbar.ax.set_title(r"{}".format(str(dic["v_n"]).replace('_', '\_')))
#
#     def plot_colobars(self):
#
#         for n_row in range(self.n_rows):
#             for n_col in range(self.n_cols):
#                 for dic in self.set_plot_dics:
#                     if n_col + 1 == int(dic['position'][1]) and n_row + 1 == int(dic['position'][0]):
#                         print("Colobar for n_row:{} n_col:{}".format(n_row, n_col))
#                         # ax  = self.sbplot_matrix[n_col][n_row]
#                         # dic = self.plot_dic_matrix[n_col][n_row]
#                         im  = self.image_matrix[n_col][n_row]
#                         if isinstance(dic, int):
#                             Printcolor.yellow("Dictionary for row:{} col:{} not set".format(n_row, n_col))
#                         else:
#                             self.plot_one_cbar(im, dic, n_row, n_col)
#
#
#         # for n_row in range(self.n_rows):
#         #     for n_col in range(self.n_cols):
#         #         print("Colobar for n_row:{} n_col:{}".format(n_row, n_col))
#         #         # ax  = self.sbplot_matrix[n_col][n_row]
#         #         dic = self.plot_dic_matrix[n_col][n_row]
#         #         im  = self.image_matrix[n_col][n_row]
#         #         if isinstance(dic, int):
#         #             Printcolor.yellow("Dictionary for row:{} col:{} not set".format(n_row, n_col))
#         #         else:
#         #             self.plot_one_cbar(im, dic, n_row, n_col)
#
#     def save_plot(self):
#
#         plt.subplots_adjust(hspace=self.gen_set["subplots_adjust_h"])
#         plt.subplots_adjust(wspace=self.gen_set["subplots_adjust_w"])
#         # plt.tight_layout()
#         plt.savefig('{}{}'.format(self.gen_set["figdir"], self.gen_set["figname"]),
#                     bbox_inches='tight', dpi=128)
#         plt.close()
#
#     def main(self):
#
#         # initializing the n_cols, n_rows
#         self.n_rows, self.n_cols = self.set_ncols_nrows()
#         # initializing the matrix of dictionaries of the
#         self.plot_dic_matrix = self.set_plot_dics_matrix()
#         # initializing the axis matrix (for all subplots) and image matrix fo colorbars
#         self.fig, self.sbplot_matrix = self.set_plot_matrix()
#         # plotting
#         self.image_matrix = self.plot_images()
#         # adding colobars
#         self.plot_colobars()
#
#
#         # saving the result
#         self.save_plot()












class PLOT_LIGHTCURVE_OLD:

    def __init__(self, sim):

        self.gen_set = {
            'figname': 'mkn.png'
        }

        self.sim = sim

        self.model_fpath = Paths.ppr_sims+self.sim+'/res_mkn/' + 'mkn_model.h5'

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
            plt.fill_between(dict_model['time'], maxs, mins, alpha=.5, color=color)

            labeled_bands = []
            for filter in dict_model.keys():
                # print("filter:{} ")
                if filter.split('_')[0] == band and filter in dict_obs_filters.keys():
                    if not band in labeled_bands:
                        ax.errorbar(dict_obs_filters[filter][:, 0], dict_obs_filters[filter][:, 1],
                                    yerr=dict_obs_filters[filter][:, 2],
                                    fmt='o', color=color, label=band)
                    else:
                        ax.errorbar(dict_obs_filters[filter][:, 0], dict_obs_filters[filter][:, 1],
                                    yerr=dict_obs_filters[filter][:, 2],
                                    fmt='o', color=color)
                    labeled_bands.append(band)


        plt.figure(figsize=(4.5, 3.6))
        ax = plt.subplot(111)

        for band in bands:
            # print('band')
            color = color_for_mkn_band(band)
            plot_for_one_band(ax, band, color)

        plt.minorticks_on()
        plt.xticks(fontsize=12)
        plt.yticks(fontsize=12)
        plt.tight_layout()
        plt.tick_params(axis='both', which='both', labelleft=True, labelright=False,
                        tick1On=True, tick2On=True, labelsize=12, direction='in')  # labeltop
        plt.gca().invert_yaxis()
        plt.xscale("log")
        plt.xlim(xmin=0.2)
        plt.ylim(ymin=25)#, ymax=2e-1)
        # plt.xlim(xmin=50, xmax=210)
        plt.ylabel(r"AB magnitude at 40 Mpc", fontsize=12)
        plt.xlabel(r"time [days]", fontsize=12)
        plt.legend(loc='best', numpoints=1, fontsize=12)


        plt.savefig(Paths.ppr_sims+self.sim+'/res_mkn/' + self.gen_set['figname'], bbox_inches='tight', dpi=128)
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



""" --- TASK SPECIFIC FUNCTIONS --- """

def plot_many_lightcuves():

    sim = "DD2_M13641364_M0_SR"

    o_data = LOAD_LIGHTCURVES(sim)
    o_task_plots = PLOT_LIGHTCURVE()
    o_plots = PLOT_MANY_TASKS(o_data, o_task_plots)

    o_plots.gen_set["figdir"] = Paths.ppr_sims + sim + "/res_mkn/"
    o_plots.gen_set["figname"]  = "mkn.png"
    o_plots.gen_set["figsize"] =  (9.1, 9.1)  # <->, |]
    o_plots.gen_set["subplots_adjust_h"] = 0.0
    o_plots.gen_set["subplots_adjust_w"] = 0.0
    o_plots.gen_set["sharex"] = True
    o_plots.gen_set["sharey"] = True


    dd2_blue_dic = {
        'class': 'lightcurve', 'name': '-', 'sim': "DD2_M13641364_M0_SR", 'fname': 'mkn_model.h5',
        'position': (1, 1), 'obs': True,
        'band': 'g', 'fillcolor': 'gray', 'obscolor': 'black', 'marker': 'o',
        'title': 'DD2', 'cbar': None,
        'v_n_x': r"time [days]", 'v_n_y': r"AB magnitude at 40 Mpc", 'v_n': None,
        'xmin': None, 'xmax': None, 'ymin': None, 'ymax': None,
        'xscale': 'log', 'yscale': "log", 'inverty': True,
        'fancy_ticks': True, 'legend': True, 'rmylbls': False
    }

    dd2_green_dic = {
        'class': 'lightcurve', 'name': '-', 'sim': "DD2_M13641364_M0_SR", 'fname': 'mkn_model.h5',
        'position': (2, 1), 'obs': True,
        'band': 'r', 'fillcolor': 'gray', 'obscolor': 'black', 'marker': 'o',
        'title': None, 'cbar': None,
        'v_n_x': r"time [days]", 'v_n_y': r"AB magnitude at 40 Mpc", 'v_n': None,
        'xmin': None, 'xmax': None, 'ymin': None, 'ymax': None,
        'xscale': 'log', 'yscale': "log", 'inverty': True,
        'fancy_ticks': True, 'legend': True, 'rmylbls': False
    }

    dd2_yellow_dic = {
        'class': 'lightcurve', 'name': '-', 'sim': "DD2_M13641364_M0_SR", 'fname': 'mkn_model.h5',
        'position': (3, 1), 'obs': True,
        'band': 'z', 'fillcolor': 'gray', 'obscolor': 'black', 'marker': 'o',
        'title': None, 'cbar': None,
        'v_n_x': r"time [days]", 'v_n_y': r"AB magnitude at 40 Mpc", 'v_n': None,
        'xmin': None, 'xmax': None, 'ymin': None, 'ymax': None,
        'xscale': 'log', 'yscale': "log", 'inverty': True,
        'fancy_ticks': True, 'legend': True, 'rmylbls': False
    }


    ls220_blue_dic = {
        'class': 'lightcurve', 'name': '-', 'sim': "LS220_M13641364_M0_SR", 'fname': 'mkn_model.h5',
        'position': (1, 2), 'obs': True,
        'band': 'g', 'fillcolor': 'gray', 'obscolor': 'black', 'marker': 'o',
        'title': 'LS220', 'cbar': None,
        'v_n_x': r"time [days]", 'v_n_y': r"AB magnitude at 40 Mpc", 'v_n': None,
        'xmin': None, 'xmax': None, 'ymin': None, 'ymax': None,
        'xscale': 'log', 'yscale': "log", 'inverty': True,
        'fancy_ticks': True, 'legend': True, 'rmylbls': True
    }

    ls220_green_dic = {
        'class': 'lightcurve', 'name': '-', 'sim': "LS220_M13641364_M0_SR", 'fname': 'mkn_model.h5',
        'position': (2, 2), 'obs': True,
        'band': 'r', 'fillcolor': 'gray', 'obscolor': 'black', 'marker': 'o',
        'title': None, 'cbar': None,
        'v_n_x': r"time [days]", 'v_n_y': r"AB magnitude at 40 Mpc", 'v_n': None,
        'xmin': None, 'xmax': None, 'ymin': None, 'ymax': None,
        'xscale': 'log', 'yscale': "log", 'inverty': True,
        'fancy_ticks': True, 'legend': True, 'rmylbls': True
    }

    ls220_yellow_dic = {
        'class': 'lightcurve', 'name': '-', 'sim': "LS220_M13641364_M0_SR", 'fname': 'mkn_model.h5',
        'position': (3, 2), 'obs': True,
        'band': 'z', 'fillcolor': 'gray', 'obscolor': 'black', 'marker': 'o',
        'title': None, 'cbar': None,
        'v_n_x': r"time [days]", 'v_n_y': r"AB magnitude at 40 Mpc", 'v_n': None,
        'xmin': None, 'xmax': None, 'ymin': None, 'ymax': None,
        'xscale': 'log', 'yscale': "log", 'inverty': True,
        'fancy_ticks': True, 'legend': True, 'rmylbls': True
    }


    sly4_blue_dic = {
        'class': 'lightcurve', 'name': '-', 'sim': "SLy4_M13641364_M0_SR", 'fname': 'mkn_model.h5',
        'position': (1, 3), 'obs': True,
        'band': 'g', 'fillcolor': 'gray', 'obscolor': 'black', 'marker': 'o',
        'title': 'SLy4', 'cbar': None,
        'v_n_x': r"time [days]", 'v_n_y': r"AB magnitude at 40 Mpc", 'v_n': None,
        'xmin': None, 'xmax': None, 'ymin': None, 'ymax': None,
        'xscale': 'log', 'yscale': "log", 'inverty': True,
        'fancy_ticks': True, 'legend': True, 'rmylbls': True
    }

    sly4_green_dic = {
        'class': 'lightcurve', 'name': '-', 'sim': "SLy4_M13641364_M0_SR", 'fname': 'mkn_model.h5',
        'position': (2, 3), 'obs': True,
        'band': 'r', 'fillcolor': 'gray', 'obscolor': 'black', 'marker': 'o',
        'title': None, 'cbar': None,
        'v_n_x': r"time [days]", 'v_n_y': r"AB magnitude at 40 Mpc", 'v_n': None,
        'xmin': None, 'xmax': None, 'ymin': None, 'ymax': None,
        'xscale': 'log', 'yscale': "log", 'inverty': True,
        'fancy_ticks': True, 'legend': True, 'rmylbls': True
    }

    sly4_yellow_dic = {
        'class': 'lightcurve', 'name': '-', 'sim': "SLy4_M13641364_M0_SR", 'fname': 'mkn_model.h5',
        'position': (3, 3), 'obs': True,
        'band': 'z', 'fillcolor': 'gray', 'obscolor': 'black', 'marker': 'o',
        'title': None, 'cbar': None,
        'v_n_x': r"time [days]", 'v_n_y': r"AB magnitude at 40 Mpc", 'v_n': None,
        'xmin': None, 'xmax': None, 'ymin': None, 'ymax': None,
        'xscale': 'log', 'yscale': "log", 'inverty': True,
        'fancy_ticks': True, 'legend': True, 'rmylbls': True
    }



    o_plots.set_plot_dics = [dd2_blue_dic,dd2_green_dic,dd2_yellow_dic,
                             ls220_blue_dic,ls220_green_dic,ls220_yellow_dic,
                             sly4_blue_dic,sly4_green_dic,sly4_yellow_dic]

    o_plots.main()

    exit(1)


def compute_for_1_simulations_fully_informed_2_comp_varying_total_bern_mass():

    from d1analysis import ADD_METHODS_1D

    # times_to_extract_total_mass
    ti = [.020, 0.1, 0.15, 0.20]
    ti = np.arange(0.015, stop=0.250, step=0.01)
    print(ti)#; exit(1)
    sim = "DD2_M13641364_M0_LK_SR_R04"
    o_data = ADD_METHODS_1D(sim)

    # sim = "LS220_M13641364_M0_LK_SR"
    # make_model = COMPUTE_LIGHTCURVE(sim)
    # make_model.output_fname = 'mkn_model1.h5'  # time in ms for a file
    # make_model.set_dyn_iso_aniso = "aniso"
    # make_model.set_psdyn_iso_aniso = ""
    # make_model.set_wind_iso_aniso = ""
    # make_model.set_secular_iso_aniso = ""
    # make_model.set_par_war()
    # make_model.compute_save_lightcurve(write_output=True)
    # exit(1)

    for t in ti:

        # times, masses = o_data.get_extrapolated_arr(v_n_x='t_tot_flux', v_n_y='mass_tot_flux',
        #                                             criterion='_0_b_w',
        #                                             x_left= None, x_right= 180, # percent
        #                                             x_start= 0.040, x_stop= None,
        #                                             depth= 1000, method= 1)

        times, masses = o_data.get_int_extra_arr(v_n_x='t_tot_flux', v_n_y='mass_tot_flux',
                                                 criterion='_0_b_w',
                                                 x_left= None, x_right= 180,  # percent
                                                 x_extr_start= 0.040, x_extr_stop= None,
                                                 depth= 1000, method= 1)

        m = masses[find_nearest_index(times, t)]
        print("\tt:{} m:{}".format(t, m))

        make_model = COMPUTE_LIGHTCURVE(sim)
        make_model.output_fname = 'mkn_model2_m{}.h5'.format(int(m*1000)) # time in ms for a file
        make_model.set_dyn_iso_aniso = "aniso"
        make_model.set_psdyn_iso_aniso = "aniso"
        make_model.set_wind_iso_aniso = ""
        make_model.set_secular_iso_aniso = ""
        make_model.set_par_war()
        make_model.ejecta_vars['psdynamics']['override_m_ej'] = True
        make_model.ejecta_vars['psdynamics']['m_ej'] =  m

        # print(make_model.ejecta_vars['psdynamics']); exit(0)

        make_model.compute_save_lightcurve(write_output=True)

    exit(1)

def compute_for_1_simulations_fully_informed_2_comp_varying_total_bern_mass_psdyn_opacity():

    import seaborn as sns
    # flights = sns.load_dataset("flights")
    # print(type(flights))
    # exit(1)


    from d1analysis import ADD_METHODS_1D

    # times_to_extract_total_mass

    kappas = [0.5, 1., 1.5, 2., 2.5, 3.]
    # kappas = [5, 10, 15, 20, 25, 30]
    ti = [.05, 0.1, 0.15, 0.20]

    sim = "DD2_M13641364_M0_LK_SR_R04"
    o_data = ADD_METHODS_1D(sim)

    log_l = [[0 for k in range(len(kappas))] for t in range(len(ti))]
    for it, t in enumerate(ti):

        times, masses = o_data.get_extrapolated_arr(v_n_x='t_tot_flux', v_n_y='mass_tot_flux',
                                                    criterion='_0_b_w',
                                                    x_left= None, x_right= 150, # percent
                                                    x_start= 0.040, x_stop= None,
                                                    depth= 1000, method= 1)
        m = masses[find_nearest_index(times, t)]

        for ik, k in enumerate(kappas):

            make_model = COMPUTE_LIGHTCURVE(sim)
            make_model.output_fname = 'mkn_model2_t{}_kh{}.h5'.format(int(t*1000), int(k*10)) # time in ms for a file
            make_model.set_dyn_iso_aniso = "aniso"
            make_model.set_psdyn_iso_aniso = "aniso"
            make_model.set_wind_iso_aniso = ""
            make_model.set_secular_iso_aniso = ""
            make_model.set_par_war()
            make_model.ejecta_vars['psdynamics']['override_m_ej'] = True
            make_model.ejecta_vars['psdynamics']['m_ej'] =  m
            make_model.ejecta_vars['psdynamics']['kappa_low'] = k
            # print(make_model.ejecta_vars['psdynamics']); exit(0)

            print("\tt:{} m:{} k:{}".format(t, m, k))
            log_l_ = make_model.compute_save_lightcurve(write_output=True)
            log_l[it][ik] = log_l_


    print("All done")
    table = []
    for it, t in enumerate(ti):
        table = np.append(table, np.array(log_l[it][:]))

    table = np.reshape(table, (len(ti), len(kappas)))

    df = pd.DataFrame(data=(np.array(table * -1, dtype=int)), index=ti, columns=kappas)#.corr()

    print(df)


    # df = pd.DataFrame({'kappas':kappas, 'ti':ti, 'logl':table})


    sns.set()
    f, ax = plt.subplots(figsize=(4.2, 3.6))
    # flights_long = sns.load_dataset("flights")
    sns.heatmap(df, annot=True, fmt='d', linewidths=.5, ax=ax, annot_kws={"rotation":45}, cbar_kws={"format": '%.1e'})
    plt.gca().invert_yaxis()
    ax.set_title(r"DD2\_M13641364\_M0\_LK\_SR\_R04")
    ax.set_xlabel(r"$\kappa_{Ye > 0.25}$")
    ax.set_ylabel(r"$t_{evol}$")
    plt.savefig('{}.png'.format(Paths.plots+'mkn_logl_evol_kappahigh'), bbox_inches='tight', dpi=128)
    plt.close()


    exit(1)

def compute_for_1_simulations_fully_informed_2_comp_varying_total_bern_mass_psdyn_theta():

    import seaborn as sns
    # flights = sns.load_dataset("flights")
    # print(type(flights))
    # exit(1)


    from d1analysis import ADD_METHODS_1D

    # times_to_extract_total_mass

    thetas = [90., 75., 45., 25., 15., 0]
    # kappas = [5, 10, 15, 20, 25, 30]
    ti = [.05, 0.1, 0.15, 0.20]

    sim = "DD2_M13641364_M0_LK_SR_R04"
    o_data = ADD_METHODS_1D(sim)

    log_l = [[0 for k in range(len(thetas))] for t in range(len(ti))]
    for it, t in enumerate(ti):

        times, masses = o_data.get_extrapolated_arr(v_n_x='t_tot_flux', v_n_y='mass_tot_flux',
                                                    criterion='_0_b_w',
                                                    x_left= None, x_right= 150, # percent
                                                    x_start= 0.040, x_stop= None,
                                                    depth= 1000, method= 1)
        m = masses[find_nearest_index(times, t)]

        for itheta, theta in enumerate(thetas):

            make_model = COMPUTE_LIGHTCURVE(sim)
            make_model.output_fname = 'mkn_model2_t{}_theta{}.h5'.format(int(t*1000), int(theta)) # time in ms for a file
            make_model.set_dyn_iso_aniso = "aniso"
            make_model.set_psdyn_iso_aniso = "aniso"
            make_model.set_wind_iso_aniso = ""
            make_model.set_secular_iso_aniso = ""
            make_model.set_par_war()
            make_model.source_name = 'AT2017gfo view_angle={}'.format(theta)
            make_model.ejecta_vars['psdynamics']['override_m_ej'] = True
            make_model.ejecta_vars['psdynamics']['m_ej'] =  m
            # print(make_model.ejecta_vars['psdynamics']); exit(0)

            print("\tt:{} m:{} theta:{}".format(t, m, theta))
            log_l_ = make_model.compute_save_lightcurve(write_output=True)
            log_l[it][itheta] = log_l_


    print("All done")
    table = []
    for it, t in enumerate(ti):
        table = np.append(table, np.array(log_l[it][:]))

    table = np.reshape(table, (len(ti), len(thetas)))

    df = pd.DataFrame(data=(np.array(table * -1, dtype=int)), index=ti, columns=thetas)#.corr()

    print(df)


    # df = pd.DataFrame({'kappas':kappas, 'ti':ti, 'logl':table})


    sns.set()
    f, ax = plt.subplots(figsize=(4.2, 3.6))
    # flights_long = sns.load_dataset("flights")
    sns.heatmap(df, annot=True, fmt='d', linewidths=.5, ax=ax, annot_kws={"rotation":45}, cbar_kws={"format": '%.1e'})
    plt.gca().invert_yaxis()
    ax.set_title(r"DD2\_M13641364\_M0\_LK\_SR\_R04")
    ax.set_xlabel(r"Angle from rotational axis")
    ax.set_ylabel(r"$t_{evol}$")
    plt.savefig('{}.png'.format(Paths.plots+'mkn_logl_evol_theta'), bbox_inches='tight', dpi=128)
    plt.close()


    exit(1)

def compute_for_1_simulations_fully_informed_2_comp_varying_total_kappa_low_theta():

    import seaborn as sns
    # flights = sns.load_dataset("flights")
    # print(type(flights))
    # exit(1)


    from d1analysis import ADD_METHODS_1D

    # times_to_extract_total_mass

    thetas = [90., 75., 45., 25., 15., 0]
    kappas = [0.5, 1., 1.5, 2., 2.5, 3]
    ti = 0.20

    sim = "DD2_M13641364_M0_LK_SR_R04"
    o_data = ADD_METHODS_1D(sim)

    times, masses = o_data.get_extrapolated_arr(v_n_x='t_tot_flux', v_n_y='mass_tot_flux',
                                                criterion='_0_b_w',
                                                x_left=None, x_right=150,  # percent
                                                x_start=0.040, x_stop=None,
                                                depth=1000, method=1)
    m = masses[find_nearest_index(times, ti)]

    log_l = [[0 for k in range(len(thetas))] for t in range(len(kappas))]
    for ik, k in enumerate(kappas):

        for itheta, theta in enumerate(thetas):

            make_model = COMPUTE_LIGHTCURVE(sim)
            make_model.output_fname = 'mkn_model2_t{}_k{}_theta{}.h5'.format(int(ti*1000), int(k*10), int(theta)) # time in ms for a file
            make_model.set_dyn_iso_aniso = "aniso"
            make_model.set_psdyn_iso_aniso = "aniso"
            make_model.set_wind_iso_aniso = ""
            make_model.set_secular_iso_aniso = ""
            make_model.set_par_war()
            make_model.source_name = 'AT2017gfo view_angle={}'.format(theta)
            make_model.ejecta_vars['psdynamics']['override_m_ej'] = True
            make_model.ejecta_vars['psdynamics']['m_ej'] =  m
            make_model.ejecta_vars['psdynamics']['kappa_low'] = k
            # print(make_model.ejecta_vars['psdynamics']); exit(0)

            print("\tt:{} m:{} theta:{}, k:{}".format(ti, m, theta, k))
            log_l_ = make_model.compute_save_lightcurve(write_output=True)
            log_l[ik][itheta] = log_l_


    print("All done")
    table = []
    for ik, k in enumerate(kappas):
        table = np.append(table, np.array(log_l[ik][:]))

    table = np.reshape(table, (len(kappas), len(thetas)))

    df = pd.DataFrame(data=(np.array(table * -1, dtype=int)), index=kappas, columns=thetas)#.corr()

    print(df)


    # df = pd.DataFrame({'kappas':kappas, 'ti':ti, 'logl':table})


    sns.set()
    f, ax = plt.subplots(figsize=(4.2, 3.6))
    # flights_long = sns.load_dataset("flights")
    sns.heatmap(df, annot=True, fmt='d', linewidths=.5, ax=ax, annot_kws={"rotation":45}, cbar_kws={"format": '%.1e'})
    plt.gca().invert_yaxis()
    ax.set_title(r"$t:{} [ms]$ $m:{:.1f}$".format(ti, m * 1e2) + " $[10^{-2}M_{\odot}]$")
    ax.set_xlabel(r"Angle from rotational axis")
    ax.set_ylabel(r"$\kappa_{Ye < 0.25}$")
    plt.savefig('{}.png'.format(Paths.plots+'mkn_logl_evol_kappa_low_theta_m{}'.format(int(ti*1000))), bbox_inches='tight', dpi=128)
    plt.close()


    exit(1)


def compute_for_1_simulations_fully_informed_2_comp():

    sim = "DD2_M13641364_M0_LK_SR_R04"

    make_model = COMPUTE_LIGHTCURVE(sim)
    make_model.output_fname = 'mkn_model2.h5'
    make_model.set_dyn_iso_aniso = "aniso"
    make_model.set_psdyn_iso_aniso = "aniso"
    make_model.set_wind_iso_aniso = ""
    make_model.set_secular_iso_aniso = ""
    make_model.set_par_war()

    print("\t Bern. Ej. Mass: {} ".format(make_model.ejecta_vars['psdynamics']['m_ej']))

    make_model.compute_save_lightcurve(write_output=True)

    exit(1)

def compute_for_1_simulations_fully_informed_4comp():

    sim = "LS220_M13641364_M0_LK_SR" # DD2_M13641364_M0_LK_SR_R04

    make_model = COMPUTE_LIGHTCURVE(sim)
    make_model.output_fname = 'mkn_model2_2.h5'
    make_model.set_dyn_iso_aniso = "aniso"
    make_model.set_psdyn_iso_aniso = "aniso"
    make_model.set_wind_iso_aniso = "aniso"
    # make_model.ejecta_vars["wind"]['xi_disk'] = 0.2
    # make_model.ejecta_vars["wind"]['central_vel'] = 0.08
    # make_model.ejecta_vars["wind"]['high_lat_op'] = 0.5
    # make_model.ejecta_vars["wind"]['low_lat_op'] = 5.
    make_model.set_secular_iso_aniso = "aniso"
    # make_model.ejecta_vars["secular"]['xi_disk'] = 0.2
    # make_model.ejecta_vars["secular"]['central_vel'] = 0.08
    # make_model.ejecta_vars["secular"]['central_op'] = 5.0
    make_model.set_par_war()

    if sim == "LS220_M13641364_M0_LK_SR": make_model.glob_vars['m_disk'] = 0.07062 # USING non LK!

    make_model.compute_save_lightcurve(write_output=True)

    exit(1)

def compute_for_1_simulations_fully_informed_2_comp_varying_eps():

    from d1analysis import ADD_METHODS_1D

    # times_to_extract_total_mass
    eps = [2e18, 6e18, 12e18, 16e18, 20e18]
    sim = "DD2_M13641364_M0_LK_SR_R04"
    o_data = ADD_METHODS_1D(sim)

    for e in eps:

        print("\teps:{}".format(str(e).replace('+', '')))
        # exit(1)
        make_model = COMPUTE_LIGHTCURVE(sim)
        make_model.output_fname = 'mkn_model2_eps{}.h5'.format(str(e).replace('+', '')) # time in ms for a file
        make_model.set_dyn_iso_aniso = "aniso"
        make_model.set_psdyn_iso_aniso = "aniso"
        make_model.set_wind_iso_aniso = ""
        make_model.set_secular_iso_aniso = ""
        make_model.set_par_war()
        make_model.glob_vars['eps0'] = float(e) # {2; 6; 12; 16; 20}
        make_model.compute_save_lightcurve(write_output=True)

    exit(1)


def compute_for_3_simulations_fully_informed_4comp():



    sim = "DD2_M13641364_M0_SR"
    make_model = COMPUTE_LIGHTCURVE(sim)
    make_model.set_wind_iso_aniso = "aniso"
    # make_model.ejecta_vars["wind"]['xi_disk'] = 0.2
    # make_model.ejecta_vars["wind"]['central_vel'] = 0.08
    # make_model.ejecta_vars["wind"]['high_lat_op'] = 0.5
    # make_model.ejecta_vars["wind"]['low_lat_op'] = 5.

    make_model.set_secular_iso_aniso = "aniso"
    # make_model.ejecta_vars["secular"]['xi_disk'] = 0.2
    # make_model.ejecta_vars["secular"]['central_vel'] = 0.08
    # make_model.ejecta_vars["secular"]['central_op'] = 5.0
    make_model.set_par_war()
    make_model.compute_save_lightcurve(write_output=True)

    " --- --- ---"

    sim = "LS220_M13641364_M0_SR"
    make_model = COMPUTE_LIGHTCURVE(sim)
    make_model.set_wind_iso_aniso = "aniso"
    # make_model.ejecta_vars["wind"]['xi_disk'] = 0.2
    # make_model.ejecta_vars["wind"]['central_vel'] = 0.08
    # make_model.ejecta_vars["wind"]['high_lat_op'] = 0.5
    # make_model.ejecta_vars["wind"]['low_lat_op'] = 5.

    make_model.set_secular_iso_aniso = "aniso"
    # make_model.ejecta_vars["secular"]['xi_disk'] = 0.2
    # make_model.ejecta_vars["secular"]['central_vel'] = 0.08
    # make_model.ejecta_vars["secular"]['central_op'] = 5.0
    make_model.set_par_war()
    make_model.compute_save_lightcurve(write_output=True)

    " --- --- ---"

    sim = "SLy4_M13641364_M0_SR"
    make_model = COMPUTE_LIGHTCURVE(sim)
    make_model.set_wind_iso_aniso = "aniso"
    # make_model.ejecta_vars["wind"]['xi_disk'] = 0.2
    # make_model.ejecta_vars["wind"]['central_vel'] = 0.08
    # make_model.ejecta_vars["wind"]['high_lat_op'] = 0.5
    # make_model.ejecta_vars["wind"]['low_lat_op'] = 5.

    make_model.set_secular_iso_aniso = "aniso"
    # make_model.ejecta_vars["secular"]['xi_disk'] = 0.2
    # make_model.ejecta_vars["secular"]['central_vel'] = 0.08
    # make_model.ejecta_vars["secular"]['central_op'] = 5.0
    make_model.set_par_war()
    make_model.compute_save_lightcurve(write_output=True)

    exit(1)


def compute_for_3_simulations_varius_component_set():


    sims = ["DD2_M13641364_M0_SR",
            "LS220_M13641364_M0_SR",
            "SLy4_M13641364_M0_SR"]

    for sim in sims:

        make_model = COMPUTE_LIGHTCURVE(sim)
        make_model.output_fname             = 'mkn_model1.h5'
        make_model.set_dyn_iso_aniso        = "aniso"
        make_model.set_psdyn_iso_aniso      = ""
        make_model.set_wind_iso_aniso       = ""
        make_model.set_secular_iso_aniso    = ""
        make_model.set_par_war()
        make_model.compute_save_lightcurve(write_output=True)

        make_model = COMPUTE_LIGHTCURVE(sim)
        make_model.output_fname             = 'mkn_model2.h5'
        make_model.set_dyn_iso_aniso        = "aniso"
        make_model.set_psdyn_iso_aniso      = "aniso"
        make_model.set_wind_iso_aniso       = ""
        make_model.set_secular_iso_aniso    = ""
        make_model.set_par_war()
        make_model.compute_save_lightcurve(write_output=True)

        make_model = COMPUTE_LIGHTCURVE(sim)
        make_model.output_fname             = 'mkn_model3.h5'
        make_model.set_dyn_iso_aniso        = "aniso"
        make_model.set_psdyn_iso_aniso      = "aniso"
        make_model.set_wind_iso_aniso       = "aniso"
        make_model.set_secular_iso_aniso    = ""
        make_model.set_par_war()
        make_model.compute_save_lightcurve(write_output=True)

        make_model = COMPUTE_LIGHTCURVE(sim)
        make_model.output_fname             = 'mkn_model4.h5'
        make_model.set_dyn_iso_aniso        = "aniso"
        make_model.set_psdyn_iso_aniso      = "aniso"
        make_model.set_wind_iso_aniso       = "aniso"
        make_model.set_secular_iso_aniso    = "aniso"
        make_model.set_par_war()
        make_model.compute_save_lightcurve(write_output=True)

    exit(1)
# compute_for_3_simulations_varius_component_set()

if __name__ == '__main__':

    ''' test '''
    # o_dat = COMBINE_LIGHTCURVES("DD2_M13641364_M0_LK_SR_R04")
    # # o_dat.get_table(); exit()
    # o_dat.get_mkn_model("mkn_model2_m20.h5")
    # # print(o_dat.get_mkn_model("mkn_model2_m20.h5")["psdynamics"])
    # print(o_dat.get_attr("dynamics", "mkn_model2_m20.h5"))
    # o_dat.get_model_min_max('g', "mkn_model2_m20.h5")
    # print(o_dat.get_table("g"))
    # exit(1)

    ''' debugging '''
    
    compute_for_1_simulations_fully_informed_2_comp_varying_total_bern_mass()
    # compute_for_1_simulations_fully_informed_2_comp_varying_total_kappa_low_theta()
    # compute_for_1_simulations_fully_informed_2_comp_varying_total_bern_mass_psdyn_theta()
    # compute_for_1_simulations_fully_informed_2_comp_varying_total_bern_mass_psdyn_opacity()
    # compute_for_1_simulations_fully_informed_2_comp()
    # compute_for_1_simulations_fully_informed_4comp()



    # o_dat = EXTRACT_LIGHTCURVE("DD2_M13641364_M0_LK_SR_R04")
    # o_dat.get_mismatch("Ks", "mkn_model2_t50.h5")
    # exit(1)


    # compute_for_1_simulations_fully_informed_2_comp_varying_total_bern_mass()
    # compute_for_1_simulations_fully_informed_2_comp_varying_eps()

    # compute_for_3_simulations_varius_component_set()


    # plot_many_lightcuves()

    #
    # make_model = COMPUTE_LIGHTCURVE("DD2_M13641364_M0_SR")
    # make_model.compute_save_lightcurve(write_output=True)
    # make_model.plot_test_smooth_profile('_0')

    # plot_model = PLOT_LIGHTCURVE_OLD("DD2_M13641364_M0_SR")
    # plot_model.plot_for_several_bands(['g', 'r', 'z'], "tst_mkn")

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