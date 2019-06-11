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
from gw import SIM_GW



class LOAD_FILES:

    def __init__(self, sim):
        self.gen_set = {
            "outflow": "outflow",
            "collated": "collated",
            "waveforms": "waveforms",
            "res_3d": "res_3d"
        }

        self.sim = sim

        self.list_criteria = ['', "_0", "_1", "_0_b", "_1_b", "_0_b_w", "_1_b_w"]

        self.list_outflowed_files = [
                           "total_flux.dat",
                           "hist_temperature.dat",
                           "hist_theta.dat",
                           "hist_ye.dat",
                           "hist_rho.dat",
                           "hist_entropy.dat",
                           "hist_vel_inf.dat",
                           "hist_vel_inf_bern.dat",

                           "ejecta_profile.dat",
                           "ejecta_profile_bern.dat",
                           ]

        self.list_collated_files = [
            "dens_unbnd.norm1.asc",
            "dens_unbnd_bernoulli.norm1.asc",
            "dens.norm1.asc",
        ]

        self.list_gw_files = [
            "waveform_l2_m2.dat",
            "tmerger.dat",
            "tcoll.dat"
        ]

        self.list_mkn_files = [
            "mkn_model.h5",
            "AT2017gfo.h5"
        ]

        self.list_3d_data_files = [
            "disk_mass.txt"
        ]

        self.list_files = self.list_outflowed_files + \
                          self.list_collated_files + \
                          self.list_gw_files + \
                          self.list_3d_data_files + \
                          self.list_mkn_files

        self.in_data_matrix = [[np.zeros(0,)
                               for x in range(len(self.list_files))]
                               for y in range(len(self.list_criteria))]

    def check_criterion(self, criterion):
        if not criterion in self.list_criteria:
            raise NameError("criteria: {} not in the list of them: {}"
                            .format(criterion, self.list_criteria))

    def check_fname(self, fname):
        if not fname in self.list_files:
            raise NameError("fname: {} not in the lit of them: {}"
                            .format(fname, self.list_files))

    def i_fn(self, fname):
        return int(self.list_files.index(fname))

    def i_cr(self, criterion):
        return int(self.list_criteria.index(criterion))

    def load_disk_mass_ev(self):

        # check if res_3d exists:
        if not os.path.isdir(Paths.ppr_sims + self.sim + '/' + self.gen_set['res_3d'] + '/'):
            return np.zeros(0,), np.zeros(0,)


        # check if any disk_mass.txt exist:
        it_dirs = os.listdir(Paths.ppr_sims + self.sim + '/' + self.gen_set['res_3d'] + '/')
        iterations, disk_masses = [], []
        if len(it_dirs) == 0:
            Printcolor.yellow("No iteration directories found in {}".
                              format(Paths.ppr_sims + self.sim + '/' + self.gen_set['res_3d'] + '/'))
        else:
            for it_dir in it_dirs:
                path = Paths.ppr_sims + self.sim + '/' + self.gen_set['res_3d'] + '/' + it_dir + '/' + "disk_mass.txt"
                # print(path)
                if os.path.isfile(path):
                    try:
                        it = int(it_dir)
                    except:
                        raise ValueError("faild to extract interation from disk_mass.txt path:\n{}"
                                         .format(path))
                    try:
                        mass = np.float(np.loadtxt(path, unpack=True))
                    except:
                        raise ValueError("failed to load disk_mass.txt from")

                    iterations.append(it)
                    disk_masses.append(mass)

        if len(iterations) == 0:
            Printcolor.yellow("While {} exists, there are no 'disk_mass.txt files in it".format(self.sim + "/res_3d/"))
            return np.zeros(0, ), np.zeros(0, )

        # get times for iterations:
        path_ittime = Paths.ppr_sims + self.sim + "/res_3d/" + "ittime.txt"
        if os.path.isfile(path_ittime):
            it_time_masses = np.zeros(3)
            all_iterations, all_times = np.loadtxt(path_ittime, usecols=(0, 1), unpack=True)
            for it, m in zip(iterations, disk_masses):
                if it not in all_iterations:
                    raise ValueError("it:{} not in all_iteration:{} (from ittime.txt)"
                                     .format(it, all_iterations))
                time_ = int(np.where(it == all_iterations)[0])
                it_time_masses = np.vstack(it_time_masses, [it, time_, m])
            it_time_masses = it_time_masses[1:-1, 1:-1, 1:-1]
        else:
            Printcolor.yellow("file ittime.txt not found in {}".format(self.sim + "/res_3d/"))
            times = list(interpoate_time_form_it(iterations, Paths.gw170817 + self.sim + '/', time_units='s',
                                            extrapolate=True))
            # print(np.array(iterations).shape, np.array(times).shape, np.array(disk_masses).shape)

            it_time_masses = np.vstack((iterations, times, disk_masses)).T

        final_times, final_masses = x_y_z_sort(it_time_masses[:, 1], it_time_masses[:, 2], np.zeros(0, ), 0)

        # print(final_masses)

        return np.vstack((final_times, final_masses)).T

    def load_file(self, criterion, fname):

        if criterion == '' and fname in self.list_3d_data_files:
            fpath = ''
            self.in_data_matrix[self.i_cr(criterion)][self.i_fn(fname)] = self.load_disk_mass_ev()

        elif criterion == '' and fname in self.list_collated_files:
            fpath = Paths.ppr_sims + self.sim + '/' + self.gen_set['collated'] + '/' + fname

        elif criterion == '' and fname in self.list_gw_files:
            fpath = Paths.ppr_sims + self.sim + '/' + self.gen_set['waveforms'] + '/' + fname

        elif criterion != '' and fname in self.list_outflowed_files:
            fpath = Paths.ppr_sims + self.sim + '/' + self.gen_set['outflow'] + criterion + '/' + fname

        else:
            raise NameError("criterion:{} and fname:{} are not interconnected"
                            .format(criterion, fname))

        if fpath != '':
            data = np.loadtxt(fpath)
            self.in_data_matrix[self.i_cr(criterion)][self.i_fn(fname)] = data

    def is_loaded(self, criterion, fname):
        data = self.in_data_matrix[self.i_cr(criterion)][self.i_fn(fname)]
        if len(data) == 0:
            self.load_file(criterion, fname)

    def get_file(self, fname, criterion=''):
        self.check_criterion(criterion)
        self.check_fname(fname)
        self.is_loaded(criterion, fname)

        return self.in_data_matrix[self.i_cr(criterion)][self.i_fn(fname)]


class COMPUTE_STORE_ARR(LOAD_FILES):

    def __init__(self, sim):
        LOAD_FILES.__init__(self, sim)

        self.list_array_names = ['t_total_mass', 'total_mass',
                                't_unb_mass', 'unb_mass',
                                't_unb_mass_bern', 'unb_mass_bern',
                                't_tot_flux', 'tot_flux', 'mass_tot_flux',
                                 't_disk_mass', 'disk_mass',

                                 'hist_ye', 'hist_ye_m',
                                 'hist_temp', 'hist_temp_m',
                                 'hist_rho', 'hist_rho_m',
                                 'hist_vel', 'hist_vel_m',
                                 'hist_vel_inf', 'hist_vel_inf_m',
                                 'hist_theta', 'hist_theta_m',
                                 'hist_entropy, hist_entropy_m',
                                 'hist_vel_inf_bern', 'hist_vel_inf_bern_m']

        self.array_matrix = [[np.zeros(0, )
                                for x in range(len(self.list_array_names))]
                               for y in range(len(self.list_criteria))]

    def check_array_name(self, v_n):
        if not v_n in self.list_array_names:
            raise NameError("Array name: {} not in the corresponding list: {}"
                            .format(v_n, self.list_array_names))

    def i_arr(self, v_n):
        return int(self.list_array_names.index(v_n))

    def compute_array(self, v_n, criterion):

        if criterion == '':

            if v_n == 't_total_mass' or v_n == 'total_mass':

                tmp = self.get_file("dens.norm1.asc")
                t_total_mass, dens = tmp[:, 1], tmp[:, 2]
                t_total_mass *= time_constant / 1000  # [s]
                m_tot = dens * volume_constant ** 3
                self.array_matrix[self.i_cr(criterion)][self.i_arr('t_total_mass')] = t_total_mass
                self.array_matrix[self.i_cr(criterion)][self.i_arr('total_mass')] = m_tot

            elif v_n == 't_unb_mass' or v_n == 'unb_mass':

                tmp2 = self.get_file("dens_unbnd.norm1.asc")
                t_unb_mass, dens_unb = tmp2[:, 1], tmp2[:, 2]
                t_unb_mass *= time_constant / 1000
                unb_mass = dens_unb * (volume_constant ** 3)
                self.array_matrix[self.i_cr(criterion)][self.i_arr('t_unb_mass')] = t_unb_mass
                self.array_matrix[self.i_cr(criterion)][self.i_arr('unb_mass')] = unb_mass

            elif v_n == 't_unb_mass_bern' or v_n == 'unb_mass_bern':

                tmp2 = self.get_file("dens_unbnd_bernoulli.norm1.asc")
                t_unb_mass, dens_unb = tmp2[:, 1], tmp2[:, 2]
                t_unb_mass *= time_constant / 1000
                unb_mass = dens_unb * (volume_constant ** 3)
                self.array_matrix[self.i_cr(criterion)][self.i_arr('t_unb_mass_bern')] = t_unb_mass
                self.array_matrix[self.i_cr(criterion)][self.i_arr('unb_mass_bern')] = unb_mass

            elif v_n == 't_disk_mass' or v_n == 'disk_mass':
                tmp = self.get_file("disk_mass.txt")
                self.array_matrix[self.i_cr(criterion)][self.i_arr('t_disk_mass')] = tmp[:, 0]
                self.array_matrix[self.i_cr(criterion)][self.i_arr('disk_mass')] = tmp[:, 1]

        else:
            if v_n.__contains__("hist"):

                if v_n.split('_')[-1] == 'm':
                    v_n_m = v_n
                    v_n_data = self.list_array_names[int(self.list_array_names.index(v_n) - 1)]
                else:
                    v_n_m = self.list_array_names[int(self.list_array_names.index(v_n) + 1)]
                    v_n_data = v_n

                if v_n_data == "hist_vel_inf" and criterion.__contains__('_b'):
                    tmp = self.get_file("hist_vel_inf_bern" + ".dat", criterion)
                else:
                    tmp = self.get_file(v_n_data + ".dat", criterion)

                # tmp = self.get_file(v_n_data+".dat", criterion)
                # print("v_n_x: {} v_n_y: {} XXX: {}".format(v_n_data, v_n_m, tmp))
                data, mass = tmp[:, 0], tmp[:, 1]
                self.array_matrix[self.i_cr(criterion)][self.i_arr(v_n_data)] = data
                self.array_matrix[self.i_cr(criterion)][self.i_arr(v_n_m)] = mass

            elif v_n == "t_tot_flux" or v_n == "tot_flux" or v_n == "mass_tot_flux":

                tmp = self.get_file("total_flux.dat", criterion)

                time, flux, mass = tmp[:,0], tmp[:,1], tmp[:,2]
                time *= time_constant / 1000

                self.array_matrix[self.i_cr(criterion)][self.i_arr('t_tot_flux')] = time
                self.array_matrix[self.i_cr(criterion)][self.i_arr('tot_flux')] = flux
                self.array_matrix[self.i_cr(criterion)][self.i_arr('mass_tot_flux')] = mass

            else:
                raise NameError("method for array v_n: {} anc criterion: {} does not listed"
                                .format(v_n, criterion))

    def is_array_computed(self, v_n, criterion = ''):

        data = self.array_matrix[self.i_cr(criterion)][self.i_arr(v_n)]
        # print(len(data))
        if len(data) == 0:
            self.compute_array(v_n, criterion)

    def get_arr(self, v_n, criterion=''):
        self.check_criterion(criterion)
        self.check_array_name(v_n)
        self.is_array_computed(v_n, criterion)

        return self.array_matrix[self.i_cr(criterion)][self.i_arr(v_n)]


class COMPUTE_STORE_PAR(COMPUTE_STORE_ARR):

    def __init__(self, sim):

        COMPUTE_STORE_ARR.__init__(self, sim)

        self.list_parameters = ["tcoll", "Mdisk", "Munb_tot", "Munb_bern_tot",
                                "Mej_tot", "Ye_ave", "s_ave", "vel_inf_ave",
                                "vel_inf_bern_ave", "theta_rms", "E_kin_ave",
                                "E_kin_bern_ave",
                                "Mdisk3D",
                                "tcoll_gw", "tmerger_gw"]

        self.parameter_matrix = [[123456789.1
                                for x in range(len(self.list_parameters))]
                               for y in range(len(self.list_criteria))]

    def check_parameter(self, par):
        if not par in self.list_parameters:
            raise NameError("Parameter: {} not in their list: {}"
                            .format(par, self.list_parameters))

    def i_par(self, par):
        return int(self.list_parameters.index(par))

    def compute_parameter(self, v_n, criterion = ''):

        if criterion == '':
            if v_n == 'tcoll' or v_n == "Mdisk":

                t, Mtot, Munb = self.get_arr('t_total_mass'), self.get_arr('total_mass'), self.get_arr('unb_mass')
                if Mtot[-1] > 1.0:
                    Mdisk = np.nan
                    tcoll = np.inf
                    Printcolor.yellow("Warning: Disk mass at tcoll not estimated")
                else:
                    i_BH = np.argmin(Mtot > 1.0)
                    tcoll = t[i_BH]  # *1000 #* utime
                    i_disk = i_BH + int(1.0 / (t[1] * 1000))  #
                    # my solution to i_disk being out of bound:
                    if i_disk > len(Mtot): i_disk = len(Mtot) - 1
                    Mdisk = Mtot[i_disk] - Munb[i_disk]
                self.parameter_matrix[self.i_cr(criterion)][self.i_par('tcoll')] = tcoll
                self.parameter_matrix[self.i_cr(criterion)][self.i_par('Mdisk')] = Mdisk

            elif v_n == 'Munb_bern_tot':
                Munb_bern_tot = self.get_arr('unb_mass_bern')
                self.parameter_matrix[self.i_cr(criterion)][self.i_par('Munb_bern_tot')] = Munb_bern_tot[-1]

            elif v_n == 'Munb_tot':
                Munb_tot = self.get_arr('unb_mass')
                self.parameter_matrix[self.i_cr(criterion)][self.i_par('Munb_tot')] = Munb_tot[-1]

            elif v_n == 'Mdisk3D':
                self.parameter_matrix[self.i_cr(criterion)][self.i_par('Mdisk3D')] = self.get_arr("disk_mass")[-1]
            elif v_n == 'tcoll_gw':
                import scivis.units as ut  # for tmerg
                time_ = np.float(self.get_file("tcoll.dat"))
                Printcolor.yellow("Warning! using defauled M_Inf=2.748, R_GW=400.0 for retardet time")
                ret_time = PHYSICS.get_retarded_time(time_, M_Inf=2.748, R_GW=400.0)
                tcoll = ut.conv_time(ut.cactus, ut.cgs, ret_time)
                self.parameter_matrix[self.i_cr(criterion)][self.i_par('tcoll_gw')] = tcoll # s

            elif v_n == 'tmerger_gw':
                import scivis.units as ut
                time_ = np.float(self.get_file("tmerger.dat"))
                Printcolor.yellow("Warning! using defauled M_Inf=2.748, R_GW=400.0 for retardet time")
                ret_time = PHYSICS.get_retarded_time(time_, M_Inf=2.748, R_GW=400.0)
                tmerg = ut.conv_time(ut.cactus, ut.cgs, ret_time)
                self.parameter_matrix[self.i_cr(criterion)][self.i_par('tmerger_gw')] = tmerg # s

            else:
                raise NameError("No method for par:{} criterion:{}".format(v_n, criterion))

        else:

            if v_n == "Mej_tot":
                M_ej_tot = self.get_arr("mass_tot_flux", criterion)
                self.parameter_matrix[self.i_cr(criterion)][self.i_par('Mej_tot')] = M_ej_tot[-1]

            elif v_n == "Ye_ave":
                Mej = self.get_par("M_ej_tot", criterion)
                ye = self.get_arr("hist_ye", criterion)
                ye_M = self.get_arr("hist_ye_m", criterion)
                ye_ave = np.sum(ye * ye_M) / Mej
                if ye_ave > 0.6: raise ValueError("Ye_ave > 0.6 crit:{}".format(criterion))
                self.parameter_matrix[self.i_cr(criterion)][self.i_par('Ye_ave')] = ye_ave

            elif v_n == 's_ave':
                Mej = self.get_par("M_ej_tot", criterion)
                s = self.get_arr("hist_entropy", criterion)
                s_m = self.get_arr("hist_entropy_m", criterion)
                s_ave = np.sum(s * s_m) / Mej
                self.parameter_matrix[self.i_cr(criterion)][self.i_par('s_ave')] = s_ave

            elif v_n == 'vel_inf_ave' and not criterion.__contains__('_b'): # NOT bernoulli v_inf
                Mej = self.get_par("M_ej_tot", criterion)
                vel_inf = self.get_arr("hist_vel_inf", criterion)
                vel_inf_m = self.get_arr("hist_vel_inf_m", criterion)
                vel_inf_ave = np.sum(vel_inf * vel_inf_m) / Mej
                self.parameter_matrix[self.i_cr(criterion)][self.i_par('vel_inf_ave')] = vel_inf_ave

            elif v_n == 'vel_inf_ave' and criterion.__contains__('_b'): # bernoulli v_inf
                Mej = self.get_par("M_ej_tot", criterion)
                vel_inf_bern = self.get_arr("hist_vel_inf_bern", criterion)
                vel_inf_bern_m = self.get_arr("hist_vel_inf_bern_m", criterion)
                vel_inf_ave_bern = np.sum(vel_inf_bern * vel_inf_bern_m) / Mej
                # print(vel_inf_ave_bern)
                self.parameter_matrix[self.i_cr(criterion)][self.i_par('vel_inf_are')] = vel_inf_ave_bern

            elif v_n == 'E_kin_ave' and not criterion.__contains__('_b'): # NOT bernoulli v_inf
                vel_inf_ave = self.get_par("vel_inf_ave", criterion)
                vel_inf_m = self.get_arr("vel_inf_m", criterion)
                E_kin_ave = np.sum(0.5 * vel_inf_ave ** 2 *vel_inf_m) * energy_constant
                self.parameter_matrix[self.i_cr(criterion)][self.i_par('E_kin_ave')] = E_kin_ave

            elif v_n == 'E_kin_bern' and criterion.__contains__('_b'): # bernoulli v_inf
                vel_inf_bern_ave = self.get_par("vel_inf_bern_ave", criterion)
                vel_inf_bern_m = self.get_arr("vel_inf_bern_m", criterion)
                E_kin_bern_ave = np.sum(0.5 * vel_inf_bern_ave ** 2 * vel_inf_bern_m) * energy_constant
                self.parameter_matrix[self.i_cr(criterion)][self.i_par('E_kin_ave')] = E_kin_bern_ave

            elif v_n == 'theta_rms':
                theta = self.get_arr("hist_theta", criterion)
                theta_m = self.get_arr("hist_theta_m", criterion)
                theta -= pi / 2
                theta_rms = 180. / pi * sqrt(np.sum(theta_m * theta ** 2) / np.sum(theta_m))
                self.parameter_matrix[self.i_cr(criterion)][self.i_par('theta_rms')] = theta_rms

            else:
                raise NameError("No method for par:{} criterion:{}".format(v_n, criterion))

    def is_parameter_computed(self, v_n, criterion=''):

        data = self.parameter_matrix[self.i_cr(criterion)][self.i_par(v_n)]
        if data == 123456789.1:
            self.compute_parameter(v_n, criterion)

    def get_par(self, v_n, criterion=''):

        self.check_parameter(v_n)
        self.check_criterion(criterion)
        self.is_parameter_computed(v_n, criterion)

        return self.parameter_matrix[self.i_cr(criterion)][self.i_par(v_n)]


class LOAD_NUCLEO:

    def __init__(self, sim):

        self.gen_set = {
            "outflow": "outflow",
            "load_file": "yields.h5",
            "solar_file": Paths.skynet + Files.solar_r
        }

        self.sim = sim

        self.list_criteria = ['', "_0", "_1", "_0_b", "_1_b", "_0_b_w", "_1_b_w", "_0 _0_b_w"]


        self.list_sims_v_ns = ["A", "Y_final", "Z"]
        self.list_sol_v_ns = ["Asun", "Ysun"]
        self.list_v_ns = self.list_sims_v_ns + self.list_sol_v_ns


        self.in_data_matrix = [[np.zeros(0,)
                               for x in range(len(self.list_sims_v_ns))]
                               for y in range(len(self.list_criteria))]

        self.sol_data_matrix = [np.zeros(0,)
                               for x in range(len(self.list_sol_v_ns))]

    def check_criterion(self, criterion):
        if not criterion in self.list_criteria:
            raise NameError("criteria: {} not in the list of them: {}"
                            .format(criterion, self.list_criteria))

    def check_sim_v_n(self, v_n):
        if not v_n in self.list_sims_v_ns:
            raise NameError("fname: {} not in the lit of sim v_ns nucleo file: {}"
                            .format(v_n, self.list_sims_v_ns))

    def check_sol_v_n(self, v_n):
        if not v_n in self.list_sol_v_ns:
            raise NameError("fname: {} not in the lit of solar file v_ns file: {}"
                            .format(v_n, self.list_sol_v_ns))

    def i_cr(self, criterion):
        return int(self.list_criteria.index(criterion))

    def i_sim_v_n(self, v_n):
        return int(self.list_sims_v_ns.index(v_n))

    def i_sol_v_n(self, v_n):
        return int(self.list_sol_v_ns.index(v_n))

    def load_and_sum_files(self, criteria):

        assert len(criteria.split(' ')) > 1
        all_As = []
        all_Ys = []
        all_Z = []
        Ys = np.zeros(1)

        for criterion in criteria.split(' '):
            As = self.get_sim_data("A", criterion)
            Ys = self.get_sim_data("Y_final", criterion)
            Z = self.get_sim_data("Z", criterion)
            all_As = As
            all_Z = Z
            all_Ys.append(Ys)

        # print(Ys.shape)

        all_Ys = np.reshape(all_Ys, (len(criteria.split(' ')), len(all_As)))
        all_Ys = np.sum(all_Ys, axis=0)

        assert Ys.shape == all_Ys.shape

        self.in_data_matrix[self.i_cr(criteria)][self.i_sim_v_n("A")] = all_As
        self.in_data_matrix[self.i_cr(criteria)][self.i_sim_v_n("Z")] = all_Z
        self.in_data_matrix[self.i_cr(criteria)][self.i_sim_v_n("Y_final")] = all_Ys



    def load_single_file(self, criterion):

        fpath = Paths.ppr_sims + self.sim + '/' + self.gen_set['outflow'] + \
                criterion + '/' + self.gen_set["load_file"]
        h5tbl = h5py.File(fpath, "r")

        for v_n in h5tbl:
            if not v_n in self.list_sims_v_ns:
                raise NameError("Found unrecognized v_n:{} in file:{}"
                                .format(v_n, fpath))
            else:
                data = np.array(h5tbl[v_n])
                self.in_data_matrix[self.i_cr(criterion)][self.i_sim_v_n(v_n)] = data

    def load_sim_file(self, criterion):

        if len(criterion.split(' ')) > 1:
            self.load_and_sum_files(criterion) # combines Yields
        else:
            self.load_single_file(criterion)






    def is_sim_loaded(self, v_n, criterion):
        data = self.in_data_matrix[self.i_cr(criterion)][self.i_sim_v_n(v_n)]
        if len(data) == 0:
            self.load_sim_file(criterion)

        data = self.in_data_matrix[self.i_cr(criterion)][self.i_sim_v_n(v_n)]
        if len(data) == 0:
            raise ValueError("failed to load simulation data: v_n:{} crit:{}"
                             .format(v_n, criterion))


    def load_sol_file(self):

        fpath = self.gen_set["solar_file"]
        Asun, Ysun = np.loadtxt(fpath, unpack=True)

        self.sol_data_matrix[self.i_sol_v_n("Asun")] = Asun
        self.sol_data_matrix[self.i_sol_v_n("Ysun")] = Ysun


    def is_sol_loaded(self, v_n):
        data = self.sol_data_matrix[self.i_sol_v_n(v_n)]
        if len(data) == 0:
            self.load_sol_file()

        data = self.sol_data_matrix[self.i_sol_v_n(v_n)]
        if len(data) == 0:
            raise ValueError("Loading solar data failed. No data has been loaded.")

    def get_sim_data(self, v_n, criterion):
        self.check_criterion(criterion)
        self.check_sim_v_n(v_n)
        self.is_sim_loaded(v_n, criterion)

        return self.in_data_matrix[self.i_cr(criterion)][self.i_sim_v_n(v_n)]


    def get_sol_data(self, v_n):

        self.check_sol_v_n(v_n)
        self.is_sol_loaded(v_n)

        return self.sol_data_matrix[self.i_sol_v_n(v_n)]


class NORMALIZE_NUCLEO(LOAD_NUCLEO):

    def __init__(self, sim):

        LOAD_NUCLEO.__init__(self, sim)

        self.list_norm_methods = ["sum", "Asol=195"]

        self.normed_sim_data_list = [[[np.zeros(0,)
                                     for x in range(len(self.list_sims_v_ns))]
                                     for y in range(len(self.list_criteria))]
                                     for z in range(len(self.list_norm_methods))]

        self.normed_sol_data_list = [[np.zeros(0, )
                                      for x in range(len(self.list_sol_v_ns))]
                                      for z in range(len(self.list_norm_methods))]

    def check_method(self, method):
        if not method in self.list_norm_methods:
            raise NameError("method:{} not in the list of normalisation methods: {}"
                            .format(method, self.list_norm_methods))



    def i_meth(self, method):
        return int(self.list_norm_methods.index(method))

    def compute_normalized_sol_yields(self, method='sum'):

        As = self.get_sol_data("Asun")
        Ys = self.get_sol_data("Ysun")

        '''Sums all Ys for a given A (for all Z)'''
        Anrm = np.arange(As.max() + 1)
        Ynrm = np.zeros(int(As.max()) + 1)
        for i in range(Ynrm.shape[0]):  # changed xrange to range
            Ynrm[i] = Ys[As == i].sum()

        if method == 'sum':
            Ynrm /= np.sum(Ynrm)
            return Anrm, Ynrm
        else:
            raise NameError("Normalisation method '{}' for the solar is not recognized. Use:{}"
                            .format(method, self.list_norm_methods))

    def is_sol_computed(self, v_n, method):

        data = self.normed_sol_data_list[self.i_meth(method)][self.i_sol_v_n(v_n)]
        if len(data) == 0:
            a_sol, y_sol = self.compute_normalized_sol_yields()
            self.normed_sol_data_list[self.i_meth(method)][self.i_sol_v_n("Asun")] = a_sol
            self.normed_sol_data_list[self.i_meth(method)][self.i_sol_v_n("Ysun")] = y_sol



    def compute_normalized_sim_yelds(self, criterion, method):

        As = self.get_sim_data("A", criterion)
        Ys = self.get_sim_data("Y_final", criterion)

        '''Sums all Ys for a given A (for all Z)'''
        Anrm = np.arange(As.max() + 1)
        Ynrm = np.zeros(int(As.max()) + 1)
        for i in range(Ynrm.shape[0]):  # changed xrange to range
            Ynrm[i] = Ys[As == i].sum()

        if method == '':
            return Anrm, Ynrm

        elif method == 'sum':
            ''' Normalizes to a sum of all A '''
            norm = Ynrm.sum()
            Ynrm /= norm
            return Anrm, Ynrm

        elif method == "Asol=195":
            ''' Normalize to the solar abundances of a given element'''
            # a_sol = self.get_sol_data("Asun")
            # y_sol = self.get_sol_data("Ysun")
            a_sol = self.get_normalized_sol_data("Asun")
            y_sol = self.get_normalized_sol_data("Ysun")

            element_a = int(method.split("=")[-1])
            if element_a not in a_sol: raise ValueError('Element: a:{} not in solar A\n{}'.format(element_a, a_sol))
            if element_a not in Anrm: raise ValueError('Element: a:{} not in a_arr\n{}'.format(element_a, Anrm))

            delta = np.float(y_sol[np.where(a_sol == element_a)] / Ynrm[np.where(Anrm == element_a)])
            Ynrm *= delta

            return Anrm, Ynrm
        else:
            raise NameError("Normalisation method '{}' for the simulation yields is not recognized. Use:{}"
                            .format(method, self.list_norm_methods))

    def is_sim_compted(self, v_n, criterion = '', method = ''):

        data = self.normed_sim_data_list[self.i_meth(method)][self.i_cr(criterion)][self.i_sim_v_n(v_n)]
        if len(data) == 0:
            a_arr, y_arr = self.compute_normalized_sim_yelds(criterion, method)
            self.normed_sim_data_list[self.i_meth(method)][self.i_cr(criterion)][self.i_sim_v_n("A")] = a_arr
            self.normed_sim_data_list[self.i_meth(method)][self.i_cr(criterion)][self.i_sim_v_n("Y_final")] = y_arr

    def get_normalized_sim_data(self, v_n, criterion='', method=''):

        self.check_sim_v_n(v_n)
        self.check_criterion(criterion)
        self.check_method(method)
        self.is_sim_compted(v_n, criterion, method)

        return self.normed_sim_data_list[self.i_meth(method)][self.i_cr(criterion)][self.i_sim_v_n(v_n)]

    def get_normalized_sol_data(self, v_n, method='sum'):

        self.check_sol_v_n(v_n)
        self.check_method(method)
        self.is_sol_computed(v_n, method)

        return self.normed_sol_data_list[self.i_meth(method)][self.i_sol_v_n(v_n)]

    def get_normalized_sim_data_miltiple(self, v_n, criteria=[], method=''):

        all_As = []
        all_Ys = []
        for criterion in criteria:
            As = self.get_sim_data("A", criterion)
            Ys = self.get_sim_data("Y_final", criterion)
            all_As = As
            all_Ys.append(Ys)

        all_Ys = np.reshape(all_Ys, (len(criteria), len(all_As)))
        all_Ys = np.sum(all_Ys, axis=0)

        assert len(all_As) == len(all_Ys)








if __name__ == '__main__':

    "--- --- TEST --- ---"

    o_nuc = NORMALIZE_NUCLEO("DD2_M13641364_M0_SR")
    print(o_nuc.get_normalized_sim_data("Y_final", "_0_b_w", "Asol=195"))


    # o_l_outflowed = LOAD_FILES("SLy4_M13641364_M0_SR")
    # print(o_l_outflowed.get_file("hist_ye.dat", "_0_b_w"))
    # print(o_l_outflowed.get_file("dens_unbnd.norm1.asc"))

    o_par = COMPUTE_STORE_PAR("DD2_M13641364_M0_SR")

    print(o_par.get_arr('hist_vel_inf_m', '_0_b_w')); exit(1)

    # print(o_par.get_arr("hist_ye", "_0_b_w")); exit(1)
    print(o_par.get_par("Munb_tot"))

    print(o_par.get_arr("hist_vel_inf_m", '_0_b_w')); exit(1)
    print(o_par.get_par("Mej_tot", "_0"))
    print(o_par.get_par("Mej_tot", "_0_b_w"))
    print(o_par.get_par("theta_rms", "_0_b_w"))
    # print(o_par.get_par("Mdisk3D", ""))
    # print(o_par.get_par("tcoll_gw", ""))
    # print(o_par.get_par("tmerger_gw", ""))
