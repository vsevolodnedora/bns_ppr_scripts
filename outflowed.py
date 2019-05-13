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

from scidata.carpet import hdf5

class SIM_DISK:

    def __init__(self, sim):

        self.sim = sim

        start_t = time.time()
        print("\t Loading Disk [{}]...".format(sim))
        self.time_disk, self.disk_mass = self.load_disk_mass_ev()
        print("\t done! (%.2f sec)" % (time.time() - start_t))

    def old_load_disk_mass_ev(self, do_remove_nans=True):

        try:
            t_disk_evol, m_disk_evol = np.loadtxt(MakePath.collated2(self.sim) + Files.disk_mass,
                                                  usecols=[1, 2], unpack=True)  # m in solar masses
        except IOError:
            Printcolor.yellow("Warning: Disk mass evolution files not found")
            return np.empty(0, ), np.empty(0, )

        t_disk_evol *= time_constant / 1000 # time in [s]

        if do_remove_nans:
            m_disk_evol = m_disk_evol[~np.isnan(m_disk_evol)]
            t_disk_evol = t_disk_evol[~np.isnan(m_disk_evol)]
            if len(m_disk_evol) != len(t_disk_evol):
                raise ValueError('Error in removing nans from disk')
        return t_disk_evol, m_disk_evol # time in [s]

    def load_disk_mass_ev(self):

        # check if res_3d exists:
        if not os.path.isdir(Paths.ppr_sims + self.sim + "/res_3d/"):
            return np.zeros(0,), np.zeros(0,)


        # check if any disk_mass.txt exist:
        it_dirs = os.listdir(Paths.ppr_sims + self.sim + "/res_3d/")
        iterations, disk_masses = [], []
        if len(it_dirs) == 0:
            Printcolor.yellow("No iteration directories found in {}".format(self.sim + "/res_3d/"))
        else:
            for it_dir in it_dirs:
                path = Paths.ppr_sims + self.sim + "/res_3d/" + it_dir + '/' + "disk_mass.txt"
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
            print(np.array(iterations).shape, np.array(times).shape, np.array(disk_masses).shape)

            it_time_masses = np.vstack((iterations, times, disk_masses)).T

        # print(it_time_masses)


        final_times, final_masses = x_y_z_sort(it_time_masses[:,1], it_time_masses[:,2],np.zeros(0,), 0)

        return final_times, final_masses

    def get_par(self, v_n):
        if v_n == 'Mdisk_last':
            if self.disk_mass.any():
                return self.disk_mass[-1]
            else:
                return np.nan

class SIM_UNBOUND:

    def __init__(self, sim):


        start_t = time.time()
        print("\t Loading Unbound [{}]...".format(sim))

        self.sim = sim

        self.time_total_mass, self.total_mass = self.load_total_mass()

        self.time_unbound_mass, self.unb_mass = \
            self.load_unbound_mass(MakePath.collated(sim) + Files.dens_unb)

        self.time_unbound_mass_bern, self.unb_mass_bern \
            = self.load_unbound_mass(MakePath.collated(sim) + Files.dens_unb_bern)

        # self.time_disk, self.disk_mass = self.load_disk_mass_ev()

        self.tcoll, self.diskmass_tcol = self.get_t_coll_and_disk_mass()

        print("\t done! (%.2f sec)" % (time.time() - start_t))

    def load_unbound_mass(self, fname):
        try:
            t_unb_mass, dens_unb = np.loadtxt(fname, usecols=[1, 2], unpack=True)
            t_unb_mass *= time_constant / 1000
            unb_mass = dens_unb * (volume_constant ** 3)
            return t_unb_mass, unb_mass
        except IOError:
            Printcolor.yellow("Warning: total unbound mass was not loaded \n"
                              "File {} not found".format(fname))
            return np.zeros(0,), np.zeros(0,) # [time in s]

    def load_total_mass(self):
        try:
            t_total_mass, dens = np.loadtxt(MakePath.collated(self.sim) + Files.dens,
                                            usecols=[1, 2], unpack=True)
            t_total_mass *= time_constant / 1000 #[s]
            m_tot = dens * volume_constant ** 3
            return t_total_mass, m_tot
        except IOError:
            Printcolor.yellow("Warning: total mass was not loaded \n"
                              "File {} not found".format(MakePath.collated(self.sim) + Files.dens))
            return np.zeros(0,), np.zeros(0,) # [time in s]

    def get_t_coll_and_disk_mass(self):
        '''
        Computes mass of the disk as a M_total - M_unbound at t = t_collapse
        '''

        t, Mtot, Munb = self.time_total_mass, self.total_mass, self.unb_mass

        if Mtot[-1] > 1.0:
            Mdisk = np.nan
            tcoll = np.inf
            Printcolor.yellow("Warning: Disk mass at tcoll not estimated")
        else:
            i_BH = np.argmin(Mtot > 1.0)
            tcoll = t[i_BH]#*1000 #* utime
            i_disk = i_BH + int(1.0 / (t[1]*1000))#
            # my solution to i_disk being out of bound:
            if i_disk > len(Mtot): i_disk = len(Mtot)-1
            Mdisk = Mtot[i_disk] - Munb[i_disk]

        return tcoll, Mdisk

    def get_par(self, v_n):

        if v_n == 'tcoll':
            return self.tcoll

        elif v_n == 'Mdisk':
            return self.diskmass_tcol

        elif v_n == 'Munb_bern_tot':
            return self.unb_mass_bern[-1]

        elif v_n == 'Munb_tot':
            return self.unb_mass[-1]
        else:
            raise NameError("Unknown v_n:{}".format(v_n))

class SIM_EJ_HIST:

    def __init__(self, sim, extension='_0', norm=True):

        self.sim = sim
        self.outflow_path = MakePath.outflow(sim, extension)

        # start_t = time.time()
        # print("\t Loading Ejecta files [{}]...".format(sim))


        self.time_total_flux, self.total_flux, self.mass_total_flux = \
            self.load_total_flux(self.outflow_path + Files.total_flux)

        #
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
        #
        # SIM.__init__(self, sim)
        #
        # '''---- SET WHAT TO LOAD ----'''
        # self.set_load_tot_flux      = True # surface integral of the material traversing a given surface and having -u_t < 1
        # self.set_load_theta_75      = False #
        # self.set_load_hist_theta    = True # for theta_rms estimation
        #
        # '''---------------------------'''
        #
        # '''---- SET FILE NAMES AND LOCATION ----'''
        # total_flux = 'total_flux.dat'
        # self.tot_flux_file = self.path_to_outflow_dir + total_flux
        #
        # hist_theta = 'hist_theta.dat'
        # self.hist_theta_file = self.path_to_outflow_dir + hist_theta
        #
        # '''-------------------------------------'''
        #
        # '''----------- LOAD SELECTED FILES -----------'''
        # if self.set_load_tot_flux:
        #     self.time_total_flux, self.total_flux, self.mass_total_flux = self.load_total_flux(self.tot_flux_file)
        #
        # if self.set_load_theta_75:
        #     self.theta_75 = np.loadtxt(self.path_to_outflow_dir + 'theta_75.dat', usecols=[0, 0], unpack=True)[0]
        #
        # if self.set_load_hist_theta:
        #     self.theta_rms = self.load_theta_rms()

        # standard histograms from outflowed.cc
        norm = False
        hist_v_ns = ['ye', 'temp', 'rho', 'vel', 'vel_inf', 'theta', 'entropy']
        self.ye, self.ye_M          = self.load_hist_v_n('ye', norm)
        self.temp, self.temp_M      = self.load_hist_v_n('temperature', norm)
        self.rho, self.rho_M        = self.load_hist_v_n('log_rho', norm)
        self.vel, self.vel_M        = self.load_hist_v_n('vel', norm)
        self.vel_inf, self.vel_inf_M= self.load_hist_v_n('vel_inf', norm)
        self.theta, self.theta_M    = self.load_hist_v_n('theta', norm)
        self.entropy, self.entropy_M= self.load_hist_v_n('entropy', norm)
        self.theta_rms = self.load_theta_rms()
        try:
            self.vel_inf_bern, self.vel_inf_M_bern = self.load_hist_v_n('vel_inf_bern', norm)
        except ValueError:
            Printcolor.yellow("vel_inf_bern not found for {} {}".format(sim, extension))

    @staticmethod
    def normalize_histogram(M):
        # M[1]  += M[0]
        # M[-1] += M[-2]
        # M[0] = M[-1] = 0
        M /= np.sum(M)
        return M

    def load_total_flux(self, fname):
        try:
            time, flux, mass = np.loadtxt(fname, usecols=[0, 1, 2], unpack=True)
            time *= time_constant / 1000 #[s]
        except IOError:
            Printcolor.yellow("Warning: total flux was not loaded \n"
                              "File {} not found".format(fname))
            return np.zeros(0,), np.zeros(0,), np.zeros(0,0) # [time in s]
        return time, flux, mass

    def load_hist_v_n(self, v_n, norm=True):

        fname =  "hist" + '_' + v_n + '.dat'
        fpath = self.outflow_path + fname

        start_t = time.time()
        print("\t Loading {} ...".format(fname)),
        try:
            s, M = np.loadtxt(fpath, usecols=(0, 1), unpack=True)
        except:
            raise ValueError("Hist File not found: {}".format(fpath))

        if norm:
            M = self.normalize_histogram(M)
        print(" done! (%.2f sec)" % (time.time() - start_t))

        return s, M

    def get_ye_average(self):

        Mej = self.mass_total_flux[-1]
        # Ye, M = np.loadtxt(self.outflow_path + Files.hist_ye, usecols=(0, 1), unpack=True)
        Ye_avg = np.sum(self.ye * self.ye_M) / Mej
        if Ye_avg > .6: print(Ye_avg); exit(1)
        return Ye_avg

    def get_entropy_average(self):

        Mej = self.mass_total_flux[-1]

        # s, M = np.loadtxt(self.outflow_path + Files.hist_entropy, usecols=(0, 1), unpack=True)
        s_avg = np.sum(self.entropy * self.entropy_M) / Mej

        return s_avg

    def get_average_vel_at_infinity(self):

        Mej = self.mass_total_flux[-1]

        # vel_inf, M = np.loadtxt(self.outflow_path + Files.hist_vel_inf, usecols=(0, 1), unpack=True)
        vel_inf_avg = np.sum(self.vel_inf * self.vel_inf_M) / Mej
        # vel_inf_avg *= 10 # [10^-1 c]

        return vel_inf_avg

    def get_average_vel_at_infinity_bern(self):

        Mej = self.mass_total_flux[-1]

        # vel_inf, M = np.loadtxt(self.outflow_path + Files.hist_vel_inf, usecols=(0, 1), unpack=True)
        vel_inf_avg = np.sum(self.vel_inf_bern * self.vel_inf_M_bern) / Mej
        # vel_inf_avg *= 10 # [10^-1 c]

        return vel_inf_avg

    def get_average_kin_energy(self):

        # _, M = np.loadtxt(self.outflow_path + Files.hist_entropy, usecols=(0, 1), unpack=True)
        vel_inf = self.get_average_vel_at_infinity()
        Ekin = np.sum(0.5 * vel_inf ** 2 * self.vel_inf_M) * energy_constant # 10^51 erg

        return Ekin

    def get_average_kin_energy_bern(self):

        # _, M = np.loadtxt(self.outflow_path + Files.hist_entropy, usecols=(0, 1), unpack=True)
        vel_inf = self.get_average_vel_at_infinity_bern()
        Ekin = np.sum(0.5 * vel_inf ** 2 * self.vel_inf_M_bern) * energy_constant # 10^51 erg

        return Ekin

    def load_theta_rms(self):

        theta_hist = np.loadtxt(self.outflow_path + Files.hist_theta, unpack=True)
        theta_hist[0][:] -= pi / 2
        theta_rms = 180. / pi * sqrt(np.sum(theta_hist[1] * theta_hist[0] ** 2) / np.sum(theta_hist[1]))

        return theta_rms

    def get_hist_s_M(self, v_n, norm):

        if v_n == 'ye' or v_n == 'Ye':              return self.ye, self.ye_M
        elif v_n == 'temp' or v_n == 'temperature': return self.temp, self.temp_M
        elif v_n == 'rho' or v_n == 'log_rho':      return self.rho, self.rho_M
        elif v_n == 'vel' or v_n == 'velocity':     return self.vel, self.vel_M
        elif v_n == 'vel_inf':                      return self.vel_inf, self.vel_inf_M
        elif v_n == 'vel_inf_bern':                 return self.load_hist_v_n('vel_inf_bern', norm)
        elif v_n == 'theta':                        return self.theta, self.theta_M
        elif v_n == 'entropy' or v_n == 's':        return self.entropy, self.entropy_M
        else: raise NameError("Variable {} is not declared in the histogram list".format(v_n))

    def get_par(self, v_n):

        if v_n == 'sim':
            return self.sim

        elif v_n == 'Mej_tot':
            return self.mass_total_flux[-1]

        elif v_n == 'dMej_tot':

            dm = np.diff(self.mass_total_flux)
            dt = np.diff(self.time_total_flux)
            dmdt = dm / dt
            # res = dmdt.max() / dmdt[-1]
            return dmdt[-1]

        elif v_n == 'Ye':
            return self.get_ye_average()

        elif v_n == 's':
            return self.get_entropy_average()

        elif v_n == 'vel_inf':
            return self.get_average_vel_at_infinity()

        elif v_n == 'vel_inf_bern':
            return self.get_average_vel_at_infinity_bern()

        elif v_n == 'theta_rms':
            return self.theta_rms

        elif v_n == 'E_kin':
            return self.get_average_kin_energy()

        elif v_n == 'E_kin_bern':
            return self.get_average_kin_energy_bern()


        else: raise NameError('v_n:{} is not available.'.format(v_n))






    # def get_hist_s_M(self, v_n, norm):
    #
    #     if v_n == 'ye' or v_n == 'Ye': return self.load_hist_v_n('ye', norm)
    #     elif v_n == 'temp' or v_n == 'temperature': return self.load_hist_v_n('temperature', norm)
    #     elif v_n == 'rho' or v_n == 'logrho' or v_n == 'log_rho': return self.load_hist_v_n('log_rho', norm)
    #     elif v_n == 'vel' or v_n == 'velocity': return self.load_hist_v_n('vel', norm)
    #     elif v_n == 'vel_inf': return self.load_hist_v_n('vel_inf', norm)
    #     elif v_n == 'vel_inf_bern': return self.load_hist_v_n('vel_inf_bern', norm)
    #     elif v_n == 'theta': return self.load_hist_v_n('theta', norm)
    #     elif v_n == 'entropy' or v_n == 's': return self.load_hist_v_n('entropy', norm)
    #     else:
    #         raise NameError("Variable {} is not declared in the histogram list".format(v_n))

class SIM_CORR:

    def __init__(self, sim, extenstion='_0'):

        self.sim = sim
        self.outflowpath = MakePath.outflow(sim, extenstion)

    def load_corr(self, v_n1, v_n2, norm=True):

        start_t = time.time()

        fname = "corr_" + v_n1 + '_' + v_n2 + '.h5'
        fname_rev = "corr_" + v_n2 + '_' + v_n1 + '.h5'
        path = self.outflowpath

        fpath = path + fname
        fpath_rev = path + fname_rev


        if os.path.isfile(fpath):
            print("\t Loading {} ...".format(fname)),
            dfile = h5py.File(fpath, "r")
            var1 = np.array(dfile[v_n1])
            var2 = np.array(dfile[v_n2])
            mass = np.array(dfile["mass"])
            print(" done! (%.2f sec)" % (time.time() - start_t))
            dfile.close()
        elif os.path.isfile(fpath_rev):
            print("\t Loading {} ...".format(fname)),
            dfile = h5py.File(fpath_rev, "r")
            var1 = np.array(dfile[v_n1])
            var2 = np.array(dfile[v_n2])
            mass = np.array(dfile["mass"]).T
            print(" done! (%.2f sec)" % (time.time() - start_t))
            dfile.close()
        else:
            raise NameError("NOT FOUND: {} or {} in {}".format(
                "corr_{}_{}.h5".format(v_n1, v_n2),
                "corr_{}_{}.h5".format(v_n2, v_n1),
                path
            ))

        # var1 = np.array(dfile[v_n1])
        # var2 = np.array(dfile[v_n2])
        # mass = np.array(dfile["mass"])

        if norm:
            mass = mass / np.sum(mass)
            mass = 2 * np.maximum(mass, 1e-15) # WHAT'S THAT?

        return var1, var2, mass

    def get_corr(self, v_n1, v_n2, norm=True):

        if v_n1 == 'vel_inf' and v_n2 == 'theta': return self.load_corr('vel_inf', 'theta', norm)
        elif v_n1 == 'ye' and v_n2 == 'entropy':  return self.load_corr('ye', 'entropy', norm)
        elif v_n1 == 'ye' and v_n2 == 'theta':    return self.load_corr('ye', 'theta', norm)
        else:
            raise NameError("Vars v_n1 {} and v_n2 {} are not declared in the correlation list".format(v_n1, v_n2))

class SIM_PROF:

    def __init__(self, sim, extension):

        self.sim = sim

        self.outflowpath = MakePath.outflow(sim, extension)


    @staticmethod
    def modify_arr(arr, v_n):
        if v_n == 'theta':
            return (180 * arr / np.pi)  # [::-1] # inverting to to get angle from the plane of the binary
        else:
            return arr

    @staticmethod
    def convert_raw_data(mdataset):
        '''
        Returns time [s], theta [radians] 2D data
        :param mdataset:
        :return:
        '''

        time_   = np.array(mdataset.time) * time_constant * 1e-3
        theta   = mdataset.frame(0).data_x
        rawdata = mdataset.data_y.reshape((mdataset.nframes, theta.shape[0]))

        return time_, theta, rawdata

    def load_profile(self, v_n):
        '''
        Returns time [s], theta [radians], DATA 2D array
        :param v_n:
        :return:
        '''

        mdataset = xg.parsefile(self.outflowpath + 'profile_{}.xg'.format(v_n))
        time, theta, arr = self.convert_raw_data(mdataset)

        # arr = np.maximum(arr, 1e-15)  # WHAT'S THAT?

        return time, theta, arr

###############################

class PLOT_PROFS:

    def __init__(self, sim):

        self.set_use_norm = True
        self.invert_theta = True
        self.sim = sim

        self.flux_dic = {'v_n': ['flux'],  #   'flux'
                         'extensions': ['_0'],  #  _0_b_w
                         'labels': ['Geo'],
                         'colors': ['black'],
                         'add_hist': True,
                         'hist_yscale': 'log',
                         'norm': LogNorm(1e-15, 1e-10),
                         'figdir': Paths.ppr_sims + sim + '/res_1d/'
                         }

        self.ye_dic = {'v_n': ['ye', 'ye'],  #   corr_vn1
                       'criteria': ['geo', 'bern wind'],  # /_b_w/corr_vn1_vn2
                       'det': [0, 0],  # /outflowed_0_b_w/corr_vn1_vn2
                       'labels': ['Geo', 'Bern.'],
                       'colors': ['black', 'red'],
                       'add_hist': True,
                       'hist_yscale': 'log',
                       'norm': Normalize(0., 0.5),
                       'figdir': Paths.ppr_sims + sim + '/res_1d/'
                       }

    @staticmethod
    def add_plot_background2(ax, arr_x, arr_y, mass, v_n_z, norm):

        im = ax.pcolormesh(arr_x, arr_y, mass, norm=norm, cmap='RdBu_r')
        im.set_rasterized(True)

        divider1 = make_axes_locatable(ax)
        cax = divider1.append_axes("right", size="5%", pad=0.05)
        cbar = plt.colorbar(im, cax=cax, norm=norm)
        cbar.set_label(LABELS.labels(v_n_z))
        cbar.ax.minorticks_off()

    def modify_arr(self, arr, v_n):
        if v_n == 'theta':
            if self.invert_theta:
                return 90 - (180 * arr / np.pi)# [::-1] # inverting to to get angle from the plane of the binary
            else:
                return (180 * arr / np.pi)  # [::-1] # inverting to to get angle from the plane of the binary
        else: return arr

    @staticmethod
    def plot_hist_theta(ax, x_arr, val_arr, color, label):

        ax.step(x_arr, np.sum(val_arr, axis=0), where='mid',  color=color, label=label)

        # mdens = np.sum(mass, axis=0)
        # fit = np.sin(theta / 180 * np.pi) ** 2
        # fit *= np.sum(mdens) / np.sum(fit)
        # ax.plot(x_arr, fit, 'k--', label=r'$\sin^2(\x_arr)$')       # fit
        # ax.step(0.5 * (x_arr[1:] + x_arr[:-1]), mdens, where='mid', color=color, label=label)

    def plot_from_dic(self, dic):

        n_rows = len(dic['v_n'])
        if dic['add_hist']: n_rows += 1

        n_cols = 1

        fig = plt.figure(figsize=(4.5, 2.5 * 3.6))  # (<->; v)
        axs = []
        for n in range(1, n_rows + 1):
            if n == 1:
                axs.append(fig.add_subplot(n_rows, n_cols, n))
            else:
                axs.append(fig.add_subplot(n_rows, n_cols, n))  # sharex=axs[n - 2]))

        time_arrs = []
        theta_arrs = []
        value_arrs = []
        # time_arrs.append(0)
        # theta_arrs.append(0)
        # value_arrs.append(0)

        # fig = plt.figure()
        # ax = fig.add_subplot(111)

        gw = SIM_GW(self.sim)
        tmerg = gw.load_tmerg()
        if np.isnan(tmerg):
            Printcolor.yellow("Warning: Merger time is not loaded")
            tmerg = 0.

        for i in range(len(dic['v_n'])):
            prof_cl = SIM_PROF(self.sim, dic['extensions'][i])
            time_arr, theta_arr, val_arr = prof_cl.load_profile(dic['v_n'][i])
            theta_arr = self.modify_arr(theta_arr, 'theta')
            # mass = np.array(mass[::-1, ...]).T
            val_arr = val_arr.transpose()
            time_arr = (time_arr-tmerg) * 1e3
            time_arrs.append(time_arr)
            theta_arrs.append(theta_arr)

            # im = ax.pcolormesh(x_arr, y_arr, mass, norm=dic['norm'])
            # plt.savefig('{}{}.png'.format(LISTS.fig_dir, 'x'), bbox_inches='tight', dpi=128)
            # plt.close()
            # exit(1)

            value_arrs.append(val_arr)
            del val_arr
            del time_arr
            del theta_arr

        tmin = np.hstack((time_arrs)).min()
        tmax = np.hstack((time_arrs)).max()

        i_plot = 0
        if dic['add_hist']:

            divider1 = make_axes_locatable(axs[i_plot])
            cax = divider1.append_axes("right", size="5%", pad=0.05)
            cax.xaxis.set_ticklabels([])
            cax.get_xaxis().set_visible(False)
            cax.get_yaxis().set_visible(False)
            # fig.patch.set_visible(False)
            cax.axis('off')
            # cb.outline.set_visible(False)

            # plot sin^2 fit
            for j in range(len(dic['v_n'])):
                self.plot_hist_theta(axs[i_plot], time_arrs[j], value_arrs[j], dic['colors'][j], dic['labels'][j])


            axs[i_plot].legend(loc='lower right')
            axs[i_plot].set_yscale(dic['hist_yscale'])
            # axs[i_plot].set_ylim(1e-3, 1e-1)
            # ax1.yaxis.tick_right()
            # ax1.yaxis.tick_left()
            axs[i_plot].set_ylabel(r'$\sum{\dot{M}}$')
            axs[i_plot].xaxis.set_ticklabels([])
            axs[i_plot].set_xlim(tmin, tmax)
            # if dic['v_n1'][0] == 'theta':
            #     axs[i_plot].get_yaxis().set_label_coords(-0.15, 0.5)
            #     axs[i_plot].set_xlim(0, 90)
                # axs[i_plot].invert_xaxis()

            i_plot += 1

        for n in range(len(dic['v_n'])):

            # if dic['v_n1'][n] == 'theta':
            #     axs[i_plot].get_yaxis().set_label_coords(-0.15, 0.5)
            #     axs[i_plot].set_xlim(0, 90)
                # axs[i_plot].invert_xaxis()
            axs[i_plot].set_ylim(0, 90)
            axs[i_plot].set_xlim(tmin, tmax)
            self.add_plot_background2(axs[i_plot], time_arrs[n], theta_arrs[n], value_arrs[n],
                                      dic['v_n'][n], dic['norm'])

            axs[i_plot].text(0.05, 0.85, dic['labels'][n], color='black', transform=axs[i_plot].transAxes)
            axs[i_plot].set_ylabel('Angle from orbital plane')
            if axs[i_plot] == axs[-1]:
                if tmerg > 0.: axs[i_plot].set_xlabel(LABELS.labels('t-tmerg'))
                else: axs[i_plot].set_xlabel(LABELS.labels('time'))
            i_plot += 1

        plt.subplots_adjust(hspace=0.1)

        # if v_n1 == 'x_arr': plt.xlim(right=90)

        if not os.path.exists(dic["figdir"]):
            os.makedirs(dic["figdir"])

        plt.savefig('{}{}.png'.format(dic["figdir"], 'prof_{}'.format(dic['v_n'][0])), bbox_inches='tight', dpi=128)
        plt.close()

    def plot_vn1_for_sim(self, vn1):

        if vn1 == 'flux':
            self.plot_from_dic(self.flux_dic)
        elif vn1 == 'ye':
            self.plot_from_dic(self.ye_dic)
        else:
            raise NameError('Not dic. for prof_{} available'.format(vn1))

class PLOT_HISTS:

    def __init__(self, sim):

        self.set_use_norm = False
        # self.sim = sim

        self.gen_set = {'figdir': Paths.ppr_sims + sim + '/res_1d/',
                        'figname': 'tst_hists',
                        'figsize':[9, 3.8]}

        self.ye_dic = {
                      'sims':[sim, sim, sim],
                      'v_ns':  ['ye', 'ye', 'ye'],               # var names to load the files
                      'extensions':   ['_0', '_0_b_w', '_0_b'],     # directories to use (_b or _b_d ... )
                      'colors':     ['black', 'red', 'orange'],       # colors for plotting
                      'labels':     ['Geo', 'Bern. Wind.', 'exclusive'], # labels to put
                      'norm':       [True, True, True]
                       }
        self.theta_dic = {
                     'sims': [sim, sim, sim],
                     'v_ns':['theta', 'theta', 'theta'],
                     'extensions':   ['_0', '_0_b_w', '_0_b'],
                     'colors':   ['black', 'red', 'orange'],
                     'labels':   ['Geo', 'Bern. Wind.', 'exclusive'],
                     'norm': [True, True, True]
                     }
        self.vel_inf_dic = {
                       'sims': [sim, sim, sim],
                       'v_ns': ['vel_inf', 'vel_inf_bern', 'vel_inf_bern'], # for bern,wind use: vel_inf_bern
                       'extensions':['_0', '_0_b_w', '_0_b'],     # directories to use (_b or _b_d ... )
                       'colors':    ['black', 'red', 'orange'],
                       'labels':    ['Geo', 'Bern. Wind.', 'exclusive'],
                       'norm': [True, True, True],
                       }
        self.entropy_dic = {
                       'sims': [sim, sim, sim],
                       'v_ns': ['entropy', 'entropy'],
                       'extensions':['_0', '_0_b_w','_0_b'],     # directories to use (_b or _b_d ... )
                       'colors':    ['black', 'red', 'orange'],
                       'labels':    ['Geo', 'Bern. Wind.', 'exclusive'],
                       'norm': [True, True, True],
                       }

    def plot_hist_ye(self, ax):

        for i in range(len(self.ye_dic['v_ns'])):
            hist_cl = SIM_EJ_HIST(self.ye_dic['sims'][i],
                                  extension=self.ye_dic['extensions'][i],
                                  norm=self.ye_dic['norm'][i])

            ye, M = hist_cl.get_hist_s_M(self.ye_dic['v_ns'][i], self.set_use_norm)
            ax.step(ye, M, color=self.ye_dic['colors'][i], where='mid', label=self.ye_dic['labels'][i])

        ax.set_xlim(xmin=0.0, xmax=0.5)
        ax.xaxis.set_major_locator(MultipleLocator(0.1))
        ax.xaxis.set_minor_locator(AutoMinorLocator(5))
        ax.set_xlabel(r"$Y_e$")
        # ax.legend(loc="upper right")

    def plot_hist_theta(self, ax):

        dtht = None
        for i in range(len(self.theta_dic['v_ns'])):
            hist_cl = SIM_EJ_HIST(self.ye_dic['sims'][i], extension=self.theta_dic['extensions'][i],
                                  norm=self.theta_dic['norm'][i])

            tht, M = hist_cl.get_hist_s_M(self.theta_dic['v_ns'][i], self.set_use_norm)
            ax.step(90. - (tht / np.pi * 180.), M,
                    color=self.theta_dic['colors'][i], where='mid', label=self.theta_dic['labels'][i])
            dtht = tht[1] - tht[0]

        ax.set_xlim(xmin=0 - dtht / np.pi * 180, xmax=90.)
        xmajorticks = np.arange(5) * 90. / 4.
        xminorticks = np.arange(17) * 90. / 16
        xmajorlabels = [r"$0^\circ$", r"$22.5^\circ$", r"$45^\circ$",
                        r"$67.5^\circ$", r"$90^\circ$"]
        ax.xaxis.set_major_locator(FixedLocator(xmajorticks))
        ax.xaxis.set_minor_locator(FixedLocator(xminorticks))
        ax.set_xticklabels(xmajorlabels)

        ax.set_xlabel(r"Angle from orbital plane")

    def plot_hist_entropy(self, ax):

        s_ = []
        for i in range(len(self.entropy_dic['v_ns'])):
            hist_cl = SIM_EJ_HIST(self.ye_dic['sims'][i], extension=self.entropy_dic['extensions'][i],
                                  norm=self.entropy_dic['norm'][i])

            s_, M_ = hist_cl.get_hist_s_M(self.entropy_dic['v_ns'][i], self.set_use_norm)
            ax.step(s_, M_, color=self.entropy_dic['colors'][i], where='mid', label=self.entropy_dic['labels'][i])

        ax.set_xlim(xmin=0, xmax=50)
        ax.xaxis.set_major_locator(FixedLocator(s_[::4]))
        ax.xaxis.set_minor_locator(FixedLocator(s_))

        ax.set_xlabel(r"$ye\ [{\rm k_B}]$")
        # ax.set_ylabel(r"$M/M_{\mathrm{ej}}$")

        # ax.set_yscale("log")
        # ax.set_ylim(ymin=1e-3, ymax=0.4)

        # ax.legend(loc="upper right")

    def plot_hist_vel_inf(self, ax):

        for i in range(len(self.vel_inf_dic['v_ns'])):
            hist_cl = SIM_EJ_HIST(self.ye_dic['sims'][i], extension=self.vel_inf_dic['extensions'][i],
                                  norm=self.vel_inf_dic['norm'][i])

            vel, M = hist_cl.get_hist_s_M(self.vel_inf_dic['v_ns'][i], self.set_use_norm)
            ax.step(vel, M, color=self.vel_inf_dic['colors'][i], where='mid', label=self.vel_inf_dic['labels'][i])

        ax.set_xlim(xmin=0., xmax=1.)
        ax.set_xticks(np.arange(0, 1.2, .2))
        # ax.xaxis.set_major_locator(FixedLocator(s_[::9]))
        # ax.xaxis.set_minor_locator(FixedLocator(s_))

        ax.set_xlabel(r"$\upsilon_{\infty}$ [c]")


    def plot_for_hists(self, tasks=[]):

        # basic_size = [4, 2]
        rows = 1
        cols = len(tasks)
        f, (ax_list) = plt.subplots(rows, cols, sharex=False, sharey=False,
                                    figsize=self.gen_set["figsize"])
                                    # figsize=[basic_size[0] * 1 * len(tasks),  # x <-->
                                    #          basic_size[1] * 1])

        # sim = self.sim

        for i, task in enumerate(tasks):
            ax = ax_list[i]
            if task == 'ye':       self.plot_hist_ye(ax)
            if task == 'vel_inf':  self.plot_hist_vel_inf(ax)
            if task == 'theta':    self.plot_hist_theta(ax)
            if task == 'entropy':  self.plot_hist_entropy(ax)
            ax.set_yscale('log')
            ax.tick_params(axis='both', which='both',
                           labelleft=True, labelright=False, tick1On=True,
                           tick2On=True, labelsize=11, direction='in')  # labeltop

            if i == 0:
                ax.legend(loc='lower center', bbox_to_anchor=(0.5, 1.02), fancybox=False, shadow=False,
                          ncol=3,
                          fontsize=8)
                ax.set_ylabel(r"$M/M_{\mathrm{ej}}$")

        plt.subplots_adjust(hspace=0.0)
        # plt.subplots_adjust(wspace=0.0)

        for a in f.axes:
            plt.setp(a.get_xticklabels()[-1], visible=False)
            plt.setp(a.get_yticklabels()[-1], visible=False)

        plt.minorticks_on()

        if not os.path.exists(self.gen_set["figdir"]):
            os.makedirs(self.gen_set["figdir"])

        plt.savefig('{}{}.png'.format(self.gen_set["figdir"], self.gen_set["figname"]), bbox_inches='tight', dpi=128)
        plt.close()

class PLOT_CORRS:

    def __init__(self, sim):

        self.set_use_norm = True
        self.invert_theta = True
        self.sim = sim

        self.theta_ye_dic = {'v_n1': ['theta', 'theta'],        #   corr_vn1
                             'v_n2': ['ye', 'ye'],              #   corr_vn1_vn2
                             'extenstions':   ['_0', '_0_b_w'],#   /_b_w/corr_vn1_vn2
                             'labels': ['Geo', 'Bern'],
                             'colors':['black', 'red'],
                             'add_hist': False,
                             'norm': LogNorm(1e-4, 1e-1),
                             'out_name': 'corr_theta_ye',
                             'figdir': Paths.ppr_sims + sim + '/res_1d/',
                             }

        self.theta_vel_inf_dic = {'v_n1': ['theta', 'theta'],
                             'v_n2': ['vel_inf', 'vel_inf_bern'],
                             'extenstions':   ['_0', '_0_b_w'],
                             'labels': ['Geo', 'Bern'],
                             'colors':['black', 'red'],
                             'add_hist': False,
                             'norm': LogNorm(1e-4, 1e-1),
                             'out_name': 'corr_theta_vel_inf',
                             'figdir': Paths.ppr_sims + sim + '/res_1d/',
                             }

    @staticmethod
    def add_plot_background2(ax, arr_x, arr_y, mass, v_n_x, v_n_y, norm):

        if v_n_x == 'theta':
            ax.set_xlim(0, 90)
        if v_n_y == 'ye':
            # ax.axhline(y=0.25, linestyle='-', color='limegreen')
            ax.set_xlim(0, 90)
            ax.set_ylim(0.05, 0.45)
            # ax.get_yaxis().set_label_coords(-0.15, 0.5)
        if v_n_y == 'vel_inf':
            ax.set_ylim(0.00, 0.95)
            # ax.get_yaxis().set_label_coords(-0.1, 1.0)

        im = ax.pcolormesh(arr_x, arr_y, mass, norm=norm, cmap='RdBu_r')
        im.set_rasterized(True)

        divider1 = make_axes_locatable(ax)
        cax = divider1.append_axes("right", size="5%", pad=0.05)
        cbar = plt.colorbar(im, cax=cax, norm=norm)
        cbar.set_label(r"Mass fraction")
        cbar.ax.minorticks_off()

    def modify_arr(self, arr, v_n):
        if v_n == 'theta':
            if self.invert_theta:
                return 90 - (180 * arr / np.pi)# [::-1] # inverting to to get angle from the plane of the binary
            else:
                return (180 * arr / np.pi)  # [::-1] # inverting to to get angle from the plane of the binary
        else: return arr

    @staticmethod
    def plot_hist_theta(ax, x_arr, mass, label, color):



        mdens = np.sum(mass, axis=0)
        # fit = np.sin(theta / 180 * np.pi) ** 2
        # fit *= np.sum(mdens) / np.sum(fit)
        # ax.plot(x_arr, fit, 'k--', label=r'$\sin^2(\x_arr)$')       # fit
        ax.step(0.5 * (x_arr[1:] + x_arr[:-1]), mdens, where='mid', color=color, label=label)

    def plot_from_dic(self, dic):

        n_rows = len(dic['v_n1'])
        if dic['add_hist']: n_rows += 1

        n_cols = 1

        fig = plt.figure(figsize=(4.5, 2.5 * 3.6))  # (<->; v)
        axs = []
        for n in range(1, n_rows + 1):
            if n == 1:
                axs.append(fig.add_subplot(n_rows, n_cols, n))
            else:
                axs.append(fig.add_subplot(n_rows, n_cols, n))  # sharex=axs[n - 2]))

        arr1 = []
        arr2 = []
        mass_arr = []

        for i in range(len(dic['v_n1'])):
            corr_cl = SIM_CORR(self.sim, dic['extenstions'][i])
            x_arr, y_arr, mass = corr_cl.load_corr(dic['v_n1'][i], dic['v_n2'][i], norm=self.set_use_norm)
            x_arr = self.modify_arr(x_arr, dic['v_n1'][i])
            mass = np.array(mass[::-1, ...]).T

            arr1.append(x_arr)
            arr2.append(y_arr)

            # im = ax.pcolormesh(x_arr, y_arr, mass, norm=dic['norm'])
            # plt.savefig('{}{}.png'.format(LISTS.fig_dir, 'x'), bbox_inches='tight', dpi=128)
            # plt.close()
            # exit(1)

            mass_arr.append(mass)
            del mass
            del x_arr
            del y_arr


        i_plot = 0
        if dic['add_hist']:

            divider1 = make_axes_locatable(axs[i_plot])
            cax = divider1.append_axes("right", size="5%", pad=0.05)
            cax.xaxis.set_ticklabels([])
            cax.get_xaxis().set_visible(False)
            cax.get_yaxis().set_visible(False)
            # fig.patch.set_visible(False)
            cax.axis('off')
            # cb.outline.set_visible(False)

            # plot sin^2 fit
            for j in range(len(dic['v_n1'])):
                self.plot_hist_theta(axs[i_plot], arr1[j], mass_arr[j], dic['labels'][j], dic['colors'][j])


            axs[i_plot].legend(loc='lower right')

            axs[i_plot].set_yscale("log")
            axs[i_plot].set_ylim(1e-3, 1e-1)
            # ax1.yaxis.tick_right()
            # ax1.yaxis.tick_left()
            axs[i_plot].set_ylabel(r"Mass fraction")
            axs[i_plot].xaxis.set_ticklabels([])

            if dic['v_n1'][0] == 'theta':
                axs[i_plot].get_yaxis().set_label_coords(-0.15, 0.5)
                axs[i_plot].set_xlim(0, 90)
                # axs[i_plot].invert_xaxis()

            i_plot += 1

        for n in range(len(dic['v_n1'])):

            if dic['v_n1'][n] == 'theta':
                axs[i_plot].get_yaxis().set_label_coords(-0.15, 0.5)
                axs[i_plot].set_xlim(0, 90)
                # axs[i_plot].invert_xaxis()

            self.add_plot_background2(axs[i_plot], arr1[n], arr2[n], mass_arr[n],
                                      dic['v_n1'][n], dic['v_n2'][n], dic['norm'])

            axs[i_plot].text(0.05, 0.85, dic['labels'][n], color='white', transform=axs[i_plot].transAxes)
            axs[i_plot].set_ylabel(LABELS.labels(dic['v_n2'][n]))
            if axs[i_plot] == axs[-1]:
                axs[i_plot].set_xlabel(r"Angle from orbital plane")
            i_plot += 1

        plt.subplots_adjust(hspace=0.1)

        if not os.path.exists(dic["figdir"]):
            os.makedirs(dic["figdir"])

        plt.savefig('{}{}.png'.format(dic["figdir"], dic['out_name']), bbox_inches='tight', dpi=128)
        plt.close()

    def plot_vn1_vn2_for_sim(self, vn1, vn2):

        if vn1 == 'theta' and vn2 == 'ye':
            self.plot_from_dic(self.theta_ye_dic)

        elif vn1 == 'theta' and vn2 == 'vel_inf':
            self.plot_from_dic(self.theta_vel_inf_dic)
        else:
            raise NameError('Not dic. for corr_{}_{} available'.format(vn1, vn2))




    # def plot_back_theta_ye(self, ax):
    #
    #     corr_cl = SIM_CORR(self.sim, self.theta_ye_dic["criteria"][], det)
    #     y_arr, x_arr, mass = corr_cl.load_corr(v_n2, v_n1, True)
    #     x_arr = modify_arr(x_arr, v_n1)
    #
    #
    #     if v_n_x == 'theta':
    #         ax.set_xlim(0, 90)
    #     if v_n_y == 'ye':
    #         ax.axhline(y=0.25, linestyle='-', color='limegreen')
    #         ax.set_xlim(0, 90)
    #         ax.set_ylim(0.05, 0.45)
    #         ax.get_yaxis().set_label_coords(-0.15, 0.5)
    #     if v_n_y == 'vel_inf':
    #         ax.set_ylim(0.00, 0.95)
    #         # ax.get_yaxis().set_label_coords(-0.1, 1.0)
    #
    #     ax.set_ylabel(LABELS.labels(v_n_y))
    #
    #
    #     im = ax.pcolormesh(arr_x, arr_y, mass, norm=norm,
    #                         cmap='RdBu_r')
    #     im.set_rasterized(True)
    #
    #     # plot histogram (sum mass for all x_arr)
    #     mdens = np.sum(mass, axis=0)
    #
    #     divider1 = make_axes_locatable(ax)
    #     cax = divider1.append_axes("right", size="5%", pad=0.05)
    #     cbar = plt.colorbar(im, cax=cax, norm=norm)
    #     cbar.set_label(r"Mass fraction")
    #     cbar.ax.minorticks_off()

# class PLOT_EJECTA:
#
#     def __init__(self, sim):
#
#         self.sim = sim
#
#         self.task_dics=[]
#
#         # task_m_tot = {
#         #     'v_n':   ['m_tot'],
#         #     'ls':    ['-'],
#         #     'color': ['black'],
#         #     'label': ['Total mass']
#         # }
#         # self.task_dics.append(task_m_tot)
#
#         task_m_unb = {
#             'v_n':      ['m_unb', 'm_unb'],
#             'ls':       ['-', '-'],
#             'color':    ['black', 'red'],
#             'label':    ['geo', 'bern']
#         }
#         self.task_dics.append(task_m_unb)
#
#         task_m_ej={
#             'v_n':      ['m_ej', 'm_ej'],
#             'extensions': ['geo', 'bern wind'],
#             'ls':       ['-', '-'],
#             'color':    ['black', 'red'],
#             'label':    ['', '']
#         }
#         self.task_dics.append(task_m_ej)
#
#         task_m_disk={
#             'v_n':      ['m_disk'],
#             'ls':       ['-'],
#             'color':    ['black'],
#             'label':    ['Disk Mass']
#         }
#         self.task_dics.append(task_m_disk)
#
#
#
#     @staticmethod
#     def plot_gw(ax, gw_cl, ls='-', color='black', label='waveform'):
#
#         tmerg = gw_cl.load_tmerg()
#
#         # gw = SIM_GW(self.sim)
#         time_gw, h_plus, _ = gw_cl.load_l2_m2()
#         ax.set_ylabel(r'$h_{+} [10^{5}]$', fontsize=12)
#         ax.plot((time_gw - tmerg) * 1000, h_plus / 1e5, ls, color=color, label=label)
#         ax.tick_params(axis='both', which='both', labelleft=True, labelright=False, tick1On=True, tick2On=True,
#                        labelsize=12, direction='in')
#         # ax.legend(loc='lower center', bbox_to_anchor=(0.5, 1.05), fancybox=False, shadow=False, ncol=2, fontsize=8)
#         ax.minorticks_on()
#
#         # print(tmerg); exit(0)
#         ax.text(0.9, 0.8, "{}".format(r'$t_{merg}$:' + '{} [ms]'.format("%.1f" % tmerg)),
#                 verticalalignment='bottom', horizontalalignment='right',
#                 transform=ax.transAxes, color='black', fontsize=12)
#     @staticmethod
#     def plot_m_tot(ax, totm_cl, t_merg=0, ls='-', color='black', label='dens.norm1'):
#         ax.set_ylabel(r'$M_{total}$ $[M_{\odot}]$', fontsize=12)
#
#         ax.plot((totm_cl.time_total_mass - t_merg) * 1000, totm_cl.total_mass, ls=ls, color=color, label=label)
#
#         ax.tick_params(axis='both', which='both', labelleft=True, labelright=False, tick1On=True, tick2On=True,
#                        labelsize=12, direction='in')
#         # ax1.legend(loc='lower center', bbox_to_anchor=(0.5, 1.05), fancybox=False, shadow=False, ncol=2, fontsize=8)
#         ax.minorticks_on()
#     @staticmethod
#     def plot_m_unb(ax, unb_cl, criteria, t_merg=0, ls='-', color='black', label='unb.norm1'):
#         ax.set_ylabel(r'$M_{unb}$ $[10^{-2}M_{\odot}]$', fontsize=12)
#
#         if criteria == 'geo':
#             ax.plot((unb_cl.time_unbound_mass - t_merg) * 1000,
#                     unb_cl.unb_mass * 1e2, ls, color=color, label=label)
#         elif criteria == 'bern':
#             ax.plot((unb_cl.time_unbound_mass_bern - t_merg) * 1000,
#                     unb_cl.unb_mass_bern * 1e2, ls, color=color, label=label)
#
#         ax.tick_params(axis='both', which='both', labelleft=True, labelright=False, tick1On=True, tick2On=True,
#                        labelsize=12, direction='in')  # labeltop
#         ax.legend(loc='lower center', bbox_to_anchor=(0.4, 0.7), fancybox=False, shadow=False, ncol=1, fontsize=8)
#         # ax.set_yscale('log')
#         ax.minorticks_on()
#     @staticmethod
#     def plot_m_ej(ax, ej_cl, t_merg=0, ls='-', color='black', label="mass_total_flux"):
#         ax.set_ylabel(r'$M_{ej}$ $[10^{-2}M_{\odot}]$', fontsize=12)
#         ax.plot((ej_cl.time_total_flux - t_merg) * 1000, ej_cl.mass_total_flux * 1e2, ls, color=color,
#                 label=label)
#         ax.tick_params(axis='both', which='both', labelleft=True, labelright=False, tick1On=True, tick2On=True,
#                        labelsize=12, direction='in')  # labeltop
#         ax.legend(loc='lower center', bbox_to_anchor=(0.4, 0.7), fancybox=False, shadow=False, ncol=1, fontsize=8)
#         ax.minorticks_on()
#     @staticmethod
#     def plot_m_disk(ax, disk_cl, t_merg=0, md_tcoll=True, ls='-', color='black', label='Disk Mass'):
#
#         ax.set_ylabel(r'$M_{disk}$ $[10^{-2}M_{\odot}]$', fontsize=12)
#         ax.plot((disk_cl.time_disk - t_merg) * 1000, disk_cl.disk_mass * 1e2, ls, color=color, label=label)
#         ax.tick_params(axis='both', which='both', labelleft=True, labelright=False, tick1On=True, tick2On=True,
#                        labelsize=12, direction='in')  # labeltop
#         ax.minorticks_on()
#         # ax4.legend(loc='lower center', bbox_to_anchor=(0.5, 1.05), fancybox=False, shadow=False, ncol=2, fontsize=8)
#         # ----------------------------------------------------------------
#         if md_tcoll:
#             tcoll, mdisk = disk_cl.get_t_coll_and_disk_mass()
#             if np.isfinite(tcoll):
#
#                 ax.plot((tcoll - t_merg) * 1000, mdisk * 1e2, 'o', color='black', label=r'$M_{d}|t_{coll}$')
#                 ax.legend(loc='lower center', bbox_to_anchor=(0.5, 0.8), fancybox=False, shadow=False, ncol=1,
#                           fontsize=8)
#
#                 ax.text(0.9, 0.8, "{}".format(r'$t_{coll}$:' + '{} [ms]'.format("%.1f" % (tcoll * 1000))),
#                         verticalalignment='bottom',
#                         horizontalalignment='right', transform=ax.transAxes, color='black', fontsize=12)
#             else:
#                 print("t_coll == inf")
#
#     def plot_task(self, ax, task_dic, tmerg=0):
#
#         if task_dic['v_n'][0] == 'm_ej':
#             for i in range(len(task_dic['v_n'])):
#                 ej_cl = SIM_EJ_HIST(self.sim, task_dic['extensions'][i], True)
#                 self.plot_m_ej(ax, ej_cl, tmerg, ls=task_dic['ls'][i], color=task_dic['color'][i],
#                                label=task_dic['label'][i])
#
#         if task_dic['v_n'][0] == 'm_unb':
#             for i in range(len(task_dic['v_n'])):
#                 unb_cl = SIM_DISK_UNBOUND(self.sim)
#                 self.plot_m_unb(ax, unb_cl, task_dic["criteria"][i], tmerg, task_dic["ls"][i],
#                                 task_dic["color"][i], task_dic["label"][i])
#
#         if task_dic['v_n'][0] == 'm_disk':
#             for i in range(len(task_dic['v_n'])):
#                 unb_cl = SIM_DISK_UNBOUND(self.sim)
#                 self.plot_m_disk(ax, unb_cl, tmerg, True, task_dic["ls"][i], task_dic["color"][i], task_dic["label"][i])
#
#         if task_dic['v_n'][0] == 'm_tot':
#             for i in range(len(task_dic['v_n'])):
#                 unb_cl = SIM_DISK_UNBOUND(self.sim)
#                 self.plot_m_tot(ax, unb_cl, tmerg, task_dic["ls"][i], task_dic["color"][i], task_dic["label"][i])
#
#         if task_dic['v_n'][0] == 'gw':
#             for i in range(len(task_dic['v_n'])):
#                 gw_cl = SIM_GW(self.sim)
#                 self.plot_gw(ax, gw_cl, task_dic["ls"][i], task_dic["color"][i], task_dic["label"][i])
#
#     def plot_from_dic_fro_1sim(self, figname='fluxes_disk'):
#
#         rows = len(self.task_dics)
#         cols = 1
#
#         f, (ax_list) = plt.subplots(rows, cols, sharex=True, sharey=False, figsize=(4.5, 2.5 * 3.6))
#         # figsize=[len(sims) * len(tasks), 2 * len(tasks)])
#
#         try:
#             gw_cl = SIM_GW(self.sim)
#             tmerg = gw_cl.load_tmerg()
#         except:
#             tmerg = 0
#             print("Warning: Failed lpading the GW (tmerger) -> time 0 is a beginning of simulation")
#
#         for ax, task_dic in zip(ax_list, self.task_dics):
#             self.plot_task(ax, task_dic, tmerg)
#             if ax == ax_list[0]:
#                 ax.legend(loc='lower left', bbox_to_anchor=(0.5, 0.85), fancybox=False, shadow=False,
#                                   ncol=2, fontsize=8)
#             if ax == ax_list[-1]:
#                 if tmerg > 0: ax.set_xlabel(r'$t-t_{merg}$ [ms]', fontsize=12)
#                 else: ax.set_xlabel(r'$t$ [ms]', fontsize=12)
#
#         plt.subplots_adjust(hspace=0.0)
#         plt.subplots_adjust(wspace=0.0)
#
#         plt.xticks(fontsize=12)
#         plt.yticks(fontsize=12)
#         for a in f.axes:
#             plt.setp(a.get_xticklabels()[-1], visible=False)
#             plt.setp(a.get_yticklabels()[-1], visible=False)
#
#         # plt.xlabel("time [ms]", fontsize=12)
#         plt.tick_params(axis='both', which='both', labelleft=True, labelright=False, tick1On=True, tick2On=True,
#                         labelsize=12, direction='in')  # labeltop
#         plt.minorticks_on()
#         plt.savefig('{}{}.png'.format(Paths.plots, figname), bbox_inches='tight', dpi=128)
#         plt.close()

class PLOT_EJECTA:

    def __init__(self, sim):

        # self.sim = sim

        self.task_dics=[]

        task_m_tot = {
            'v_n':   ['m_tot'],
            'ls':    ['-'],
            'color': ['black'],
            'label': ['Total mass']
        }
        # self.task_dics.append(task_m_tot)

        task_m_unb = {
            'sims':     ['DD2_M13641364_M0_SR',
                         'DD2_M13641364_M0_SR'],
            'v_n':      ['m_unb', 'm_unb'],
            'criteria': ['geo', 'bern'],
            'ls':       ['-', '-'],
            'color':    ['black', 'red'],
            'label':    ['geo', 'bern']
        }
        # self.task_dics.append(task_m_unb)

        task_m_ej={
            'sims':     [sim, sim, sim],
            'v_n':      ['m_ej', 'm_ej', 'm_ej'],
            'extensions':['_0', '_0_b_w', '_0_b'],
            'ls':       ['-', '-', '-'],
            'color':    ['black', 'red', 'orange'],
            'label':    ['', '', ''],
            'tcoll':    True,
            'sim_tcoll': sim,
            'yscale': None,
            'ymin': 1e-3,
            'ymax': None
        }
        self.task_dics.append(task_m_ej)

        task_m_disk={
            'sims':     ['DD2_M13641364_M0_SR',
                         'DD2_M13641364_M0_SR'],
            'v_n':      ['m_disk'],
            'ls':       ['-'],
            'color':    ['black'],
            'label':    ['Disk Mass']
        }
        # self.task_dics.append(task_m_disk) # Paths.ppr_sims + sim + '/postprocess/', 'figname': 'fluxes.png'

        self.gen_set = {'figdir': Paths.ppr_sims + sim + '/res_1d/',
                        'figname': 'fluxes',
                        'legend_loc': 'lower right',
                        'legend_box': (0.5, 0.85),
                        'legend_ncols': 2}

        # self.figname = 'fluxes_DD2'

    @staticmethod
    def plot_gw(ax, gw_cl, ls='-', color='black', label='waveform'):

        tmerg = gw_cl.load_tmerg()

        # gw = SIM_GW(self.sim)
        time_gw, h_plus, _ = gw_cl.load_l2_m2()
        ax.set_ylabel(r'$h_{+} [10^{5}]$', fontsize=12)
        ax.plot((time_gw - tmerg) * 1000, h_plus / 1e5, ls, color=color, label=label)
        ax.tick_params(axis='both', which='both', labelleft=True, labelright=False, tick1On=True, tick2On=True,
                       labelsize=12, direction='in')
        # ax.legend(loc='lower center', bbox_to_anchor=(0.5, 1.05), fancybox=False, shadow=False, ncol=2, fontsize=8)
        ax.minorticks_on()

        # print(tmerg); exit(0)
        ax.text(0.9, 0.8, "{}".format(r'$t_{merg}$:' + '{} [ms]'.format("%.1f" % tmerg)),
                verticalalignment='bottom', horizontalalignment='right',
                transform=ax.transAxes, color='black', fontsize=12)
    @staticmethod
    def plot_m_tot(ax, totm_cl, t_merg=0, ls='-', color='black', label='dens.norm1'):
        ax.set_ylabel(r'$M_{total}$ $[M_{\odot}]$', fontsize=12)

        ax.plot((totm_cl.time_total_mass - t_merg) * 1000, totm_cl.total_mass, ls=ls, color=color, label=label)

        ax.tick_params(axis='both', which='both', labelleft=True, labelright=False, tick1On=True, tick2On=True,
                       labelsize=12, direction='in')
        # ax1.legend(loc='lower center', bbox_to_anchor=(0.5, 1.05), fancybox=False, shadow=False, ncol=2, fontsize=8)
        ax.minorticks_on()
    @staticmethod
    def plot_m_unb(ax, unb_cl, criteria, t_merg=0, ls='-', color='black', label='unb.norm1'):
        ax.set_ylabel(r'$M_{unb}$ $[10^{-2}M_{\odot}]$', fontsize=12)

        if criteria == 'geo':
            ax.plot((unb_cl.time_unbound_mass - t_merg) * 1000,
                    unb_cl.unb_mass * 1e2, ls, color=color, label=label)
        elif criteria == 'bern':
            ax.plot((unb_cl.time_unbound_mass_bern - t_merg) * 1000,
                    unb_cl.unb_mass_bern * 1e2, ls, color=color, label=label)

        ax.tick_params(axis='both', which='both', labelleft=True, labelright=False, tick1On=True, tick2On=True,
                       labelsize=12, direction='in')  # labeltop
        ax.legend(loc='lower center', bbox_to_anchor=(0.4, 0.7), fancybox=False, shadow=False, ncol=1, fontsize=8)
        # ax.set_yscale('log')
        ax.minorticks_on()
    @staticmethod
    def plot_m_ej(ax, ej_cl, t_merg=0, ls='-', color='black', label="mass_total_flux"):
        ax.set_ylabel(r'$M_{ej}$ $[10^{-2}M_{\odot}]$', fontsize=12)
        ax.plot((ej_cl.time_total_flux - t_merg) * 1000, ej_cl.mass_total_flux * 1e2, ls, color=color,
                label=label)
        ax.tick_params(axis='both', which='both', labelleft=True, labelright=False, tick1On=True, tick2On=True,
                       labelsize=12, direction='in')  # labeltop
        ax.legend(loc='lower center', bbox_to_anchor=(0.4, 0.7), fancybox=False, shadow=False, ncol=1, fontsize=8)
        ax.minorticks_on()
    @staticmethod
    def plot_m_disk(ax, disk_cl, t_merg=0, ls='-', color='black', label='Disk Mass'):

        ax.set_ylabel(r'$M_{disk}$ $[10^{-2}M_{\odot}]$', fontsize=12)
        ax.plot((disk_cl.time_disk - t_merg) * 1000, disk_cl.disk_mass * 1e2, ls, color=color, label=label)
        ax.tick_params(axis='both', which='both', labelleft=True, labelright=False, tick1On=True, tick2On=True,
                       labelsize=12, direction='in')  # labeltop
        ax.minorticks_on()
        # ax4.legend(loc='lower center', bbox_to_anchor=(0.5, 1.05), fancybox=False, shadow=False, ncol=2, fontsize=8)
    @staticmethod
    def plot_m_disk_tcoll(ax, unb_cl, t_merg=0):
        tcoll, mdisk = unb_cl.tcoll, unb_cl.diskmass_tcol
        if np.isfinite(tcoll):
            ax.plot((tcoll - t_merg) * 1000, mdisk * 1e2, 'o', color='black', label=r'$M_{d}|t_{coll}$')
            ax.legend(loc='lower center', bbox_to_anchor=(0.5, 0.8), fancybox=False, shadow=False, ncol=1,
                      fontsize=8)

            ax.text(0.9, 0.8, "{}".format(r'$t_{coll}$:' + '{} [ms]'.format("%.1f" % (tcoll * 1000))),
                    verticalalignment='bottom',
                    horizontalalignment='right', transform=ax.transAxes, color='black', fontsize=12)
        else:
            Printcolor.yellow("t_coll == inf")

    def plot_task(self, ax, task_dic, tmerg=0):

        if task_dic['v_n'][0] == 'm_ej':
            if task_dic['tcoll']:
                unb_cl = SIM_UNBOUND(task_dic['sim_tcoll'])
                tcoll, mdisk = unb_cl.tcoll, unb_cl.diskmass_tcol
                if np.isfinite(tcoll):
                    ax.axvline(x=((tcoll - tmerg) * 1000), linestyle='--', linewidth=0.5, color='black')

            print(tmerg)

                # print(ej_cl.mass_total_flux)
        if task_dic['v_n'][0] == 'm_ej':
            for i in range(len(task_dic['v_n'])):
                print('plotting sim: {}'.format(task_dic['sims'][i]))
                ej_cl = SIM_EJ_HIST(task_dic['sims'][i], task_dic['extensions'][i], True)  # True for norm
                self.plot_m_ej(ax, ej_cl, 0, ls=task_dic['ls'][i], color=task_dic['color'][i],
                               label=task_dic['label'][i])

        if task_dic['v_n'][0] == 'm_unb':
            for i in range(len(task_dic['v_n'])):
                unb_cl = SIM_UNBOUND(task_dic['sims'][i])
                self.plot_m_unb(ax, unb_cl, task_dic["criteria"][i], tmerg, task_dic["ls"][i],
                                task_dic["color"][i], task_dic["label"][i])

        if task_dic['v_n'][0] == 'm_disk':
            for i in range(len(task_dic['v_n'])):
                disk_cl = SIM_DISK(task_dic['sims'][i])
                self.plot_m_disk(ax, disk_cl, tmerg, task_dic["ls"][i],
                                 task_dic["color"][i], task_dic["label"][i])

        if task_dic['v_n'][0] == 'm_disk':
            for i in range(len(task_dic['v_n'])):
                unb_cl = SIM_UNBOUND(task_dic['sims'][i])
                self.plot_m_disk_tcoll(ax, unb_cl, tmerg)

        if task_dic['v_n'][0] == 'm_tot':
            for i in range(len(task_dic['v_n'])):
                unb_cl = SIM_UNBOUND(task_dic['sims'][i])
                self.plot_m_tot(ax, unb_cl, tmerg, task_dic["ls"][i], task_dic["color"][i], task_dic["label"][i])

        if task_dic['v_n'][0] == 'gw':
            for i in range(len(task_dic['v_n'])):
                gw_cl = SIM_GW(task_dic['sims'][i])
                self.plot_gw(ax, gw_cl, task_dic["ls"][i], task_dic["color"][i], task_dic["label"][i])

        if task_dic['yscale'] == 'log':
            ax.set_yscale("log")
        if task_dic["ymin"] != None:
            ax.set_ylim(ymin=task_dic["ymin"])
        if task_dic["ymax"] != None:
            ax.set_ymax(ymin=task_dic["ymax"])


    def plot_from_task_dics(self):

        rows = len(self.task_dics)
        cols = 1

        f, (ax_list) = plt.subplots(rows, cols, sharex=True, sharey=False, figsize=(4.5, 3.6)) # 2.5 *
        # figsize=[len(sims) * len(tasks), 2 * len(tasks)])

        # gw_cl = SIM_GW(self.sim)
        # tmerg = gw_cl.load_tmerg()
        # if np.isnan(tmerg):
        #     print("Warning: Failed lpading the GW (tmerger) -> time 0 is a beginning of simulation")

        tmerg = 0

        if rows == 1:
            ax_list = [ax_list]

        for ax, task_dic in zip(ax_list, self.task_dics):
            self.plot_task(ax, task_dic, tmerg)
            if ax == ax_list[0]:
                ax.legend(loc=self.gen_set["legend_loc"],
                          # bbox_to_anchor=self.gen_set["legend_box"],
                          fancybox=False, shadow=False,
                          ncol=self.gen_set["legend_ncols"], fontsize=8)

            if ax == ax_list[-1]:
                if tmerg > 0: ax.set_xlabel(r'$t-t_{merg}$ [ms]', fontsize=12)
                else: ax.set_xlabel(r'$t$ [ms]', fontsize=12)

        plt.subplots_adjust(hspace=0.0)
        plt.subplots_adjust(wspace=0.0)

        plt.xticks(fontsize=12)
        plt.yticks(fontsize=12)
        for a in f.axes:
            plt.setp(a.get_xticklabels()[-1], visible=False)
            plt.setp(a.get_yticklabels()[-1], visible=False)

        # plt.xlabel("time [ms]", fontsize=12)
        plt.tick_params(axis='both', which='both', labelleft=True, labelright=False, tick1On=True, tick2On=True,
                        labelsize=12, direction='in')  # labeltop
        plt.minorticks_on()

        if not os.path.exists(self.gen_set["figdir"]):
            os.makedirs(self.gen_set["figdir"])

        plt.savefig('{}{}.png'.format(self.gen_set["figdir"], self.gen_set["figname"]), bbox_inches='tight', dpi=128)
        plt.close()


''' TASK SPECIFIC FUNCTIONS '''

def plot_DD2_M13641364_M0_LK_resolution_comparison_ejecta():

    fluxs = PLOT_EJECTA('-')

    fluxs.gen_set = {'figdir': Paths.plots + 'comparison_ejecta' + '/',
                     'figname': 'DD2_M13641364_M0_LK_R04_fluxes_convergence',
                     'legend_loc': 'upper left',
                     'legend_box': (0.5, 0.85),
                     'legend_ncols': 3}
    task_dics = []

    task_m_ej = {
        'sims': ["DD2_M13641364_M0_LK_SR_R04", "DD2_M13641364_M0_LK_SR_R04", "DD2_M13641364_M0_LK_SR_R04",
                 "DD2_M13641364_M0_LK_LR_R04", "DD2_M13641364_M0_LK_LR_R04", "DD2_M13641364_M0_LK_LR_R04",
                 "DD2_M13641364_M0_LK_HR_R04", "DD2_M13641364_M0_LK_HR_R04", "DD2_M13641364_M0_LK_HR_R04"],
        'v_n': ['m_ej', 'm_ej', 'm_ej',
                'm_ej', 'm_ej', 'm_ej',
                'm_ej', 'm_ej', 'm_ej'],
        'extensions': ['_0', '_0_b_w', '_0_b',
                       '_0', '_0_b_w', '_0_b',
                       '_0', '_0_b_w', '_0_b'],
        'ls': ['-', '-.', '--',
               '-', '-.', '--',
               '-', '-.', '--'],
        'color': ['black', 'black', 'black',
                  'red', 'red', 'red',
                  'blue', 'blue', 'blue'],
        'label': ['SR G', 'SR B W', 'SR B',
                  'LR G', 'LR B W', 'LR B',
                  'HR G', 'HR B W', 'HR B'],
        'tcoll': False,
        'sim_tcoll': None,
        'yscale': '',
        'ymin': 1e-3,
        'ymax': None
    }
    task_dics.append(task_m_ej)
    fluxs.task_dics = task_dics

    fluxs.plot_from_task_dics()

def plot_DD2_M13641364_M0_LK_resolution_comparison_hists():

    hists = PLOT_HISTS('-')

    hists.gen_set = {'figdir': Paths.plots + 'comparison_hists' + '/',
                     'figname': 'DD2_M13641364_M0_LK_R04_convergance_hists',
                     'figsize':[9, 3.2],
                     'legend_loc': 'lower right',
                     'legend_box': (0.5, 0.85)}

    hists.ye_dic = {
        'sims': ["DD2_M13641364_M0_LK_SR_R04",
                 "DD2_M13641364_M0_LK_LR_R04",
                 "DD2_M13641364_M0_LK_HR_R04"],
        'v_ns': ['ye', 'ye', 'ye'],  # var names to load the files
        'extensions': ['_0_b_w', '_0_b_w', '_0_b_w'],  # directories to use (_b or _b_d ... )
        'colors': ['black', 'red', 'blue'],  # colors for plotting
        'labels': ['SR', 'LR', 'HR'],  # labels to put
        'norm': [True, True, True]
    }
    hists.theta_dic = {
        'sims': ["DD2_M13641364_M0_LK_SR_R04",
                 "DD2_M13641364_M0_LK_LR_R04",
                 "DD2_M13641364_M0_LK_HR_R04"],
        'v_ns': ['theta', 'theta', 'theta'],
        'extensions': ['_0_b_w', '_0_b_w', '_0_b_w'],
        'colors': ['black', 'red', 'blue'],
        'labels': ['SR', 'LR', 'HR'],
        'norm': [True, True, True]
    }
    hists.vel_inf_dic = {
        'sims': ["DD2_M13641364_M0_LK_SR_R04",
                 "DD2_M13641364_M0_LK_LR_R04",
                 "DD2_M13641364_M0_LK_HR_R04"],
        'v_ns': ['vel_inf_bern', 'vel_inf_bern', 'vel_inf_bern'],  # for bern,wind use: vel_inf_bern
        'extensions': ['_0_b_w', '_0_b_w', '_0_b_w'],  # directories to use (_b or _b_d ... )
        'colors': ['black', 'red', 'blue'],
        'labels': ['SR', 'LR', 'HR'],
        'norm': [True, True, True],
    }


    hists.plot_for_hists(['theta', 'vel_inf', 'ye'])

def plot_DD2_M13641364_M0_SR_effect_of_LK_ejecta():

    fluxs = PLOT_EJECTA('-')

    fluxs.gen_set = {'figdir': Paths.plots + 'comparison_ejecta' + '/',
                     'figname': 'DD2_M13641364_M0_SR_R04_comparison_LK_fluxes_log',
                     'legend_loc': 'lower right',
                     'legend_box': (0.5, 0.85),
                     'legend_ncols': 2}
    task_dics = []

    task_m_ej = {

        'sims': ["DD2_M13641364_M0_SR_R04", "DD2_M13641364_M0_SR_R04", "DD2_M13641364_M0_SR_R04",
                 "DD2_M13641364_M0_LK_SR_R04", "DD2_M13641364_M0_LK_SR_R04", "DD2_M13641364_M0_LK_SR_R04"],
        'v_n': ['m_ej', 'm_ej', 'm_ej',
                'm_ej', 'm_ej', 'm_ej'],
        'extensions': ['_0', '_0_b_w', '_0_b',
                       '_0', '_0_b_w', '_0_b'],
        'ls': ['-', '-.', '--',
               '-', '-.', '--'],
        'color': ['black', 'black', 'black',
                  'red', 'red', 'red'],
        'label': ['no LK G', 'no LK B W', 'no LK B',
                  'LK G', 'LK B W', 'LK B'],
        'tcoll': False,
        'sim_tcoll': None,
        'yscale': 'log',
        'ymin': 1e-3,
        'ymax': None
    }
    task_dics.append(task_m_ej)
    fluxs.task_dics = task_dics

    fluxs.plot_from_task_dics()

def plot_DD2_M13641364_M0_SR_effect_of_LK_hists():

    hists = PLOT_HISTS('-')

    hists.gen_set = {'figdir': Paths.plots + 'comparison_hists' + '/',
                     'figname': 'DD2_M13641364_M0_SR_R04_effect_of_LK_hists_B',
                     'figsize': [9, 3.2],
                     'legend_loc': 'lower right',
                     'legend_box': (0.5, 0.85)}

    hists.ye_dic = {
        'sims': ["DD2_M13641364_M0_SR_R04",
                 "DD2_M13641364_M0_LK_LR_R04"],
        'v_ns': ['ye', 'ye'],  # var names to load the files
        'extensions': ['_0_b', '_0_b'],  # directories to use (_b or _b_d ... )
        'colors': ['black', 'red'],  # colors for plotting
        'labels': ['no LK', 'LK'],  # labels to put
        'norm': [True, True, True]
    }
    hists.theta_dic = {
        'sims': ["DD2_M13641364_M0_LK_SR_R04",
                 "DD2_M13641364_M0_LK_LR_R04"],
        'v_ns': ['theta', 'theta'],
        'extensions': ['_0_b', '_0_b'],
        'colors': ['black', 'red'],
        'labels': ['no LK', 'LK'],
        'norm': [True, True, True]
    }
    hists.vel_inf_dic = {
        'sims': ["DD2_M13641364_M0_LK_SR_R04",
                 "DD2_M13641364_M0_LK_LR_R04"],
        'v_ns': ['vel_inf_bern', 'vel_inf_bern'],  # for bern,wind use: vel_inf_bern
        'extensions': ['_0_b', '_0_b',],  # directories to use (_b or _b_d ... )
        'colors': ['black', 'red'],
        'labels': ['no LK', 'LK'],
        'norm': [True, True, True],
    }

    hists.plot_for_hists(['theta', 'vel_inf', 'ye'])


def plot_LS220_M13641364_M0_SR_effect_of_LK_ejecta():

    fluxs = PLOT_EJECTA('-')

    fluxs.gen_set = {'figdir': Paths.plots + 'comparison_ejecta' + '/',
                     'figname': 'LS220_M13641364_M0_SR_R04_comparison_LK_fluxes_log',
                     'legend_loc': 'upper left',
                     'legend_box': (0.5, 0.85),
                     'legend_ncols': 2}
    task_dics = []

    task_m_ej = {

        'sims': ["LS220_M13641364_M0_SR", "LS220_M13641364_M0_SR", "LS220_M13641364_M0_SR",
                 "LS220_M13641364_M0_LK_SR", "LS220_M13641364_M0_LK_SR", "LS220_M13641364_M0_LK_SR"],
        'v_n': ['m_ej', 'm_ej', 'm_ej',
                'm_ej', 'm_ej', 'm_ej'],
        'extensions': ['_0', '_0_b_w', '_0_b',
                       '_0', '_0_b_w', '_0_b'],
        'ls': ['-', '-.', '--',
               '-', '-.', '--'],
        'color': ['black', 'black', 'black',
                  'red', 'red', 'red'],
        'label': ['no LK G', 'no LK B W', 'no LK B',
                  'LK G', 'LK B W', 'LK B'],
        'tcoll': False,
        'sim_tcoll': None,
        'yscale': 'log',
        'ymin': 1e-3,
        'ymax': None
    }
    task_dics.append(task_m_ej)
    fluxs.task_dics = task_dics

    fluxs.plot_from_task_dics()

def plot_LS220_M13641364_M0_SR_effect_of_LK_hists():

    hists = PLOT_HISTS('-')

    hists.gen_set = {'figdir': Paths.plots + 'comparison_hists' + '/',
                     'figname': 'LS220_M13641364_M0_SR_effect_of_LK_hists_B',
                     'figsize': [9, 3.2],
                     'legend_loc': 'lower right',
                     'legend_box': (0.5, 0.85),
                     'legend_ncols': 2}

    _0_b = "_0_b"

    hists.ye_dic = {
        'sims': ["LS220_M13641364_M0_SR",
                 "LS220_M13641364_M0_LK_SR"],
        'v_ns': ['ye', 'ye'],  # var names to load the files
        'extensions': [_0_b, _0_b],  # directories to use (_b or _b_d ... )
        'colors': ['black', 'red'],  # colors for plotting
        'labels': ['no LK', 'LK'],  # labels to put
        'norm': [True, True, True]
    }
    hists.theta_dic = {
        'sims': ["LS220_M13641364_M0_SR",
                 "LS220_M13641364_M0_LK_SR"],
        'v_ns': ['theta', 'theta'],
        'extensions': [_0_b, _0_b],
        'colors': ['black', 'red'],
        'labels': ['no LK', 'LK'],
        'norm': [True, True, True]
    }
    hists.vel_inf_dic = {
        'sims': ["LS220_M13641364_M0_SR",
                 "LS220_M13641364_M0_LK_SR"],
        'v_ns': ['vel_inf_bern', 'vel_inf_bern'],  # for bern,wind use: vel_inf_bern
        'extensions': [_0_b, _0_b],  # directories to use (_b or _b_d ... )
        'colors': ['black', 'red'],
        'labels': ['no LK', 'LK'],
        'norm': [True, True, True],
    }

    hists.plot_for_hists(['theta', 'vel_inf', 'ye'])

def plot_LS220_M13641364_M0_resolution_comparison_ejecta():

    fluxs = PLOT_EJECTA('-')

    fluxs.gen_set = {'figdir': Paths.plots + 'comparison_ejecta' + '/',
                     'figname': 'LS220_M13641364_M0_fluxes_convergence',
                     'legend_loc': 'upper left',
                     'legend_box': (0.5, 0.85),
                     'legend_ncols': 3}
    task_dics = []

    task_m_ej = {
        'sims': ["LS220_M13641364_M0_SR", "LS220_M13641364_M0_SR", "LS220_M13641364_M0_SR",
                 "LS220_M13641364_M0_LR", "LS220_M13641364_M0_LR", "LS220_M13641364_M0_LR",
                 "LS220_M13641364_M0_HR", "LS220_M13641364_M0_HR", "LS220_M13641364_M0_HR"],
        'v_n': ['m_ej', 'm_ej', 'm_ej',
                'm_ej', 'm_ej', 'm_ej',
                'm_ej', 'm_ej', 'm_ej'],
        'extensions': ['_0', '_0_b_w', '_0_b',
                       '_0', '_0_b_w', '_0_b',
                       '_0', '_0_b_w', '_0_b'],
        'ls': ['-', '-.', '--',
               '-', '-.', '--',
               '-', '-.', '--'],
        'color': ['black', 'black', 'black',
                  'red', 'red', 'red',
                  'blue', 'blue', 'blue'],
        'label': ['SR G', 'SR B W', 'SR B',
                  'LR G', 'LR B W', 'LR B',
                  'HR G', 'HR B W', 'HR B'],
        'tcoll': False,
        'sim_tcoll': None,
        'yscale': '',
        'ymin': 1e-3,
        'ymax': None
    }
    task_dics.append(task_m_ej)
    fluxs.task_dics = task_dics

    fluxs.plot_from_task_dics()

def plot_LS220_M13641364_M0_resolution_comparison_hists():

    hists = PLOT_HISTS('-')

    hists.gen_set = {'figdir': Paths.plots + 'comparison_hists' + '/',
                     'figname': 'LS220_M13641364_M0_convergance_hists_B',
                     'figsize':[9, 3.2],
                     'legend_loc': 'lower right',
                     'legend_box': (0.5, 0.85),
                     'legend_ncols': 3}

    _b_w = '_0_b'

    hists.ye_dic = {
        'sims': ["LS220_M13641364_M0_SR",
                 "LS220_M13641364_M0_LR"],
        'v_ns': ['ye', 'ye'],  # var names to load the files
        'extensions': [_b_w, _b_w],  # directories to use (_b or _b_d ... )
        'colors': ['black', 'red'],  # colors for plotting
        'labels': ['SR', 'LR'],  # labels to put
        'norm': [True, True]
    }
    hists.theta_dic = {
        'sims': ["LS220_M13641364_M0_SR",
                 "LS220_M13641364_M0_LR"],
        'v_ns': ['theta', 'theta'],
        'extensions': [_b_w, _b_w],
        'colors': ['black', 'red'],
        'labels': ['SR', 'LR'],
        'norm': [True, True]
    }
    hists.vel_inf_dic = {
        'sims': ["LS220_M13641364_M0_SR",
                 "LS220_M13641364_M0_LR"],
        'v_ns': ['vel_inf_bern', 'vel_inf_bern'],  # for bern,wind use: vel_inf_bern
        'extensions': [_b_w, _b_w],  # directories to use (_b or _b_d ... )
        'colors': ['black', 'red'],
        'labels': ['SR', 'LR'],
        'norm': [True, True],
    }


    hists.plot_for_hists(['theta', 'vel_inf', 'ye'])

if __name__ == '__main__':

    # plot_DD2_M13641364_M0_LK_resolution_comparison_hists(); exit(1)
    # plot_DD2_M13641364_M0_LK_resolution_comparison_ejecta(); exit(1)
    # plot_DD2_M13641364_M0_SR_effect_of_LK_ejecta(); exit(1)
    # plot_DD2_M13641364_M0_SR_effect_of_LK_hists(); exit(1)


    # plot_LS220_M13641364_M0_SR_effect_of_LK_ejecta(); exit(1)
    # plot_LS220_M13641364_M0_SR_effect_of_LK_hists(); exit(1)
    # plot_LS220_M13641364_M0_resolution_comparison_ejecta(); exit(1)
    # plot_LS220_M13641364_M0_resolution_comparison_hists(); exit(1)

    # TODO For Average you are using only NOT normalized. SO you have to set a choce to get a normalized separately

    sim = "SFHo_M14521283_M0_LR"

    SIM_UNBOUND(sim)
    SIM_DISK(sim)
    SIM_EJ_HIST(sim)
    profs = PLOT_PROFS(sim);  profs.plot_vn1_for_sim('flux')
    hists = PLOT_HISTS(sim);  hists.plot_for_hists(['theta', 'vel_inf', 'ye'])
    corrs = PLOT_CORRS(sim);  corrs.plot_vn1_vn2_for_sim('theta', 'vel_inf')
    corrs.plot_vn1_vn2_for_sim('theta', 'ye')
    fluxs = PLOT_EJECTA(sim); fluxs.plot_from_task_dics()

