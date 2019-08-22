###################################################################################
#                                                                                 #
# This is a comprehensive set of analysis methods and tools                       #
# for the standart output of the Neutron Star Merger simulations                  #
# done with WhiskyTHC code.                                                       #
# For D1 analysis:                                                                #
#   :required: output of the Hydro Thorn of Cactus (for extraction spheres)       #
# For D2 analysis:                                                                #
#   :required: variable.xy.h5 files containing the xy ad xz slices of variables   #
# For D3 analysis:                                                                #
#   :required: iteration.h5 - profiles, containing 3D data for multiple variables #
#                             already put on one grid (for each refinment level)  #
# To use: run the analysis.py                                                     #
###################################################################################
from __future__ import division
from sys import path

from py.path import local

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
import click
import time
import copy
import h5py
import csv
import os
import gc
import re
warnings.filterwarnings("ignore",category=matplotlib.mplDeprecation)
# from visit_utils import *
from scidata.utils import locate
import scidata.carpet.hdf5 as h5
import scidata.xgraph as xg
from scidata.carpet.interp import Interpolator
import scivis.data.carpet2d
from scipy import interpolate

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
import sys
import os
from os import path
import multiprocessing as mp
from argparse import ArgumentParser
from functools import partial
import json
import h5py
# import matplotlib.pyplot as plt
# from matplotlib.colors import LogNorm
import time
import datetime

import scivis.units as units

from preanalysis import SIM_STATUS, PRINT_SIM_STATUS, LOAD_ITTIME
from d1analysis import ADD_METHODS_1D, NORMALIZE_NUCLEO
from plotting_nexus import PLOT_MANY_TASKS
from d2analysis import COMPUTE_STORE
from d3analysis import MAINMETHODS_STORE, LOAD_PROFILE_XYXZ

''' ==============================================| SETTINGS |====================================================== '''
task_list = ["init_data", "berntime", "tmerg_tcoll", "update_status", "print_status", "print_ppr_status", "d1plots",
             "analyze.sh", "d3all", "d3allplots", "d3dmplot", "d2movies", "d2movie", "d2plots"]
general_plots=["berntime"]
__d1plots__ =  ["correlations_vinf", "correlations_ye",
                "ejecta_profile_dyn", "ejecta_profile_dyn_bern", "histograms", "yields"]
__d2movievns__ = ["rho", "Y_e", "entropy", "temperature"]
__d3slicesvns__ = ["x", "y", "z", "rho", "w_lorentz", "vol", "press", "eps", "lapse", "velx", "vely", "velz",
             "gxx", "gxy", "gxz", "gyy", "gyz", "gzz", "betax", "betay", "betaz", 'temp', 'Ye'] + \
           ["density",  "enthalpy", "vphi", "vr", "dens_unb_geo", "dens_unb_bern", "dens_unb_garch",
            "ang_mom", "ang_mom_flux", "theta", "r", "phi" ]
__d3slicesplanes__ = ["xy", "xz"]
__d3dmass__ = ["disk_mass"]
__d3corrs__ = ["rho_r", "rho_Ye", "temp_Ye", "rho_theta", "velz_theta", "rho_ang_mom",
                 "rho_ang_mom_flux", "rho_dens_unb_bern", "ang_mom_flux_theta",
                 "ang_mom_flux_dens_unb_bern", "inv_ang_mom_flux_dens_unb_bern"]
__d3sliceplotvns__ = ["ang_mom_flux","ang_mom","dens_unb_garch","dens_unb_bern","dens_unb_geo","vr","vphi","enthalpy",
         "density","Ye","temp","velz","vely","velx","lapse","eps","press","vol","w_lorentz","rho"]
__d3dm__ = ["density_modes_lap15"]
scrdic = "/data01/numrel/vsevolod.nedora/scripts_server/"
reslist= ["*.done", "output-0003/parfile.par", "berntime.*", "WARNING*",
          "outflow*", "collated*", "waveforms", "ittime.h5"]
nlevels = 7
''' ===========================================| PreProcessing |==================================================== '''

class PRINT_SIM_PPR_STATUS(LOAD_ITTIME):

    def __init__(self, sim):

        LOAD_ITTIME.__init__(self, sim)

        self.sim = sim

        self.path_in_data = Paths.gw170817 + sim + '/'
        self.prof_in_data = Paths.gw170817 + sim + '/profiles/3d/'
        self.path_out_data = Paths.ppr_sims + sim + '/'
        self.file_for_gw_time = "data/dens.norm1.asc"
        self.file_for_ppr_time = "collated/dens.norm1.asc"

        self.print_data_from_parfile(self.path_out_data + "parfile.par")

        status = self.assert_initial_data(dir_="initial_data/")
        if not status:
            from preanalysis import INIT_DATA
            INIT_DATA(sim)

        if not os.path.isdir(self.path_out_data):
            Printcolor.red("\tNo postprocessed data")
            exit(0)

        print(" ")
        print("\tSimulation parameters:")

        is_data = self.print_init_data()
        if not is_data:
            from preanalysis import INIT_DATA
            INIT_DATA(sim)

        print("\tAsserting postprocessing:")

        # check if init_data.csv exists
        # if not os.path.isfile(Paths.ppr_sims + sim + '/init_data.csv'):
        #     try:
        #         from preanalysis import INIT_DATA
        #         INIT_DATA(sim)
        #     except:
        #         Printcolor.red("\tFailed to compute/save 'init_data.csv'")
        #     else:
        #         Printcolor.green("\tinit_data.csv is computed and saved")
        # else:
        #     Printcolor.green("\tinit_data.csv is present")

        self.tgw, self.tppr = self.print_gw_ppr_time()
        self.print_assert_collated_data()
        self.print_assert_gw_data()
        self.print_assert_outflow_data(criterion="_0")
        self.print_assert_outflow_data(criterion="_0_b_w")

        self.print_assert_d1_plots(dir_="res_1d/")
        self.assert_d2_movies(dir_="res_2d/")
        self.assert_d3_data(dir_= "res_3d/")

        print('=' * 100)
        print('\n')

    def assert_initial_data(self, dir_, fname="calcul.d"):

        status = False

        if os.path.isdir(self.path_out_data+dir_) and os.path.isfile(self.path_out_data+dir_+fname):
            Printcolor.green("\tInitial Data is copied", comma=False)
            status = True
        elif os.path.isdir(self.path_out_data+dir_) and not os.path.isfile(self.path_out_data+dir_+fname):
            Printcolor.red("\tInitial Data file {} is absent".format(fname), comma=False)
        else:
            Printcolor.red("\tInitial Data is not found", comma=False)

        return status

    def print_init_data(self, fname = '/init_data.csv'):
        # Printcolor.blue("\tPars: ", comma=True)
        fpath = Paths.ppr_sims + sim + fname
        is_file = False
        if os.path.isfile(fpath):
            is_file = True
            from preanalysis import LOAD_INIT_DATA
            data = LOAD_INIT_DATA(self.sim)
            counter = 1
            for key, val in data.par_dic.iteritems():
                if counter % 10 == 0:
                    comma = False
                else:
                    comma = True

                if int(counter-1) % 10 == 0:
                    print("\t"),

                Printcolor.blue(key, comma=True)
                try:
                    Printcolor.green("%.2f"%float(val), comma=comma)
                except:
                    Printcolor.green(val, comma=comma)
                counter = counter + 1

        else:
            Printcolor.red("\t[ {} is absent ]".format(fname), comma=False)

        print('\n')
        return is_file

    def print_data_from_parfile(self, fpath_parfile):

        parlist_to_print = [
            "PizzaIDBase::eos_file",
            "LoreneID::lorene_bns_file",
            "EOS_Thermal_Table3d::eos_filename",
            "WeakRates::table_filename"

        ]

        if not os.path.isfile(fpath_parfile):
            Printcolor.red("\tParfile is absent")
        else:
            flines = open(fpath_parfile, "r").readlines()
            for fname in parlist_to_print:
                found = False
                for fline in flines:
                    if fline.__contains__(fname):
                        Printcolor.blue("\t{}".format(fline), comma=True)
                        found = True
                if not found:
                    Printcolor.red("\t{} not found in parfile".format(fname))


    def print_gw_ppr_time(self, comma=True):

        # --- --- gw ---
        # last_output = list(self.get_outputdirs())[-1]

        # it_time_i = np.loadtxt(self.path_in_data+last_output+self.file_for_gw_time, usecols=(0, 1))
        # it_time_i[:, 1] *= 0.004925794970773136 * 1e-3  # ms
        # gw_iterations = np.array(it_time_i[:, 0], dtype=int)
        # gw_timesteps = np.array(it_time_i[:, 1], dtype=float)

        if not os.path.isfile(self.path_out_data + self.file_for_ppr_time):
            Printcolor.red("\tNo collated files; exiting...")
            exit(0)

        isdata, gw_iterations, gw_timesteps = self.get_ittime("overall", d1d2d3prof="d1")

        # --- --- ppr ----
        it_time_i = np.loadtxt(self.path_out_data + self.file_for_ppr_time, usecols=(0, 1))
        it_time_i[:, 1] *= 0.004925794970773136 * 1e-3  # ms
        ppr_iterations = np.array(it_time_i[:, 0], dtype=int)
        ppr_timesteps = np.array(it_time_i[:, 1], dtype=float)

        print('\t'),

        if int(gw_iterations[-1]) > int(ppr_iterations[-1]):
            Printcolor.red("t(ititme): {:.1f} > {:.1f} t(collated) [ms] ppr required"
                           .format(gw_timesteps[-1] * 1e3, ppr_timesteps[-1] * 1e3))
        elif int(gw_iterations[-1]) < int(ppr_iterations[-1]):
            Printcolor.red("t(ititme): {:.1f} < {:.1f} t(collated) [ms] ppr's longer"
                           .format(gw_timesteps[-1] * 1e3, ppr_timesteps[-1] * 1e3))
        else:
            Printcolor.green("t(ititme): {:.1f} = {:.1f} t(collated) [ms] - up to date"
                             .format(gw_timesteps[-1] * 1e3, ppr_timesteps[-1] * 1e3))

        if comma:
            print(' '),
        else:
            print(' ')

        return gw_timesteps[-1], ppr_timesteps[-1]

    def is_outflow_ppr_required(self):
        if self.tgw == self.tppr:
            return True
        else:
            return False

    @staticmethod
    def print_assert_content(dir, expected_files, marker1='.', marker2='x'):
        """
        If all files are found:  return "full", []
        else:                    return "partial", [missing files]
        or  :                    return "empty",   [missing files]
        :param expected_files:
        :param dir:
        :return:
        """
        status = "full"
        missing_files = []

        assert os.path.isdir(dir)
        print('['),
        for file_ in expected_files:
            if os.path.isfile(dir + file_):
                Printcolor.green(marker1, comma=True)
            else:
                Printcolor.red(marker2, comma=True)
                status = "partial"
                missing_files.append(file_)
        print(']'),
        if len(missing_files) == len(expected_files):
            status = "empty"

        return status, missing_files

    @staticmethod
    def print_assert_contents(dir, expected_files1, expected_files2,
                              marker1=':', marker2='.', marker3='x'):
        """
        If all files are found:  return "full", []
        else:                    return "partial", [missing files]
        or  :                    return "empty",   [missing files]
        :param expected_files:
        :param dir:
        :return:
        """
        status = "full"
        missing_files1 = []
        missing_files2 = []

        assert len(expected_files1) == len(expected_files2)

        assert os.path.isdir(dir)
        print('['),
        for file1, file2 in zip(expected_files1, expected_files2):
            if os.path.isfile(dir + file1) and os.path.isfile(dir + file2):
                Printcolor.green(marker1, comma=True)
            elif os.path.isfile(dir + file1) and not os.path.isfile(dir + file2):
                Printcolor.green(marker2, comma=True)
                missing_files2.append(file2)
                status = "partial"
            else:
                Printcolor.red(marker3, comma=True)
                status = "partial"
                missing_files1.append(file1)
                missing_files2.append(file2)
        print(']'),
        if len(missing_files1) == len(expected_files1):
            status = "empty"

        return status, missing_files1, missing_files2

    ''' --- POSTPROCESSING --- '''

    def print_assert_collated_data(self, comma=False):
        dir = self.path_out_data + '/collated/'
        flist = copy.deepcopy(Lists.tarball)
        flist.remove("outflow_surface_det_0_fluxdens.asc")
        flist.remove("outflow_surface_det_1_fluxdens.asc")
        Printcolor.blue("\tCollated files ", comma=True)
        status, missing = \
            self.print_assert_content(dir, flist, marker1='.', marker2='x')

        if status == "full":
            Printcolor.green(" complete", comma=True)
        elif status == "partial":
            Printcolor.yellow(" partial, ({}) missing".format(len(missing)), comma=True)
        else:
            Printcolor.red(" absent", comma=True)

        if comma:
            print(' '),
        else:
            print(' ')


    def print_assert_outflow_data(self, criterion='_0', comma=False):
        dir = self.path_out_data + '/outflow{}/'.format(criterion)
        flist = copy.deepcopy(Lists.outflow)
        if not criterion.__contains__("_b"):
            # if the criterion is not Bernoulli
            flist.remove("corr_vel_inf_bern_theta.h5", )
            flist.remove("hist_vel_inf_bern.dat")
            flist.remove("ejecta_profile_bern.dat")

        if criterion == "_0":
            name = "geodiesic 0"
        elif criterion == "_0_b_w":
            name = "bernoilli 0"
        else:
            name = criterion

        Printcolor.blue("\tOutflow {}".format(name), comma=True)
        status, missing = \
            self.print_assert_content(dir, flist, marker1='.', marker2='x')

        if status == "full":
            Printcolor.green(" complete", comma=True)
        elif status == "partial":
            Printcolor.yellow(" partial, ({}) missing".format(len(missing)), comma=True)
        else:
            Printcolor.red(" absent", comma=True)

        if comma:
            print(' '),
        else:
            print(' ')


    def print_assert_gw_data(self, comma=False):

        dir = self.path_out_data + '/waveforms/'
        flist = copy.deepcopy(Lists.gw)
        Printcolor.blue("\tGW data files ", comma=True)
        status, missing = \
            self.print_assert_content(dir, flist, marker1='.', marker2='x')

        if status == "full":
            Printcolor.green(" complete", comma=True)
        elif status == "partial":
            Printcolor.yellow(" partial, ({}) missing".format(len(missing)), comma=True)
        else:
            Printcolor.red(" absent", comma=True)

        if comma:
            print(' '),
        else:
            print(' ')


    def print_assert_d1_plots(self, dir_="res_1d/", res=".png", comma=False):
        dir = self.path_out_data + dir_
        if os.path.isdir(dir):
            flist = [plot_ + res for plot_ in __d1plots__]
            # if not criterion.__contains__("_b"):
            #     # if the criterion is not Bernoulli
            #     flist.remove("corr_vel_inf_bern_theta.h5")
            Printcolor.blue("\tD1 plots", comma=True)
            status, missing = \
                self.print_assert_content(dir, flist)

            if status == "full":
                Printcolor.green(" complete", comma=True)
            elif status == "partial":
                Printcolor.yellow(" partial, ({}) missing".format(len(missing)), comma=True)
            else:
                Printcolor.red(" absent", comma=True)

            if comma:
                print(' '),
            else:
                print(' ')
        else:
            Printcolor.red("[] absent", comma=False)

    def assert_d2_movies(self, dir_="res_2d/", extension=".mp4", comma=False):

        isd2, itd2, td2 = self.get_ittime("overall", d1d2d3prof="d2")

        Printcolor.blue("\tD2 movies, rl:", comma=True)
        v_ns = __d2movievns__

        missing = []
        outdated = []
        for rl in [0, 1, 2, 3, 4, 5, 6]:
            Printcolor.blue("{}".format(rl), comma=True),
            print('['),
            for v_n in v_ns:
                moviedir = v_n + "_rl" + str(rl) + "_movie/"
                moviename = v_n + "_rl" + str(rl) + extension

                # moviedir = str(v_n) + "_movie" + '/'
                # print(self.path_out_data+dir_+moviedir)
                # assert os.path.isdir(self.path_out_data + dir_ + moviedir)
                if os.path.isdir(self.path_out_data + dir_ + moviedir):
                    max_file = max(glob(self.path_out_data + dir_ + moviedir + "???????.png"))
                    int_max_file = int(str(max_file.split('/')[-1]).split('.png')[0])
                    is_movie = os.path.isfile(self.path_out_data + dir_ + moviedir + moviename)
                    if is_movie and itd2[-1] == int_max_file:
                        Printcolor.green(".", comma=True)
                    elif is_movie and itd2 != int_max_file:
                        Printcolor.yellow(".", comma=True)
                        outdated.append(v_n)
                    else:
                        Printcolor.red("x", comma=True)
                        missing.append(v_n)
                else:
                    Printcolor.red("x", comma=True)
                    missing.append(v_n)
            print("]"),

        if len(missing) == 0:
            Printcolor.green(" complete", comma=True)
        elif len(missing) != len(v_ns):
            Printcolor.yellow(" partial, ({}) missing".format(len(missing)), comma=True)
        else:
            Printcolor.red(" absent", comma=True)

        if comma:
            print(' '),
        else:
            print(' ')


    def assert_d3_data(self, dir_="res_3d/", comma=False):

        isprof, itprof, tprof = self.get_ittime("profiles", d1d2d3prof="prof")

        Printcolor.blue("\tD3 (", comma=True)
        Printcolor.green("%d" % len(tprof), comma=True)
        Printcolor.blue(") profiles: ", comma=False)
        for it, t in zip(itprof, tprof):
            data_dir = self.path_out_data + dir_ + str(it)
            if os.path.isdir(data_dir):
                Printcolor.green("\t {} ({:.1f}[ms])".format(it, t * 1e3), comma=True)

                # disk mass:
                flist = [v_n + ".txt" for v_n in __d3dmass__]
                data_dir = self.path_out_data + dir_ + str(it) + "/"
                assert os.path.isdir(data_dir)
                status, missing_files = self.print_assert_content(data_dir, flist)
                # profile.xy.h5
                flist = ["profile." + v_n + ".h5" for v_n in __d3slicesplanes__]
                data_dir = self.path_out_data + dir_ + str(it) + "/"
                assert os.path.isdir(data_dir)
                status, missing_files = self.print_assert_content(data_dir, flist)
                # corr_vn1_vn2.h5
                flist1 = ["corr_" + v_n + ".h5" for v_n in __d3corrs__]
                flist2 = [v_n + ".png" for v_n in __d3corrs__]
                data_dir = self.path_out_data + dir_ + str(it) + "/"
                assert os.path.isdir(data_dir)
                status, missing_files1, missing_files2 = self.print_assert_contents(data_dir, flist1, flist2)
                # corr_vn1_vn2.h5
                flist = [v_n + "_rl0.png" for v_n in __d3sliceplotvns__]
                data_dir = self.path_out_data + dir_ + str(it) + "/slices/"
                if not os.path.isdir(data_dir):
                    os.mkdir(data_dir)
                status, missing_files = self.print_assert_content(data_dir, flist)
                print(' ')
            else:
                Printcolor.red("\t {} ({:.1f}[ms]) [] absent".format(it, t * 1e3), comma=False)

        # # disk mass:
        # flist = [v_n + ".txt" for v_n in __d3dmass__]
        # for it, t in zip(itprof, tprof):
        #     data_dir = self.path_out_data + dir_ + str(it) + "/"
        #     status, missing_files = self.print_assert_content(data_dir, flist)
        #
        # # profile.xy.h5
        # flist = ["profile." + v_n + ".h5" for v_n in __d3slicesplanes__]
        # for it, t in zip(itprof, tprof):
        #     data_dir = self.path_out_data + dir_ + str(it) + "/"
        #     status, missing_files = self.print_assert_content(data_dir, flist)
        #
        # # corr_vn1_vn2.h5
        # flist = ["corr_" + v_n + ".h5" for v_n in __d3corrs__]
        # for it, t in zip(itprof, tprof):
        #     data_dir = self.path_out_data + dir_ + str(it) + "/"
        #     status, missing_files = self.print_assert_content(data_dir, flist)
        #
        # flist = [v_n + "rl0.png" for v_n in __d3sliceplotvns__]
        # for it, t in zip(itprof, tprof):
        #     data_dir = self.path_out_data + dir_ + str(it) + "/"
        #     status, missing_files = self.print_assert_content(data_dir, flist)

        if comma:
            print(' '),
        else:
            print(' ')

    #
    #
    #
    # def print_assert_outflowed_data(self, criterion):
    #     dir = Paths.ppr_sims + self.sim + "/outflow{}/".format(criterion)
    #     flist = copy.deepcopy(LOAD_FILES.list_outflowed_files)
    #     if not criterion.__contains__("_b"):
    #         # if the criterion is not Bernoulli
    #         flist.remove("hist_vel_inf_bern.dat")
    #         flist.remove("ejecta_profile_bern.dat")
    #     status, missing = \
    #         self.print_assert_data_status("outflow {} flist"
    #                                       .format(criterion), dir, flist, comma=False)
    #
    # def print_assert_outflowed_corr_data(self, criterion):
    #     dir = Paths.ppr_sims + self.sim + "/outflow{}/".format(criterion)
    #     flist = copy.deepcopy(LOAD_FILES.list_of_outflowed_h5_files)
    #     if not criterion.__contains__("_b"):
    #         # if the criterion is not Bernoulli
    #         flist.remove("corr_vel_inf_bern_theta.h5")
    #     status, missing = \
    #         self.print_assert_data_status("outflow corr {} flist"
    #                                       .format(criterion), dir, flist, comma=False)
    #
    # def print_assert_gw_data(self):
    #     dir = self.path_out_data + '/waveforms/'
    #     flist = copy.deepcopy(LOAD_FILES.list_gw_files)
    #     status, missing = \
    #         self.print_assert_data_status("waveforms", dir, flist, comma=False)
    #
    # def print_assert_mkn_data(self, criterion):
    #     dir = Paths.ppr_sims + self.sim + "/outflow{}/".format(criterion)
    #     flist = copy.deepcopy(LOAD_FILES.list_mkn_files)
    #     # if not criterion.__contains__("_b"):
    #     #     # if the criterion is not Bernoulli
    #     #     flist.remove("corr_vel_inf_bern_theta.h5")
    #     status, missing = \
    #         self.print_assert_data_status("mkn ej. data {} flist"
    #                                       .format(criterion), dir, flist, comma=False)
    #
    # ''' --- POSTPOSTPROCESSING --- '''
    #
    # def print_assert_d1_plots(self,  dir="res_1d/", res=".png"):
    #     dir = Paths.ppr_sims+self.sim+'/'+dir
    #     flist = [plot_+ res for plot_ in __d1plots__]
    #     # if not criterion.__contains__("_b"):
    #     #     # if the criterion is not Bernoulli
    #     #     flist.remove("corr_vel_inf_bern_theta.h5")
    #     status, missing = \
    #         self.print_assert_data_status("plots D1", dir, flist, comma=False)
    #
    # def print_assert_d2_movies(self, dir="res_2d/", res=".mp4"):
    #     dir = Paths.ppr_sims+self.sim+'/'+dir
    #     flist = [v_n + "_movie/" + v_n + res for v_n in __d2movievns__]
    #     # if not criterion.__contains__("_b"):
    #     #     # if the criterion is not Bernoulli
    #     #     flist.remove("corr_vel_inf_bern_theta.h5")
    #     status, missing = \
    #         self.print_assert_data_status("movies D2", dir, flist, comma=False)

def get_fname_for_init_data_from_parfile(sim):

    assert os.path.isfile(Paths.ppr_sims+sim+'/parfile.par')

    initial_data = ""
    lines = open(Paths.ppr_sims+sim+'/parfile.par', "r").readlines()
    for line in lines:
        if line.__contains__("LoreneID::lorene_bns_file"):
            initial_data = line
            break

    assert initial_data != ""

    run = initial_data.split("/")[-3]
    initial_data_dirname = initial_data.split("/")[-2]

    return run, initial_data_dirname

def copy_extract_initial_data(sim, finaldir = 'initial_data'):

    run, dirnmame = get_fname_for_init_data_from_parfile(sim)

    assert os.path.isdir(Paths.lorene+run+'/')

    assert os.path.isfile(Paths.lorene+run+'/'+dirnmame+'.tar.gz')

    if not os.path.isdir(Paths.ppr_sims+sim+'/'+ finaldir + '/'):
        os.mkdir(Paths.ppr_sims+sim+'/'+ finaldir + '/')

    # Andrea's fucking new approach
    if run == "R05":
        os.system("tar -xzf {} --directory {}".format(
            Paths.lorene+run+'/'+dirnmame+".tar.gz",
            Paths.ppr_sims+sim+"/"+finaldir+'/'
        ))
    else:
        os.system("rm -r {}".format(Paths.ppr_sims+sim+'/'+ finaldir + '/'))
        os.system("tar -xzf {} --directory {}".format(
            Paths.lorene+run+'/'+dirnmame+".tar.gz",
            Paths.ppr_sims+sim+"/"
        ))
        os.system("mv {} {}".format(
            Paths.ppr_sims + sim + '/' + dirnmame,
            Paths.ppr_sims + sim + '/' + finaldir
        ))

    Printcolor.blue("Initial data {} ({}) has been extracted and moved"
                    .format(dirnmame, run))


def get_it_time_cum_mass(sim, asci_file="outflow_det_0.asc"):

    simdir = Paths.gw170817 + sim + '/'
    resdir = Paths.ppr_sims + sim + '/'

    fpath = simdir + "output-*" + '/data/' + asci_file
    # assert os.path.isfile(fpath)

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

def get_bern_time(sim, fraction = 0.98, berntimefile = "berntime.txt", plot = "berntime.png"):

    simdir = Paths.gw170817 + sim + '/'
    resdir = Paths.ppr_sims + sim + '/'

    times, masses1, masses2 = get_it_time_cum_mass(sim)

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

def plot_ittime_status2(toffset=0., ax=None):

    simdir = Paths.gw170817 + sim + '/'
    resdir = Paths.ppr_sims + sim + '/'

    plotname = "ittime.png"
    color_d1_data = "lightgreen"
    color_d2_data = "green"
    color_d3_data = "blue"
    textcolor='black'

    import matplotlib.pyplot as plt
    # exit(1)
    if ax == None:
        Printcolor.blue("\tPlotting {}".format('{}'.format(plotname)))
        fig = plt.figure(figsize=(17., 3.0)) # <-> |
        ax = fig.add_subplot(111)

    o_ittime = LOAD_ITTIME(sim)

    alld1data, allitd1, alltd1 = o_ittime.get_ittime(output="overall", d1d2d3prof="d1")

    # setting up a rising sequence for annotations
    n_outputs = len(o_ittime.get_list_outputs())
    ylim1, ylim2 = ax.get_ylim()
    y_coords = np.mgrid[0:ylim2:n_outputs*1j]


    if alld1data: alltd1 -= toffset
    for i_o, output in enumerate(o_ittime.get_list_outputs()):

        d1data, itd1, td1 = o_ittime.get_ittime(output=output, d1d2d3prof="d1")
        d2data, itd2, td2 = o_ittime.get_ittime(output=output, d1d2d3prof="d2")
        d3data, itd3, td3 = o_ittime.get_ittime(output=output, d1d2d3prof="d3")
        # if d3data: print(output); exit(1)
        band_color = "gray"
        if d1data: td1 -= toffset; band_color = color_d1_data
        if d2data: td2 -= toffset; band_color = color_d2_data
        if d3data: td3 -= toffset; band_color = color_d3_data
        intoutput = int(output.split("output-")[-1])
        tcoord = td1[0] + 0.5 * (td1[-1] - td1[0])

        # print output-xxxx annotation
        ax.annotate(str(intoutput), xy=(tcoord*1e3, y_coords[i_o]), textcoords='data', color=textcolor)
        # print output-xxxx band
        ax.axvspan(td1[0]*1e3, td1[-1]*1e3, ymin=0.0, ymax=1.0, alpha=0.5, color=band_color)

        if d3data:
            for t in td3:
                y_coord = (t-td1[0])/(td3[-1]-td3[0])
                # print(y_coord)
                ax.plot(t*1e3, y_coord, marker='.', color=textcolor)
                ax.annotate("{:.1f}".format(t*1e3), xy=(t, y_coord), textcoords='data', color=textcolor)


    profdata, profit, proftime = o_ittime.get_ittime("profiles", d1d2d3prof="prof")
    if profdata: proftime -= toffset
    if profdata:
        for t in proftime:
            ax.axvline(x=t*1e3, color=textcolor, linestyle='dashed', linewidth=.5)
    # else:
    #     print("No profiles found")

    # import matplotlib.pyplot as plt



    # ax.axvline(x=times[idx_] * 0.004925794970773136 * 1e-3, color='black', lw='0.5', label='Bern time')
    # ax.plot(times * 0.004925794970773136 * 1e-3, masses2 * 1e2, '-', color='black', label='Bernoulli wind')
    # ax.plot(times * 0.004925794970773136 * 1e-3, masses1 * 1e2, '-', color='gray', label='Geodesic')
    ax.set_xlabel(r'time $[ms]$', fontsize=12)
    # ax.set_ylabel(r'$M_{ej}$ $[10^{-2}M_{\odot}]$', fontsize=12)
    # ax.minorticks_on()
    # ax.set_title('Outflow', fontsize=20)

    # ax.legend(loc='upper right', numpoints=1)
    # ax.tick_params(
    #     axis='both', which='both', labelleft=True,
    #     labelright=False, tick1On=True, tick2On=True,
    #     labelsize=12,
    #     direction='in',
    #     bottom=True, top=True, left=True, right=True
    # )
    plt.savefig(resdir + plotname, bbox_inches='tight', dpi=128)
    Printcolor.blue("\tfinished plotting {}".format('{}.png'.format(plotname)))
    plt.close()

def plot_ittime_with_flux_and_berntime(toffset=0., ax=None):

    simdir = Paths.gw170817 + sim + '/'
    resdir = Paths.ppr_sims + sim + '/'

    plotname = "flux_ittime_tbern.png"
    d1times_0, d1mflux_0, d1cumflux_0 = get_it_time_cum_mass(simdir,asci_file="outflow_det_0.asc")
    # d1times_b, d1mflux_b, d1cumflux_b = get_it_time_cum_mass(simdir, asci_file="outflow_det_0_b_w.asc")
    tbern = get_bern_time(simdir)

    d1times_0 = d1times_0 * 0.004925794970773136 * 1e-3
    # d1times_b = d1times_b * 0.004925794970773136 * 1e-3
    tbern = tbern * 0.004925794970773136 * 1e-3


    if ax == None:
        Printcolor.blue("\tPlotting {}".format('{}'.format(resdir+plotname)))
        fig = plt.figure(figsize=(17.5, 3.0)) # <-> |
        ax = fig.add_subplot(111)

    ax.plot(d1times_0*1e3, d1cumflux_0 * 1e2, lw=1.7, color='blue', marker='.', label='Dyn. Ejecta')
    # ax.plot(d1times_b, d1cumflux_b, lw=1.7, color='red', ls='.')

    ax.axvline(x=tbern*1e3, linewidth=1.2, linestyle='dashed', color='black')
    ax.set_ylabel(r'$M_{ej}$ $[10^{-2}M_{\odot}]$', fontsize=12)


    try:
        d1class = ADD_METHODS_1D(sim)
        d1times_b = d1class.get_arr(v_n="t_tot_flux", criterion="_0_b_w")
        d1cumflux_b = d1class.get_arr(v_n="mass_tot_flux", criterion="_0_b_w")
        ax.plot(d1times_b*1e3, d1cumflux_b * 1e2, lw=1.7, color='red', marker='.', label="Bern. Wind")

        tmerg = d1class.get_par(v_n="tmerger_gw")
        ax.axvline(x=tmerg*1e3, linewidth=2.0, linestyle='solid', color='pink', label="tmerg")

        # new_tick_locations = np.array(np.arange(start=d1times_0[0] * 1e3, stop=d1times_0[-1] * 1e3, step=5.))
        # new_tick_locations = np.array(np.arange(start=d1times_0[0] * 1e3, stop=d1times_0[-1] * 1e3, step=5.))
        o_ittime = LOAD_ITTIME(sim)
        t_ticks = []
        for output in o_ittime.get_list_outputs():
            isdata, itd1, td1 = o_ittime.get_ittime(output=output, d1d2d3prof="d1")
            t_ticks.append(td1[-1])

        new_tick_locations= np.array(t_ticks) * 1e3
        ax.set_xticks(new_tick_locations)
        ax.set_xticklabels(["%.1f" % x for x in new_tick_locations], rotation=90)

        ax2 = ax.twiny()

        _, itd1, td1 = o_ittime.get_ittime("overall", "d1")
        ax2.set_xticks(new_tick_locations)
        ax2.set_xticklabels(["%.1f" % (x-tmerg*1e3) for x in new_tick_locations], rotation=90)
        ax2.set_xlabel(r"$t-t_{merg}$ [ms]")
        # ax.set_xticklabels(xticklabels, rotation=45, ha="right")


        ax.set_xlim(td1[0]*1e3, td1[-1]*1e3)
        ax2.set_xlim(td1[0]*1e3, td1[-1]*1e3)
    except:
        Printcolor.yellow("Failed to use outfllow data for ittime plitting")

    plt.xticks(rotation=90)
    ax.legend(loc='upper right', bbox_to_anchor=(0.1, 0.95), fancybox=False)
    plot_ittime_status2(toffset, ax=ax)

''' =================================================| D1 |========================================================= '''

def tmerg_tcoll(sim, fraction=0.01):
    """

    :param sim:
    :param fraction: Fraction of the maximum strain amplitude to be used to set a BH formation.
    :return:
    """

    strain_fname = "strain_l2_m2.dat"
    outfile_tmerg = "tmerger2.dat"
    outfile_tcoll = "tcoll.dat"
    plot_name = "tmergtcoll.png"

    if os.path.isdir(Paths.ppr_sims + sim + '/' + "waveforms/") and \
            os.path.isfile(Paths.ppr_sims + sim + '/' + "waveforms/" + strain_fname):

        w_dir = Paths.ppr_sims + sim + '/' + "waveforms/"
        strain_file = w_dir + strain_fname

        time, reh, imh = np.genfromtxt(strain_file, unpack=True, usecols=[0, 1, 2])
        h = reh + 1j * imh
        amp = np.abs(h)
        tmerg = time[np.argmax(amp)]
        time -= tmerg

        maxA = np.max(amp)
        indx_collapse = np.min(np.where((amp < fraction * maxA) & (time >= 0)))
        tcoll = time[indx_collapse]
        tmerg = 0
        print "tmerg: {:.4f}".format(tmerg)
        print "tcoll: {:.4f}".format(tcoll)

        Printcolor.blue("\tSaving t merger:   {}".format(w_dir + outfile_tmerg))
        open(w_dir+ outfile_tmerg, "w").write("{}\n".format(float(tmerg)))
        Printcolor.blue("\tSaving t collapse: {}".format(w_dir + outfile_tmerg))
        open(w_dir+ outfile_tcoll, "w").write("{}\n".format(float(tcoll)))

        f = plt.figure()
        plt.plot(time, reh, c='k')
        plt.plot(time, amp, c='r', label="amplitude")
        plt.axvline(tmerg, ls='-.', c='b', label="tmerg")
        plt.axvline(tcoll, ls='-.', c='c', label="tcoll")

        plt.ylabel(r'strain [Re]', fontsize=12)
        plt.xlabel(r'time [s]', fontsize=12)
        plt.minorticks_on()
        plt.xticks(fontsize=12)
        plt.yticks(fontsize=12)
        plt.title('GW analysis', fontsize=20)
        plt.legend(loc='upper right', numpoints=1)
        plt.savefig(w_dir + plot_name, bbox_inches='tight', dpi=128)
        plt.close()
    else:
        Printcolor.red("\tError. No waferom dir or/and no strann file found")

def ejecta_profiles(sim, d1class, plot_bern = True, fig_dir ='res_1d/', figname = 'ejecta_profile.png'):
    simdir = Paths.gw170817 + sim + '/'
    resdir = Paths.ppr_sims + sim + '/'

    Printcolor.blue("\tInitializing ejecta profile plotting.")

    if not os.path.isdir(resdir):
        raise IOError("No simulation results frolder found: {}".format(resdir))
    if not os.path.isdir(resdir + fig_dir):
        Printcolor.yellow("\tCreating {} ".format(fig_dir))
        os.mkdir(resdir+fig_dir)

    o_plot = PLOT_MANY_TASKS()
    o_plot.gen_set["figdir"] = resdir + fig_dir
    o_plot.gen_set["type"] = "cartesian"
    o_plot.gen_set["dpi"] = 512
    o_plot.gen_set["figsize"] = (4.2, 3.5)  # <->, |]
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
        'xmin':None, 'xmax':None, 'ymin': 0, 'ymax': 1.4,
        'xlabel': r"$t-t_{merg}$ [ms]", 'ylabel': r'$M_{ej}$ $[10^{-2}M_{\odot}]$',
        'xunits': 'ms', '-t': d1class.get_par("tmerger_gw"), 'yunits': '1e-2Msun',
        'xscale':None, 'yscale': None,
        'label': "Dyn.", 'title': {},
        'fancyticks': True, 'minorticks': True,
        'fontsize': 16,
        'labelsize': 16,
        'legend': {'loc': 'upper left', 'ncol': 1, 'fontsize': 14} # 'loc': 'best',
    }
    o_plot.set_plot_dics.append(dic_ej_prof)
    # bern
    if plot_bern:
        dic_ej_prof = {
            'task': 'ejprof', 'ptype': 'cartesian',
            'position': (1, 1),
            'data': d1class, 'criterion': '_0_b_w',
            'v_n_x': 't_tot_flux', 'v_n_y': 'mass_tot_flux',
            'color': 'blue', 'ls': '-', 'lw': 1., 'ds': 'default','alpha':1.,
            'xmin':None, 'xmax':None, 'ymin': 0, 'ymax': 1.4,
            'xlabel': r"$t-t_{merg}$ [ms]", 'ylabel': r'$M_{ej}$ $[10^{-2}M_{\odot}]$',
            'xunits': 'ms', '-t': d1class.get_par("tmerger_gw"), 'yunits': '1e-2Msun',
            'label': "Wind ", 'yscale': None, 'title': {},
            'fancyticks': True, 'minorticks': True,
            'fontsize': 16,
            'labelsize': 16,
            'legend': {'loc': 'upper left', 'ncol': 1, 'fontsize': 14} # 'loc': 'best',
        }
        o_plot.set_plot_dics.append(dic_ej_prof)
    Printcolor.blue("\tPlotting ejecta profile")
    o_plot.main()
    Printcolor.blue("\tPlotting ejecta profile is done.")
    print('\n')

def ejecta_properties(sim, d1class, plot_bern = True, fig_dir ='res_1d/', figname = 'histograms.png'):
    simdir = Paths.gw170817 + sim + '/'
    resdir = Paths.ppr_sims + sim + '/'

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

def ejecta_correlations_vinf(sim, d1class, plot_bern = True, fig_dir ='res_1d/', figname = 'correlations_vinf.png'):
    simdir = Paths.gw170817 + sim + '/'
    resdir = Paths.ppr_sims + sim + '/'

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

def ejecta_correlations_ye(sim, d1class, plot_bern = True, fig_dir ='res_1d/', figname = 'correlations_ye.png'):
    simdir = Paths.gw170817 + sim + '/'
    resdir = Paths.ppr_sims + sim + '/'

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


    coll_time_vert = {
        'task': 'vertline', 'dtype': '-', 'ptype': 'cartesian',
        'value': 0.25,
        'position': (1, 1),
        'ls': '-.', 'color': 'black', 'lw': 0.8,
    }
    o_plot.set_plot_dics.append(coll_time_vert)


    Printcolor.blue("\tPlotting ejecta correlations")
    o_plot.main()
    Printcolor.blue("\tPlotting ejecta correlations is done.")
    print('\n')

def plot_nucleo_yields(sim, nuc_class, plot_bern = True, fig_dir ='res_1d/', figname = 'yields.png'):

    simdir = Paths.gw170817 + sim + '/'
    resdir = Paths.ppr_sims + sim + '/'

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

''' =================================================| D2 |========================================================= '''

def plots_d2(d1class, d2class, v_n, iterations, times, rl=3, fin_dir_end='plots', rewritefigs=False):
    simdir = Paths.gw170817 + d2class.sim + '/'
    resdir = Paths.ppr_sims + d2class.sim + '/'

    def_dic_xz = {'task': 'slice', 'dtype': '2d rl', 'ptype': 'cartesian', 'aspect': 1.,
                  'data': d2class, 'it': 00000, 'plane': 'xz', 'rl': rl,
                  'position': (1, 1),  # 'title': '[{:.1f} ms]'.format(time_),
                  'cbar': {'location': 'right .03 -0.125', 'label': r'$\rho$ [geo]',  # 'fmt': '%.1e',
                           'labelsize': 14,
                           'fontsize': 14},
                  'v_n_x': 'x', 'v_n_y': 'z', 'v_n': 'rho',
                  'xmin': None, 'xmax': None, 'ymin': None, 'ymax': None, 'vmin': 1e-10, 'vmax': 1e-4,
                  'fill_vmin': False,  # fills the x < vmin with vmin
                  'xscale': None, 'yscale': None,
                  'mask': None, 'cmap': 'inferno_r', 'norm': "log",  # 'inferno_r'
                  'fancyticks': True,
                  'title': {"text": r'$t-t_{merg}:$' + r'${:.1f}$'.format(0), 'fontsize': 14},
                  'sharex': True,  # removes angular citkscitks
                  'fontsize': 14,
                  'labelsize': 14
                  }

    def_dic_xy = {'task': 'slice', 'dtype': '2d rl', 'ptype': 'cartesian', 'aspect': 1.,
                  'data': d2class, 'it': 00000, 'plane': 'xy', 'rl': rl,
                  'position': (2, 1),  # 'title': '[{:.1f} ms]'.format(time_),
                  'cbar': {},
                  'v_n_x': 'x', 'v_n_y': 'y', 'v_n': 'rho',
                  'xmin': None, 'xmax': None, 'ymin': None, 'ymax': None, 'vmin': 1e-10, 'vmax': 1e-4,
                  'fill_vmin': False,  # fills the x < vmin with vmin
                  'xscale': None, 'yscale': None,
                  'mask': None, 'cmap': 'inferno_r', 'norm': "log",
                  'fancyticks': True,
                  'title': {},
                  'sharex': False,  # removes angular citkscitks
                  'fontsize': 14,
                  'labelsize': 14
                  }

    if v_n == "rho":
        def_dic_xz['v_n'] = 'rho'
        def_dic_xz['vmin'] = 1e-10
        def_dic_xz['vmax'] = 1e-4
        def_dic_xz['cbar']['label'] = r'$\rho$ [geo]'

        def_dic_xy['v_n'] = 'rho'
        def_dic_xy['vmin'] = 1e-10
        def_dic_xy['vmax'] = 1e-4
    elif v_n == "Y_e":
        def_dic_xz['v_n'] = 'Y_e'
        def_dic_xz['vmin'] = 0.05
        def_dic_xz['vmax'] = 0.5
        def_dic_xz['cbar']['label'] = r'$Y_e$ [geo]'
        def_dic_xz['norm'] = "linear"
        def_dic_xz['cmap'] = 'inferno'

        def_dic_xy['v_n'] = 'Y_e'
        def_dic_xy['vmin'] = 0.05
        def_dic_xy['vmax'] = 0.5
        def_dic_xy['norm'] = "linear"
        def_dic_xy['cmap'] = 'inferno'
    elif v_n == "temp" or v_n == "temperature":
        def_dic_xz['v_n'] = 'temperature'
        def_dic_xz['vmin'] = 1e-2
        def_dic_xz['vmax'] = 1e2
        def_dic_xz['cbar']['label'] = r'$Temperature$ [geo]'

        def_dic_xy['v_n'] = 'temperature'
        def_dic_xy['vmin'] = 1e-2
        def_dic_xy['vmax'] = 1e2
    elif v_n == 'entropy':
        def_dic_xz['v_n'] = 'entropy'
        def_dic_xz['vmin'] = 1e-1
        def_dic_xz['vmax'] = 1e2
        def_dic_xz['cbar']['label'] = r'$Entropy$ [geo]'

        def_dic_xy['v_n'] = 'entropy'
        def_dic_xy['vmin'] = 1e-1
        def_dic_xy['vmax'] = 1e2
    else:
        raise NameError("v_n:{} not recogmized".format(v_n))
        # for rl_ in [5,6]:
    #     print("rl:{}".format(rl))
    #     for it, t in zip(d2class.iterations, d2class.times):
    #         x = d2class.get_grid_v_n_rl(it=it,rl=rl_, plane='xy',v_n='x')[:, 0]
    #         print("it:{} rl:{} x:[{:.1f}, {:.1f}]".format(it, rl_, float(x[0]), float(x[-1])))

    # exit(1)

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

    def_dic_xy['xmin'], def_dic_xy['xmax'] = xmin, xmax
    def_dic_xy['ymin'], def_dic_xy['ymax'] = ymin, ymax
    def_dic_xz['xmin'], def_dic_xz['xmax'] = xmin, xmax
    def_dic_xz['ymin'], def_dic_xz['ymax'] = zmin, zmax

    # rho_ave_dic_xz = {'task': 'slice', 'dtype': '2d rl', 'ptype': 'cartesian',
    #     'data': d2class, 'it': 00000, 'plane': 'xz', 'rl':rl,
    #     'position': (1, 1),  # 'title': '[{:.1f} ms]'.format(time_),
    #     'cbar': {'location': 'right .03 .0', 'label': r'$\rho$ [geo]',  # 'fmt': '%.1e',
    #              'labelsize': 14,
    #              'fontsize': 14},
    #     'v_n_x': 'x', 'v_n_y': 'z', 'v_n': 'rho',
    #     'xmin': None, 'xmax': None, 'ymin': None, 'ymax': None, 'vmin': 1e-10, 'vmax': 1e-4,
    #     'fill_vmin': False,  # fills the x < vmin with vmin
    #     'xscale': None, 'yscale': None,
    #     'mask': None, 'cmap': 'inferno_r', 'norm': "log",
    #     'fancyticks': False,
    #     'title': {"text": r'$t-t_{merg}:$' + r'${:.1f}$'.format(0), 'fontsize': 14},
    #     'sharex': True,  # removes angular citkscitks
    #     'fontsize': 14,
    #     'labelsize': 14
    #     }
    # rho_ave_dic_xy = {'task': 'slice', 'dtype': '2d rl', 'ptype': 'cartesian',
    #     'data': d2class, 'it': 00000, 'plane': 'xy', 'rl':rl,
    #     'position': (2, 1),  # 'title': '[{:.1f} ms]'.format(time_),
    #     'cbar': {},
    #     'v_n_x': 'x', 'v_n_y': 'y', 'v_n': 'rho',
    #     'xmin': None, 'xmax': None, 'ymin': None, 'ymax': None, 'vmin': 1e-10, 'vmax': 1e-4,
    #     'fill_vmin': False,  # fills the x < vmin with vmin
    #     'xscale': None, 'yscale': None,
    #     'mask': None, 'cmap': 'inferno_r', 'norm': "log",
    #     'fancyticks': False,
    #     'title': {},
    #     'sharex': False,  # removes angular citkscitks
    #     'fontsize': 14,
    #     'labelsize': 14
    #     }

    # ye_dic_xz = {'task': 'slice', 'dtype': '2d rl', 'ptype': 'cartesian',
    #     'data': d2class, 'it': 00000, 'plane': 'xz', 'rl':rl,
    #     'position': (1, 1),  # 'title': '[{:.1f} ms]'.format(time_),
    #     'cbar': {'location': 'right .03 .0', 'label': r'$Y_e$ [geo]',  # 'fmt': '%.1e',
    #              'labelsize': 14,
    #              'fontsize': 14},
    #     'v_n_x': 'x', 'v_n_y': 'z', 'v_n': 'Y_e',
    #     'xmin': None, 'xmax': None, 'ymin': None, 'ymax': None, 'vmin': 0.05, 'vmax': 0.45,
    #     'fill_vmin': False,  # fills the x < vmin with vmin
    #     'xscale': None, 'yscale': None,
    #     'mask': None, 'cmap': 'inferno', 'norm': "linear",
    #     'fancyticks': False,
    #     'title': {"text": r'$t-t_{merg}:$' + r'${:.1f}$'.format(0), 'fontsize': 14},
    #     'sharex': True,  # removes angular citkscitks
    #     'fontsize': 14,
    #     'labelsize': 14
    #     }
    # ye_dic_xy = {'task': 'slice', 'dtype': '2d rl', 'ptype': 'cartesian',
    #     'data': d2class, 'it': 00000, 'plane': 'xy', 'rl':rl,
    #     'position': (2, 1),  # 'title': '[{:.1f} ms]'.format(time_),
    #     'cbar': {},
    #     'v_n_x': 'x', 'v_n_y': 'y', 'v_n': 'Y_e',
    #     'xmin': None, 'xmax': None, 'ymin': None, 'ymax': None, 'vmin': 0.05, 'vmax': 0.45,
    #     'fill_vmin': False,  # fills the x < vmin with vmin
    #     'xscale': None, 'yscale': None,
    #     'mask': None, 'cmap': 'inferno', 'norm': "linear",
    #     'fancyticks': False,
    #     'title': {},
    #     'sharex': False,  # removes angular citkscitks
    #     'fontsize': 14,
    #     'labelsize': 14
    #     }

    # temp_dic_xz = {'task': 'slice', 'dtype': '2d rl', 'ptype': 'cartesian',
    #     'data': d2class, 'it': 00000, 'plane': 'xz', 'rl':rl,
    #     'position': (1, 1),  # 'title': '[{:.1f} ms]'.format(time_),
    #     'cbar': {'location': 'right .03 .0', 'label': r'$Temperature$ [geo]',  # 'fmt': '%.1e',
    #              'labelsize': 14,
    #              'fontsize': 14},
    #     'v_n_x': 'x', 'v_n_y': 'z', 'v_n': 'temperature',
    #     'xmin': None, 'xmax': None, 'ymin': None, 'ymax': None, 'vmin': 1e-2, 'vmax': 1e2,
    #     'fill_vmin': False,  # fills the x < vmin with vmin
    #     'xscale': None, 'yscale': None,
    #     'mask': None, 'cmap': 'inferno_r', 'norm': "log",
    #     'fancyticks': False,
    #     'title': {"text": r'$t-t_{merg}:$' + r'${:.1f}$'.format(0), 'fontsize': 14},
    #     'sharex': True,  # removes angular citkscitks
    #     'fontsize': 14,
    #     'labelsize': 14
    #     }
    # temp_dic_xy = {'task': 'slice', 'dtype': '2d rl', 'ptype': 'cartesian',
    #     'data': d2class, 'it': 00000, 'plane': 'xy', 'rl':rl,
    #     'position': (2, 1),  # 'title': '[{:.1f} ms]'.format(time_),
    #     'cbar': {},
    #     'v_n_x': 'x', 'v_n_y': 'y', 'v_n': 'temperature',
    #     'xmin': None, 'xmax': None, 'ymin': None, 'ymax': None, 'vmin': 1e-2, 'vmax': 1e2,
    #     'fill_vmin': False,  # fills the x < vmin with vmin
    #     'xscale': None, 'yscale': None,
    #     'mask': None, 'cmap': 'inferno_r', 'norm': "log",
    #     'fancyticks': False,
    #     'title': {},
    #     'sharex': False,  # removes angular citkscitks
    #     'fontsize': 14,
    #     'labelsize': 14
    #     }

    # s_dic_xz = {'task': 'slice', 'dtype': '2d rl', 'ptype': 'cartesian',
    #     'data': d2class, 'it': 00000, 'plane': 'xz', 'rl':rl,
    #     'position': (1, 1),  # 'title': '[{:.1f} ms]'.format(time_),
    #     'cbar': {'location': 'right .03 .0', 'label': r'$Entropy$ [geo]',  # 'fmt': '%.1e',
    #              'labelsize': 14,
    #              'fontsize': 14},
    #     'v_n_x': 'x', 'v_n_y': 'z', 'v_n': 'entropy',
    #     'xmin': None, 'xmax': None, 'ymin': None, 'ymax': None, 'vmin': 1e-1, 'vmax': 1e2,
    #     'fill_vmin': False,  # fills the x < vmin with vmin
    #     'xscale': None, 'yscale': None,
    #     'mask': None, 'cmap': 'inferno_r', 'norm': "log",
    #     'fancyticks': False,
    #     'title': {"text": r'$t-t_{merg}:$' + r'${:.1f}$'.format(0), 'fontsize': 14},
    #     'sharex': True,  # removes angular citkscitks
    #     'fontsize': 14,
    #     'labelsize': 14
    #     }
    # s_dic_xy = {'task': 'slice', 'dtype': '2d rl', 'ptype': 'cartesian',
    #     'data': d2class, 'it': 00000, 'plane': 'xy', 'rl':rl,
    #     'position': (2, 1),  # 'title': '[{:.1f} ms]'.format(time_),
    #     'cbar': {},
    #     'v_n_x': 'x', 'v_n_y': 'y', 'v_n': 'entropy',
    #     'xmin': None, 'xmax': None, 'ymin': None, 'ymax': None, 'vmin': 1e-1, 'vmax': 1e2,
    #     'fill_vmin': False,  # fills the x < vmin with vmin
    #     'xscale': None, 'yscale': None,
    #     'mask': None, 'cmap': 'inferno_r', 'norm': "log",
    #     'fancyticks': False,
    #     'title': {},
    #     'sharex': False,  # removes angular citkscitks
    #     'fontsize': 14,
    #     'labelsize': 14
    #     }

    # if v_n == 'rho':
    #     dic_xy = rho_dic_xy
    #     dic_xz = rho_dic_xz
    # elif v_n == 'Y_e':
    #     dic_xy = ye_dic_xy
    #     dic_xz = ye_dic_xz
    # elif v_n == 'temperature':
    #     dic_xy = temp_dic_xy
    #     dic_xz = temp_dic_xz
    # elif v_n == 'entropy':
    #     dic_xy = s_dic_xy
    #     dic_xz = s_dic_xz
    # else:
    #     raise NameError("v_n:{} not recogmized".format(v_n))

    """ --- --- --- """

    # Printcolor.blue("\tInitializing profile for plotting the movie.")

    if not os.path.isdir(resdir):
        raise IOError("No simulation results frolder found: {}".format(resdir))

    if not os.path.isdir(resdir + 'res_2d/'):
        Printcolor.blue("\tCreating {} ".format(resdir + 'res_2d/'))
        os.mkdir(resdir + 'res_2d/')

    full_path = resdir + 'res_2d/' + str(v_n) + '_rl{}'.format(rl) + '_' + fin_dir_end + '/'
    if not os.path.isdir(full_path):
        Printcolor.blue("\tCreating {} ".format(full_path))
        os.mkdir(full_path)


    # d2class.set_all_it_times_from_outputs(rho_dic_xy["plane"], rho_dic_xy["v_n"])
    # d2class.load_all(dic_xy["plane"], dic_xy["v_n"])
    # d2class.load_all(dic_xz["plane"], dic_xz["v_n"])
    # d2class.set_use_new_output_if_duplicated = True

    o_plot = PLOT_MANY_TASKS()
    o_plot.gen_set["figdir"] = full_path
    o_plot.gen_set["type"] = "cartesian"
    o_plot.gen_set["dpi"] = 128
    o_plot.gen_set["figsize"] = (4.2, 8.0)  # <->, |] # to match hists with (8.5, 2.7)
    o_plot.gen_set["figname"] = "0000000.png"
    o_plot.gen_set["sharex"] = False
    o_plot.gen_set["sharey"] = False
    o_plot.gen_set["subplots_adjust_h"] = -0.35
    o_plot.gen_set["subplots_adjust_w"] = 0.2
    o_plot.set_plot_dics = []

    i = 1
    tmerg = d1class.get_par("tmerger_gw")
    for it, t in zip(iterations, times):  # zip([346112],[0.020]):# #

        figname = "{0:07d}.png".format(int(it))
        tr = (t - tmerg) * 1e3  # ms
        # print("\nProcessing v_n:{} it:{} ".format(v_n, it) + "{:.2f}ms".format(tr))
        # if os.path.isfile(full_path + figname) and not rewritefigs:
        #     Printcolor.blue("Skipping v_n:{} it:{} [{}/{}] t:"
        #                     .format(v_n, it, i, len(d2class.iterations)) + "{:.2f}ms".format(tr))
        #     o_plot.set_plot_dics = []
        #     continue
        # else:
        #     Printcolor.blue("\nProcessing v_n:{} rl:{} it:{} [{}/{}] t:"
        #                     .format(v_n, rl, it, i, len(d2class.iterations)) + "{:.2f}ms".format(tr))
        #     dic_xz["it"] = int(it)
        #     dic_xz["title"]["text"] = r'$t-t_{merg}:$' + r'${:.2f}ms$'.format(float(tr))
        #     o_plot.gen_set["figname"] = figname  # 7 digit output
        #     o_plot.set_plot_dics.append(dic_xz)
        #
        #     dic_xy["it"] = int(it)
        #     # rho_dic_xy["title"]["text"] = r'$t-t_{merg}:$' + r'${:.2f}ms$'.format(float(tr))
        #     o_plot.gen_set["figname"] = figname  # 7 digit output
        #     o_plot.set_plot_dics.append(dic_xy)
        #
        #     o_plot.main()
        #     o_plot.set_plot_dics = []
        try:
            if os.path.isfile(full_path + figname) and not rewritefigs:
                Printcolor.blue("Skipping v_n:{} it:{} [{}/{}] t:"
                                .format(v_n, it, i, len(d2class.iterations)) + "{:.2f}ms".format(tr))
                o_plot.set_plot_dics = []
                continue
            else:
                Printcolor.blue("\nProcessing v_n:{} rl:{} it:{} [{}/{}] t:"
                                .format(v_n, rl, it, i, len(d2class.iterations)) + "{:.2f}ms".format(tr))
                def_dic_xz["it"] = int(it)
                def_dic_xz["title"]["text"] = r'$t-t_{merg}:$' + r'${:.2f}ms$'.format(float(tr))
                # dic_xz["cbar"]['location'] = 'right .03 -.2',
                o_plot.gen_set["figname"] = figname  # 7 digit output
                o_plot.set_plot_dics.append(def_dic_xz)

                def_dic_xy["it"] = int(it)
                # rho_dic_xy["title"]["text"] = r'$t-t_{merg}:$' + r'${:.2f}ms$'.format(float(tr))
                o_plot.gen_set["figname"] = figname  # 7 digit output
                o_plot.set_plot_dics.append(def_dic_xy)

                o_plot.main()
                o_plot.set_plot_dics = []
        except KeyboardInterrupt:
            exit(0)
        except:
            Printcolor.red("\nFailed v_n:{} rl:{} it:{} [{}/{}] t:"
                           .format(v_n, rl, it, i, len(d2class.iterations)) + "{:.2f}ms".format(tr))
        i = i + 1
        o_plot.set_plot_dics = []
        # raise ValueError("done?")

def movie_2d(d1class, d2class, v_n='rho', rl=3, fin_dir_end ='movie/', rewritefigs=False, rewritemovie=False):

    all_iterations = d2class.iterations
    all_times = d2class.times

    Printcolor.blue("\tInitializing profile for plotting the movie.")

    if not os.path.isdir(resdir):
        raise IOError("No simulation results frolder found: {}".format(resdir))


    full_path = resdir + 'res_2d/' + str(v_n) + '_rl{}'.format(rl) + '_' + fin_dir_end + '/'
    moviename = v_n + "_rl" + str(rl) + ".mp4"
    if os.path.isfile(full_path + moviename) and not rewritemovie:
        Printcolor.blue("File:{} already exists. Skipping."
                        .format(moviename))
        return 0
    if os.path.isfile(full_path + moviename) and rewritemovie:
        Printcolor.blue("File:{} already exists. Rewriting."
                        .format(moviename))

    # --- ---

    plots_d2(d1class, d2class, v_n, all_iterations, all_times, rl=rl, fin_dir_end=fin_dir_end, rewritefigs=rewritefigs)

    # --- ---
    if os.path.isfile(full_path + moviename):
        Printcolor.blue("Rewriting movie: {}".format(full_path + moviename))
        os.remove(full_path + moviename)
    else:
        Printcolor.blue("Saving movie: {}".format(full_path + moviename))
    # ffmpeg -framerate 10 -pattern_type glob -i "*.png" -s:v 1280x720 -c:v libx264 -profile:v high -crf 20 -pix_fmt yuv420p dt.mp4
    os.system("ffmpeg -framerate 10 -pattern_type glob -i '{}*.png' -s:v 1280x720 "
              "-c:v libx264 -profile:v high -crf 20 -pix_fmt yuv420p {}"
              .format(full_path, full_path + moviename))

''' =============================================| D3 [SLICES] |==================================================== '''

def do_save_xy_xz_slices(d3class, plane, v_ns, save=True, overwrite=False):

    isprofs, itprofs, tprofs = d3class.get_ittime(output="profiles", d1d2d3prof='prof')

    i = 1
    if isprofs:
        for it, t in zip(itprofs, tprofs):
            Printcolor.blue("\tSaving slice:", comma=True)
            Printcolor.green("{}".format(plane), comma=True)
            Printcolor.blue("it:", comma=True)
            Printcolor.green("{}".format(it), comma=True)
            Printcolor.blue("({}/{})".format(i,len(itprofs)), comma=True)
            Printcolor.blue("t:", comma=True)
            Printcolor.green("{:.1f}".format(t*1e3), comma=True)
            Printcolor.blue("ms", comma=True)
            print('')

            #
            #
            # Printcolor.blue("\tit:{} ({}/{}) t:{:.1f}ms extracting slice {} [d3class: {}]"
            #                 .format(it,i,len(itprofs), t * 1e3, plane, sys.getsizeof(d3class)))

            d3class.get_slice(it, plane, v_ns, save=save, overwrite=overwrite, description=None)
            i = i+1
    else:
        Printcolor.yellow("\nNo profiles found. No histograms computed.")

def plot_slice_from_3D(d1class, d3class, rl, v_n, figdir='slices/', rewritefigs=False):

    def_dic_xz = {'task': 'slice', 'dtype': '3d rl', 'ptype': 'cartesian',
                  'data': None, 'it': 00000, 'plane': 'xz', 'rl': rl,
                  'position': (1, 1),  # 'title': '[{:.1f} ms]'.format(time_),
                  'cbar': {'location': 'right .03 .0', 'label': r'$\rho$ [geo]',  # 'fmt': '%.1e',
                           'labelsize': 14,
                           'fontsize': 14},
                  'v_n_x': 'x', 'v_n_y': 'z', 'v_n': 'rho',
                  'xmin': None, 'xmax': None, 'ymin': None, 'ymax': None, 'vmin': 1e-10, 'vmax': 1e-4,
                  'fill_vmin': False,  # fills the x < vmin with vmin
                  'xscale': None, 'yscale': None,
                  'mask': None, 'cmap': 'inferno_r', 'norm': "log",
                  'fancyticks': True,
                  'title': {"text": r'$t-t_{merg}:$' + r'${:.1f}$'.format(0), 'fontsize': 14},
                  'sharex': True,  # removes angular citkscitks
                  'fontsize': 14,
                  'labelsize': 14
                  }
    def_dic_xy = {'task': 'slice', 'dtype': '3d rl', 'ptype': 'cartesian',
                  'data': None, 'it': 00000, 'plane': 'xy', 'rl': rl,
                  'position': (2, 1),  # 'title': '[{:.1f} ms]'.format(time_),
                  'cbar': {},
                  'v_n_x': 'x', 'v_n_y': 'y', 'v_n': 'rho',
                  'xmin': None, 'xmax': None, 'ymin': None, 'ymax': None, 'vmin': 1e-10, 'vmax': 1e-4,
                  'fill_vmin': False,  # fills the x < vmin with vmin
                  'xscale': None, 'yscale': None,
                  'mask': None, 'cmap': 'inferno_r', 'norm': "log",
                  'fancyticks': True,
                  'title': {},
                  'sharex': False,  # removes angular citkscitks
                  'fontsize': 14,
                  'labelsize': 14
                  }

    # tmp = d3class.get_data(688128, 3, "xy", "ang_mom_flux")
    # print(tmp.min(), tmp.max())
    # print(tmp)
    # exit(1) # dens_unb_geo

    if v_n == 'rho':
        pass
    elif v_n == 'w_lorentz':
        def_dic_xy['v_n'] = 'w_lorentz'
        def_dic_xy['vmin'] = 1
        def_dic_xy['vmax'] = 1.3
        def_dic_xy['norm'] = None
        def_dic_xz['v_n'] = 'w_lorentz'
        def_dic_xz['vmin'] = 1
        def_dic_xz['vmax'] = 1.3
        def_dic_xz['norm'] = None
    elif v_n == 'vol':
        def_dic_xy['v_n'] = 'vol'
        def_dic_xy['vmin'] = 1
        def_dic_xy['vmax'] = 10
        # def_dic_xy['norm'] = None
        def_dic_xz['v_n'] = 'vol'
        def_dic_xz['vmin'] = 1
        def_dic_xz['vmax'] = 10
        # def_dic_xz['norm'] = None
    elif v_n == 'press':
        def_dic_xy['v_n'] = 'press'
        def_dic_xy['vmin'] = 1e-12
        def_dic_xy['vmax'] = 1e-6

        def_dic_xz['v_n'] = 'press'
        def_dic_xz['vmin'] = 1e-12
        def_dic_xz['vmax'] = 1e-6
    elif v_n == 'eps':
        def_dic_xy['v_n'] = 'eps'
        def_dic_xy['vmin'] = 5e-3
        def_dic_xy['vmax'] = 5e-1
        def_dic_xz['v_n'] = 'eps'
        def_dic_xz['vmin'] = 5e-3
        def_dic_xz['vmax'] = 5e-1
    elif v_n == 'lapse':
        def_dic_xy['v_n'] = 'lapse'
        def_dic_xy['vmin'] = 0.15
        def_dic_xy['vmax'] = 1
        def_dic_xy['norm'] = None
        def_dic_xz['v_n'] = 'lapse'
        def_dic_xz['vmin'] = 0.15
        def_dic_xz['vmax'] = 1
        def_dic_xz['norm'] = None
    elif v_n == 'velx':
        def_dic_xy['v_n'] = 'velx'
        def_dic_xy['vmin'] = 0.01
        def_dic_xy['vmax'] = 1.
        # def_dic_xy['norm'] = None
        def_dic_xz['v_n'] = 'velx'
        def_dic_xz['vmin'] = 0.01
        def_dic_xz['vmax'] = 1.
        # def_dic_xz['norm'] = None
    elif v_n == 'vely':
        def_dic_xy['v_n'] = 'vely'
        def_dic_xy['vmin'] = 0.01
        def_dic_xy['vmax'] = 1.
        # def_dic_xy['norm'] = None
        def_dic_xz['v_n'] = 'vely'
        def_dic_xz['vmin'] = 0.01
        def_dic_xz['vmax'] = 1.
        # def_dic_xz['norm'] = None
    elif v_n == 'velz':
        def_dic_xy['v_n'] = 'velz'
        def_dic_xy['vmin'] = 0.01
        def_dic_xy['vmax'] = 1.
        # def_dic_xy['norm'] = None
        def_dic_xz['v_n'] = 'velz'
        def_dic_xz['vmin'] = 0.01
        def_dic_xz['vmax'] = 1.
        # def_dic_xz['norm'] = None
    elif v_n == 'temp':
        def_dic_xy['v_n'] = 'temp'
        def_dic_xy['vmin'] =  1e-2
        def_dic_xy['vmax'] = 1e2

        def_dic_xz['v_n'] = 'temp'
        def_dic_xz['vmin'] =  1e-2
        def_dic_xz['vmax'] = 1e2
    elif v_n == 'Ye':
        def_dic_xy['v_n'] = 'Ye'
        def_dic_xy['vmin'] = 0.05
        def_dic_xy['vmax'] = 0.5
        def_dic_xy['norm'] = None
        def_dic_xy['cmap'] = 'inferno'

        def_dic_xz['v_n'] = 'Ye'
        def_dic_xz['vmin'] = 0.05
        def_dic_xz['vmax'] = 0.5
        def_dic_xz['norm'] = None
        def_dic_xz['cmap'] = 'inferno'
    elif v_n == 'density':
        def_dic_xy['v_n'] = 'density'
        def_dic_xy['vmin'] = 1e-9
        def_dic_xy['vmax'] = 1e-5
        # def_dic_xy['norm'] = None

        def_dic_xz['v_n'] = 'density'
        def_dic_xz['vmin'] = 1e-9
        def_dic_xz['vmax'] = 1e-5
        # def_dic_xz['norm'] = None
    elif v_n == 'enthalpy':
        def_dic_xy['v_n'] = 'enthalpy'
        def_dic_xy['vmin'] = 1.
        def_dic_xy['vmax'] = 1.5
        def_dic_xy['norm'] = None

        def_dic_xz['v_n'] = 'enthalpy'
        def_dic_xz['vmin'] = 1.
        def_dic_xz['vmax'] = 1.5
        def_dic_xz['norm'] = None
    elif v_n == 'vphi':
        def_dic_xy['v_n'] = 'vphi'
        def_dic_xy['vmin'] = 0.01
        def_dic_xy['vmax'] = 10.
        # def_dic_xy['norm'] = None
        def_dic_xz['v_n'] = 'vphi'
        def_dic_xz['vmin'] = 0.01
        def_dic_xz['vmax'] = 10.
        # def_dic_xz['norm'] = None
    elif v_n == 'vr':
        def_dic_xy['v_n'] = 'vr'
        def_dic_xy['vmin'] = 0.01
        def_dic_xy['vmax'] = 0.5
        # def_dic_xy['norm'] = None
        def_dic_xz['v_n'] = 'vr'
        def_dic_xz['vmin'] = 0.01
        def_dic_xz['vmax'] = 0.5
        # def_dic_xz['norm'] = None
    elif v_n == 'dens_unb_geo':
        def_dic_xy['v_n'] = 'dens_unb_geo'
        def_dic_xy['vmin'] = 1e-10
        def_dic_xy['vmax'] = 1e-5
        # def_dic_xy['norm'] = None
        def_dic_xz['v_n'] = 'dens_unb_geo'
        def_dic_xz['vmin'] = 1e-10
        def_dic_xz['vmax'] = 1e-5
        # def_dic_xz['norm'] = None
    elif v_n == 'dens_unb_bern':
        def_dic_xy['v_n'] = 'dens_unb_bern'
        def_dic_xy['vmin'] = 1e-10
        def_dic_xy['vmax'] = 1e-5
        # def_dic_xy['norm'] = None
        def_dic_xz['v_n'] = 'dens_unb_bern'
        def_dic_xz['vmin'] = 1e-10
        def_dic_xz['vmax'] = 1e-5
        # def_dic_xz['norm'] = None
    elif v_n == 'dens_unb_garch':
        def_dic_xy['v_n'] = 'dens_unb_garch'
        def_dic_xy['vmin'] = 1e-10
        def_dic_xy['vmax'] = 1e-6
        # def_dic_xy['norm'] = None
        def_dic_xz['v_n'] = 'dens_unb_garch'
        def_dic_xz['vmin'] = 1e-10
        def_dic_xz['vmax'] = 1e-6
        # def_dic_xz['norm'] = None
    elif v_n == 'ang_mom':
        def_dic_xy['v_n'] = 'ang_mom'
        def_dic_xy['vmin'] = 1e-8
        def_dic_xy['vmax'] = 1e-3
        # def_dic_xy['norm'] = None
        def_dic_xz['v_n'] = 'ang_mom'
        def_dic_xz['vmin'] = 1e-8
        def_dic_xz['vmax'] = 1e-3
        # def_dic_xz['norm'] = None
    elif v_n == 'ang_mom_flux':
        def_dic_xy['v_n'] = 'ang_mom_flux'
        def_dic_xy['vmin'] = 1e-9
        def_dic_xy['vmax'] = 1e-5
        # def_dic_xy['norm'] = None
        def_dic_xz['v_n'] = 'ang_mom_flux'
        def_dic_xz['vmin'] = 1e-9
        def_dic_xz['vmax'] = 1e-5
        # def_dic_xz['norm'] = None
    else:
        raise NameError("v_n:{} not recogmized".format(v_n))

    """ --- --- --- """
    datafpath = Paths.ppr_sims + d3class.sim + '/res_3d/'
    # outfpath = Paths.ppr_sims + d3_slices.sim + '/res_3d/'
    # if not os.path.isdir(datafpath):
    #     os.mkdir(datafpath)
    figname = "{}_rl{}.png".format(v_n, rl)

    o_plot = PLOT_MANY_TASKS()
    o_plot.gen_set["figdir"] = datafpath
    o_plot.gen_set["type"] = "cartesian"
    o_plot.gen_set["figsize"] = (4.2, 8.0)  # <->, |] # to match hists with (8.5, 2.7)
    o_plot.gen_set["figname"] = figname
    o_plot.gen_set["sharex"] = False
    o_plot.gen_set["sharey"] = False
    o_plot.gen_set["subplots_adjust_h"] = 0.0
    o_plot.gen_set["subplots_adjust_w"] = 0.2
    o_plot.set_plot_dics = []

    i = 1
    tmerg = d1class.get_par("tmerger_gw")
    for it, t in zip(d3class.list_iterations, d3class.times):  # zip([346112],[0.020]):# #

        # figname = "{0:07d}.png".format(int(it))
        tr = (t - tmerg) * 1e3  # ms
        # print("\nProcessing v_n:{} it:{} ".format(v_n, it) + "{:.2f}ms".format(tr))

        if not os.path.isfile(datafpath + str(it) + '/' + "profile.xy.h5") \
                or not os.path.isfile(datafpath + str(it) + '/' + "profile.xz.h5"):
            Printcolor.yellow("Required data ia missing: {}".format(datafpath + str(it) + '/' + "profile.xy(or yz).h5"))
            continue

        if os.path.isfile(datafpath + str(it) + '/' + figdir + figname) and not rewritefigs:
            Printcolor.blue("Skipping v_n:", True)
            Printcolor.green("{}".format(v_n), True)
            Printcolor.blue("it:", True)
            Printcolor.green("{}".format(it), True)
            Printcolor.blue("rl:", True),
            Printcolor.green("{}".format(rl), True),
            Printcolor.blue("[{}/{}] t:".format(i, len(d3class.list_iterations)), True)
            Printcolor.green("{:.2f}".format(tr), True)
            Printcolor.blue("ms", False)
            continue
        else:
            Printcolor.blue("Processing v_n:", True)
            Printcolor.green("{}".format(v_n), True)
            Printcolor.blue("it:", True)
            Printcolor.green("{}".format(it), True)
            Printcolor.blue("rl:", True),
            Printcolor.green("{}".format(rl), True),
            Printcolor.blue("[{}/{}] t:".format(i, len(d3class.list_iterations)), True)
            Printcolor.green("{:.2f}".format(tr), True)
            Printcolor.blue("ms", False)

            if not os.path.isdir(datafpath + str(it) + '/' + figdir):
                os.mkdir(datafpath + str(it) + '/' + figdir)

            if v_n in ["velx", "vely", "velz", "vphi", "vr", "ang_mom_flux"]:
                # make separate plotting >0 and <0 with log scales
                o_plot.gen_set["figdir"] = datafpath + str(it) + '/' + figdir

                def_dic_xz['cmap'] = 'Reds'
                def_dic_xz["mask"] = "negative"
                def_dic_xz['cbar'] = {'location': 'right .03 .0', 'label': v_n.replace('_', '\_') + r"$<0$",  # 'fmt': '%.1e',
                                      'labelsize': 14,
                                      'fontsize': 14}
                def_dic_xz["it"] = int(it)
                def_dic_xz["title"]["text"] = r'$t-t_{merg}:$' + r'${:.2f}ms$'.format(float(tr))

                n_def_dic_xz = def_dic_xz.copy() #copy.deepcopy(def_dic_xz)
                def_dic_xz['data'] = d3class
                o_plot.set_plot_dics.append(def_dic_xz)

                n_def_dic_xz['data'] = d3class
                n_def_dic_xz['cmap'] = 'Blues'
                n_def_dic_xz["mask"] = "positive"
                n_def_dic_xz['cbar'] = {}
                n_def_dic_xz["it"] = int(it)
                n_def_dic_xz["title"]["text"] = r'$t-t_{merg}:$' + r'${:.2f}ms$'.format(float(tr))

                o_plot.set_plot_dics.append(n_def_dic_xz)

                # --- ---
                def_dic_xy["it"] = int(it)
                def_dic_xy['cmap'] = 'Blues'
                def_dic_xy['mask'] = "positive"
                def_dic_xy['cbar'] = {'location': 'right .03 .0', 'label': v_n.replace('_', '\_') + r"$>0$",  # 'fmt': '%.1e',
                                      'labelsize': 14,
                                      'fontsize': 14}
                # n_def_dic_xy = copy.deepcopy(def_dic_xy)
                n_def_dic_xy = def_dic_xy.copy()
                def_dic_xy['data'] = d3class
                o_plot.set_plot_dics.append(def_dic_xy)

                n_def_dic_xy['data'] = d3class
                n_def_dic_xy['cbar'] = {}
                n_def_dic_xy['cmap'] = 'Reds'
                n_def_dic_xy['mask'] = "negative"
                o_plot.set_plot_dics.append(n_def_dic_xy)

                for dic in o_plot.set_plot_dics:
                    if not 'cbar' in dic.keys():
                        raise IOError("dic:{} no cbar".format(dic))

                # ---- ----
                o_plot.main()
                o_plot.set_plot_dics = []

                n_def_dic_xz = {}
                n_def_dic_xy = {}
            else:
                def_dic_xz['data'] = d3class
                def_dic_xz['cbar']['label'] = v_n.replace('_', '\_')
                def_dic_xz["it"] = int(it)
                def_dic_xz["title"]["text"] = r'$t-t_{merg}:$' + r'${:.2f}ms$'.format(float(tr))
                o_plot.gen_set["figdir"] = datafpath + str(it) + '/' + figdir
                o_plot.set_plot_dics.append(def_dic_xz)

                def_dic_xy['data'] = d3class
                def_dic_xy["it"] = int(it)
                # rho_dic_xy["title"]["text"] = r'$t-t_{merg}:$' + r'${:.2f}ms$'.format(float(tr))
                # o_plot.gen_set["figname"] =   # 7 digit output
                o_plot.set_plot_dics.append(def_dic_xy)

                o_plot.main()
                o_plot.set_plot_dics = []

        i = i + 1
    # exit(1)

''' ============================================| D3 [CORRELATIONS] |=============================================== '''

def get_disk_mass(d3corrclass, tmin=None, tmax=None, fname="disk_mass.txt", save=True, overwrite=False):
    isprofs, itprofs, tprofs = d3corrclass.get_ittime(output="profiles", d1d2d3prof='prof')

    if tmin != None:
        if tmin < tprofs.min():
            raise ValueError("tmin:{} < time from profiles minimum : {}"
                             .format(tmin, tprofs.min()))
        if tmin > tprofs.max():
            raise ValueError("tmin:{} > time from profiles maximum : {}"
                             .format(tmin, tprofs.min()))
        tprofs = tprofs[tprofs >= tmin]
        itprofs = itprofs[tprofs >= tmin]
        Printcolor.blue("{} profiles selected for time interval: [{} {}]"
                          .format(len(itprofs), tprofs.min(), tprofs.max()))

    if tmax != None:
        if tmax > tprofs.max():
            raise ValueError("tmax:{} > time from profiles maximum : {}"
                             .format(tmax, tprofs.min()))
        if tmax < tprofs.min():
            raise ValueError("tmax:{} < time from profiles minimum : {}"
                             .format(tmax, tprofs.min()))

        tprofs = tprofs[tprofs <= tmax]
        itprofs = itprofs[tprofs <= tmax]
        Printcolor.blue("{} profiles selected for time interval: [{} {}]"
                          .format(len(itprofs), tprofs.min(), tprofs.max()))

    if isprofs:
        for it, t in zip(itprofs, tprofs):
            iscomputed, data = d3corrclass.get_total_mass(it, fname=fname, save=save, overwrite=overwrite)
            if iscomputed:
                Printcolor.blue("it:{} t:{:.1f}ms Mass disk: {}".format(it, t*1e3, data))
            else:
                pass
    else:
        Printcolor.yellow("\nNo profiles found. No histograms computed.")

def do_histogram_processing_of_iterations(d3corrclass, corr_task, tmin=None, tmax=None, save=True, overwrite=False):

    isprofs, itprofs, tprofs = d3corrclass.get_ittime(output="profiles", d1d2d3prof='prof')

    if corr_task == "rho_r":
        corr_task_dic = d3corrclass.corr_task_dic_rho_r
    elif corr_task == "rho_Ye":
        corr_task_dic = d3corrclass.corr_task_dic_rho_ye
    elif corr_task == "temp_Ye":
        corr_task_dic = d3corrclass.corr_task_dic_temp_ye
    elif corr_task == "rho_theta":
        corr_task_dic = d3corrclass.corr_task_dic_rho_theta
    elif corr_task == "velz_theta":
        corr_task_dic = d3corrclass.corr_task_dic_velz_theta
    elif corr_task == "rho_ang_mom":
        corr_task_dic = d3corrclass.corr_task_dic_rho_ang_mom
    elif corr_task == "rho_ang_mom_flux":
        corr_task_dic = d3corrclass.corr_task_dic_rho_ang_mom_flux
    elif corr_task == "rho_dens_unb_bern":
        corr_task_dic = d3corrclass.corr_task_dic_rho_dens_unb_bern
    elif corr_task == "ang_mom_flux_theta":
        corr_task_dic = d3corrclass.corr_task_dic_ang_mom_flux_theta
    elif corr_task == "ang_mom_flux_dens_unb_bern":
        corr_task_dic = d3corrclass.corr_task_dic_ang_mom_flux_dens_unb_bern
    elif corr_task == "inv_ang_mom_flux_dens_unb_bern":
        corr_task_dic = d3corrclass.corr_task_dic_inv_ang_mom_flux_dens_unb_bern
    else:
        raise NameError("unknown task for correlation computation: {}"
                        .format(corr_task))

    if tmin != None:
        if tmin < tprofs.min():
            raise ValueError("tmin:{} < time from profiles minimum : {}"
                             .format(tmin, tprofs.min()))
        if tmin > tprofs.max():
            raise ValueError("tmin:{} > time from profiles maximum : {}"
                             .format(tmin, tprofs.min()))
        tprofs = tprofs[tprofs >= tmin]
        itprofs = itprofs[tprofs >= tmin]
        Printcolor.yellow("\t{} profiles selected for time interval: [{} {}]"
                          .format(len(itprofs), tprofs.min(), tprofs.max()))

    if tmax != None:
        if tmax > tprofs.max():
            raise ValueError("tmax:{} > time from profiles maximum : {}"
                             .format(tmax, tprofs.min()))
        if tmax < tprofs.min():
            raise ValueError("tmax:{} < time from profiles minimum : {}"
                             .format(tmax, tprofs.min()))

        tprofs = tprofs[tprofs <= tmax]
        itprofs = itprofs[tprofs <= tmax]
        Printcolor.yellow("\t{} profiles selected for time interval: [{} {}]"
                          .format(len(itprofs), tprofs.min(), tprofs.max()))

    if isprofs:
        for it, t in zip(itprofs, tprofs):
            try:
                iscomputed, data = d3corrclass.get_correlation(it, corr_task_dic, save=save, overwrite=overwrite)
                if iscomputed:
                    Printcolor.blue("\tit:{} t:{:.1f}ms {}-{} correlation is done. M:{}"
                                    .format(it, t * 1e3, corr_task_dic[0]["v_n"], corr_task_dic[1]["v_n"], np.sum(data)))
                else:
                    pass
            except ValueError:
                Printcolor.red("\tit:{} t:{:.1f}ms {}-{} correlation has Failed."
                               .format(it, t * 1e3, corr_task_dic[0]["v_n"], corr_task_dic[1]["v_n"]))
    else:
        Printcolor.yellow("\nNo profiles found. No histograms computed.")

def plot_hostograms_from_3D(d1class, d3histclass, vn1vn2 = "ang_mom_flux_dens_unb_bern", fig_dir ='res_3d/', rewrite=False):

    iterations = d3histclass.list_iterations
    times = d3histclass.times

    default_dic = {  # relies on the "get_res_corr(self, it, v_n): " method of data object
        'task': '2d colormesh', 'dtype': 'corr', 'ptype': 'cartesian',
        'data': d3histclass, 'it': 0,
        'position': (1, 1),
        'v_n_x': 'ang_mom_flux', 'v_n_y': 'dens_unb_bern', 'v_n': 'mass',
        'xmin': 1e-11, 'xmax': 1e-7, 'ymin': 1e-11, 'ymax': 1e-7, 'vmin': 1e-7, 'vmax': 1e-3,
        'xscale': 'log', 'yscale': 'log',
        'mask_below': None, 'mask_above': None, 'cmap': 'inferno_r', 'norm': 'log', 'todo': None,
        'cbar': {'location': 'right .03 .0', 'label': r'mass',
                 'labelsize': 14,
                 'fontsize': 14},
        'title': {"text": r'$t-t_{merg}:$' + r'${:.1f}$'.format(0), 'fontsize': 14},
        'fontsize': 14,
        'labelsize': 14,
        'minorticks': True,
        'fancyticks': True,
        'sharey': False,
        'sharex': False,
    }

    if vn1vn2 == "rho_r":
        default_dic['v_n_x'] = 'rho'
        default_dic['v_n_y'] = 'r'
        default_dic['xmin'] = 1e-9
        default_dic['xmax'] = 2e-5
        default_dic['ymin'] = 0
        default_dic['ymax'] = 500
        default_dic['yscale'] = None
    elif vn1vn2 == "rho_Ye":
        default_dic['v_n_x'] = 'rho'
        default_dic['v_n_y'] = 'Ye'
        default_dic['xmin'] = 1e-9
        default_dic['xmax'] = 2e-5
        default_dic['ymin'] = 0.01
        default_dic['ymax'] = 0.5
        default_dic['yscale'] = None
    elif vn1vn2 == "temp_Ye":
        default_dic['v_n_x'] = 'temp'
        default_dic['v_n_y'] = 'Ye'
        default_dic['xmin'] = 1e-2
        default_dic['xmax'] = 1e2
        default_dic['ymin'] = 0.01
        default_dic['ymax'] = 0.5
        default_dic['yscale'] = None
    elif vn1vn2 == "rho_theta":
        default_dic['v_n_x'] = 'rho'
        default_dic['v_n_y'] = 'theta'
        default_dic['xmin'] = 1e-9
        default_dic['xmax'] = 2e-5
        default_dic['ymin'] = 0
        default_dic['ymax'] = 1.7
        default_dic['yscale'] = None
    elif vn1vn2 == "velz_theta":
        default_dic['v_n_x'] = 'velz'
        default_dic['v_n_y'] = 'theta'
        default_dic['xmin'] = .5
        default_dic['xmax'] = -.5
        default_dic['ymin'] = 1.
        default_dic['ymax'] = 1.7
        default_dic['yscale'] = None
        default_dic['xscale'] = None
    elif vn1vn2 == "rho_ang_mom":
        default_dic['v_n_x'] = 'rho'
        default_dic['v_n_y'] = 'ang_mom'
        default_dic['xmin'] = 1e-9
        default_dic['xmax'] = 2e-5
        default_dic['ymin'] = 1e-9
        default_dic['ymax'] = 1e-3
    elif vn1vn2 == "rho_ang_mom_flux":
        default_dic['v_n_x'] = 'rho'
        default_dic['v_n_y'] = 'ang_mom_flux'
        default_dic['xmin'] = 1e-9
        default_dic['xmax'] = 2e-5
        default_dic['ymin'] = 1e-9
        default_dic['ymax'] = 8e-5
    elif vn1vn2 == "rho_dens_unb_bern":
        default_dic['v_n_x'] = 'rho'
        default_dic['v_n_y'] = 'dens_unb_bern'
        default_dic['xmin'] = 1e-9
        default_dic['xmax'] = 2e-5
        default_dic['ymin'] = 1e-9
        default_dic['ymax'] = 2e-6
    elif vn1vn2 == "ang_mom_flux_theta":
        default_dic['v_n_x'] = 'ang_mom_flux'
        default_dic['v_n_y'] = 'theta'
        default_dic['xmin'] = 1e-9
        default_dic['xmax'] = 8e-5
        default_dic['ymin'] = 0
        default_dic['ymax'] = 1.7
        default_dic['yscale'] = None
    elif vn1vn2 == "ang_mom_flux_dens_unb_bern":
        pass
    elif vn1vn2 == "inv_ang_mom_flux_dens_unb_bern":
        default_dic['v_n_x'] = 'inv_ang_mom_flux'
    else:
        raise NameError("vn1vn2:{} is not recognized"
                        .format(vn1vn2))

    outfpath = resdir + fig_dir

    o_plot = PLOT_MANY_TASKS()
    o_plot.gen_set["figdir"] = outfpath
    o_plot.gen_set["type"] = "cartesian"
    o_plot.gen_set["figsize"] = (4.2, 3.8)  # <->, |] # to match hists with (8.5, 2.7)
    o_plot.gen_set["figname"] = "{}.png".format(vn1vn2)
    o_plot.gen_set["sharex"] = False
    o_plot.gen_set["sharey"] = False
    o_plot.gen_set["subplots_adjust_h"] = 0.0
    o_plot.gen_set["subplots_adjust_w"] = 0.2
    o_plot.set_plot_dics = []

    tmerg = d1class.get_par("tmerger_gw")

    i = 0
    for it, t in zip(iterations, times):

        tr = (t - tmerg) * 1e3  # ms

        o_plot.gen_set["figdir"] = outfpath + str(it) + '/'
        if not os.path.isdir(o_plot.gen_set["figdir"]):
            raise IOError("dir with data for the plot ({}) does not seem to exist"
                          .format(o_plot.gen_set["figdir"]))
        fpath = o_plot.gen_set["figdir"] + o_plot.gen_set["figname"]

        if os.path.isfile(fpath) and not rewrite:
            Printcolor.blue("Skipping vn1vn2:", True)
            Printcolor.green("{}".format(vn1vn2), True)
            Printcolor.blue("it:", True)
            Printcolor.green("{}".format(it), True)
            Printcolor.blue("[{}/{}] t:".format(i, len(iterations)), True)
            Printcolor.green("{:.2f}".format(tr), True)
            Printcolor.blue("ms", False)
            continue
            #
            #
            #
            # Printcolor.blue("Skipping v_n:{} it:{} [{}/{}] t:"
            #                 .format(vn1vn2, it, i, len(iterations)) + "{:.2f}ms".format(tr))
            # continue
        else:
            data_fname = outfpath + str(it) + '/' + "corr_" + vn1vn2 + ".h5"
            if os.path.isfile(data_fname):
                Printcolor.blue("Skipping vn1vn2:", True)
                Printcolor.green("{}".format(vn1vn2), True)
                Printcolor.blue("it:", True)
                Printcolor.green("{}".format(it), True)
                Printcolor.blue("[{}/{}] t:".format(i, len(iterations)), True)
                Printcolor.green("{:.2f}".format(tr*1e3), True)
                Printcolor.blue("ms", False)

                default_dic["it"] = it
                default_dic["title"]["text"] = r'$t-t_{merg}:$' + r'${:.2f}ms$'.format(float(tr))
                o_plot.set_plot_dics.append(default_dic)

                o_plot.main()
                o_plot.set_plot_dics = []
            else:
                Printcolor.yellow("Warning. Data:{} is not available for it:{}"
                                  .format("corr_" + vn1vn2 + ".h5", it))
            i = i + 1

    # Printcolor.blue("\nHistogram processing is done.")

''' ============================================| D3 [DENSITY MODES] |============================================== '''

def do_compute_density_modes(d3class, rl, mmax, tmin=None, tmax=None, save=True, overwrite=False):

    fname = "density_modes_lap15.h5"
    isprofs, itprofs, tprofs = d3class.get_ittime(output="profiles", d1d2d3prof='prof')

    if tmin != None:
        if tmin < tprofs.min():
            raise ValueError("tmin:{} < time from profiles minimum : {}"
                             .format(tmin, tprofs.min()))
        if tmin > tprofs.max():
            raise ValueError("tmin:{} > time from profiles maximum : {}"
                             .format(tmin, tprofs.min()))
        tprofs = tprofs[tprofs >= tmin]
        itprofs = itprofs[tprofs >= tmin]
        Printcolor.yellow("{} profiles selected for time interval: [{} {}]"
                          .format(len(itprofs), tprofs.min(), tprofs.max()))

    if tmax != None:
        if tmax > tprofs.max():
            raise ValueError("tmax:{} > time from profiles maximum : {}"
                             .format(tmax, tprofs.min()))
        if tmax < tprofs.min():
            raise ValueError("tmax:{} < time from profiles minimum : {}"
                             .format(tmax, tprofs.min()))

        tprofs = tprofs[tprofs <= tmax]
        itprofs = itprofs[tprofs <= tmax]
        Printcolor.yellow("{} profiles selected for time interval: [{} {}]"
                          .format(len(itprofs), tprofs.min(), tprofs.max()))

    if not os.path.isdir(Paths.ppr_sims+d3class.sim+"/res_3d/"):
        os.mkdir(Paths.ppr_sims+d3class.sim+"/res_3d/")



    if isprofs:
        if os.path.isfile(Paths.ppr_sims+d3class.sim+"/res_3d/"+fname) and not overwrite:
            Printcolor.blue("Skipping:", True)
            Printcolor.green("{}".format("density modes"), True)
            Printcolor.blue("rl:", True)
            Printcolor.green("{}".format(rl), False)
        else:
            Printcolor.blue("Computing:", True)
            Printcolor.green("{}".format("density modes"), True)
            Printcolor.blue("rl:", True)
            Printcolor.green("{}".format(rl), False)

            d3class.get_dens_modes_for_rl(rl=rl, mmax = mmax, tmin=tmin, tmax=tmax, fname=fname,
                                  save=save, overwrite=overwrite)
    else:
        Printcolor.yellow("\nNo profiles found. No density modes computed computed.")

def plot_density_modes(d1class, dmclass, fname="density_modes", overwrite=False):

    path = Paths.ppr_sims + dm_class.sim + "/res_3d/"
    # fpath = path + fname
    dmclass.gen_set['fname'] = path + fname + ".h5"#"density_modes_lap15.h5"

    o_plot = PLOT_MANY_TASKS()
    o_plot.gen_set["figdir"] = path
    o_plot.gen_set["type"] = "cartesian"
    o_plot.gen_set["figsize"] = (9.0, 2.7)  # <->, |]
    o_plot.gen_set["figname"] = fname + ".png"
    o_plot.gen_set["sharex"] = True
    o_plot.gen_set["sharey"] = False
    o_plot.set_plot_dics = []

    densmode_m1 = {
        'task': 'line', 'dtype': 'dm', 'ptype': 'cartesian',
        'data': dmclass,
        'position': (1, 1),
        'v_n_x': 'times', 'v_n_y': 'int_phi_r abs',
        'xunits': 'ms', '-t': d1class.get_par("tmerger_gw"),
        'mode': 1, 'norm_to_m': 0,
        'ls': '-', 'color': 'black', 'lw': 1., 'ds': 'default', 'alpha': 1.,
        'label': r'$m=1$', 'ylabel': r'$C_m/C_0$ Magnitude', 'xlabel': r'time [ms]',
        'xmin': None, 'xmax': None, 'ymin': None, 'ymax': None,
        'xscale': None, 'yscale': 'log', 'legend': {},
        'fancyticks': True, 'minorticks': True,
        'fontsize': 14,
        'labelsize': 14,
    }
    densmode_m2 = {
        'task': 'line', 'dtype': 'dm', 'ptype': 'cartesian',
        'data': dmclass,
        'position': (1, 1),
        'v_n_x': 'times', 'v_n_y': 'int_phi_r abs',
        'xunits': 'ms', '-t': d1class.get_par("tmerger_gw"),
        'mode': 2, 'norm_to_m': 0,
        'ls': ':', 'color': 'black', 'lw': 0.8, 'ds': 'default', 'alpha': 1.,
        'label': r'$m=2$', 'ylabel': r'$C_m/C_0$ Magnitude', 'xlabel': r'time [ms]',
        'xmin': None, 'xmax': None, 'ymin': None, 'ymax': None,
        'xscale': None, 'yscale': 'log',
        'fancyticks': True, 'minorticks': True,
        'legend': {'loc': 'best', 'ncol': 3, 'fontsize': 14},
        'fontsize': 14,
        'labelsize': 14,
    }
    # o_plot.set_plot_dics.append(densmode_m0)

    if os.path.isfile(o_plot.gen_set["figname"]) and not overwrite:
        Printcolor.blue("Skipping:", True)
        Printcolor.green("{}".format("density modes"), True)
        # Printcolor.blue("rl:", True)
        # Printcolor.green("{}".format(rl), False)
    else:
        Printcolor.blue("Plotting:", True)
        Printcolor.green("{}".format("density modes"), True)
        # Printcolor.blue("rl:", True)
        # Printcolor.green("{}".format(rl), False)

        o_plot.set_plot_dics.append(densmode_m1)
        o_plot.set_plot_dics.append(densmode_m2)

        o_plot.main()


if __name__ == '__main__':
    parser = ArgumentParser(description="postprocessing pipeline")
    parser.add_argument("-s", dest="sim", required=True, help="task to perform")
    parser.add_argument("-t", dest="task", required=False, default='no', help="task to perform")
    parser.add_argument("--v_n", dest="v_n", required=False, default='no', help="variable (or group) name")
    parser.add_argument("--rl", dest="rl", required=False, default=-1, help="reflevel")
    parser.add_argument("--it", dest="it", required=False, default=-1, help="iteration")
    parser.add_argument('--times', nargs='+', help='Timesteps to use', required=False)
    parser.add_argument("-o", dest="overwrite", required=False, default=[], help="overwrite if exists")
    parser.add_argument("-p", dest="plot", required=False, default='no', help="task to perform")
    parser.add_argument("-i", dest="simdir", required=False, default=Paths.gw170817, help="path to sim dir")
    # parser.add_argument("-i", dest="filenames",
    #                     nargs='+', required=True,
    #                     help="inputfiles", metavar="FILE")
    # parser.add_argument("-p", dest="procs",
    #                     nargs='?', default=1,
    #                     help="processes in the pool")
    # parser.add_argument("-o", dest="outputdir",
    #                     nargs='?', default=os.getcwd(),
    #                     help="output directory")
    # parser.add_argument("-s", dest="setup", type=json.loads,
    #                    nargs='?', default=default_setup(),
    #                    help="JSON file for setup/parameters")

    # Info/checks
    args = parser.parse_args()
    task = args.task
    plot = args.plot
    sim = args.sim
    v_n = args.v_n
    overwrite = args.overwrite
    Paths.gw170817 = args.simdir # '/data1/numrel/WhiskyTHC/Backup/2018/SLy4_M130130_SR_physics/'
    simdir = Paths.gw170817 + sim + '/'
    resdir = Paths.ppr_sims + sim + '/'
    rl = int(args.rl)
    it = int(args.it)
    times = np.array(args.times, dtype=np.float32)
    # assert that given input is correct
    assert os.path.isdir(args.simdir)
    assert os.path.isdir(args.simdir+sim+'/')

    if task == "no" and plot == "no":
        print("No task/plot given. Exiting...")
        exit(0)
    elif task != "no" and plot == "no":
        pass
    elif task == "no" and plot != "no":
        pass
    elif task != "no" and plot != "no":
        raise IOError("Chose either -t your_task for task or -p your_plot for plotting")


    if task != "no" and not task in task_list:
        raise NameError("Task: {} is not available. See available list:\n{}"
                        .format(task, task_list))

    if task in ["d2movie"] and v_n == "no":
        raise IOError("task: {} requires option '-v' variable name"
                      .format(task))


    """ ------------------------------------------------------------------------------------------------------------ """

    if task == "init_data":
        copy_extract_initial_data(sim)


    if task == "update_status":
        SIM_STATUS(sim)

    if task == "print_status":
        PRINT_SIM_STATUS(sim)

    if task == "tmerg_tcoll":
        tmerg_tcoll(sim)

    if task == "berntime":
        get_bern_time(sim, fraction = 0.98, berntimefile = "berntime.txt", plot = "berntime.png")
        os.system("cp {} {}".format(Paths.gw170817+sim+'/'+"berntime.*", Paths.ppr_sims+sim+'/'))

    if task == "analyze.sh":
        os.system(scrdic + 'analyze.sh' + ' ' + simdir)
        if not os.path.isdir(resdir):
            os.mkdir(resdir)
        for item in reslist:
            fname = simdir + item
            os.system("cp -r " + fname + ' ' + resdir)

    if task == "print_ppr_status":
        PRINT_SIM_PPR_STATUS(sim)

    if task == 'plot_status':
        plot_ittime_with_flux_and_berntime()

    if task == "d1plots":
        # perform the plotting of the outflowed.cpp
        d1class = ADD_METHODS_1D(sim)
        ejecta_profiles(sim, d1class, plot_bern=False, fig_dir='res_1d/', figname='ejecta_profile_dyn.png')
        ejecta_profiles(sim, d1class, plot_bern=True, fig_dir='res_1d/', figname='ejecta_profile_dyn_bern.png')
        ejecta_properties(sim, d1class, plot_bern=True, fig_dir='res_1d/', figname='histograms.png')
        ejecta_correlations_vinf(sim, d1class, plot_bern=True, fig_dir='res_1d/', figname='correlations_vinf.png')
        ejecta_correlations_ye(sim, d1class, plot_bern=True, fig_dir='res_1d/', figname='correlations_ye.png')
        nucclass = NORMALIZE_NUCLEO(sim)
        plot_nucleo_yields(sim, nucclass, plot_bern=True, fig_dir='res_1d/', figname='yields.png')

    if task == "d3corr":
        d3class = MAINMETHODS_STORE(sim)
        assert v_n != "no"
        assert v_n in __d3corrs__
        do_histogram_processing_of_iterations(d3class, corr_task=v_n,
                                              tmin=None, tmax=None,
                                              save=True, overwrite=True)

    if task == "d3corrs":
        d3class = MAINMETHODS_STORE(sim)
        do_histogram_processing_of_iterations(d3class, corr_task=task,
                                              tmin=None, tmax=None, save=True, overwrite=overwrite)

    if task == "d3all":
        d3class = MAINMETHODS_STORE(sim)

        for plane in __d3slicesplanes__:
            do_save_xy_xz_slices(d3class, plane=plane, v_ns=__d3slicesvns__, save=True, overwrite=overwrite)

        ''' --- --- --- PROCESSING 3D [DENSITY MODES] --- --- --- '''

        do_compute_density_modes(d3class, rl=6, mmax=8, tmin=None, tmax=None, save=True, overwrite=overwrite)

        ''' --- --- --- PROCESSING 3D [CORRELATIONS] --- --- --- '''
        d3class.mask_setup = {'rm_rl': True,  # REMOVE previouse ref. level from the next
                              'rho': [6.e4 / 6.176e+17, 1.e13 / 6.176e+17],  # REMOVE atmo and NS
                              'lapse': [0.15, 1.]}  # remove apparent horizon

        get_disk_mass(d3class, tmin=None, tmax=None, fname=__d3dmass__[0]+".txt", save=True, overwrite=overwrite)

        for task in __d3corrs__:
            do_histogram_processing_of_iterations(d3class, corr_task=task,
                                                  tmin=None, tmax=None, save=True, overwrite=overwrite)

    if task == "d3dmplot":
        from d3analysis import LOAD_DENSITY_MODES
        d1class = ADD_METHODS_1D(sim)
        dm_class = LOAD_DENSITY_MODES(sim)
        plot_density_modes(d1class, dm_class, fname=__d3dm__[0], overwrite=True)

    if task == "d3slices":
        pass

    if task == "d3allplots":
        d1class = ADD_METHODS_1D(sim)

        ''' --- --- --- PLOTTING 3D [SLICES] --- --- --- '''
        rewrite = False
        d3_slices = LOAD_PROFILE_XYXZ(sim)
        for rl in [0, 1, 2, 3, 4, 5, 6]:
            for v_n in __d3sliceplotvns__:
                plot_slice_from_3D(d1class, d3_slices, rl=rl, v_n=v_n, rewritefigs=rewrite)

        ''' --- --- --- PLOTTING 3D [CORRELATIONS] --- --- --- '''

        from d3analysis import LOAD_RES_CORR
        d3_corr = LOAD_RES_CORR(sim)
        for task in __d3corrs__:
            plot_hostograms_from_3D(d1class, d3_corr, vn1vn2=task, rewrite=True)

        ''' --- --- --- PLOTTING 3D [DENSITY MODES] --- --- --- '''

        from d3analysis import LOAD_DENSITY_MODES
        dm_class = LOAD_DENSITY_MODES(sim)
        plot_density_modes(d1class, dm_class, fname=__d3dm__[0], overwrite=True)

    if task == "d2plots":
        assert len(times) > 0
        assert rl != -1
        d1class = ADD_METHODS_1D(sim)
        d2class = COMPUTE_STORE(sim)
        _iterations = []
        _times = []
        tmerg = d1class.get_par("tmerger_gw")
        for t in times:
            _t = (t / 1e3 + tmerg)
            _it = d2class.get_it_for_time(_t, d1d2d3="d2")
            _iterations.append(_it)
            _times.append(d2class.get_time_for_it(_it))
        plots_d2(d1class,d2class,v_n,_iterations, _times, rl, 'plots', overwrite)


    if task == "d2movies":
        d1class = ADD_METHODS_1D(sim)
        d2class = COMPUTE_STORE(sim)
        if click.confirm(Printcolor.yellow(
                'Produce ALL movies for {} iterations x {} '
                'reflevels ({} plots) Overwrite is set to: {}'
                .format(len(d2class.iterations), nlevels,
                        len(d2class.iterations) * nlevels, overwrite),
                comma=False), default=True):
            for rl_ in range(nlevels):
                for task in __d2movievns__:
                    movie_2d(d1class, d2class, rl=rl_, v_n=task,
                             rewritefigs=overwrite, rewritemovie=overwrite)

    if task == "d2movie":
        d1class = ADD_METHODS_1D(sim)
        d2class = COMPUTE_STORE(sim)
        if rl == -1:
            rl = 3
            Printcolor.blue("selected defauled reflevel -rl {}".format(rl))
        for task in __d2movievns__:
            movie_2d(d1class, d2class, rl=rl, v_n=v_n, rewritefigs=overwrite, rewritemovie=overwrite)

    """ ------------------------------------------------------------------------------------------------------------ """

    # if plot != "no" and not plot in __d1plots__:
    #     raise NameError("Plot: {} is not available. See available list:\n{}"
    #                     .format(task, __d1plots__))

    if plot != "no" and plot in __d1plots__:
        assert os.path.isfile(resdir+"res_1d/"+plot+".png")
        os.system("display {}".format(resdir+"res_1d/"+plot+".png"))
        exit(0)

    if plot != "no" and plot in general_plots:
        assert os.path.isfile(resdir+plot+".png")
        os.system("display {}".format(resdir+plot+".png"))
        exit(0)

    if plot != "no":
        raise NameError("Plot: {} is not available. See available list:\n{}"
                        .format(task, __d1plots__))