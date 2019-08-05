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
from glob import glob
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

from d1analysis import LOAD_ITTIME, LOAD_FILES
from analyze import __d1plots__, __d2movievns__, __d3dmass__, __d3dm__, \
    __d3sliceplotvns__, __d3corrs__, __d3slicesplanes__, __d3slicesvns__


class PRINT_SIM_STATUS(LOAD_ITTIME):

    def __init__(self, sim):

        LOAD_ITTIME.__init__(self, sim)

        self.sim = sim

        self.path_in_data = Paths.gw170817 + sim + '/'
        self.prof_in_data = Paths.gw170817 + sim + '/profiles/3d/'
        self.path_out_data = Paths.ppr_sims + sim + '/'
        self.file_for_gw_time = "/data/dens.norm1.asc"
        self.file_for_ppr_time = "/collated/dens.norm1.asc"


        ''' --- '''

        tstart = 0.
        tend = 130.
        tstep = 1.
        prec = 0.5

        ''' --- PRINTING ---  '''
        print('='*100)
        print("<<< {} >>>".format(sim))
        # assert that the ittime.h5 file is upt to date
        isgood = self.assert_ittime()
        if not isgood:
            from analyze import SIM_STATUS
            SIM_STATUS(sim, save=True, clean=True)
            Printcolor.green("\tittime.h5 is updated")

        self.print_what_output_tarbal_dattar_present(comma=False)
        print("\tAsserting output contnet:")
        self.print_assert_tarball_content()
        print("\tAsserting data availability: ")
        
        tstart, tend = self.get_overall_tstart_tend()
        Printcolor.green("\tOverall Data span: {:.1f} to {:.1f} [ms]"
                         .format(tstart-1, tend-1))
        
        self.print_timemarks_output(start=tstart, stop=tend, tstep=tstep, precision=0.5)
        self.print_timemarks(start=tstart, stop=tend, tstep=tstep, tmark=10., comma=False)
        self.print_ititme_status("overall", d1d2d3prof="d1", start=tstart, stop=tend, tstep=tstep, precision=prec)
        self.print_ititme_status("overall", d1d2d3prof="d2", start=tstart, stop=tend, tstep=tstep, precision=prec)
        self.print_ititme_status("overall", d1d2d3prof="d3", start=tstart, stop=tend, tstep=tstep, precision=prec)
        self.print_ititme_status("profiles", d1d2d3prof="prof", start=tstart, stop=tend, tstep=tstep, precision=prec)

        # self.print_gw_ppr_time(comma=True)
        # self.print_assert_collated_data()
        #
        # self.print_assert_outflowed_data(criterion="_0")
        # self.print_assert_outflowed_data(criterion="_0_b_w")
        # self.print_assert_outflowed_corr_data(criterion="_0")
        # self.print_assert_outflowed_corr_data(criterion="_0_b_w")
        # self.print_assert_gw_data()
        # self.print_assert_mkn_data("_0")
        # self.print_assert_mkn_data("_0_b_w")
        #
        # self.print_assert_d1_plots()
        # self.print_assert_d2_movies()


    def get_tars(self):
        tars = glob(self.path_in_data + 'output-????.tar')
        tars = [str(tar.split('/')[-1]).split('.tar')[0] for tar in tars]
        return tars

    def get_dattars(self):
        dattars = glob(self.path_in_data + 'output-????.dat.tar')
        dattars = [str(dattar.split('/')[-1]).split('.dat.tar')[0] for dattar in dattars]
        return dattars

    def get_outputdirs(self):

        def get_number(output_dir):
            return int(str(output_dir.split('/')[-1]).split("output-")[-1])

        dirs = os.listdir(self.path_in_data)
        output_dirs = []
        for dir_ in dirs:
            if str(dir_).__contains__("output-") and \
                    not str(dir_).__contains__('.tar') and \
                    not str(dir_).__contains__('.dat.tar'):
                output_dirs.append(dir_)
        output_dirs.sort(key=get_number)
        return output_dirs

    def get_outputs(self):
        return [str(output_dir.split('/')[-1]) for output_dir in self.get_outputdirs()]

    def get_profiles(self):
        if not os.path.isdir(self.prof_in_data):
            return 0, np.empty(0,), np.empty(0,)
        profiles = glob(self.prof_in_data + "*.h5")
        return profiles

    def assert_ittime(self):

        new_output_dirs = self.get_outputdirs()
        new_outputs = [str(output) for output in new_output_dirs]
        old_outputs = self.get_list_outputs()

        if sorted(old_outputs) == sorted(new_outputs):
            last_output = list(new_output_dirs)[-1]
            it_time_i = np.loadtxt(self.path_in_data + last_output + self.file_for_gw_time, usecols=(0, 1))
            new_it_end = int(it_time_i[-1, 0])
            _, itd1, _= self.get_ittime("overall", d1d2d3prof="d1")
            old_it_end = itd1[-1]

            if int(new_it_end) == int(old_it_end):
                is_up_to_data = True

                new_profiles = glob(self.prof_in_data + "*.h5")
                _, itprofs, _ = self.get_ittime("profiles", d1d2d3prof="prof")
                if len(new_profiles) == len(itprofs):
                    Printcolor.green("\tittime.h5 is up to date")
                else:
                    Printcolor.red("\tittime.h5 is NOT up to date: profiles (old{:d} != new{:d})"
                                   .format(len(itprofs), len(new_profiles)))
            else:
                is_up_to_data = False
                Printcolor.red("\tittime.h5 is NOT up to date: d1 iterations (old{:d} != new{:d})"
                               .format(old_it_end, new_it_end))
        else:
            Printcolor.red("\tittime.h5 is NOT up to date: outputs: (old{} != new{})"
                           .format(old_outputs[-1], new_outputs[-1]))
            return False

        return is_up_to_data

    def get_overall_tstart_tend(self):

        t1, t2 = [], []
        isd1, itd1, td1 = self.get_ittime("overall", d1d2d3prof="d1")
        isd2, itd2, td2 = self.get_ittime("overall", d1d2d3prof="d2")
        isd3, itd3, td3 = self.get_ittime("overall", d1d2d3prof="d3")
        isprof, itprof, tprof = self.get_ittime("profiles", d1d2d3prof="prof")

        if len(td1)>0:
            t1.append(td1[0])
            t2.append(td1[-1])
        if len(td2)>0:
            t1.append(td2[0])
            t2.append(td2[-1])
        if len(td3)>0:
            t1.append(td3[0])
            t2.append(td3[-1])
        if len(tprof)>0:
            t1.append(tprof[0])
            t2.append(tprof[-1])

        return np.array(t1).min() * 1e3 + 1, np.array(t2).max() * 1e3 + 1

    ''' --- '''

    def print_what_output_tarbal_dattar_present(self, comma=False):
        print("\toutputs: "),
        Printcolor.blue(str(len(self.get_outputs())), comma=True)
        print("\ttars: "),
        Printcolor.yellow(str(len(self.get_tars())), comma=True)
        print("\tdattars: "),
        Printcolor.red(str(len(self.get_dattars())), comma=True)
        print("\tprofiles: "),
        Printcolor.blue(str(len(self.get_profiles())), comma=True)
        if comma: print(' '),
        else: print(' ')


    ''' --- '''

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

    def print_assert_data_status(self, name, path, flist, comma=True):

        Printcolor.blue("\t{}: ".format(name), comma=True)
        # flist = copy.deepcopy(LOAD_FILES.list_collated_files)

        status, missing = self.print_assert_content(path, flist)

        if status == "full":
            Printcolor.green(" complete", comma=True)
        elif status == "partial":
            Printcolor.yellow(" partial, ({}) missing".format(len(missing)), comma=True)
        else:
            Printcolor.red(" absent", comma=True)

        if comma: print(' '),
        else: print(' ')

        return status, missing

    def print_assert_tarball_content(self, comma=False):

        outputs = self.get_outputdirs()
        for output in outputs:

            output = self.path_in_data + output
            assert os.path.isdir(output)
            output_n = int(str(output.split('/')[-1]).split('output-')[-1])
            n_files = len([name for name in os.listdir(output + '/data/')])
            Printcolor.blue("\toutput: {0:03d}".format(output_n), comma=True)
            # print('('),
            if n_files == 259 or n_files == 258:
                Printcolor.green("{0:05d} files".format(n_files), comma=True)
            else:
                Printcolor.yellow("{0:05d} files".format(n_files), comma=True)
            # print(')'),
            status, missing = self.print_assert_content(output + '/data/', Lists.tarball)
            if status == "full":
                Printcolor.green(" complete", comma=True)
            elif status == "partial":
                Printcolor.yellow(" partial, ({}) missing".format(len(missing)), comma=True)
            else:
                Printcolor.red(" absent", comma=True)
            print('')

        if comma:
            print(' '),
        else:
            print(' ')

    def print_timemarks(self, start=0., stop=30., tstep = 1., tmark=10., comma=False):

        trange = np.arange(start=start, stop=stop, step=tstep)
        Printcolor.blue("\tTimesteps {}ms ".format(tmark, tstep), comma=True)
        print('['),
        for t in trange:
            if t % tmark == 0:
                print("{:d}".format(int(t / tmark))),
            else:
                print(' '),
        print(']'),
        if comma:
            print(' '),
        else:
            print(' ')

    def print_timemarks_output(self, start=0., stop=30., tstep = 1., comma=False, precision=0.5):

        tstart = []
        tend = []
        dic_outend = {}
        for output in self.get_outputs():
            isdata, itd1, td1 = self.get_ittime(output=output, d1d2d3prof="d1")
            tstart.append(td1[0] * 1e3)
            tend.append(td1[-1] * 1e3)
            dic_outend["%.3f"%(td1[-1]*1e3)] = output.split("output-")[-1]



        for digit, letter, in zip(range(4), ['o', 'u', 't', '-']):
            print("\t       {}         ".format(letter)),
            # Printcolor.blue("\tOutputs end [ms] ", comma=True)
            trange = np.arange(start=start, stop=stop, step=tstep)
            print('['),
            for t in trange:
                tnear = tend[self.find_nearest_index(tend, t)]
                if abs(tnear - t) < precision:  # (tnear - t) >= 0
                    output = dic_outend["%.3f"%tnear]
                    numbers = []
                    for i in [0,1,2,3]:
                        numbers.append(str(output[i]))

                    if digit != 3 and int(output[digit]) == 0:
                        print(' '),
                        # Printcolor.blue(output[digit], comma=True)
                    else:
                        Printcolor.blue(output[digit], comma=True)


                    # for i in range(len(numbers)-1):
                    #     if numbers[i] == "0" and numbers[i+1] != "0":
                    #         Printcolor.blue(numbers[i], comma=True)
                    #     else:
                    #         Printcolor.yellow(numbers[i], comma=True)
                    # print("%.2f"%tnear, t)
                else:
                    print(' '),
            print(']')


    def print_ititme_status(self, output, d1d2d3prof, start=0., stop=30., tstep = 1., precision=0.5):

        isdi1, itd1, td = self.get_ittime(output, d1d2d3prof=d1d2d3prof)
        td = td * 1e3  # ms
        # print(td); exit(1)
        # trange = np.arange(start=td[0], stop=td[-1], step=tstep)
        trange = np.arange(start=start, stop=stop, step=tstep)

        _name_ = ''
        if d1d2d3prof == 'd1': _name_ = "D1  "
        elif d1d2d3prof == "d2": _name_ = "D2  "
        elif d1d2d3prof == "d3": _name_ = "D3  "
        elif d1d2d3prof == "prof": _name_ = "prof"

        # print(td)

        if len(td) > 0:
            Printcolor.blue("\tTime {} [{}ms]".format(_name_, tstep), comma=True)
            print('['),
            for t in trange:
                tnear = td[self.find_nearest_index(td, t)]
                if abs(tnear - t) < precision: # (tnear - t) >= 0
                    Printcolor.green('.', comma=True)
                    # print("%.2f"%tnear, t)
                else:
                    print(' '),
                    # print("%.2f"%tnear, t)

            print(']'),
            Printcolor.green("{:.1f}ms".format(td[-1]), comma=False)
        else:
            Printcolor.red("\tTime {} No Data".format(_name_), comma=False)

        # ---

        # isdi2, itd2, td2 = self.get_ittime("overall", d1d2d3prof="d2")
        # td2 = td2 * 1e3  # ms
        # trange = np.arange(start=td2[0], stop=td2[-1], step=tstep)
        #
        # Printcolor.blue("\tTime 2D [1ms]", comma=True)
        # print('['),
        # for t in trange:
        #     tnear = td2[self.find_nearest_index(td2, t)]
        #     if abs(tnear - t) < tstep:
        #         Printcolor.green('.', comma=True)
        # print(']'),
        # Printcolor.green("{:.1f}ms".format(td2[-1]), comma=False)
        #
        #
        # exit(1)
        #
        # isdi1, itd1, td = self.get_ittime("overall", d1d2d3prof="d1")
        # td = td * 1e3 # ms
        # # print(td); exit(1)
        # Printcolor.blue("\tTime 1D [1ms]", comma=True)
        # n=1
        # print('['),
        # for it, t in enumerate(td[1:]):
        #     # tcum = tcum + td[it]
        #     # print(tcum, tstart + n*tstep)
        #     if td[it] > n*tstep:
        #         Printcolor.green('.', comma=True)
        #         n = n+1
        # print(']'),
        # Printcolor.green("{:.1f}ms".format(td[-1]), comma=False)
        #
        # isd2, itd2, td2 = self.get_ittime("overall", d1d2d3prof="d2")
        # td2 = td2 * 1e3 # ms
        # # print(td); exit(1)
        # Printcolor.blue("\tTime 2D [1ms]", comma=True)
        # n=1
        # print('['),
        # for it, t in enumerate(td2[1:]):
        #     # tcum = tcum + td[it]
        #     # print(tcum, tstart + n*tstep)
        #     if td2[it] > n*tstep:
        #         Printcolor.green('.', comma=True)
        #         n = n+1
        # print(']'),
        # Printcolor.green("{:.1f}ms".format(td2[-1]), comma=False)

    def print_ititme_status_(self, tstep = 1.):

        isdi1, itd1, td1 = self.get_ittime("overall", d1d2d3prof="d1")
        td1 = td1 * 1e3 # ms
        # print(td1); exit(1)
        Printcolor.blue("\tTime 1D [1ms]", comma=True)
        n=1
        print('['),
        for it, t in enumerate(td1[1:]):
            # tcum = tcum + td1[it]
            # print(tcum, tstart + n*tstep)
            if td1[it] > n*tstep:
                Printcolor.green('.', comma=True)
                n = n+1
        print(']'),
        Printcolor.green("{:.1f}ms".format(td1[-1]), comma=False)

        isd2, itd2, td2 = self.get_ittime("overall", d1d2d3prof="d2")
        td2 = td2 * 1e3 # ms
        # print(td1); exit(1)
        Printcolor.blue("\tTime 2D [1ms]", comma=True)
        n=1
        print('['),
        for it, t in enumerate(td2[1:]):
            # tcum = tcum + td1[it]
            # print(tcum, tstart + n*tstep)
            if td2[it] > n*tstep:
                Printcolor.green('.', comma=True)
                n = n+1
        print(']'),
        Printcolor.green("{:.1f}ms".format(td2[-1]), comma=False)





    # def print_assert_outflowed_data(self, criterion):
    #
    #     flist = copy.deepcopy(LOAD_FILES.list_outflowed_files)
    #     if not criterion.__contains__("_b"):
    #         # if the criterion is not Bernoulli
    #         flist.remove("hist_vel_inf_bern.dat")
    #         flist.remove("ejecta_profile_bern.dat")
    #
    #     outflow_status, outflow_missing = \
    #         self.__assert_content(Paths.ppr_sims + self.sim + "/outflow{}/".format(criterion),
    #                               flist)
    #
    #     return outflow_status, outflow_missing

class PRINT_SIM_PPR_STATUS(LOAD_ITTIME):

    def __init__(self, sim):

        LOAD_ITTIME.__init__(self, sim)

        self.sim = sim

        self.path_in_data = Paths.gw170817 + sim + '/'
        self.prof_in_data = Paths.gw170817 + sim + '/profiles/3d/'
        self.path_out_data = Paths.ppr_sims + sim + '/'
        self.file_for_gw_time = "/data/dens.norm1.asc"
        self.file_for_ppr_time = "/collated/dens.norm1.asc"

        print(" ")
        print("\tAsserting postprocessing:")

        self.tgw, self.tppr = self.print_gw_ppr_time()
        self.print_assert_collated_data()
        self.print_assert_gw_data()
        self.print_assert_outflow_data(criterion="_0")
        self.print_assert_outflow_data(criterion="_0_b_w")


        self.print_assert_d1_plots(dir_="res_1d/")
        self.assert_d2_movies(dir_="res_2d/")
        self.assert_d3_data()

        print('='*100)
        print('\n')

    def print_gw_ppr_time(self, comma=True):

        # --- --- gw ---
        # last_output = list(self.get_outputdirs())[-1]


        # it_time_i = np.loadtxt(self.path_in_data+last_output+self.file_for_gw_time, usecols=(0, 1))
        # it_time_i[:, 1] *= 0.004925794970773136 * 1e-3  # ms
        # gw_iterations = np.array(it_time_i[:, 0], dtype=int)
        # gw_timesteps = np.array(it_time_i[:, 1], dtype=float)

        isdata, gw_iterations, gw_timesteps = self.get_ittime("overall", d1d2d3prof="d1")


        # --- --- ppr ----
        it_time_i = np.loadtxt(self.path_out_data+self.file_for_ppr_time, usecols=(0, 1))
        it_time_i[:, 1] *= 0.004925794970773136 * 1e-3  # ms
        ppr_iterations = np.array(it_time_i[:, 0], dtype=int)
        ppr_timesteps = np.array(it_time_i[:, 1], dtype=float)

        print('\t'),

        if int(gw_iterations[-1]) > int(ppr_iterations[-1]):
            Printcolor.red(  "t(ititme): {:.1f} > {:.1f} t(collated) [ms] ppr required"
                           .format(gw_timesteps[-1]*1e3, ppr_timesteps[-1]*1e3))
        elif int(gw_iterations[-1]) < int(ppr_iterations[-1]):
            Printcolor.red(  "t(ititme): {:.1f} < {:.1f} t(collated) [ms] ppr's longer"
                           .format(gw_timesteps[-1]*1e3, ppr_timesteps[-1]*1e3))
        else:
            Printcolor.green("t(ititme): {:.1f} = {:.1f} t(collated) [ms] - up to date"
                           .format(gw_timesteps[-1]*1e3, ppr_timesteps[-1]*1e3))

        if comma: print(' '),
        else: print(' ')

        return gw_timesteps[-1],  ppr_timesteps[-1]

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
                missing_files.append(file2)
                status = "partial"
            else:
                Printcolor.red(marker3, comma=True)
                status = "partial"
                missing_files.append(file1)
                missing_files.append(file2)
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
            flist.remove("corr_vel_inf_bern_theta.h5",)
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

    def print_assert_d1_plots(self,  dir_="res_1d/", res=".png", comma=False):
        dir = self.path_out_data+dir_
        flist = [plot_+ res for plot_ in __d1plots__]
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

    def assert_d2_movies(self, dir_="res_2d/", extension=".mp4", comma=False):

        isd2, itd2, td2 = self.get_ittime("overall",d1d2d3prof="d2")

        Printcolor.blue("\tD2 movies ", comma=True)
        v_ns = __d2movievns__
        print("["),
        missing = []
        outdated = []
        for v_n in v_ns:
            moviedir = str(v_n) + "_movie" + '/'
            # print(self.path_out_data+dir_+moviedir)
            assert os.path.isdir(self.path_out_data+dir_+moviedir)
            max_file = max(glob(self.path_out_data+dir_+moviedir+"???????.png"))
            int_max_file = int(str(max_file.split('/')[-1]).split('.png')[0])
            is_movie = os.path.isfile(self.path_out_data+dir_+moviedir+v_n+extension)
            if is_movie and itd2[-1] == int_max_file:
                Printcolor.green(".", comma=True)
            elif is_movie and itd2 != int_max_file:
                Printcolor.yellow(".", comma=True)
                outdated.append(v_n)
            else:
                Printcolor.red("x", comma=True)
                missing.append(v_n)
        print("]")

    def assert_d3_data(self, dir_="res_3d/", comma=False):

        isprof, itprof, tprof = self.get_ittime("profiles", d1d2d3prof="prof")

        Printcolor.blue("\tD3 ({}) profiles: ".format(len(tprof)))
        for it, t in zip(itprof, tprof):
            Printcolor.green("\t {} ({:.1f}[ms])".format(it, t*1e3), comma=True)
            # disk mass:
            flist = [v_n + ".txt" for v_n in __d3dmass__]
            data_dir = self.path_out_data + dir_ + str(it) + "/"
            status, missing_files = self.print_assert_content(data_dir, flist)
            # profile.xy.h5
            flist = ["profile." + v_n + ".h5" for v_n in __d3slicesplanes__]
            data_dir = self.path_out_data + dir_ + str(it) + "/"
            status, missing_files = self.print_assert_content(data_dir, flist)
            # corr_vn1_vn2.h5
            flist1 = ["corr_" + v_n + ".h5" for v_n in __d3corrs__]
            flist2 = [v_n + ".png" for v_n in __d3corrs__]
            data_dir = self.path_out_data + dir_ + str(it) + "/"
            status, missing_files1, missing_files2 = self.print_assert_contents(data_dir, flist1, flist2)
            # corr_vn1_vn2.h5
            flist = [v_n + "_rl0.png" for v_n in __d3sliceplotvns__]
            data_dir = self.path_out_data + dir_ + str(it) + "/slices/"
            status, missing_files = self.print_assert_content(data_dir, flist)
            print(' ')

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


class SIM_PPR_STATUS:

    def __init__(self, sim):
        self.sim = sim
        # LOAD_ITTIME.__init__(self, sim)
        self.res_3d_itdirs = get_list_iterationsfrom_res_3d(self.sim)

    def __assert_content(self, dir, expected_files):
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

        for file_ in expected_files:
            if os.path.isfile(dir + file_):
                pass
            else:
                status = "partial"
                missing_files.append(file_)
        if len(missing_files) == len(expected_files):
            status = "empty"

        return status, missing_files

    def assert_collated_data(self):

        flist = copy.deepcopy(LOAD_FILES.list_collated_files)

        collated_status, collated_missing = \
            self.__assert_content(Paths.ppr_sims + self.sim + '/collated/',
                                  flist)

        return collated_status, collated_missing

    def assert_outflowed_data(self, criterion):

        flist = copy.deepcopy(LOAD_FILES.list_outflowed_files)
        if not criterion.__contains__("_b"):
            # if the criterion is not Bernoulli
            flist.remove("hist_vel_inf_bern.dat")
            flist.remove("ejecta_profile_bern.dat")

        outflow_status, outflow_missing = \
            self.__assert_content(Paths.ppr_sims + self.sim + "/outflow{}/".format(criterion),
                                  flist)

        return outflow_status, outflow_missing

    def assert_outflowed_corr_data(self, criterion):

        flist = copy.deepcopy(LOAD_FILES.list_of_outflowed_h5_files)
        if not criterion.__contains__("_b"):
            # if the criterion is not Bernoulli
            flist.remove("corr_vel_inf_bern_theta.h5")

        outflow_status, outflow_missing = \
            self.__assert_content(Paths.ppr_sims + self.sim + "/outflow{}/".format(criterion),
                                  flist)

        return outflow_status, outflow_missing

    def assert_gw_data(self):

        flist = copy.deepcopy(LOAD_FILES.list_gw_files)

        gw_status, gw_missing = \
            self.__assert_content(Paths.ppr_sims + self.sim + '/waveforms/',
                                  flist)

        return gw_status, gw_missing

    def assert_mkn_data(self, criterion):

        flist = copy.deepcopy(LOAD_FILES.list_mkn_files)
        # flist.remove("AT2017gfo.h5")

        mkn_status, mkn_missing = \
            self.__assert_content(Paths.ppr_sims + self.sim + "/outflow{}/".format(criterion),
                                  flist)

        return mkn_status, mkn_missing

    # ---

    def assert_d1_plots(self, criterion="_0", dir="res_1d/", res=".png"):

        flist = [plot_+res for plot_ in __d1plots__]

        plot_status, missing_plots = \
            self.__assert_content(Paths.ppr_sims+self.sim+'/'+dir,
                                  flist)

        return plot_status, missing_plots

    def assert_d2_movies(self, dir="res_2d/", res=".mp4"):

        flist = [v_n + "_movie/" + v_n + res for v_n in __d2movievns__]

        movie_status, missing_movies = \
            self.__assert_content(Paths.ppr_sims+self.sim+'/'+dir,
                                  flist)

        return movie_status, missing_movies

    def assert_d3_dmass(self, dir="res_3d/", res=".txt"):

        flist = [str(int(itdir)) + '/' + __d3dmass__[0] + res
                 for itdir in self.res_3d_itdirs]

        dmass_status, dmass_missing = self.__assert_content(Paths.ppr_sims + self.sim + '/' + dir,
                                                            flist)

        return dmass_status, dmass_missing

    def assert_d3_corrs(self, dir="res_3d/", intro="corr_", res=".h5"):

        flist = []
        for it in self.res_3d_itdirs:
            for v_n in __d3corrs__:
                flist.append(str(int(it))+'/'+intro+v_n+res)

        corr_status, corr_missing = self.__assert_content(Paths.ppr_sims + self.sim + '/' + dir,
                                                            flist)

        return corr_status, corr_missing

    def assert_d3_corr_plots(self, dir="res_3d/", res=".png"):

        flist = []
        for it in self.res_3d_itdirs:
            for v_n in __d3corrs__:
                flist.append(str(int(it)) + '/' + v_n + res)

        corr_plot_status, corr_plot_missing = self.__assert_content(Paths.ppr_sims + self.sim + '/' + dir,
                                                          flist)

        return corr_plot_status, corr_plot_missing

    def assert_d3_slices(self, dir="res_3d/", intro="profile.", res=".h5"):

        flist = []
        for it in self.res_3d_itdirs:
            for plane in __d3slicesplanes__:
                flist.append(str(int(it)) + '/' + intro + plane + res)

        slices_status, slices_missing = self.__assert_content(Paths.ppr_sims + self.sim + '/' + dir,
                                                                    flist)

        return slices_status, slices_missing

    def assert_d3_densmodes(self, dir="res_3d/", res=".h5"):

        flist = [file_+res for file_ in __d3dm__]

        dm_status, dm_missing = self.__assert_content(Paths.ppr_sims+self.sim+'/'+dir,
                                                      flist)

        return dm_status, dm_missing

    def assert_d3_slice_plots(self, dir="res_3d", ddir = "slices/", rl=0, res=".png"):

        flist = []
        for it in self.res_3d_itdirs:
            for v_n in __d3sliceplotvns__:
                flist.append(str(int(it)) + '/' + ddir + v_n + "_rl" + str(rl) + res)

        slices_plots_status, slices_plots_missing = self.__assert_content(Paths.ppr_sims + self.sim + '/' + dir,
                                                              flist)

        return slices_plots_status, slices_plots_missing















''' --- --- --- '''

# asci_file = "dens.norm1.asc"
# gws = os.listdir(Paths.gw170817)
# pprs = os.listdir(Paths.ppr_sims)
#
# print("N sims in GW170817: {}".format(len(gws)))
# print("N sims in postprocessed: {}".format(len(pprs)))
#
# def get_outputs(sim):
#     _outputs = os.listdir(Paths.gw170817 + sim + '/')
#     outputs = []
#     tars = []
#     dattar = []
#     for output in _outputs:
#
#         if output.__contains__("output-") and output.__contains__(".dat.tar"):
#             dattar.append(output)
#         elif output.__contains__("output-") and output.__contains__(".tar") \
#                 and not output.__contains__(".dat.tar"):
#             tars.append(output)
#         elif output.__contains__("output-") \
#                 and not output.__contains__(".tar") \
#                 and not output.__contains__(".dat.tar"):
#             outputs.append(output)
#         else:
#             pass
#
#     return outputs, tars, dattar
#
# def get_last_it_gw(sim, fname=asci_file):
#     _, it_time = set_it_output_map(Paths.gw170817 + sim + '/', fname)
#     return it_time[:,0].max(), it_time[:,1].max() *  0.004925794970773136 * 1e-3
#
# def get_last_it_ppr(sim, fname="dens.norm1.asc"):
#     path = MakePath.collated(sim)
#     if not os.listdir(path):
#         raise ImportError("directory: {} not found".format(path))
#     if not os.path.isfile(path+fname):
#         raise ImportError("file: {} does not exist".format(path + fname))
#
#     it_time = np.loadtxt(path+fname, usecols=(0, 1))
#
#     return it_time[:,0].max(), it_time[:,1].max() *  0.004925794970773136 * 1e-3
#
# def analyze_sim(sim, time_limit=0.030):
#
#     if not sim in pprs:
#         Printcolor.red(", NOT in ppr", True)
#
#         it, t = get_last_it_gw(sim)
#
#         if t < time_limit:
#             Printcolor.red(", tend {:.3f}s".format(t), True)
#         else:
#             Printcolor.green(", tend {:.3f}s".format(t), True)
#         #
#         # print(", GW [it:{:d} time:{:.3f}]".format(int(it), t)),
#
#     else:
#         Printcolor.blue(", in ppr", True)
#
#         it, t = get_last_it_gw(sim)
#         # print(", GW [it:{:d} time:{:.3f}]".format(int(it), t)),
#         it_, t_ = get_last_it_ppr(sim)
#         # print(", PPR: [it:{:d} time:{:.3f}]".format(int(it_), t_)),
#         if int(it) != int(it_):
#             Printcolor.red(", ppr REQUIRED [t_gw:{:.3f} t_ppr:{:.3f}]".format(t, t_), True),
#         else:
#             Printcolor.green(", ppr is complete", ',')
#
#         if t < time_limit:
#             Printcolor.red(", tend {:.3f}s".format(t), True)
#         else:
#             Printcolor.green(", tend {:.3f}s".format(t), True)
#
#         if not os.path.isfile(Paths.ppr_sims+sim+"/parfile.par"):
#             Printcolor.red(", parfile not found |", True)
#
#         outflowed = os.path.isdir(MakePath.outflow(sim, "_0"))
#         if outflowed:
#             Printcolor.green(" _0 done", True)
#         else:
#             Printcolor.red(" _0 NOT done", True)
#
#         if not sim in Lists.dyn_not_pas:
#             Printcolor.green(" Dyn pass", True)
#         else:
#             Printcolor.red(" Dyn NOT pass", True)
#
#         if sim in Lists.bern_pass:
#             Printcolor.green(" Bern pass", True)
#         else:
#             Printcolor.yellow(" Bern NOT pass", True)
#
#
#         outflowed_b_w = os.path.isdir(MakePath.outflow(sim, "_0_b_w"))
#         if outflowed_b_w:
#             Printcolor.green(", _0_b_w done", True)
#         else:
#             Printcolor.red(", _0_b_w NOT done", True)
#
# for sim in gws:
#     # print("\n------------------------------------------")
#     name = str(sim)
#     def_length = 30
#     name += ' ' * (def_length-len(name))
#     print("{} ".format(name)),
#     outputs, tars, dattars = get_outputs(sim)
#     # print("\toutput({}) .tar({}) .dat.tar({})".format(len(outputs), len(tars), len(dattars)))
#
#     if len(outputs) > 0:
#         Printcolor.blue("\t{} outputs".format(len(outputs)), True)
#
#         # analyze_sim(sim)
#
#         try:
#             analyze_sim(sim)
#         except ValueError:
#             Printcolor.red("Error. FAILD with ValueError", True)
#     else:
#         Printcolor.red("\tNO outputs", True)
#
#     print(" |")
#
# if len(gws) < len(pprs):
#     Printcolor.red("NOT in gw170717:")
#     for sim in pprs:
#         if not sim in gws:
#             Printcolor.red("\t{}".format(sim))


''' --- '''
print('\n')
''' --- '''





# if len(gws) > len(pprs):
#     Printcolor.yellow("NOT in postprocessed:")
#     for sim in gws:
#         if not sim in pprs:
#             print("\t\t{}".format(sim)),
#             _outputs = os.listdir(Paths.gw170817 + sim + '/')
#             outputs = []
#             for output in _outputs:
#                 if output.__contains__("output-"):
#                     outputs.append(output)
#             print(" {} outputs, ".format(len(outputs)))

if __name__ == '__main__':
    "--- --- TEST --- --- "

    PRINT_SIM_STATUS("DD2_M15091235_M0_LK_HR")
    PRINT_SIM_PPR_STATUS("DD2_M15091235_M0_LK_HR")
    exit(1)

    ppr_status = SIM_PPR_STATUS("SFHo_M14521283_M0_LK_HR")
    print("Collated:           {}".format(ppr_status.assert_collated_data()))
    print("Outflow _0          {}".format(ppr_status.assert_outflowed_data('_0')))
    print("Outflow _0_b_w      {}".format(ppr_status.assert_outflowed_data('_0_b_w')))
    print("Outflow corr _0     {}".format(ppr_status.assert_outflowed_corr_data('_0')))
    print("Outflow corr _0_b_w {}".format(ppr_status.assert_outflowed_corr_data('_0_b_w')))
    print("GW                  {}".format(ppr_status.assert_gw_data()))
    print("mkn _0              {}".format(ppr_status.assert_mkn_data("_0")))
    print("mkn _0_b_w          {}".format(ppr_status.assert_mkn_data("_0_b_w")))
    print('-'*50)
    print("D1 plots            {}".format(ppr_status.assert_d1_plots("_0", "res_1d/", ".png")))
    print("D2 movies           {}".format(ppr_status.assert_d2_movies("res_2d/", ".mp4")))
    print("D3 disk mass        {}".format(ppr_status.assert_d3_dmass("res_3d/", ".txt")))
    print("D3 correlations     {}".format(ppr_status.assert_d3_corrs("res_3d/", "corr_", ".h5")))
    print("D3 corr. plots      {}".format(ppr_status.assert_d3_corr_plots("res_3d/", ".png")))
    print("D3 prof.slices      {}".format(ppr_status.assert_d3_slices("res_3d/", "profile.", ".h5")))
    print("D3 dens.modes       {}".format(ppr_status.assert_d3_densmodes("res_3d/", ".h5")))
    print("D3 dens.modes plot  {}".format(ppr_status.assert_d3_densmodes("res_3d/", ".png")))
    print("D3 prof.slice.plots {}".format(ppr_status.assert_d3_slice_plots("res_3d/", "slices/",0, ".png")))