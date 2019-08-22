from __future__ import division
from sys import path
path.append('modules/')
from _curses import raw
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib import ticker
import matplotlib.pyplot as plt
# from matplotlib import rc
# plt.rc('text', usetex=True)
# plt.rc('font', family='serif')
# import scivis.units as ut # for tmerg
# import statsmodels.formula.api as smf
# import scipy.optimize as opt
# from math import pi, sqrt
# import matplotlib as mpl
# from glob import glob
# import pandas as pd
from glob import glob
import numpy as np
import itertools
import os.path
import cPickle
import time
import copy
import h5py
import csv
import os
import gc

from scidata.utils import locate
import scidata.carpet.hdf5 as h5
import scidata.xgraph as xg
from scidata.carpet.interp import Interpolator

from scipy import interpolate
# cmap = plt.get_cmap("viridis")
# from sklearn.linear_model import LinearRegression-
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


""" --- scidata problem --- """

# fpath = "/data1/numrel/WhiskyTHC/Backup/2018/SLy4_M130130_SR_physics/SLy4_M130130_M0_SR/tmp/output-0012/data/"
# v_n = "thc_M0_abs_energy.file_*.h5"
# v_nn = "/THC_LEAKAGEM0::thc_M0_abs_energy"
# it=475136
# flist = glob(fpath + v_n)
# print("\tfiles found: {}".format(len(flist)))
# # grid = h5.
# dset = h5.dataset(flist)
# reflevel = dset.get_reflevel(variable=v_nn, iteration=it, reflevel=0, timelevel=0)
# data = dset.get_reflevel_data(reflevel=reflevel, variable="/THC_LEAKAGEM0::thc_M0_ndens_nua", iteration=it, timelevel=0)
# print(data)
# exit(0)

""" --- COPY TARS FROM LRZ TO TULLIO """

list_lrz_sim = [
             ""
             "BHBlp_13641364_M0_LK_HR",
             "DD2_13641364_M0_LK_LR_PI",
             "DD2_15091235_M0_LK_HR",
             "LS220_13641364_M0_LK_HR",
             "LS220_14691268_M0_LK_HR",
             "SFHo_13641364_M0_LK_HR",
             "SFHo_14521283_M0_LK_HR_2",
             "SLy4_13641364_M0_LK_HR"]

dic_to_tul={
             "BHBlp_13641364_M0_LK_HR": "BHBlp_M13641364_M0_LK_HR",
             "DD2_13641364_M0_LK_LR_PI":"DD2_M13641364_M0_LK_LR_PI",
             "DD2_15091235_M0_LK_HR":   "DD2_M15091235_M0_LK_HR",
             "LS220_13641364_M0_LK_HR": "LS220_M13641364_M0_LK_HR",
             "LS220_14691268_M0_LK_HR": "LS220_M14691268_M0_LK_HR",
             "SFHo_13641364_M0_LK_HR":  "SFHo_M13641364_M0_LK_HR",
             "SFHo_14521283_M0_LK_HR_2":"SFHo_M14521283_M0_LK_HR",
             "SLy4_13641364_M0_LK_HR":  "SLy4_M13641364_M0_LK_HR"
}

lrz_sim_list = [
    "DD2_13641364_LK_SR_PI",
    "DD2_13641364_M0_LK_LR_PI",
    "DD2_13641364_M0_LK_SR_PI",
    "DD2_15091235_M0_LK_HR",
    "LS220_13641364_M0_LK_HR",
    "LS220_13641364_M0_LK_HR",
    "LS220_14691268_M0_LK_HR",
    "LS220_14691268_M0_LK_HR",
    "SFHo_14521283_M0_LK_HR_2",
    "SLy4_13641364_M0_LK_HR",
    "SLy4_13641364_M0_LK_LR"
]

# for sim in lrz_sim_list:
#     print("mkdir $WORK/old_simulations/{}".format(sim))


for sim in lrz_sim_list:
    print("rsync -arvP --append $SCRATCH/simulations/deep_simulations/{}/*.tar $WORK/old_simulations/{}/".format(sim, sim))
print("echo 'done'")
# for sim in lrz_sim_list:
#     print("find {}".format(sim)),
#     print(" -exec touch '{}' \;")
# print("echo 'done'")









print(4.925 * 10**(-6) * 299792458)
print(-26.4717 / 1.476)
print(18.5283 / 1.476)

os.system("display ~/Tullio/postprocessed3/LS220_M14691268_M0_LK_SR/res_3d/1490944/slices/Ye_rl0.png")
exit(1)





for i in [14, 20, 24, 28, 32, 38, 44, 50, 53]:
    print("tar -xvf output-{0:04d}.dat.tar --directory tmp/;".format(i))
print("echo 'done'")

string = "rsync -arvP --append di52pub@hw.supermuc.lrz.de"
for i in range(1, 11):
    string = string + ":/gss/scratch/pn56zo/di52pub/simulations/deep_simulations/SFHo_14521283_M0_LK_HR/output-{0:04d}.tar ".format(i)
string = string + "./"
print string



for i in range(56):
    print("tar -xvf output-{0:04d}.tar".format(i) +
          " --directory /data1/numrel/WhiskyTHC/Backup/2018/SLy4_M130130_SR_physics/SLy4_M130130_VM0_SR/;\n")

exit(1)
# 25641369fmjks
for sim in list_lrz_sim:
    try:
        os.system("rsync -arvP --append "
                  "di52pub@hw.supermuc.lrz.de:/gss/scratch/pn56zo/di52pub/simulations/deep_simulations/" + \
                  "{}/output-????.tar ".format(sim) + \
                  "/data1/numrel/WhiskyTHC/Backup/2018/GW170817/{}/".format(dic_to_tul[sim]))
    except:
        print("Error for {}".format(sim))

print("Finished copying tars")
exit(1)






# exit(1)

strain_fname = "strain_l2_m2.dat"
outfile_tmerg = "tmerger2.dat"
outfile_tcoll = "tcoll.dat"
plot_name = "tmergtcoll.png"



flist = ["hists_ej_profiles.png", "map_Jflux_m1Phase_3sim.png", "corr_Dunb_Jflux_3sims.png", "dens_modes_3sim.png"]
fpath = "/home/vsnedo/Tullio/figs/combine_test/"
for file_ in flist:
    print("cp {}{} ./ ; git add {}; git commit -m 'update'; git push;"
          .format(fpath, file_, file_))


for sim in os.listdir(Paths.ppr_sims):
    print(" --- {} --- ".format(sim))
    os.chdir(Paths.ppr_sims + sim + '/outflow_0/')
    os.system("/data01/numrel/vsevolod.nedora/scripts_server/ejecta/profile_for_mkn.py")

exit(1)

for sim in os.listdir(Paths.ppr_sims):
    print("Processing: {}".format(sim))

    if os.path.isdir(Paths.ppr_sims + sim + '/' + "waveforms/") and \
            os.path.isfile(Paths.ppr_sims + sim + '/' + "waveforms/" + strain_fname):

        w_dir = Paths.ppr_sims + sim + '/' + "waveforms/"
        strain_file = w_dir + strain_fname

        time, reh, imh = np.genfromtxt(strain_file, unpack=True, usecols=[0, 1, 2])
        h = reh + 1j * imh
        amp = np.abs(h)
        tmerg = time[np.argmax(amp)]
        time -= tmerg
        tmerg = 0
        maxA = np.max(amp)
        indx_collapse = np.min(np.where(amp < 0.01 * maxA))
        tcoll = time[indx_collapse]
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
        Printcolor.red("Error: no {} in wafeforms/ found".format(strain_file))





