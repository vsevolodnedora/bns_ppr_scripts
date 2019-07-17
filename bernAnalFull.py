from __future__ import division
from sys import path
path.append('modules/')
import os
import numpy as np
from glob import glob
import click

from general import *
from lists import *
from filework import *
from units import time_constant, volume_constant, energy_constant

from overall import MODELS
models = MODELS()


loc_in_data = "/data1/numrel/WhiskyTHC/Backup/2018/GW170817/"
sims = ["DD2_M15091235_M0_LK_HR"]


fraction = 0.98
outfile="berntime.txt"
asci_file = "outflow_det_0.asc"
gws = os.listdir(Paths.gw170817) # list of sims to postrocess
pprs_selected = models.get_selected_models({"dyn_phase":"passed"})



# files to opy into the ppr tre
flist = ["*.done", "output-0000/parfile.par", outfile, "WARNING*"]
dirlist = ["outflow*", "collated*", "waveforms"]

def get_it_time_cum_mass(sim, asci_file="outflow_det_0.asc"):
    fpath = Paths.gw170817 + sim + '/' + "output-*" + '/data/' + asci_file
    files = glob(fpath)
    print("files: {}".format(len(files)))
    if len(files) == 0:
        raise ValueError("No files found for {} found in outputs"
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

def get_bern_time(simdir, fraction = 0.98, fnameout = Paths.gw170817 + isim + '/' + outfile):

    times, masses1, masses2 = get_it_time_cum_mass(simdir)

    idx_ = find_nearest_index(masses2, fraction * masses2.max())

    if idx_ == len(times) - 1:
        raise ValueError("{} of total mass is at the last timestep. Cannot do Bernoulli with that")

    if times[idx_] > 0.95 * times.max():
        Printcolor.yellow("Warning! berntime is {} total time".format(times[idx_] * 100 / times.max()))


    if save:
        Printcolor.blue("\tSaving {}".format(fnameout))
        if os.path.isfile(fnameout):
            Printcolor.yellow("Warning. File:{} already exist".format(fnameout))
        open(fnameout, "w").write("{}\n".format(float(times[idx_])))

    return float(times[idx_])




if __name__ == "__main__":
    import matplotlib.pyplot as plt

    for isim, sim in pprs_selected.iterrows():
        # print("isim:{}".format(isim))
        # print(gws)
        if isim in gws:

            # if not os.path.isdir(Paths.ppr_sims + isim + '/'):

            if not os.path.isdir(Paths.ppr_sims + isim + '/'):
                if click.confirm("\nFind berntime for {} ?".format(isim), default=True):

                    if not os.path.isdir(Paths.ppr_sims + isim + '/'):
                        Printcolor.blue("\tcreating {}".format(Paths.ppr_sims + isim + '/'))
                        os.mkdir(Paths.ppr_sims + isim + '/')

                    print("sim: {}".format(isim))
                    times, masses1, masses2 = get_it_time_cum_mass(isim)

                    idx_ = find_nearest_index(masses2, fraction * masses2.max())

                    if idx_ == len(times)-1:
                        raise ValueError("{} of total mass is at the last timestep. Cannot do Bernoulli with that")

                    if times[idx_] > 0.95 * times.max():
                        Printcolor.yellow("Warning! berntime is {} total time".format(times[idx_]*100 / times.max()))

                    Printcolor.blue("\tSaving {}".format(Paths.gw170817 + isim + '/' + outfile))
                    open(Paths.gw170817 + isim + '/' + outfile, "w").write("{}\n".format(float(times[idx_])))

                    Printcolor.blue("\tPlotting {}".format(Paths.ppr_sims+isim+'/{}.png'.format(str(outfile.split('.')[0]))))
                    plt.axvline(x=times[idx_] * 0.004925794970773136 * 1e-3, color='black', lw='0.5')
                    plt.plot(times  * 0.004925794970773136 * 1e-3, masses2, '-', color='black')
                    plt.plot(times  * 0.004925794970773136 * 1e-3, masses1, '-', color='gray')
                    plt.ylabel(r'M $[M_{\odot}$]', fontsize=12)
                    plt.xlabel(r'time $[M_{\odot}$]', fontsize=12)
                    plt.minorticks_on()
                    plt.xticks(fontsize=12)
                    plt.yticks(fontsize=12)
                    plt.title('Outflow', fontsize=20)
                    plt.legend(loc='upper right', numpoints=1)
                    plt.savefig(Paths.ppr_sims+isim+'/{}.png'.format(str(outfile.split('.')[0])), bbox_inches='tight', dpi=128)
                    plt.close()

                    if click.confirm("\nDo postprocessing?", default=True):
                        command = Paths.scripts + "analyze.sh" + ' ' + Paths.gw170817 + isim + '/'
                        Printcolor.blue("Initializing postprocessing...\n {}"
                              .format("# {}".format(command)))
                        os.system(command)

                    if click.confirm("\nCopy results?", default=True):
                        Printcolor.blue("Copying results from \n{} to {}".format(Paths.gw170817+isim+'/', Paths.ppr_sims + isim + '/'))

                        flist = ["*.done", "output-0000/parfile.par", outfile, "WARNING*"]



                        for file_ in flist:
                            os.system("cp " + Paths.gw170817+isim+'/'+file_ + ' ' + Paths.ppr_sims + isim + '/')
                        for dir_ in dirlist:
                            os.system("cp -r " + Paths.gw170817 + isim + '/' + dir_ + ' ' + Paths.ppr_sims + isim + '/')

                else:
                    Printcolor.yellow("Skipping {}".format(isim))
            else:
                Printcolor.blue("\tAlready done for {}".format(isim))

# print("N sims in GW170817: {}".format(len(gws)))
# print("N sims in postprocessed: {}".format(len(pprs)))


#
#
#
#
#
#
# def set_it_output_map(path_to_sim, file_for_it="dens.norm1.asc", clean=True):
#     """
#     Loads set of files that have '1:it 2:time ...' structure to get a map
#     of what output-xxxx contains what iteration (and time)
#     """
#
#
#
#
#
#
#
#
#
#
#
#
#     gen_set = {"indir": path_to_sim,
#                "file_for_it": file_for_it}
#
#     output_it_map = {}
#
#     fpath = gen_set["indir"] + "output-*" + "/data/" + gen_set['file_for_it']
#     files = glob(fpath)
#
#     if not clean:
#         print('-' * 25 + 'LOADING it list ({})'
#               .format(gen_set['file_for_it']) + '-' * 25)
#         print("\t loading from: {}".format(gen_set['indir']))
#
#     if len(files) == 0:
#         raise ValueError("No files found for it_time mapping searched:\n{}"
#                          .format(fpath))
#     if not clean:
#         print("\t   files found: {}".format(len(files)))
#
#
#
#     it_time = np.zeros(2)
#     for file in files:
#         o_name = file.split('/')
#         o_dir = ''
#         for o_part in o_name:
#             if o_part.__contains__('output-'):
#                 o_dir = o_part
#         if o_dir == '':
#             raise NameError("Did not find output-xxxx in {}".format(o_name))
#         it_time_i = np.loadtxt(file, usecols=(0, 1))
#         output_it_map[o_dir] = it_time_i
#         it_time = np.vstack((it_time, it_time_i))
#     it_time = np.delete(it_time, 0, 0)
#     if not clean:
#         print('outputs:{} iterations:{} [{}->{}]'.format(len(files),
#                                                          len(it_time[:, 0]),
#                                                          int(it_time[:, 0].min()),
#                                                          int(it_time[:, 0].max())))
#     if len(it_time[:, 0]) != len(set(it_time[:, 0])):
#         if not clean:
#             Printcolor.yellow("Warning: repetitions found in the "
#                               "loaded iterations")
#         iterations = np.unique(it_time[:, 0])
#         timestpes = np.unique(it_time[:, 1])
#         if not len(iterations) == len(timestpes):
#             raise ValueError("Failed attmept to remove repetitions from "
#                              "\t it and time lists. Wrong lengths: {} {}"
#                              .format(len(iterations), len(timestpes)))
#     else:
#         if not clean:
#             Printcolor.blue("No repetitions found in loaded it list")
#         iterations = np.unique(it_time[:, 0])
#         timestpes = np.unique(it_time[:, 1])
#
#     if not clean:
#         print('-' * 30 + '------DONE-----' + '-' * 30)
#
#     return output_it_map, np.vstack((iterations, timestpes)).T
#
#
#
# if __name__ == "__main__":
#
#     percentage = 98
#     berntimefname = "berntime.txt"