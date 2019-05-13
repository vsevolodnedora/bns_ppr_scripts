
from __future__ import division
from sys import path
path.append('modules/')
import os
import numpy as np
from glob import glob

from general import *
from lists import *
from filework import *
from units import time_constant, volume_constant, energy_constant

asci_file = "dens.norm1.asc"
gws = os.listdir(Paths.gw170817)
pprs = os.listdir(Paths.ppr_sims)

print("N sims in GW170817: {}".format(len(gws)))
print("N sims in postprocessed: {}".format(len(pprs)))

def get_outputs(sim):
    _outputs = os.listdir(Paths.gw170817 + sim + '/')
    outputs = []
    tars = []
    dattar = []
    for output in _outputs:

        if output.__contains__("output-") and output.__contains__(".dat.tar"):
            dattar.append(output)
        elif output.__contains__("output-") and output.__contains__(".tar") \
                and not output.__contains__(".dat.tar"):
            tars.append(output)
        elif output.__contains__("output-") \
                and not output.__contains__(".tar") \
                and not output.__contains__(".dat.tar"):
            outputs.append(output)
        else:
            pass

    return outputs, tars, dattar

def get_last_it_gw(sim, fname=asci_file):
    _, it_time = set_it_output_map(Paths.gw170817 + sim + '/', fname)
    return it_time[:,0].max(), it_time[:,1].max() *  0.004925794970773136 * 1e-3

def get_last_it_ppr(sim, fname="dens.norm1.asc"):
    path = MakePath.collated(sim)
    if not os.listdir(path):
        raise ImportError("directory: {} not found".format(path))
    if not os.path.isfile(path+fname):
        raise ImportError("file: {} does not exist".format(path + fname))

    it_time = np.loadtxt(path+fname, usecols=(0, 1))

    return it_time[:,0].max(), it_time[:,0].max() *  0.004925794970773136 * 1e-3

def analyze_sim(sim):

    if not sim in pprs:
        Printcolor.red("\t{} NOT in postprocessed2".format(sim))
        print("\t\t{} outputs".format(len(outputs)))
        it, t = get_last_it_gw(sim)
        print("\t\tGW: it:{} time:{}".format(it, t))

    else:
        Printcolor.blue("\t{} in postprocessed".format(sim))
        print("\t\t{} outputs".format(len(outputs)))
        it, t = get_last_it_gw(sim)
        print("\t\tGW: it:{} time:{}".format(it, t))
        it_, t_ = get_last_it_ppr(sim)
        print("\t\tPPR: it:{} time:{}".format(it_, t_))
        if it != it_:
            Printcolor.red("\t\tit({}) != ppr it({}): Postrocessing required".format(it, it_))

        if not os.path.isfile(Paths.ppr_sims+sim+"/parfile.par"):
            Printcolor.red("\tparfile not found")

        outflowed = os.path.isdir(MakePath.outflow(sim, "_0"))
        if outflowed:
            Printcolor.blue("\t\t\tGEO outflow is done")
        else:
            Printcolor.red("\t\t\tGEO outflow is NOT done")

        outflowed_b_w = os.path.isdir(MakePath.outflow(sim, "_0_b_w"))
        if outflowed_b_w:
            Printcolor.blue("\t\t\tBernoulli Wind outflow is done")
        else:
            Printcolor.red("\t\t\tBernoulli Wind outflow is NOT done")

for sim in gws:
    print("\n------------------------------------------")
    print("{}".format(sim))
    outputs, tars, dattars = get_outputs(sim)

    print("\toutput({}) .tar({}) .dat.tar({})".format(len(outputs), len(tars), len(dattars)))

    if len(outputs) > 0:
        Printcolor.blue("\t{} has {} outputs".format(sim, len(outputs)))

        try:
            analyze_sim(sim)
        except ValueError:
            Printcolor.red("Error. {} files for it_time not found".format(sim))

if len(gws) < len(pprs):
    Printcolor.red("NOT in gw170717:")
    for sim in pprs:
        if not sim in gws:
            Printcolor.red("\t{}".format(sim))


''' --- '''
print('\n')
''' --- '''

if len(gws) > len(pprs):
    Printcolor.yellow("NOT in postprocessed:")
    for sim in gws:
        if not sim in pprs:
            print("\t\t{}".format(sim)),
            _outputs = os.listdir(Paths.gw170817 + sim + '/')
            outputs = []
            for output in _outputs:
                if output.__contains__("output-"):
                    outputs.append(output)
            print(" {} outputs, ".format(len(outputs)))