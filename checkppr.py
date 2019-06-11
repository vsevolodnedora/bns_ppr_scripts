
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

    return it_time[:,0].max(), it_time[:,1].max() *  0.004925794970773136 * 1e-3

def analyze_sim(sim, time_limit=0.030):

    if not sim in pprs:
        Printcolor.red(", NOT in ppr"),

        it, t = get_last_it_gw(sim)

        if t < time_limit:
            Printcolor.red(", tend {:.3f}s".format(t), ','),
        else:
            Printcolor.green(", tend {:.3f}s".format(t), ','),
        #
        # print(", GW [it:{:d} time:{:.3f}]".format(int(it), t)),

    else:
        Printcolor.blue(", in ppr"),

        it, t = get_last_it_gw(sim)
        # print(", GW [it:{:d} time:{:.3f}]".format(int(it), t)),
        it_, t_ = get_last_it_ppr(sim)
        # print(", PPR: [it:{:d} time:{:.3f}]".format(int(it_), t_)),
        if int(it) != int(it_):
            Printcolor.red(", ppr REQUIRED [t_gw:{:.3f} t_ppr:{:3f}]".format(t, t_), ','),
        else:
            Printcolor.green(", ppr is complete", ',')

        if t < time_limit:
            Printcolor.red(", tend {:.3f}s".format(t), ','),
        else:
            Printcolor.green(", tend {:.3f}s".format(t), ','),

        if not os.path.isfile(Paths.ppr_sims+sim+"/parfile.par"):
            Printcolor.red(", parfile not found |"),

        outflowed = os.path.isdir(MakePath.outflow(sim, "_0"))
        if outflowed:
            Printcolor.green(" _0 done", ',')
        else:
            Printcolor.red(" _0 NOT done", ',')

        if not sim in Lists.dyn_not_pas:
            Printcolor.green(" Dyn pass", ',')
        else:
            Printcolor.red(" Dyn NOT pass", ',')

        if sim in Lists.bern_pass:
            Printcolor.green(" Bern pass", ',')
        else:
            Printcolor.yellow(" Bern NOT pass", ',')


        outflowed_b_w = os.path.isdir(MakePath.outflow(sim, "_0_b_w"))
        if outflowed_b_w:
            Printcolor.green(", _0_b_w done", ',')
        else:
            Printcolor.red(", _0_b_w NOT done", ',')

for sim in gws:
    # print("\n------------------------------------------")
    print("{}\t ".format(sim)),
    outputs, tars, dattars = get_outputs(sim)
    # print("\toutput({}) .tar({}) .dat.tar({})".format(len(outputs), len(tars), len(dattars)))

    if len(outputs) > 0:
        Printcolor.blue("\t{} outputs".format(len(outputs))),

        # analyze_sim(sim)

        try:
            analyze_sim(sim)
        except ValueError:
            Printcolor.red("Error. FAILD with ValueError" ',')
    else:
        Printcolor.red("\tNO outputs")

    print(" |")

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