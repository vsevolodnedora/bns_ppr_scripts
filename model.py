# bring together all methods like 'model.disk, model.ejecta, model.gw, ... "
from __future__ import division
from sys import path
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

path.append('modules/')
from general import *
from lists import *
from filework import *
from units import time_constant, volume_constant, energy_constant


# print simulations


marker_list = ["o", "s", "^", "d"]

def eos_color(eos):
    if eos == 'DD2':
        return 'blue'
    elif eos == 'BHBlp':
        return 'purple'
    elif eos == 'LS220':
        return 'orange'
    elif eos == 'SFHo':
        return 'red'
    elif eos == 'SLy4':
        return 'green'
    else:
        Printcolor.yellow("Unknown eos:{} [color->black]".format(eos))
        return 'black'

def eos_colors(eoss):
   colors = []
   for eos in eoss:
       colors.append(eos_color(eos))
   return colors

if __name__ == '__main__':
    print(simulations.index.values)