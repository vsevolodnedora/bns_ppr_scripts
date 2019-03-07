# bring together all methods like 'model.disk, model.ejecta, model.gw, ... "
from sys import path
path.append('modules/')

from general import Printcolor


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