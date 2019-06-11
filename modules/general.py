# contains methods that are very general
import numpy as np
from lists import *
from units import time_constant
from glob import glob
import os
import re
from math import log

class Printcolor:

    HEADER  = '\033[95m'
    OKBLUE  = '\033[94m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL    = '\033[91m'
    ENDC    = '\033[0m'
    BOLD    = '\033[1m'
    UNDERLINE='\033[4m'

    def __init__(self):
        pass

    @staticmethod
    def red(text, comma=','):
        if comma == ',':
            print(Printcolor.FAIL + text + Printcolor.ENDC),
        else:
            print(Printcolor.FAIL + text + Printcolor.ENDC)

    @staticmethod
    def yellow(text, comma=','):
        if comma == ',':
            print(Printcolor.WARNING + text + Printcolor.ENDC),
        else:
            print(Printcolor.WARNING + text + Printcolor.ENDC)

    @staticmethod
    def blue(text, comma=','):
        if comma == ',':
            print(Printcolor.OKBLUE + text + Printcolor.ENDC),
        else:
            print(Printcolor.OKBLUE + text + Printcolor.ENDC)

    @staticmethod
    def green(text, comma=','):
        if comma == ',':
            print(Printcolor.OKGREEN + text + Printcolor.ENDC),
        else:
            print(Printcolor.OKGREEN + text + Printcolor.ENDC)
            

    @staticmethod
    def bold(text, comma=','):
        if comma == ',':
            print(Printcolor.BOLD + text + Printcolor.ENDC),
        else:
            print(Printcolor.BOLD + text + Printcolor.ENDC)

class LABELS:
    def __init__(self):
        pass

    @staticmethod
    def labels(v_n):
        # solar

        if v_n == 'flux':
            return r'$\dot{M}$'

        if v_n == 'time':
            return r'$t$ [ms]'

        if v_n == 't-tmerg':
            return r'$t-t_{merg}$ [ms]'

        if v_n == 'vel_inf' or v_n == 'vel_inf_bern':
            return r'$\upsilon$ $[c]$'

        if v_n == 'l':
            return r'$\log(L/L_{\odot})$'  # (L_{\odot})

        if v_n == 'r':
            return r'$R(R_{\odot})$'

        if v_n == 'theta':
            return r'$\theta$'

        if v_n == 'm' or v_n == 'xm':
            return r'$M(M_{\odot})$'

        if v_n == 'mass':
            return r'$M/M_{ej}$'

        if v_n == 'v' or v_n == 'u' or v_n == 'vel':
            return 'v (km/s)'

        if v_n == 'rho':
            return r'$\log(\rho)$'

        if v_n == 'k' or v_n == 'kappa':
            return r'$\kappa$'

        if v_n == 't':
            return r'log(T/K)'

        if v_n == 'lm':
            return r'$\log(L/M)$'

        if v_n == 'mdot':
            return r'$\log(\dot{M}$)'

        if v_n == 'Ye' or v_n == 'ye':
            return r'$Ye$'

        if v_n == 'entropy' or v_n == 's':
            return r'$s$'

        if v_n == 't_eff' or v_n == 'T_eff':
            return r'$\log($T$_{eff}/K)$'

        if v_n == 't_*' or v_n == 'T_*':
            return r'$\log($T$_{*}$/K$)$'

        if v_n == 'r_eff' or v_n == 'R_eff':
            return r'$\log($R$_{eff}/R_{\odot})$'

        if v_n == 'rho':
            return r'$\log(\rho)$'

        if v_n == 'tau':
            return r'$\tau$'

        if v_n == 'mfp':
            return r'$\log(\lambda)$'

class FromName:
    def __init__(self):
        pass

    @staticmethod
    def get_eos(self):
        eos = self.sim.split('_')[0]
        if eos in Lists.eos: return eos
        else:raise NameError("Unrecognised EOS: {} is sim: {}".format(eos, self.sim))

    @staticmethod
    def get_m1(self):
        m_str = self.sim.split('_')[1]
        if m_str[0] != 'M': raise NameError("Unrecognized masses naming (M is missing) in sim: {}".format(self.sim))
        if len(m_str) != 9: raise NameError("Unrecognized masses naming (length != 9) in sim: {}".format(self.sim))
        m1 = float(self.sim.split('_')[1][1:5]) / 1000
        return m1

    @staticmethod
    def get_m2(self):
        m_str = self.sim.split('_')[1]
        if m_str[0] != 'M': raise NameError("Unrecognized masses naming (M is missing) in sim: {}".format(self.sim))
        if len(m_str) != 9: raise NameError("Unrecognized masses naming (length != 9) in sim: {}".format(self.sim))
        m2 = float(self.sim.split('_')[1][5:]) / 1000
        return m2

    @staticmethod
    def get_q(self):
        m_str = self.sim.split('_')[1]
        if m_str[0] != 'M': raise NameError("Unrecognized masses naming (M is missing) in sim: {}".format(self.sim))
        if len(m_str) != 9: raise NameError("Unrecognized masses naming (length != 9) in sim: {}".format(self.sim))
        q = self.get_m1() / self.get_m2()
        return q

    @staticmethod
    def get_res(self):

        for res in Lists.res:
            if res in self.sim.split("_"):
                return res
        else: raise NameError("Unrecognized resolution in sim {}".format(self.sim))

class PHYSICS:

    def __init__(self):
        pass

    @staticmethod
    def get_dens_decomp_2d(dens_2d, phi_2d, dphi_2d, dr_2d, m=1):
        '''
        Uses a 2d slice at z=0
        Returns complex arrays [\int(d\phi)] and [\int(dr d\phi)]
        '''
        # dens_2d = dens_3d[:, :, 0]
        # phi_2d  = phi_3d[:, :, 0]
        # dr_2d   = dr_3d[:, :, 0]
        # dphi_2d = dphi_3d[:, :, 0]

        integ_over_phi = np.sum(dens_2d * np.exp(1j * m * phi_2d) * dphi_2d, axis=1)

        integ_over_phi_r = np.sum(integ_over_phi * dr_2d[:, 0])

        return integ_over_phi, integ_over_phi_r

    @staticmethod
    def get_dens_decomp_3d(dens_3d, r, phi_3d, dphi_3d, dr_3d, dz_3d, m=1):
        '''
        Integrates density over 'z'
        Returns complex arrays [\int(d\phi)] and [\int(dr d\phi)]
        '''

        integ_dens_z = np.sum(dens_3d * dz_3d[:, :, :], axis=2) # -> 2d array

        integ_over_z_phi = np.sum(integ_dens_z * np.exp(1j * m * phi_3d[:, :, 0]) * dphi_3d[:, :, 0], axis=1) # -> 1d array

        integ_over_z_phi_r = np.sum(integ_over_z_phi * dr_3d[:, 0, 0] * r[:, 0, 0]) # -> number

        return integ_over_z_phi, integ_over_z_phi_r

    @staticmethod
    def get_retarded_time(t, M_Inf=2.7, R_GW=400.0):
        R = R_GW * (1 + M_Inf / (2 * R_GW)) ** 2
        rstar = R + 2 * M_Inf * log(R / (2 * M_Inf) - 1)
        return t - rstar

def find_nearest_index(array, value):
        ''' Finds index of the value in the array that is the closest to the provided one '''
        idx = (np.abs(array - value)).argmin()
        return idx

def combine(x, y, xy, corner_val=None):
    '''creates a 2d array  1st raw    [0, 1:]-- x -- density      (log)
                           1st column [1:, 0] -- y -- lemperature (log)
                           Matrix     [1:,1:] -- xy --Opacity     (log)
       0th element in 1st raw (column) - can be used a corner value

    '''
    x = np.array(x)
    y = np.array(y)
    xy = np.array((xy))

    if len(x) != len(y):
        print('\t__Warning. x({}) != y({}) (combine)'.format(len(x), len(y)))
    if len(x) != len(xy[0, :]):
        raise ValueError('\t__Warning. x({}) != xy[0, :]({}) (combine)'.format(len(x), len(xy[0, :])))
    if len(y) != len(xy[:, 0]):
        raise ValueError('\t__Warning. y({}) != xy[:, 0]({}) (combine)'.format(len(y), len(xy[:, 0])))

    res = np.insert(xy, 0, x, axis=0)
    new_y = np.insert(y, 0, 0, axis=0)  # inserting a 0 to a first column of a
    res = np.insert(res, 0, new_y, axis=1)

    if corner_val != None:
        res[0, 0] = corner_val

    return res

def get_uniqe_eoss(eoss):

    uinique_eos = []
    uinique_eos.append(eoss[0])
    for eos in eoss:
        if not eos in uinique_eos:
            uinique_eos.append(eos)

    return uinique_eos

def color_for_mkn_band(band):

    if band == 'K' or band == 'Ks':
        return 'darkred'
    elif band == 'H':
        return 'red'
    elif band == 'J':
        return 'darkorangered'
    elif band == 'y':
        return 'orange'
    elif band == 'z':
        return 'gold'
    elif band == 'i':
        return 'greenyellow'
    elif band == 'r':
        return 'lawngreen'
    elif band == 'V':
        return 'lightskyblue'
    elif band == 'g':
        return 'blue'
    else:
        Printcolor.yellow("Warning: no color is set for band:{} returning black".format(band))
        return 'black'

def clrs(arr):
    colors = []
    for i in range(len(arr)):
        colors.append('C'+str(i))
    return colors

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

def get_color_for_q(q):
    return Lists.colors_q[round(q, 3)]

def set_it_output_map(path_to_sim, file_for_it="dens.norm1.asc", clean=True):
    """
    Loads set of files that have '1:it 2:time ...' structure to get a map
    of what output-xxxx contains what iteration (and time)
    """

    gen_set = {"indir": path_to_sim,
               "file_for_it": file_for_it}

    output_it_map = {}

    fpath = gen_set["indir"] + "output-*" + "/data/" + gen_set['file_for_it']
    files = glob(fpath)

    if not clean:
        print('-' * 25 + 'LOADING it list ({})'
              .format(gen_set['file_for_it']) + '-' * 25)
        print("\t loading from: {}".format(gen_set['indir']))

    if len(files) == 0:
        raise ValueError("No files found for it_time mapping searched:\n{}"
                         .format(fpath))
    if not clean:
        print("\t   files found: {}".format(len(files)))



    it_time = np.zeros(2)
    for file in files:
        o_name = file.split('/')
        o_dir = ''
        for o_part in o_name:
            if o_part.__contains__('output-'):
                o_dir = o_part
        if o_dir == '':
            raise NameError("Did not find output-xxxx in {}".format(o_name))
        it_time_i = np.loadtxt(file, usecols=(0, 1))
        output_it_map[o_dir] = it_time_i
        it_time = np.vstack((it_time, it_time_i))
    it_time = np.delete(it_time, 0, 0)
    if not clean:
        print('outputs:{} iterations:{} [{}->{}]'.format(len(files),
                                                         len(it_time[:, 0]),
                                                         int(it_time[:, 0].min()),
                                                         int(it_time[:, 0].max())))
    if len(it_time[:, 0]) != len(set(it_time[:, 0])):
        if not clean:
            Printcolor.yellow("Warning: repetitions found in the "
                              "loaded iterations")
        iterations = np.unique(it_time[:, 0])
        timestpes = np.unique(it_time[:, 1])
        if not len(iterations) == len(timestpes):
            raise ValueError("Failed attmept to remove repetitions from "
                             "\t it and time lists. Wrong lengths: {} {}"
                             .format(len(iterations), len(timestpes)))
    else:
        if not clean:
            Printcolor.blue("No repetitions found in loaded it list")
        iterations = np.unique(it_time[:, 0])
        timestpes = np.unique(it_time[:, 1])

    if not clean:
        print('-' * 30 + '------DONE-----' + '-' * 30)

    return output_it_map, np.vstack((iterations, timestpes)).T

def interpoate_time_form_it(it_list, path_to_sim, time_units='s', extrapolate=True):
    '''
    From list of iterations, returns a list of timesteps (in seconds)
    '''
    from scipy import interpolate
    # start_t = time.time()
    print("\t Interpolating time for interations...")

    # path_to_template = file_it_use # "/collated/outflow_det_0.asc"

    _, it_time = set_it_output_map(path_to_sim)# np.loadtxt(path_to_template, unpack=False, usecols=[0, 1])  # it t

    if not any(it_list):
        raise ValueError("Passes empty iteration list")
    if len(it_time[:,0]) == 0:
        raise ValueError("Failed to get it_time map")

    if np.array(it_list).min() < it_time[:, 0].min():
        raise ValueError("list it:{} < min it:{} in sim: {}"
                         .format(np.array(it_list).min(), it_time[:, 0].min(), path_to_sim))
    if np.array(it_list).max() > it_time[:, 0].max():
        if extrapolate:
            Printcolor.yellow("list it:{} > max it:{} in sim: {}"
                              .format(np.array(it_list).max(), it_time[:, 0].max(), path_to_sim))
        else:
            raise ValueError("list it:{} > max it:{} in sim: {}"
                             .format(np.array(it_list).max(), it_time[:, 0].max(), path_to_sim))

    if extrapolate:
        f = interpolate.interp1d(it_time[:, 0], it_time[:, 1], kind='linear', bounds_error=False,
                                 fill_value="extrapolate")
    else:
        f = interpolate.interp1d(it_time[:, 0], it_time[:, 1], kind='linear')
    time_list = f(it_list)
    if time_units == 's': time_list *= (time_constant * 1e-3)
    elif time_units == 'ms': time_list *= (time_constant)
    elif time_units == 'msol': time_list *= 1.
    else: raise NameError("Time units:{} are not recognized. Useon of [s, ms, msol]")

    # print(" done! (%.2f sec)" % (time.time() - start_t))

    return time_list

def interpoate_it_form_time(time_list, path, time_units='s'):
    '''
    From list of iterations, returns a list of timesteps (in seconds)
    '''
    from scipy import interpolate
    # start_t = time.time()
    print("\t Interpolating iterations for given times...")

    time_list = np.array(time_list, dtype=float)
    if time_units == 's': time_list /= (time_constant * 1e-3)
    elif time_units == 'ms': time_list /= (time_constant)
    elif time_units == 'msol': time_list /= 1.
    else: raise NameError("Time units:{} are not recognized. Useon of [s, ms, msol]")

    # path_to_template = file_it_use # "/collated/outflow_det_0.asc"

    _, it_time = set_it_output_map(path) # np.loadtxt(path_to_template, unpack=False, usecols=[0, 1])  # it t

    if np.array(time_list).min() < it_time[:, 1].min():
        raise ValueError("list [sol.mass] time:{} < min time:{} in file: {}"
                         .format(np.array(time_list).min(), it_time[:, 1].min(), path))
    if np.array(time_list).max() > it_time[:, 1].max():
        raise ValueError("list [sol.mass] time:{} > max time:{} in file: {}"
                         .format(np.array(time_list).max(), it_time[:, 1].max(), path))

    f = interpolate.interp1d(it_time[:, 1], it_time[:, 0], kind='linear')
    it_list = f(time_list)

    # print(" done! (%.2f sec)" % (time.time() - start_t))

    return np.array(it_list, dtype=int)

def get_list_iterationsfrom_res_3d(sim):
    """
    Checks the /res_3d/ for 12345 folders, (iterations) retunrs their sorted list
    :param sim:
    :return:
    """

    if not os.path.isdir(Paths.ppr_sims + sim + '/res_3d/'):
        raise IOError("no {} directory found".format(Paths.ppr_sims + sim + '/res_3d/'))
    itdirs = os.listdir(Paths.ppr_sims + sim + '/res_3d/')
    if len(itdirs) == 0:
        raise NameError("No iteration-folders found in the {}".format(Paths.ppr_sims + sim + '/res_3d/'))


    # this is a f*cking masterpiece of programming)))
    list_iterations = np.array(np.sort(np.array(list([int(itdir) for itdir in itdirs if re.match("^[-+]?[0-9]+$", itdir)]))))
    if len(list_iterations) == 0:
        raise ValueError("Error extracting the iterations")

    return list(list_iterations)

def get_list_iterationsfrom_profiles_3d(sim, in_sim_dir="/profiles/3d/"):
    """
    Checks the /res_3d/ for 12345 folders, (iterations) retunrs their sorted list
    :param sim:
    :return:
    """

    if not os.path.isdir(Paths.gw170817 + sim + in_sim_dir):
        raise IOError("no {} directory found".format(Paths.gw170817 + sim + in_sim_dir))
    profiles = glob(Paths.gw170817 + sim + in_sim_dir + "*.h5")
    if len(profiles) == 0:
        raise NameError("No iteration-folders found in the {}".format(Paths.ppr_sims + sim + in_sim_dir))

    # print(profiles)

    # this is a f*cking masterpiece of programming)))
    list_iterations = np.array(np.sort(np.array(list([int(profile.split('/')[-1].split(".h5")[0])
                                                      for profile in profiles
                                                      if re.match("^[-+]?[0-9]+$",
                                                                  profile.split('/')[-1].split(".h5")[0])]))))
    if len(list_iterations) == 0:
        raise ValueError("Error extracting the iterations")

    return list(list_iterations)

def find_itdir_with_grid(sim, gridfname='cyl_grid.h5'):

    # for dir_ in os.listdir()




    path = Paths.ppr_sims + sim + '/res_3d/' + '*' + '/' + gridfname
    files = glob(path)

    # print(files)

    if len(files) == 0:
        raise ValueError("No grid file ({}) found for {}".format(gridfname, sim))

    if len(files) > 1:
        Printcolor.yellow("More than 1({}) grid file ({}) found for {}"
                          .format(len(files), gridfname, sim))
        return files[0]

    return str(files)


def get_it_from_itdir(itdir):
    return int(itdir.split('/')[-2])


def get_list_ppr_simulations():

    sims = os.listdir(Paths.ppr_sims)
    return sims




def x_y_z_sort(x_arr, y_arr, z_arr=np.empty(0, ), sort_by_012=0):
    '''
    RETURNS x_arr, y_arr, (z_arr) sorted as a matrix by a row, given 'sort_by_012'
    :param x_arr:
    :param y_arr:
    :param z_arr:
    :param sort_by_012:
    :return:
    '''

    if len(z_arr) == 0 and sort_by_012 < 2:
        if len(x_arr) != len(y_arr):
            raise ValueError('len(x)[{}]!= len(y)[{}]'.format(len(x_arr), len(y_arr)))

        x_y_arr = []
        for i in range(len(x_arr)):
            x_y_arr = np.append(x_y_arr, [x_arr[i], y_arr[i]])

        x_y_sort = np.sort(x_y_arr.view('float64, float64'), order=['f{}'.format(sort_by_012)], axis=0).view(
            np.float)
        x_y_arr_shaped = np.reshape(x_y_sort, (int(len(x_y_sort) / 2), 2))
        return x_y_arr_shaped[:, 0], x_y_arr_shaped[:, 1]

    if len(z_arr) > 0 and len(z_arr) == len(y_arr):
        if len(x_arr) != len(y_arr) or len(x_arr) != len(z_arr):
            raise ValueError('len(x)[{}]!= len(y)[{}]!=len(z_arr)[{}]'.format(len(x_arr), len(y_arr), len(z_arr)))

        x_y_z_arr = []
        for i in range(len(x_arr)):
            x_y_z_arr = np.append(x_y_z_arr, [x_arr[i], y_arr[i], z_arr[i]])

        x_y_z_sort = np.sort(x_y_z_arr.view('float64, float64, float64'), order=['f{}'.format(sort_by_012)],
                             axis=0).view(np.float)
        x_y_z_arr_shaped = np.reshape(x_y_z_sort, (int(len(x_y_z_sort) / 3), 3))
        return x_y_z_arr_shaped[:, 0], x_y_z_arr_shaped[:, 1], x_y_z_arr_shaped[:, 2]