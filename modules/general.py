# contains methods that are very general
import numpy as np
from lists import *
from units import time_constant

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
    def red(text):
        print(Printcolor.FAIL + text + Printcolor.ENDC)

    @staticmethod
    def yellow(text):
        print(Printcolor.WARNING + text + Printcolor.ENDC)

    @staticmethod
    def blue(text):
        print(Printcolor.OKBLUE + text + Printcolor.ENDC)

    @staticmethod
    def green(text):
        print(Printcolor.OKGREEN + text + Printcolor.ENDC)

    @staticmethod
    def bold(text):
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
    def get_dens_decomp_3d(dens_3d, phi_3d, dphi_3d, dr_3d, dz_3d, m=1):
        '''
        Integrates density over 'z'
        Returns complex arrays [\int(d\phi)] and [\int(dr d\phi)]
        '''

        integ_dens_z = np.sum(dens_3d * dz_3d[:, :, 0], axis=2) # -> 2d array

        integ_over_z_phi = np.sum(integ_dens_z * np.exp(1j * m * phi_3d[:, :, 0]) * dphi_3d[:, :, 0], axis=1) # -> 1d array

        integ_over_z_phi_r = np.sum(integ_over_z_phi * dr_3d[:, 0, 0]) # -> number

        return integ_over_z_phi, integ_over_z_phi_r

def find_nearest_index(array, value):
        ''' Finds index of the value in the array that is the closest to the provided one '''
        idx = (np.abs(array - value)).argmin()
        return idx

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


def interpoate_time_form_it(it_list, file_it_use, time_units='s'):
    '''
    From list of iterations, returns a list of timesteps (in seconds)
    '''
    from scipy import interpolate
    # start_t = time.time()
    print("\t Interpolating time for interations...")

    path_to_template = file_it_use # "/collated/outflow_det_0.asc"

    it_time = np.loadtxt(path_to_template, unpack=False, usecols=[0, 1])  # it t

    if np.array(it_list).min() < it_time[:, 0].min():
        raise ValueError("list it:{} < min it:{} in file: {}"
                         .format(np.array(it_list).min(), it_time[:, 0].min(), file_it_use))
    if np.array(it_list).max() > it_time[:, 0].max():
        raise ValueError("list it:{} > max it:{} in file: {}"
                         .format(np.array(it_list).max(), it_time[:, 0].max(), file_it_use))

    f = interpolate.interp1d(it_time[:, 0], it_time[:, 1], kind='linear')
    time_list = f(it_list)
    if time_units == 's': time_list *= (time_constant * 1e-3)
    elif time_units == 'ms': time_list *= (time_constant)
    elif time_units == 'msol': time_list *= 1.
    else: raise NameError("Time units:{} are not recognized. Useon of [s, ms, msol]")

    # print(" done! (%.2f sec)" % (time.time() - start_t))

    return time_list

def interpoate_it_form_time(time_list, file_it_use, time_units='s'):
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

    path_to_template = file_it_use # "/collated/outflow_det_0.asc"

    it_time = np.loadtxt(path_to_template, unpack=False, usecols=[0, 1])  # it t

    if np.array(time_list).min() < it_time[:, 1].min():
        raise ValueError("list [sol.mass] time:{} < min time:{} in file: {}"
                         .format(np.array(time_list).min(), it_time[:, 1].min(), file_it_use))
    if np.array(time_list).max() > it_time[:, 1].max():
        raise ValueError("list [sol.mass] time:{} > max time:{} in file: {}"
                         .format(np.array(time_list).max(), it_time[:, 1].max(), file_it_use))

    f = interpolate.interp1d(it_time[:, 1], it_time[:, 0], kind='linear')
    it_list = f(time_list)

    # print(" done! (%.2f sec)" % (time.time() - start_t))

    return np.array(it_list, dtype=int)

