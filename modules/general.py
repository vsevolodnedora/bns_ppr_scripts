# contains methods that are very general

from lists import *

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