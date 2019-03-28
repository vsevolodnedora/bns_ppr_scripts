from __future__ import division
from sys import path
path.append('modules/')

from scipy import interpolate
import numpy as np
import csv
import os

from utils import get_retarded_time
import units as ut
from general import *
from lists import *

# table is a global here
table = []
header = "global"

# load save .csv file
def load_table(file_name):

    global header
    with open(file_name, "r") as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            table.append(row)
            # print(row)
        header= reader.fieldnames
def save_table(file_name):

    if os.path.exists(file_name):
        os.remove(file_name)

    with open(file_name, "w") as csvfile:
        writer = csv.DictWriter(csvfile, header)
        writer.writeheader()
        for run in table:
            writer.writerow(run)

# from the name of the simulation, fill the mass, eos,
def fill_self_parameters():
    for run in table:
        eos = run["name"].split('_')[0]
        m1m2 = run["name"].split('_')[1]
        if m1m2[0] != 'M':
            Printcolor.yellow(
                "Warning. m1m2 is not [1] component of name. Using [2] (run:{})".format(run["name"]))
            # print("Warning. m1m2 is not [1] component of name. Using [2] (run:{})".format(run["name"]))
            m1m2 = run["name"].split('_')[2]
        else:
            m1m2 = ''.join(m1m2[1:])
        m1 = float(''.join(m1m2[:4])) / 1000
        m2 = float(''.join(m1m2[4:])) / 1000
        print(run["name"])

        run["M1"] = m1
        run["M2"] = m2
        run["EOS"] = eos

def fill_self_dyn_phase_status():

    Printcolor.blue("...adding dynamical phase info from static list")

    for run in table:
        if run["name"] in Lists.dyn_not_pas:
            run["dyn_phase"] = "not passed"
        else:
            run["dyn_phase"] = "passed"

# Parfile

def get_par_form_parfile():
    dics_lorene = []

    for run in table:
        basepath = Paths.ppr_sims + run["name"] + '/parfile.par'

        dic_lorene = {}
        dic_lorene["name"] = run["name"]
        try:
            lorene_run_dir = ''
            with open(basepath, "r") as parfile:
                for line in parfile:
                    if line.__contains__('LoreneID::lorene_bns_file'):
                        lorene_run_dir = line.split()[-1]
                if lorene_run_dir == '':
                    raise IOError('Line with *LoreneID::lorene_bns_file* in not found in parfile: {}.par'
                                  .format(run["name"]))
            dic_lorene["lorene_file"] = lorene_run_dir
        except:
            Printcolor.yellow('Warning Parfile {}.par is not found'.format(run["name"]))
            # print('Warning Parfile {}.par is not found'.format(run["name"]))
            dic_lorene["lorene_file"] = ''

        # fix for those runs that have mess in the lorene name
        if run["name"] in ["DD2_T05_15491205_45km_M0corot", "DD2_T05_15491205_45km_M0LKcorot"]:
            dic_lorene["lorene_file"] = 'x/R04_corot/x/'

        dics_lorene.append(dic_lorene)
    return dics_lorene

# Lorene

def get_lorene_data(run, lorene_file='properties_R04.txt'):
    '''
    Searches for the corresponding run inside the Lorene summary file, based on *EOS* and *masses* (in the name)
    and returns the following variables if the entry is found:
    omega, orb_sep, m_adm, j_adm, mb1, mb2
    '''

    spec_list = ["DD2_T05_15491205_45km_M0corot", "DD2_T05_15491205_45km_M0LKcorot"]

    # extract run's mass and eos for matching with Lorene entry
    run_parts = run.split('_')
    run_eos = run_parts[0]
    if run in spec_list:
        Printcolor.yellow('Run {} in spec. list. Possible error with matching with Lorene Data.'.format(run))
        # print('Run {} in spec. list. Possible error with matching with Lorene Data.'.format(run))
        # return None, None, None, None, None, None
        run_masses = run_parts[2]
    else:
        run_masses = ''.join(run_parts[1][1:])

    with open(Paths.lorene + lorene_file, 'r') as lorene:
        next(lorene)
        for line in lorene:
            line_parts = line.split()
            if len(line_parts) != 7:
                raise ValueError('Lorene table line contains <> 7 entries: {}, {}'.format(line_parts, len(line_parts)))

            name = line_parts[0]
            if len(name.split('_')) != 4: raise NameError('Run: {} in Lorene has a too long name (<>4)'.format(name))
            eos = name.split('_')[0]
            temp = name.split('_')[1]
            masses = name.split('_')[2]

            separ = name.split('_')[3]

            if run_eos == eos and run_masses == masses:
                omega = line.split()[1]  # Omega [rad/s]
                orb_sep = line.split()[2]  # km
                m_adm = str(2 * float(line.split()[3]))  # MADM [Msun]
                j_adm = line.split()[4]  # JADM [GMsun^2/c]
                mb1 = line.split()[5]  # MbA [Msun]
                mb2 = line.split()[6]  # MbB [Msun]

                return omega, orb_sep, m_adm, j_adm, mb1, mb2

            # else:
        Printcolor.yellow("Warning: Lorene entry is not found for "
                                 "eos:{} m:{} (lorene file:{})"
                                 .format(run_eos, run_masses, lorene_file))
        return None, None, None, None, None, None

    # if run == 'DD2_M14971246_M0_LR': print("Error: {}".format("DD2_M14971246_M0_LR"))
    # exit(1)
    # lorene.close()

    return None, None, None, None, None, None

def fill_table_value(run_dic, v_n, value, replace=False):
    '''
    Simply takes the dictionary of the run and places the variable value.
    May be forced not to replace already existing value
    '''

    if v_n not in run_dic.keys(): raise NameError('Variable {} not in run dictionary {}'.format(v_n, run_dic.keys()))

    old_value = run_dic[v_n]
    if old_value == '':
        run_dic[v_n] = value
    if old_value != '' and replace:
        run_dic[v_n] = value

def fill_from_lorene_data(lorene_file='properties_R02.txt', replace_value=False):
    '''
    Reads the *lorene_file*,
    for every run in GLOBAL table finds an entry in this file based on the name eos_masses,
    that has variables:
    [omega, orb_sep, m_adm, j_adm, mb1, mb2]
    adds some of them:
    m_adm, j_adm, mb1, mb2
    into the GLOBAL table
    '''

    for run in table:
        omega, orb_sep, m_adm, j_adm, mb1, mb2 = get_lorene_data(run['name'], lorene_file)
        if omega != None:
            mb = float(mb1) + float(mb2)

            fill_table_value(run, 'Mb1', mb1, replace_value)
            fill_table_value(run, 'Mb2', mb2, replace_value)
            fill_table_value(run, 'Mb', mb, replace_value)
            fill_table_value(run, 'MADM', m_adm, replace_value)
            fill_table_value(run, 'JADM', j_adm, replace_value)

            fill_table_value(run, 'comment', '{}'.format(lorene_file.split('_')[-1].split('.')[0]), replace_value)

def fill_from_lorene(dics_lorene, replace_value=False):
    for run, dic_lorene in zip(table, dics_lorene):

        if dic_lorene["lorene_file"] != '':
            lorene_run_fold = dic_lorene["lorene_file"].split('/')[-3]
            lorene_run_summary_file = 'properties_{}.txt'.format(lorene_run_fold)
            print(lorene_run_summary_file)

            omega, orb_sep, m_adm, j_adm, mb1, mb2 = get_lorene_data(run['name'], lorene_run_summary_file)
            if omega != None:
                mb = float(mb1) + float(mb2)

                fill_table_value(run, 'f0', float(omega) / (2 * np.pi), replace_value)
                fill_table_value(run, 'Mb1', mb1, replace_value)
                fill_table_value(run, 'Mb2', mb2, replace_value)
                fill_table_value(run, 'Mb', mb, replace_value)
                fill_table_value(run, 'MADM', m_adm, replace_value)
                fill_table_value(run, 'JADM', j_adm, replace_value)

                fill_table_value(run, 'comment', '{}'.format(lorene_run_fold), replace_value)

# TOVs

def fill_from_TOVs():
    sfho_tov = np.loadtxt(Paths.TOVs + 'SFHo_love.dat')
    dd2_tov = np.loadtxt(Paths.TOVs + 'DD2_love.dat')
    ls220_tov = np.loadtxt(Paths.TOVs + 'LS220_love.dat')
    sly4_tov = np.loadtxt(Paths.TOVs + 'SLy4_love.dat')
    bhblp_tov = np.loadtxt(Paths.TOVs + 'BHBlp_love.dat')

    for run in table:

        if run["EOS"] == 'SFHo':
            tov_table = sfho_tov
        elif run["EOS"] == 'DD2':
            tov_table = dd2_tov
        elif run["EOS"] == 'LS220':
            tov_table = ls220_tov
        elif run["EOS"] == 'SLy4':
            tov_table = sly4_tov
        elif run["EOS"] == 'BHBlp' or run["EOS"] == 'BHB':
            tov_table = bhblp_tov
        else:
            raise NameError("TOV sequences are not found for EOS:{} ".format(run["EOS"]))

        m_grav = tov_table[:, 1]
        m_bary = tov_table[:, 2]
        r = tov_table[:, 3]
        comp = tov_table[:, 4]  # compactness
        kl = tov_table[:, 5]
        lamb = tov_table[:, 6]  # lam

        interp_grav_bary = interpolate.interp1d(m_bary, m_grav, kind='cubic')
        interp_lamb_bary = interpolate.interp1d(m_bary, lamb, kind='cubic')
        interp_comp_bary = interpolate.interp1d(m_bary, comp, kind='cubic')
        interp_k_bary = interpolate.interp1d(m_bary, kl, kind='cubic')
        interp_r_bary = interpolate.interp1d(m_bary, r, kind='cubic')

        if run["Mb1"] != '':
            run["lam21"] = interp_lamb_bary(float(run["Mb1"]))  # lam21
            run["Mg1"] = interp_grav_bary(float(run["Mb1"]))
            run["C1"] = interp_comp_bary(float(run["Mb1"]))  # C1
            run["k21"] = interp_k_bary(float(run["Mb1"]))
            run["R1"] = interp_r_bary(float(run["Mb1"]))
            # run["R1"] = run["M1"] / run["C1"]

        if run["Mb2"] != '':
            run["lam22"] = interp_lamb_bary(float(run["Mb2"]))  # lam22
            run["Mg2"] = interp_grav_bary(float(run["Mb2"]))
            run["C2"] = interp_comp_bary(float(run["Mb2"]))  # C2
            run["k22"] = interp_k_bary(float(run["Mb2"]))
            run["R2"] = interp_r_bary(float(run["Mb2"]))
            # run["R2"] = run["M2"] / run["C2"]

        if run["Mg1"] != '' and run["Mg2"] != '':
            mg1 = float(run["Mg1"])
            mg2 = float(run["Mg2"])
            mg_tot = mg1 + mg2
            k21 = float(run["k21"])
            k22 = float(run["k22"])
            c1 = float(run["C1"])
            c2 = float(run["C2"])
            lam1 = float(run["lam21"])
            lam2 = float(run["lam22"])

            kappa21 = 2 * ((mg1 / mg_tot) ** 5) * (mg2 / mg1) * (k21 / (c1 ** 5))

            kappa22 = 2 * ((mg2 / mg_tot) ** 5) * (mg1 / mg2) * (k22 / (c2 ** 5))

            run["k2T"] = kappa21 + kappa22

            tmp1 = (mg1 + (12 * mg2)) * (mg1 ** 4) * lam1
            tmp2 = (mg2 + (12 * mg1)) * (mg2 ** 4) * lam2
            run["Lambda"] = (16. / 13.) * (tmp1 + tmp2) / (mg_tot ** 5.)

# GW

def fill_fpeak():
    for run in table:
        try:
            f, hr, hi = np.loadtxt(Paths.ppr_sims + run["name"] + "/waveforms/postmerger_psd_l2_m2.dat",
                                   usecols=(0, 1, 2), unpack=True)
            idx = f > 0.0
            f, h = f[idx], np.sqrt(f[idx]) * np.sqrt(hr[idx] ** 2 + hi[idx] ** 2)
            ipeak = np.argmax(np.abs(h))

            run['fpeak'] = f[ipeak]
        except:
            print('Failed *fpeak* extraction for run:{}'.format(run["name"]))
            run['fpeak'] = "nan"

def fill_tmerg():
    dic_tmergs = []
    for run in table:
        dic_tmerg = {}

        try:
            tmrgr = float(open(Paths.ppr_sims + run["name"]
                               + "/waveforms/tmerger.dat", "r").readline())

            tmrgr_r = get_retarded_time(tmrgr, M_Inf=float(run['M1'])
                                                     + float(run['M2']), R_GW=400)
            tmrgr = ut.conv_time(ut.cactus, ut.cgs, tmrgr)  # * 1e3
            tmrgr_r = ut.conv_time(ut.cactus, ut.cgs, tmrgr_r)  # * 1e3

        except:
            print('Failed *tmerg* extraction for run:{}'.format(run["name"]))
            tmrgr = np.nan
            tmrgr_r = np.nan

        dic_tmerg["name"] = run["name"]
        dic_tmerg["tmrgr"] = tmrgr
        dic_tmerg["tmrgr_r"] = tmrgr_r

        dic_tmergs.append(dic_tmerg)

    return dic_tmergs

def fill_jgw_egw(tmrgr_dic):
    for run, tmrgrs in zip(table, tmrgr_dic):

        basepath = '../sim_anal_tst/' + run["name"]

        tmrgr = tmrgrs["tmrgr"]

        try:
            t, EGW, JGW = np.loadtxt(Paths.ppr_sims +
                                     run["name"] + "/waveforms/EJ.dat",
                                     usecols=(0, 2, 4), unpack=True)
            t = 1e3 * ut.conv_time(ut.cactus, ut.cgs, t)

            idx = np.argmin(np.abs(t - tmrgr - 20.0))

            # if np.isfinite(tcoll) or t[-1] - tmrgr > 21.0:
            #     run["finished"] = "yes"
            # else:
            #     run["finished"] = "no ({} ms to go)".format(
            #         round(21.0 - t[-1] + tmrgr, 2))

            # print('{}, EGW[idx]{} JGW[idx]{}'.format(run["name"], EGW[idx], JGW[idx]))

            run['EGW'] = EGW[idx]
            run['JGW'] = JGW[idx]
        except:
            Printcolor.yellow("Warning: No /waveforms/EJ.dat found")
            run['EGW'] = 'nan'
            run['JGW'] = 'nan'

# Ejecta & Disk

def fill_ej_times(tmrgr_dics, extension='_0'):
    from outflowed import SIM_EJ_HIST

    for run, tmrgr_dic in zip(table, tmrgr_dics):

        try:
            ej = SIM_EJ_HIST(run["name"], extension, True) # geodesic, at det.0 by default
            run['Mej']      = ej.get_par('Mej_tot')
            run["dMej"]     = ej.get_par('dMej_tot')
            run['ThetaRMS'] = ej.get_par('theta_rms')
            run['Yeej']     = ej.get_par('Ye')
            run['Sej']      = ej.get_par('s')
            run['vej']      = ej.get_par('vel_inf')
            run['Ekej']     = ej.get_par('E_kin')
            run['tend']     = ej.time_total_flux[-1]
        except IOError:
            Printcolor.yellow("Warning: Ej IOError for {}".format(run["name"]))

def fill_unb_disk_times(tmrgr_dics):
    from outflowed import SIM_UNBOUND, SIM_DISK

    for run, tmrgr_dic in zip(table, tmrgr_dics):

        tmrgr_r = tmrgr_dic['tmrgr_r']
        run['tmerg_r'] = tmrgr_r

        try:
            unb = SIM_UNBOUND(run["name"])
            run['Munb']     = unb.get_par('Munb_tot')
            run['Munb_bern']= unb.get_par('Munb_bern_tot')
            run['tcoll']    = unb.get_par('tcoll') - tmrgr_r
            run['Mdisk']    = unb.get_par('Mdisk')
        except IOError:
            Printcolor.yellow("Warning: Unbound not IOError for {}".format(run["name"]))

        try:
            disk = SIM_DISK(run["name"])
            run['Mdisk3D'] = disk.get_par('Mdisk_last')
        except IOError:
            Printcolor.yellow("Warning: Disk evol. IOError for {}".format(run["name"]))

if __name__ == '__main__':
    '''-------------------------------| FILLING DATA FROM SELF NAME |------------------------------'''
    load_table(Paths.output + Files.models_empty)

    # get sim names from the source dir
    # initial_data = fill_from_source_dir("/data1/numrel/WhiskyTHC/Backup/2018/GW170817/")

    # for run in table:
    #     print('cp {}/output-0002/parfile.par parfiles/{}.par'.format(run["name"], run["name"]))
    # fill_from_log('./simulations.log')
    # exit(1)

    dics_lorene = get_par_form_parfile()
    # print([dic_lorene["lorene_file"] for dic_lorene in dics_lorene if dic_lorene["name"] == 'DD2_T05_15491205_45km_M0LKcorot'])
    # exit(1)
    fill_self_parameters()

    fill_self_dyn_phase_status()

    # save_table('models.csv')

    '''---------------------------| FILLING DATA FROM LORENE SUMMARY FILES |------------------------'''

    # load_table('models.csv')

    fill_from_lorene(dics_lorene, False)
    # print(table)
    # fill_from_lorene_data('properties_R04.txt', False)
    # fill_from_lorene_data('properties_R03.txt', False)
    # fill_from_lorene_data('properties_R02.txt', False)

    # save_table('models.csv')

    '''------------------------------| FILLING DATA FROM TOV SEQUENCES |----------------------------'''

    fill_from_TOVs()

    # TODO place 'fate' column, Study Waveforms for it, total mass and TOVs

    '''-------------------------------| FILLING DATA FROM GW ANALYSIS |-----------------------------'''

    # load_table('models.csv')

    fill_fpeak()

    dic_tmergs = fill_tmerg()
    fill_jgw_egw(dic_tmergs)

    # save_table('models.csv')

    '''-----------------------------| FILLING DATA FROM EJECTA ANALYSIS |---------------------------'''

    # load_table('models.csv')

    fill_unb_disk_times(dic_tmergs)
    fill_ej_times(dic_tmergs, extension='_0') # geodesic


    save_table(Paths.output + Files.models)