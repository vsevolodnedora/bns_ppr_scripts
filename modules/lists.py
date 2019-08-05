# contains the list of paths and other constant information

#############################################################3
class Paths:

    scripts  =  '/data01/numrel/vsevolod.nedora/scripts_server/'
    ppr_sims =  '/data01/numrel/vsevolod.nedora/postprocessed3/'
    lorene =    '/data01/numrel/vsevolod.nedora/Data/Lorene/'
    TOVs =      '/data01/numrel/vsevolod.nedora/Data/TOVs/'
    gw170817 =   '/data1/numrel/WhiskyTHC/Backup/2018/SLy4_M130130_SR_physics/' # '/data1/numrel/WhiskyTHC/Backup/2018/GW170817/' # "/data01/numrel/vsevolod.nedora/tmp/" # '/data1/numrel/WhiskyTHC/Backup/2018/SLy4_M130130_SR_physics/'
    skynet =    '/data01/numrel/vsevolod.nedora/Data/skynet/'
    output =    '/data01/numrel/vsevolod.nedora/output/'
    plots =     '/data01/numrel/vsevolod.nedora/figs/'
    mkn =       '/data01/numrel/vsevolod.nedora/macrokilonova_bayes/mk_source/source/'
    home =      '/data01/numrel/vsevolod.nedora/bns_ppr_scripts/'
    SLy4_hydo=  '/data01/numrel/vsevolod.nedora/Data/EOS/SLy4/SLy4_hydro_14-Dec-2017.h5'

class Lists:

    dd2_sim_list = [
        "DD2_M13641364_M0_LK_SR_R04",
        "DD2_M13641364_M0_SR",
    ]


    eos = ["DD2", "LS220", "SFHo", "SLy4"]
    res = ["VLR", "LR", "SR", "HR"]
    neut= ["M0", "M1"]
    visc= ["LK"]
    # q   = [1.0, 1.05, 1.1, 1.11, 1.13, 1.16, 1.19, 1.2, 1.22, 1.29]
    # q   = [1.0, 1.053, 1.102, 1.106, 1.132, 1.159, 1.185, 1.201, 1.222, 1.285]
    colors_q={1.000:"firebrick", 1.053:"red", 1.102:"gold",
              1.106:"darkkhaki", 1.132:"olivedrab", 1.159:"darkgreen",
              1.185:"lightseagreen", 1.201:"darkgreen", 1.222:"royalblue", 1.285:"navy"}

    dyn_not_pas = [
                  "DD2_M13641364_M0_HR_R04", # not in sim2
                  "DD2_M13641364_M0_HR", # not in sim2
                  "DD2_M14861254_M0_HR", # not in sim2
                  "LS220_M14001330_M0_HR", # not in sim2
                  "LS220_M14351298_M0_HR", # in sim2
                  "SLy4_M13641364_M0_HR", # not in sim2
                  "SLy4_M13641364_M0_LR", # in sim2 BUT problem with rocketing dynimical ejecta
                  ]

    bern_pass=[
                  "LS220_M13641364_M0_LK_SR",
                  "LS220_M13641364_M0_SR",
                  "SFHo_M14521283_M0_LR",
                  "SFHo_M14521283_M0_LK_SR",
                  "SFHo_M13641364_M0_LK_SR_2019pizza",
                  "SFHo_M13641364_M0_LK_SR",
                  "SFHo_M14521283_M0_HR"
    ]

    tarball = [
        # "bnstrackergen - bns_positions..asc",
        "dens.norm1.asc",
        "dens_unbnd.norm1.asc",
        "dens_unbnd_bernoulli.norm1.asc",
        "dens_unbnd_garching.norm1.asc",
        "H.norm2.asc",
        "luminosity_nua.norm1.asc",
        "luminosity_nue.norm1.asc",
        "luminosity_nux.norm1.asc",
        "mp_Psi4_l0_m0_r400.00.asc",
        "mp_Psi4_l1_m0_r400.00.asc",
        "mp_Psi4_l1_m1_r400.00.asc",
        "mp_Psi4_l2_m0_r400.00.asc",
        "mp_Psi4_l2_m1_r400.00.asc",
        "mp_Psi4_l2_m2_r400.00.asc",
        "mp_Psi4_l3_m0_r400.00.asc",
        "mp_Psi4_l3_m1_r400.00.asc",
        "mp_Psi4_l3_m2_r400.00.asc",
        "mp_Psi4_l3_m3_r400.00.asc",
        "mp_Psi4_l4_m0_r400.00.asc",
        "mp_Psi4_l4_m1_r400.00.asc",
        "mp_Psi4_l4_m2_r400.00.asc",
        "mp_Psi4_l4_m3_r400.00.asc",
        "mp_Psi4_l4_m4_r400.00.asc",
        "outflow_det_0.asc",
        "outflow_det_1.asc",
        "outflow_det_2.asc",
        "outflow_det_3.asc",
        "rho.maximum.asc",
        "temperature.maximum.asc",
        # "thc_leakagem0 - thc_leakage_m0_flux..asc",
        "outflow_surface_det_0_fluxdens.asc",
        "outflow_surface_det_1_fluxdens.asc",
        # "BH_diagnostics.ah1.gp"
    ]

    outflow = [
        "corr_vel_inf_bern_theta.h5", #
        "corr_vel_inf_theta.h5",
        "corr_ye_entropy.h5",
        "corr_ye_entropy.png",
        "corr_ye_theta.h5",
        "ejecta.h5",
        "ejecta_profile_bern.dat", #
        "ejecta_profile.dat",
        "hist_entropy.dat",
        "hist_entropy.xg",
        "hist_log_rho.dat",
        "hist_log_rho.xg",
        "hist_vel_inf_bern.dat", #
        "hist_temperature.dat",
        "hist_temperature.xg",
        "hist_theta.dat",
        "hist_theta.xg",
        "hist_vel.dat",
        "hist_vel_inf.dat",
        "hist_vel_inf.xg",
        "hist_vel.xg",
        "hist_ye.dat",
        "hist_ye.xg",
        "mass_averages.dat",
        "profile_entropy.xg",
        "profile_flux.xg",
        "profile_rho.xg",
        "profile_temperature.xg",
        "profile_vel_inf.xg",
        "profile_vel.xg",
        "profile_ye.xg",
        "theta_75.dat",
        "total_flux.dat",
        "yields.h5"
    ]

    gw = [
        "E_GW_dot_l2.png",
        "E_GW_dot.png",
        "E_GW_l2.png",
        "E_GW.png",
        "EJ.dat",
        "fpeak.dat",
        "postmerger_psd_l2_m2.dat",
        "postmerger_strain_l2_m2.dat",
        # "psd_l0_m0.dat",
        # "psd_l1_m0.dat",
        # "psd_l1_m1.dat",
        # "psd_l2_m0.dat",
        # "psd_l2_m1.dat",
        "psd_l2_m2.dat",
        # "psd_l3_m0.dat",
        # "psd_l3_m1.dat",
        # "psd_l3_m2.dat",
        # "psd_l3_m3.dat",
        # "psd_l4_m0.dat",
        # "psd_l4_m1.dat",
        # "psd_l4_m2.dat",
        # "psd_l4_m3.dat",
        # "psd_l4_m4.dat",
        # "psi4_l0_m0.dat",
        # "psi4_l1_m0.dat",
        # "psi4_l1_m1.dat",
        # "psi4_l2_m0.dat",
        # "psi4_l2_m1.dat",
        "psi4_l2_m2.dat",
        # "psi4_l3_m0.dat",
        # "psi4_l3_m1.dat",
        # "psi4_l3_m2.dat",
        # "psi4_l3_m3.dat",
        # "psi4_l4_m0.dat",
        # "psi4_l4_m1.dat",
        # "psi4_l4_m2.dat",
        # "psi4_l4_m3.dat",
        # "psi4_l4_m4.dat",
        # "strain_l0_m0.dat",
        # "strain_l1_m0.dat",
        # "strain_l1_m1.dat",
        # "strain_l2_m0.dat",
        # "strain_l2_m1.dat",
        "strain_l2_m2.dat",
        # "strain_l3_m0.dat",
        # "strain_l3_m1.dat",
        # "strain_l3_m2.dat",
        # "strain_l3_m3.dat",
        # "strain_l4_m0.dat",
        # "strain_l4_m1.dat",
        # "strain_l4_m2.dat",
        # "strain_l4_m3.dat",
        # "strain_l4_m4.dat",
        "waveform_l2_m2.dat",
        "tmerger.dat",
        "tcoll.dat"
    ]

    h5_disk = [
        "rho.file *.h5",
        "temperature.file *.h5",
        "Y_e.file *.h5",
        "volform.file *.h5",
        "w_lorentz.file *.h5"
    ]

class Constants:

    ns_rho = 1.6191004634e-5

class Files:
    it_time     = 'collated/rho.maximum.asc'
    models      = 'models.csv'
    models_empty= 'models_tmp2.csv'

    # collated2
    disk_mass   = 'disk_mass.asc'

    # collated
    dens_unb    = 'dens_unbnd.norm1.asc'
    dens_unb_bern = 'dens_unbnd_bernoulli.norm1.asc'
    dens        = 'dens.norm1.asc'

    # ejecta
    total_flux  = 'total_flux.dat'
    hist_theta  = 'hist_theta.dat'
    hist_ye     = 'hist_ye.dat'
    hist_entropy= 'hist_entropy.dat'
    hist_vel_inf= 'hist_vel_inf.dat'

    # nucle
    yields      = 'yields.h5'
    solar_r     = 'solar_r.dat'

    # waveforms
    l2_m2       = 'waveform_l2_m2.dat'
    tmerg       = 'tmerger.dat'

    # ejecta profile:
    ejecta_profile="ejecta_profile.dat"
    ejecta_profile_bern= "ejecta_profile_bern.dat"

    # mkn
    mkn_model   = 'mkn_model.h5'
    filt_at2017gfo= 'AT2017gfo.h5'