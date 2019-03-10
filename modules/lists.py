# contains the list of paths and other constant information

#############################################################3
class Paths:

    ppr_sims =  '/data01/numrel/vsevolod.nedora/postprocessed/'
    lorene =    '/data01/numrel/vsevolod.nedora/Data/Lorene/'
    TOVs =      '/data01/numrel/vsevolod.nedora/Data/TOVs/'
    gw170817 =  '/data1/numrel/WhiskyTHC/Backup/2018/GW170817/'
    skynet =    '/data01/numrel/vsevolod.nedora/Data/skynet/'
    output =    '/data01/numrel/vsevolod.nedora/output/'
    plots =     '/data01/numrel/vsevolod.nedora/figs/'
    mkn =       '/data01/numrel/vsevolod.nedora/macrokilonova_bayes/mk_source/source/'
    home =      '/data01/numrel/vsevolod.nedora/bns_ppr_scripts/'

class Lists:

    eos = ["DD2", "LS220", "SFHo", "SLy4"]
    res = ["VLR", "LR", "SR", "HR"]
    neut= ["M0", "M1"]

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

    h5_disk = [
    "rho.file *.h5",
    "temperature.file *.h5",
    "Y_e.file *.h5",
    "volform.file *.h5",
    "w_lorentz.file *.h5"
]

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

    # mkn
    mkn_model   = 'mkn_model.h5'
    filt_at2017gfo= 'AT2017gfo.h5'