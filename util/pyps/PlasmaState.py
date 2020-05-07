import pyplasmastate as pyps
import collections
import numpy as np
import os
from functools import wraps

def flush_decorate(func):
   @wraps(func)
   def func_wrapper(*args, **kwargs):
       args[0]._flush_cache()
       return func(*args, **kwargs)
   return func_wrapper


class PlasmaState(collections.MutableMapping):
    oneD_double = {
        "ang_href_coil", "ang_vref_coil", "anom_evar", "ap2_halfheight",
        "ap2_halfwidth", "ap2_horiz_offset", "ap2_vert_offset",
        "ap_halfheight", "ap_halfwidth", "ap_horiz_offset", "ap_vert_offset",
        "area", "b_halfheight", "b_halfwidth", "b_hdivergence",
        "b_hfocal_length", "b_surfmax", "b_surfmin", "b_vdivergence",
        "b_vfocal_length", "balstab", "chie_anom", "chie_drbm", "chie_etg",
        "chie_glf", "chie_gtcneo", "chie_mmm71", "chie_neo", "chie_paleo",
        "chie_tglf", "chie_trans", "chie_w19", "chii_anom", "chii_drbm",
        "chii_glf", "chii_gtcneo", "chii_mmm71", "chii_neo", "chii_tglf",
        "chii_trans", "chii_w19", "chimo_trans", "chiphi_anom", "chiphi_drbm",
        "chiphi_glf", "chiphi_gtcneo", "chiphi_mmm71", "chiphi_neo",
        "chiphi_tglf", "chiphi_w19", "coil_apt", "coil_resispt", "cur_rw",
        "curbeam", "curech", "curfusn", "curich", "curlh", "curmino",
        "curr_bootstrap", "curr_ohmic", "curt", "d_einj_standard",
        "d_ffull_standard", "d_fhalf_standard", "dcur_rw_dvloop", "difb_fusi",
        "difb_nbi", "difb_rfmi", "dr0_momeq_drho", "dx_fshield",
        "dz0_momeq_drho", "e0_av", "e_anom", "e_anom2", "ec_beam_elongation",
        "ec_half_power_angle", "ec_omode_fraction", "ec_phi_aim",
        "ec_theta_aim", "ee_mobility_factor", "einj_max", "einj_min", "elong",
        "eperp_rw", "epll_rw", "epot", "eta_parallel", "fi_depletion",
        "frac_full", "frac_half", "fracmin", "freq_drbm", "freq_ec",
        "freq_glf", "freq_ic", "freq_lh", "g_eq", "gamma_nc", "gb1", "gb2",
        "gb2i", "gbr2", "gncfb2h", "gncfh", "gr1", "gr2", "gr2i", "gr2rho2",
        "gr3i", "grho1", "grho2", "grho2b2i", "grho2r2i", "grho2r3i", "gri",
        "grirhoi", "growthrate_drbm", "growthrate_exb", "growthrate_glf",
        "h_einj_standard", "h_ffull_standard", "h_fhalf_standard",
        "hsize_coil", "iota", "jdotb", "kvolt_nbi", "ky_anomth", "ky_mmm",
        "lbscap", "lbscap2", "lbsctan", "lpol", "m_all", "m_alla", "m_rfmin",
        "m_s", "m_sa", "m_sfus", "m_sgas", "m_simp0", "m_simpi", "m_snbi",
        "nbeami_bdy", "nfusi_bdy", "ngradb2_av", "ni", "nmini_bdy", "nmodel",
        "nrw", "ns_bdy", "omegat", "p0_reco", "p_eq", "pb0_halo", "pbe", "pbi",
        "pbth", "pcx_halo", "pcx_reco", "pe_anom", "pe_trans", "peech", "pelh",
        "pfuse", "pfusi", "pfusth", "phi_ec_launch", "phibsc", "phit",
        "pi_anom", "pi_trans", "picrf_abs", "picrf_ext", "picth", "pilh",
        "pmine", "pmini", "pminth", "pohme", "power_ec", "power_ic",
        "power_lh", "power_nbi", "power_nbi_trace", "prad", "prad_br",
        "prad_cy", "prad_li", "psc_halo", "psc_reco", "psipol", "psmom_errck",
        "q_all", "q_alla", "q_eq", "q_rfmin", "q_s", "q_sa", "q_sfus",
        "q_sgas", "q_simp0", "q_simpi", "q_snbi", "qatom_all", "qatom_alla",
        "qatom_rfmin", "qatom_s", "qatom_sa", "qatom_sfus", "qatom_sgas",
        "qatom_simp0", "qatom_simpi", "qatom_snbi", "qcomp_e", "qcomp_i",
        "qie", "qrot_conv", "qrot_diff", "r0_momeq", "r_ant", "r_ec_launch",
        "r_grid", "r_midp_in", "r_midp_out", "r_ripple", "r_surfmax",
        "r_surfmin", "res_te", "res_ti", "res_tq", "rho", "rho_anom",
        "rho_anom2", "rho_anomth", "rho_bdy_ns", "rho_ecrf", "rho_eq",
        "rho_eq_geo", "rho_fus", "rho_gas", "rho_icrf", "rho_lhrf", "rho_lmhd",
        "rho_mist", "rho_nbi", "rho_rad", "rho_rw", "rlim", "rloc_coil",
        "rmajor_mean", "rminor_mean", "s0reco_e", "sbedep", "sbehalo", "sc0",
        "squarelo", "squareuo", "srtcen", "surf", "t_einj_standard",
        "t_ffull_standard", "t_fhalf_standard", "th_eq", "ti", "tmodel",
        "tq0_reco", "tq_anom", "tq_trans", "tqb0_halo", "tqbe", "tqbi",
        "tqbjxb", "tqbth", "tqcx_halo", "tqcx_reco", "tqsc_halo", "tqsc_reco",
        "triang", "triang_miller_l", "triang_miller_u", "triangl", "triangu",
        "upwind_pfrac_omega", "upwind_pfrac_te", "upwind_pfrac_ti", "v_loop",
        "vee_trans", "velb_fusi", "velb_nbi", "velb_rfmi", "velni_drbm",
        "velpphi_anom", "velpphi_glf", "velpphi_gtcneo", "velpphi_mmm71",
        "velpphi_neo", "velpphi_tglf", "velte_anom", "velte_glf",
        "velte_gtcneo", "velte_mmm71", "velte_neo", "velte_paleo",
        "velte_tglf", "velti_anom", "velti_glf", "velti_gtcneo", "velti_mmm71",
        "velti_neo", "velti_tglf", "vie_trans", "vmo_trans", "vol", "vphi0_av",
        "vrep_coil_spacing", "vsize_coil", "xs", "z0_momeq", "z_ant",
        "z_ec_launch", "z_grid", "z_mid_ant", "z_midp", "z_ripple",
        "z_surfmax", "z_surfmin", "zbap", "zbsc", "zeff", "zeff_fi", "zeff_th",
        "zlim", "zloc_coil"
    }

    twoD_double = {
        "bphirz", "brrz", "bzrz", "cdicrf", "chie_transcat", "chii_transcat",
        "chimo_transcat", "curech_src", "curlh_src", "d2r_geo_drhodth",
        "d2z_geo_drhodth", "diff_trans", "difns_anom", "difns_drbm",
        "difns_glf", "difns_gtcneo", "difns_mmm71", "difns_neo", "difns_tglf",
        "difns_w19", "dr_geo_drho", "dr_geo_dth", "dxrjcos_momeq_drho",
        "dxrjsin_momeq_drho", "dxzjcos_momeq_drho", "dxzjsin_momeq_drho",
        "dz_geo_drho", "dz_geo_dth", "eperp_beami", "eperp_fusi", "eperp_mini",
        "epll_beami", "epll_fusi", "epll_mini", "freq_mmm71", "freq_tglf",
        "gamma_lmhd", "growthrate_mmm71", "growthrate_tglf", "growthrate_vpar",
        "imag_ant_coef", "kyspectrum_tglf", "n0_halo", "n0_reco", "nbeami",
        "nfusi", "nmini", "ns", "omeg0_halo", "omeg0_reco", "omeg0cx",
        "pe_transcat", "peech_src", "pelh_src", "pi_transcat", "picrf_totals",
        "pilh_src", "psirz", "psmom_nc", "qbeami", "qioniz", "qqcx",
        "r_antgeo", "r_geo", "rate_sinb0i", "rate_sinb0x", "rate_sinf0i",
        "rate_sinf0x", "real_ant_coef", "res_sn", "s0reco", "s0reco_recap",
        "sb0halo", "sb0halo_recap", "sbsce", "sbtherm", "sc_anom", "sfsce",
        "sftherm", "sn_trans", "sprof0e", "t0_halo", "t0_reco", "t0cx",
        "tfrip2_log", "tfrip_log", "tq_transcat", "tqioniz", "tqqcx", "ts",
        "upwind_pfrac_ns", "v_pars", "v_pers", "vee_transcat", "velns_anom",
        "velns_glf", "velns_gtcneo", "velns_mmm71", "velns_neo", "velns_tglf",
        "vie_transcat", "vmo_transcat", "vn_trans", "vpol_inmp", "vpol_omp",
        "vtor_inmp", "vtor_omp", "wt_nphi", "wt_nphi_abs", "wt_nphi_ext",
        "xrjcos_momeq", "xrjsin_momeq", "xzjcos_momeq", "xzjsin_momeq",
        "z_antgeo", "z_geo"
    }

    oneD_int = {
        "all_index", "all_type", "alla_index", "alla_type", "balloon_stable",
        "balloon_status", "deltaw_status", "icircuit_coil",
        "id_bidiff_bpass_co", "id_bidiff_bpass_ctr", "id_bidiff_btrap",
        "id_bidiff_dpass_co", "id_bidiff_dpass_ctr", "id_bidiff_dtrap",
        "id_cdicrf", "id_chie_transcat", "id_chii_transcat",
        "id_chimo_transcat", "id_curech_src", "id_curlh_src", "id_diff_trans",
        "id_difns_anom", "id_difns_drbm", "id_difns_glf", "id_difns_gtcneo",
        "id_difns_mmm71", "id_difns_neo", "id_difns_tglf", "id_difns_w19",
        "id_eperp_beami", "id_eperp_fusi", "id_eperp_mini", "id_epll_beami",
        "id_epll_fusi", "id_epll_mini", "id_fidiff_bpass_co",
        "id_fidiff_bpass_ctr", "id_fidiff_btrap", "id_fidiff_dpass_co",
        "id_fidiff_dpass_ctr", "id_fidiff_dtrap", "id_growthrate_vpar",
        "id_n0_halo", "id_n0_reco", "id_nbeami", "id_nfusi", "id_nmini",
        "id_ns", "id_omeg0_halo", "id_omeg0_reco", "id_omeg0cx",
        "id_pe_transcat", "id_peech_src", "id_pelh_src", "id_pi_transcat",
        "id_picrf_totals", "id_pilh_src", "id_psmom_nc", "id_qbeami",
        "id_qioniz", "id_qqcx", "id_rate_sinb0i", "id_rate_sinb0x",
        "id_rate_sinf0i", "id_rate_sinf0x", "id_res_sn", "id_s0reco",
        "id_s0reco_recap", "id_sb0halo", "id_sb0halo_recap", "id_sbsce",
        "id_sbtherm", "id_sc_anom", "id_sfsce", "id_sftherm", "id_sn_trans",
        "id_sprof0e", "id_t0_halo", "id_t0_reco", "id_t0cx", "id_tq_transcat",
        "id_tqioniz", "id_tqqcx", "id_ts", "id_upwind_pfrac_ns", "id_v_pars",
        "id_v_pers", "id_vee_transcat", "id_velns_anom", "id_velns_glf",
        "id_velns_gtcneo", "id_velns_mmm71", "id_velns_neo", "id_velns_tglf",
        "id_vie_transcat", "id_vmo_transcat", "id_vn_trans", "id_vpol_inmp",
        "id_vpol_omp", "id_vtor_inmp", "id_vtor_omp", "id_xrjcos_momeq",
        "id_xrjsin_momeq", "id_xzjcos_momeq", "id_xzjsin_momeq",
        "is_recycling", "isthermal", "mercier_stable", "n_straps",
        "newcomb_status", "nrz_antgeo", "ns_is_input", "nturns", "num_nphi",
        "num_nphi_vac", "rfmin_to_all", "rfmin_to_alla", "rfmin_type",
        "s_to_sa", "s_type", "sa_index", "sa_type", "sc0_to_sgas",
        "sfus_to_all", "sfus_to_alla", "sfus_type", "sgas_to_s", "sgas_type",
        "simp0_to_s", "simp0_type", "simpi_to_s", "simpi_type", "snbi_to_all",
        "snbi_to_alla", "snbi_type", "tor_mode_no", "tornum_balcrit",
        "trace_flag", "ts_is_input", "vpol_is_input", "vrep_coil_count",
        "vtor_is_input"
    }

    twoD_int = {
        "deltaw_stable", "id_cdicrf_nphi", "id_diff_transcat", "id_n0norm",
        "id_omeg0sc0", "id_picrf_srcs", "id_rate_sinb0xs", "id_rate_sinf0xs",
        "id_sii_xs", "id_sn_transcat", "id_sprof0", "id_t0sc0",
        "id_vn_transcat", "newcomb_stable", "nphi"
    }

    scalar_int = {
        "id_anom_evar", "id_area", "id_b_surfmax", "id_b_surfmin",
        "id_balstab", "id_bphirz", "id_brrz", "id_bzrz", "id_chie_anom",
        "id_chie_drbm", "id_chie_etg", "id_chie_glf", "id_chie_gtcneo",
        "id_chie_mmm71", "id_chie_neo", "id_chie_paleo", "id_chie_tglf",
        "id_chie_trans", "id_chie_w19", "id_chii_anom", "id_chii_drbm",
        "id_chii_glf", "id_chii_gtcneo", "id_chii_mmm71", "id_chii_neo",
        "id_chii_tglf", "id_chii_trans", "id_chii_w19", "id_chimo_trans",
        "id_chiphi_anom", "id_chiphi_drbm", "id_chiphi_glf",
        "id_chiphi_gtcneo", "id_chiphi_mmm71", "id_chiphi_neo",
        "id_chiphi_tglf", "id_chiphi_w19", "id_cur_rw", "id_curbeam",
        "id_curech", "id_curfusn", "id_curich", "id_curlh", "id_curmino",
        "id_curr_bootstrap", "id_curr_ohmic", "id_curt", "id_dcur_rw_dvloop",
        "id_difb_fusi", "id_difb_nbi", "id_difb_rfmi", "id_ee_mobility_factor",
        "id_elong", "id_eperp_rw", "id_epll_rw", "id_epot", "id_eta_parallel",
        "id_fi_depletion", "id_freq_drbm", "id_freq_glf", "id_freq_mmm71",
        "id_freq_tglf", "id_g_eq", "id_gamma_nc", "id_gb1", "id_gb2",
        "id_gb2i", "id_gbr2", "id_gncfb2h", "id_gncfh", "id_gr1", "id_gr2",
        "id_gr2i", "id_gr2rho2", "id_gr3i", "id_grho1", "id_grho2",
        "id_grho2b2i", "id_grho2r2i", "id_grho2r3i", "id_gri", "id_grirhoi",
        "id_growthrate_drbm", "id_growthrate_exb", "id_growthrate_glf",
        "id_growthrate_mmm71", "id_growthrate_tglf", "id_iota", "id_jdotb",
        "id_kyspectrum_tglf", "id_lpol", "id_ngradb2_av", "id_ni", "id_nmodel",
        "id_nrw", "id_omegat", "id_p0_reco", "id_p_eq", "id_pb0_halo",
        "id_pbe", "id_pbi", "id_pbth", "id_pcx_halo", "id_pcx_reco",
        "id_pe_anom", "id_pe_trans", "id_peech", "id_pelh", "id_pfuse",
        "id_pfusi", "id_pfusth", "id_phit", "id_pi_anom", "id_pi_trans",
        "id_picth", "id_pilh", "id_pmine", "id_pmini", "id_pminth", "id_pohme",
        "id_prad", "id_prad_br", "id_prad_cy", "id_prad_li", "id_psc_halo",
        "id_psc_reco", "id_psipol", "id_psirz", "id_psmom_errck", "id_q_eq",
        "id_qcomp_e", "id_qcomp_i", "id_qie", "id_qrot_conv", "id_qrot_diff",
        "id_r0_momeq", "id_r_geo", "id_r_midp_in", "id_r_midp_out",
        "id_r_surfmax", "id_r_surfmin", "id_res_te", "id_res_ti", "id_res_tq",
        "id_rmajor_mean", "id_rminor_mean", "id_s0reco_e", "id_sbedep",
        "id_sbehalo", "id_squarelo", "id_squareuo", "id_surf", "id_tfrip2_log",
        "id_tfrip_log", "id_ti", "id_tmodel", "id_tq0_reco", "id_tq_anom",
        "id_tq_trans", "id_tqb0_halo", "id_tqbe", "id_tqbi", "id_tqbjxb",
        "id_tqbth", "id_tqcx_halo", "id_tqcx_reco", "id_tqsc_halo",
        "id_tqsc_reco", "id_triang", "id_triang_miller_l",
        "id_triang_miller_u", "id_triangl", "id_triangu",
        "id_upwind_pfrac_omega", "id_upwind_pfrac_te", "id_upwind_pfrac_ti",
        "id_v_loop", "id_vee_trans", "id_velb_fusi", "id_velb_nbi",
        "id_velb_rfmi", "id_velni_drbm", "id_velpphi_anom", "id_velpphi_glf",
        "id_velpphi_gtcneo", "id_velpphi_mmm71", "id_velpphi_neo",
        "id_velpphi_tglf", "id_velte_anom", "id_velte_glf", "id_velte_gtcneo",
        "id_velte_mmm71", "id_velte_neo", "id_velte_paleo", "id_velte_tglf",
        "id_velti_anom", "id_velti_glf", "id_velti_gtcneo", "id_velti_mmm71",
        "id_velti_neo", "id_velti_tglf", "id_vie_trans", "id_vmo_trans",
        "id_vol", "id_z0_momeq", "id_z_geo", "id_z_midp", "id_z_surfmax",
        "id_z_surfmin", "id_zeff", "id_zeff_fi", "id_zeff_th", "kccw_bphi",
        "kccw_jphi", "lmhd_status", "lock_machine_descr", "lock_shot_config",
        "lock_sim_init", "max_n_straps", "max_nrz_antgeo", "max_num_nphi",
        "max_num_nphi_vac", "mhd_eq_status", "mmodes", "n_d_einj_standard",
        "n_h_einj_standard", "n_t_einj_standard", "nballoon", "nbeam",
        "ncircuits", "ncoils", "ndelta_w", "necrf_src", "nefi_anom",
        "nefi_anom2", "neqmom", "ngsc0", "nicrf_src", "nky_anomth", "nky_mmm",
        "nlhrf_src", "nmdifb", "nmom", "nnewcomb", "npsmom", "nr", "nr_rip",
        "nrho", "nrho_anom", "nrho_anom2", "nrho_anomth", "nrho_ecrf",
        "nrho_eq", "nrho_eq_geo", "nrho_fus", "nrho_gas", "nrho_icrf",
        "nrho_lhrf", "nrho_lmhd", "nrho_mist", "nrho_nbi", "nrho_rad",
        "nrho_rw", "nspec_all", "nspec_alla", "nspec_beam", "nspec_fusion",
        "nspec_gas", "nspec_imp0", "nspec_impi", "nspec_rfmin", "nspec_th",
        "nspec_tha", "ntf_coil2", "ntf_coils", "nth_eq", "ntrcat", "num_rzlim",
        "nxs", "nz", "nz_rip", "shot_number", "z0max", "zimp1", "zimp2"
    }

    scalar_double = {
        "b_axis", "b_axis_vac", "b_max_lcfs", "b_min_lcfs", "dn0out", "eldot",
        "nebar", "nmodel_bdy", "omegat_bdy", "psi_to_machine_axis", "r_axis",
        "r_max_box", "r_max_lcfs", "r_min_box", "r_min_lcfs", "rho_bdy_omegat",
        "rho_bdy_te", "rho_bdy_ti", "t0", "t1", "te_bdy", "tfinal", "ti_bdy",
        "tinit", "tmodel_bdy", "vsur", "z_axis", "z_max_box", "z_max_lcfs",
        "z_min_box", "z_min_lcfs"
    }

    oneD_string = {
        "all_name", "alla_name", "ant_model", "balloon_code", "balloon_mapper",
        "beam_type", "circuit_name", "coil_in_circuit", "coil_name",
        "deltaw_code", "deltaw_mapper", "dist_fun", "ecrf_src_name",
        "eqmom_num", "gas_atom", "gs_name", "icrf_src_name", "lhrf_src_name",
        "nbap2_shape", "nbap_shape", "nbi_src_name", "nbion", "nbion_trace",
        "nbshape", "newcomb_code", "newcomb_mapper", "psmom_num",
        "ql_operator", "rfmin_name", "s_name", "sa_name", "sfus_name",
        "sgas_name", "simp0_name", "simpi_name", "snbi_name", "tor_mode_label",
        "trcat", "trmodels", "xs_name"
    }

    scalar_string = {
        "anom_code_info", "cur_data_info", "ec_code_info", "ec_data_info",
        "eq_code_info", "eq_data_info", "eqdsk_file", "fus_code_info",
        "gas_code_info", "geometry", "global_label", "ic_code_info",
        "ic_data_info", "kdens_rfmin", "lh_code_info", "lh_data_info",
        "lmhd_eq_code", "nbi_code_info", "nbi_data_info", "ns_data_info",
        "plasma_code_info", "rad_code_info", "rad_data_info",
        "runaway_code_info", "runid", "tf_data_info", "tokamak_id",
        "ts_data_info", "vpol_data_info", "vsur_data_info", "vtor_data_info",
        "zeff_data_info"
    }

    threeD_int = {"id_picrf_nphi_srcs"}

    threeD_double = {
        "bidiff_bpass_co", "bidiff_bpass_ctr", "bidiff_btrap",
        "bidiff_dpass_co", "bidiff_dpass_ctr", "bidiff_dtrap", "cdicrf_nphi",
        "diff_transcat", "fidiff_bpass_co", "fidiff_bpass_ctr", "fidiff_btrap",
        "fidiff_dpass_co", "fidiff_dpass_ctr", "fidiff_dtrap", "n0norm",
        "omeg0sc0", "picrf_srcs", "rate_sinb0xs", "rate_sinf0xs", "sii_xs",
        "sn_transcat", "sprof0", "t0sc0", "vn_transcat"
    }

    scalar_fun_map = {}
    P = pyps.PlasmaState
    scalar_fun_map.update({k: pyps.PlasmaState.getDataInt for k in scalar_int})
    scalar_fun_map.update(
        {k: pyps.PlasmaState.getDataDouble
         for k in scalar_double})
    scalar_fun_map.update(
        {k: pyps.PlasmaState.getDataString
         for k in scalar_string})

    oneD_fun_map = {}
    oneD_fun_map.update({
        k: (pyps.PlasmaState.getDataIntVec, pyps.IntVector, int)
        for k in oneD_int
    })
    oneD_fun_map.update({
        k: (pyps.PlasmaState.getDataDoubleVec, pyps.DoubleVector, float)
        for k in oneD_double
    })
    oneD_fun_map.update({
        k: (pyps.PlasmaState.getDataStringVec, pyps.StringVector, str)
        for k in oneD_string
    })

    twoD_fun_map = {}
    twoD_fun_map.update({
        k: (pyps.PlasmaState.getDataIntVec, pyps.MultiArrayInt, pyps.IntVector,
            int)
        for k in twoD_int
    })
    twoD_fun_map.update({
        k: (pyps.PlasmaState.getDataDoubleVec, pyps.MultiArrayDouble,
            pyps.DoubleVector, float)
        for k in twoD_double
    })

    threeD_fun_map = {}
    threeD_fun_map.update({
        k: (pyps.PlasmaState.getDataIntVec, pyps.MultiArrayInt, pyps.IntVector,
            int)
        for k in threeD_int
    })
    threeD_fun_map.update({
        k: (pyps.PlasmaState.getDataDoubleVec, pyps.MultiArrayDouble,
            pyps.DoubleVector, float)
        for k in threeD_double
    })

    def __init__(self, label, debug_level):
        self._ps = pyps.PlasmaState(label, debug_level)
        self.cache = {}

    def __getitem__(self, item):

        key = item.lower()
        if key in list(self.cache.keys()):
            return self.cache[key]
        else:
            retval = self._get_item_internal(key)
            self.cache[key] = retval
            return retval


    def __delitem__(self, key):
        pass

    def _get_item_internal(self, name):
        key = name.lower()
        try:
            retval = self.scalar_fun_map[key](self._ps, key)
            if type(retval) == str:
                return retval.strip()
            else:
                return retval
        except KeyError:
            pass
            #print "%s is not a Plasma State scalar" % key

        try:
            return np.array(self.oneD_fun_map[key][0](self._ps, key))
        except KeyError:
            pass
            #print "%s is not a Plasma State vector" % key

        try:
            vec = np.array(self.twoD_fun_map[key][0](self._ps, key))
            shape = self._ps.getDims(key)
            m = vec.reshape(shape[::-1], order='C')
            return m
        except KeyError:
            pass
            #print "%s is not a Plasma State 2-D array" % key

        try:
            vec = np.array(self.threeD_fun_map[key][0](self._ps, key))
            shape = self._ps.getDims(key)
            m = vec.reshape(shape[::-1], order='C')
            return m
        except KeyError:
            pass
            #print "%s is not a Plasma State 3-D array" % key

        raise KeyError("Unknown Plasma state variable : %d" % key)

    def _flush_cache(self):
        for (k,v) in self.cache.items():
            self._set_item_internal(k, v)

    def __setitem__(self, name, value):
        key = name.lower()
        try:
            del self.cache[key]
        except KeyError:
            pass
        try:
            self._set_item_internal(key, value)
        except Exception :
            raise

    def _set_item_internal(self, key, value):
        try:
            getf = self.scalar_fun_map[key]
            return self._ps.setData(key, value)
        except KeyError:
            #print "%s is not a Plasma State scalar" % key
            pass
        try:
            (_, _, d_type) = self.oneD_fun_map[key]
            # Convert to native Python from possible numpy types
            value_list = [d_type(v) for v in value]
            return self._ps.setData(key, value_list)
        except KeyError:
            pass
            #print "%s is not a Plasma State one-D vector" % key

        try:
            (_, ma_class, vec_class, d_type) = self.twoD_fun_map[key]
            if type(value) != np.ndarray:
                raise TypeError(
                    "Only Numpy arrays can be used to set 2-D plasma state variables")
            if len(value.shape) != 2:
                raise ValueError(
                    "Only 2-D Numpy arrays can be used to set 2-D plasma state variables")
            # Convert to native Python from possible numpy types
            array_as_list = [d_type(v) for v in value.flatten(order='C')]
            ma = ma_class()
            ma.setDims(value.shape[::-1])
            ma.setData(vec_class(array_as_list))
            return self._ps.setData(key, ma)
        except KeyError:
            pass
            #print "%s is not a Plasma State two-D array" % key

        try:
            (_, ma_class, vec_class, d_type) = self.threeD_fun_map[key]
            if type(value) != np.ndarray:
                raise TypeError(
                    "Only Numpy arrays can be used to set 3-D plasma state variables")
            if len(value.shape) != 3:
                raise ValueError(
                    "Only 3-D Numpy arrays can be used to set 3-D plasma state variables")
            # Convert to native Python from possible numpy types
            array_as_list = [d_type(v) for v in value.flatten(order='C')]
            ma = ma_class()
            ma.setDims(value.shape[::-1])
            ma.setData(vec_class(array_as_list))
            return self._ps.setData(key, ma)
        except KeyError:
            pass
            #print "%s is not a Plasma State three-D array" % key
        raise KeyError("Unknown Plasma state variable : %d" % key)

    def __iter__(self):
        pass

    def __len__(self):
        pass

    def read(self, filename):
        if not os.path.isfile(filename):
            raise IOError('Cannot access file %s' % filename)
        self._ps.getState(filename)
        self.cache = {}
        return

    @flush_decorate
    def store(self, filename):
        self._ps.store(filename)
        return

    @flush_decorate
    def setFusionSpecies(self, zAtom, zCharge, amu):
        self._ps.setFusionSpecies(zAtom, zCharge, amu)
        return

    @flush_decorate
    def alloc(self):
        self._ps.alloc()
        return

    @flush_decorate
    def deriveMhdEq(self, action, rhowk=None):
        if rhowk:
            vec = pyps.DoubleVector([float(v) for v in rhowk])
        else:
            vec = pyps.DoubleVector()
        return self._ps.deriveMhdEq(action, vec)

    @flush_decorate
    def finishSpecies(self):
        return self._ps.finishSpecies()

    @flush_decorate
    def getVersionId(self):
        return self._ps.getVersionId()

    @flush_decorate
    def hashDiff(self, lbl1, obj2, lbl2, sstr):
        return self._ps.hashDiff(lbl1, obj2._ps, lbl2, sstr)

    @flush_decorate
    def interp1d(self, nm, deriv, xs, index = None):
        vec = pyps.DoubleVector([float(v) for v in xs])
        if index != None:
            ret_vec = self._ps.interp1d(nm, deriv, vec, index)
        else:
            ret_vec = self._ps.interp1d(nm, deriv, vec)
        return np.array(ret_vec)

    @flush_decorate
    def interp2d(self, nm, deriv1, deriv2, x1s, x2s, index = None):
        x1s_vec = pyps.DoubleVector([float(v) for v in x1s])
        x2s_vec = pyps.DoubleVector([float(v) for v in x2s])
        if index:
            ret_vec = self._ps.interp2d(nm, deriv1, deriv2, x1s_vec, x2s_vec,
                                        index)
        else:
            ret_vec = self._ps.interp2d(nm, deriv1, deriv2, x1s_vec, x2s_vec)
        return np.array(ret_vec)

    @flush_decorate
    def readMachDescr(self, filepath, action, g_filepath=""):
        if (len(filepath) >= 256):
            raise ValueError('filepath length exceeds 256 characters')
        if len(g_filepath) >= 256:
            raise ValueError('g_filepath length exceeds 256 characters')
        if action not in ['NEW', 'INIT']:
            raise ValueError('action can only be \'NEW\' or \'INIT\'')
        if (g_filepath == ""):
            self._ps.readMachDescr(filepath, action)
        else:
            self._ps.readMachDescr(filepath, action, g_filepath)
        return

    @flush_decorate
    def readShotConfig(self, filepath):
        if not os.path.isfile(filepath):
            raise IOError('Cannot access file %s' % filepath)
        self._ps.readShotConfig(filepath)
        return

    @flush_decorate
    def setFusionSpecies(self, zAtom, zCharge, amu):
        self._ps.setFusionSpecies(zAtom, zCharge, amu)
        return

    @flush_decorate
    def setImpuritySpecies(self, zAtom, zCharge, amu):
        self._ps.setImpuritySpecies(zAtom, zCharge, amu)
        return

    @flush_decorate
    def setNeutralBeamSpecies(self, zAtom, zCharge, amu):
        self._ps.setNeutralBeamSpecies(zAtom, zCharge, amu)
        return

    @flush_decorate
    def setRFMinoritySpecies(self, zAtom, zCharge, amu):
        self._ps.setRFMinoritySpecies(zAtom, zCharge, amu)
        return

    @flush_decorate
    def setThermalSpecies(self, zAtom, zCharge, amu):
        self._ps.setThermalSpecies(zAtom, zCharge, amu)
        return

    @flush_decorate
    def updateEquilibrium(self, g_filepath, bdy_crat, kcur_option, rho_curbrk):
        if not os.path.isfile(g_filepath):
            raise IOError('Cannot access file %s' % g_filepath)
        self._ps.updateEquilibrium(g_filepath, bdy_crat, kcur_option,
                                   rho_curbrk)
        return

    @flush_decorate
    def updateHashCodes(self, mdescr_flag, sconfig_flag, simInit_flag):
        self._ps.updateHashCodes(mdescr_flag, sconfig_flag, simInit_flag)
        return

    @flush_decorate
    def verifyMhdEq(self,
                    rhowk=None,
                    reltol=0.02,
                    update_phit=False,
                    update_psi=False,
                    update_q=False,
                    check_phit=True,
                    check_q=True):

        if rhowk:
            rhowk_vec = pyps.DoubleVector(rhowk)
        else:
            rhowk_vec = pyps.DoubleVector()
        self._ps.verifyMhdEq(rhowk_vec, reltol, update_phit, update_psi,
                             update_q, check_phit, check_q)
        return

    @flush_decorate
    def writeGeqdsk(self, filename):
        self._ps.writeGeqdsk(filename)
        return

    @flush_decorate
    def writeMachDescr(self, filename):
        self._ps.writeMachDescr(filename)
        return

    @flush_decorate
    def writeShotConfig(self, filename):
        self._ps.writeShotConfig(filename)
        return

    @flush_decorate
    def writeUpdateFile(self, filename, hashflag):
        self._ps.writeUpdateFile(filename, hashflag)
        return
