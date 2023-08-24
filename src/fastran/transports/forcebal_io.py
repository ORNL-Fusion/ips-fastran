import os
import sys
import numpy   as npy
import netCDF4 as ncdf

from fastran.equilibrium import cheasefiles

def get_forcebal_vars():
    #>============================================
    #>Notes on meaning of subscripts:
    #>============================================
    #_r             Radial grid (all variables presently on same grid)
    #_p             Poloidal
    #_t             Toroidal
    #_in            INside of torus in axial plane
    #_out           OUTside of torus in axial plane
    #_par           PARallel to B
    #_e             Electron value
    #_i             Ion (main) value
    #_im            IMpurity value
    #_con           CONduction
    #_bs            BootStrap
    #_wp            Ware Pinch
    #_g             <Grad(r_t)**2> normalization
    #_n             deNsity term

    forcebal = {}

    #>============================================
    #> INPUT PARAMETERS
    #>============================================

    forcebal['l_banana']         = {}
    forcebal['l_banana']['data'] = None
    forcebal['l_banana']['unit'] = None
    forcebal['l_banana']['info'] = "Include banana viscosity (1=yes)"

    forcebal['l_pfirsch']         = {}
    forcebal['l_pfirsch']['data'] = None
    forcebal['l_pfirsch']['unit'] = None
    forcebal['l_pfirsch']['info'] = "Include Pfirsch-Schlueter viscosity (1=yes)"

    forcebal['l_classical']         = {}
    forcebal['l_classical']['data'] = None
    forcebal['l_classical']['unit'] = None
    forcebal['l_classical']['info'] = "Include classical transport (1=yes)"

    forcebal['l_potato']         = {}
    forcebal['l_potato']['data'] = None
    forcebal['l_potato']['unit'] = None
    forcebal['l_potato']['info'] = "Include potato viscosity (1=yes)"

    forcebal['l_squeeze']         = {}
    forcebal['l_squeeze']['data'] = None
    forcebal['l_squeeze']['unit'] = None
    forcebal['l_squeeze']['info'] = "Include orbit squeezing (1=yes)"

    forcebal['diiid']         = {}
    forcebal['diiid']['data'] = None
    forcebal['diiid']['unit'] = None
    forcebal['diiid']['info'] = "Device Name"

    forcebal['id_shot']         = {}
    forcebal['id_shot']['data'] = None
    forcebal['id_shot']['unit'] = None
    forcebal['id_shot']['info'] = "Discharge Shot Number"

    forcebal['time']         = {}
    forcebal['time']['data'] = None
    forcebal['time']['unit'] = "s"
    forcebal['time']['info'] = "Time"

    forcebal['nr_r']         = {}
    forcebal['nr_r']['data'] = None
    forcebal['nr_r']['unit'] = None
    forcebal['nr_r']['info'] = "Number of Radial Grid"

    #>============================================
    #> PLASMA PARAMETERS
    #>============================================

    forcebal['a0']         = {}
    forcebal['a0']['data'] = None
    forcebal['a0']['unit'] = "m"
    forcebal['a0']['info'] = "Plama Minor Radius at Midplane"

    forcebal['a1']         = {}
    forcebal['a1']['data'] = None
    forcebal['a1']['unit'] = "m"
    forcebal['a1']['info'] = "Reference Plamsa Minor Radius"

    forcebal['R0']         = {}
    forcebal['R0']['data'] = None
    forcebal['R0']['unit'] = "m"
    forcebal['R0']['info'] = "Reference Plasma Major Radius at Midplane"

    forcebal['Bt0']         = {}
    forcebal['Bt0']['data'] = None
    forcebal['Bt0']['unit'] = "T"
    forcebal['Bt0']['info'] = "Toroidal Field at R0"

    forcebal['elongation']         = {}
    forcebal['elongation']['data'] = None
    forcebal['elongation']['unit'] = None
    forcebal['elongation']['info'] = "Plasma Elongation at Edge"

    #>============================================
    #> RADIAL GRID PARAMETERS
    #>============================================

    forcebal['rho_p']         = {}
    forcebal['rho_p']['data'] = None
    forcebal['rho_p']['unit'] = None
    forcebal['rho_p']['info'] = "Normalized Poloidal Flux Grid - Proportional to Poloidal Flux"

    forcebal['rho_t']         = {}
    forcebal['rho_t']['data'] = None
    forcebal['rho_t']['unit'] = None
    forcebal['rho_t']['info'] = "Normalized Toroidal Flux Grid - Proportional to Square-Root Toroidal Flux"

    forcebal['r_p']         = {}
    forcebal['r_p']['data'] = None
    forcebal['r_p']['unit'] = "m"
    forcebal['r_p']['info'] = "a1 * Normalized Poloidal Flux Grid"

    forcebal['r_t']         = {}
    forcebal['r_t']['data'] = None
    forcebal['r_t']['unit'] = "m"
    forcebal['r_t']['info'] = "a1 * Normalized Toroidal Flux Grid"

    forcebal['R_in']         = {}
    forcebal['R_in']['data'] = None
    forcebal['R_in']['unit'] = "m"
    forcebal['R_in']['info'] = "Major Radius at Inboard Midplane"

    forcebal['R_o']         = {}
    forcebal['R_o']['data'] = None
    forcebal['R_o']['unit'] = "m"
    forcebal['R_o']['info'] = "Major Radius at Outboard Midplane"

    #>============================================
    #> MAGNETIC FIELDS AND FLUXES
    #>============================================

    forcebal['B_p_o']         = {}
    forcebal['B_p_o']['data'] = None
    forcebal['B_p_o']['unit'] = "T"
    forcebal['B_p_o']['info'] = "Poloidal Magnetic Field Outboard Midplane"

    forcebal['B_t_o']         = {}
    forcebal['B_t_o']['data'] = None
    forcebal['B_t_o']['unit'] = "T"
    forcebal['B_t_o']['info'] = "Toroidal Magnetic Field Outboard Midplane"

    forcebal['Phi_t']         = {}
    forcebal['Phi_t']['data'] = None
    forcebal['Phi_t']['unit'] = "Wb"
    forcebal['Phi_t']['info'] = "Toroidal magnetic flux"

    forcebal['Psi']         = {}
    forcebal['Psi']['data'] = None
    forcebal['Psi']['unit'] = "Wb/rad"
    forcebal['Psi']['info'] = "Poloidal Magnetic Flux Divided by (2*pi)"

    forcebal['q']         = {}
    forcebal['q']['data'] = None
    forcebal['q']['unit'] = None
    forcebal['q']['info'] = "Safety Factor"

    #>============================================
    #> FLUX FUNCTIONS AND METRICS
    #>============================================

    forcebal['F']         = {}
    forcebal['F']['data'] = None
    forcebal['F']['unit'] = "A"
    forcebal['F']['info'] = "Poloidal Current External to Surface"

    forcebal['f_trap']         = {}
    forcebal['f_trap']['data'] = None
    forcebal['f_trap']['unit'] = None
    forcebal['f_trap']['info'] = "Trapped Particles Fraction"

    forcebal['grad_rho_sq']         = {}
    forcebal['grad_rho_sq']['data'] = None
    forcebal['grad_rho_sq']['unit'] = None
    forcebal['grad_rho_sq']['info'] = "Metric, <grad(rho)**2>"

    forcebal['grad_rho']         = {}
    forcebal['grad_rho']['data'] = None
    forcebal['grad_rho']['unit'] = None
    forcebal['grad_rho']['info'] = "Metric, <|grad(rho)|>"

    forcebal['dV7dr_t']         = {}
    forcebal['dV7dr_t']['data'] = None
    forcebal['dV7dr_t']['unit'] = "m^2"
    forcebal['dV7dr_t']['info'] = "Radial derivative of volume, dV/drho/a1"

    forcebal['Vol']         = {}
    forcebal['Vol']['data'] = None
    forcebal['Vol']['unit'] = "m^3"
    forcebal['Vol']['info'] = "Enclosed Plasma Volume"

    #>============================================
    #> PLASMA PROFILES
    #>============================================

    forcebal['den_e1']         = {}
    forcebal['den_e1']['data'] = None
    forcebal['den_e1']['unit'] = "/m^3"
    forcebal['den_e1']['info'] = "Density of e1"

    forcebal['T_e1']         = {}
    forcebal['T_e1']['data'] = None
    forcebal['T_e1']['unit'] = "keV"
    forcebal['T_e1']['info'] = "Temperature of e1"

    forcebal['prs_e1']         = {}
    forcebal['prs_e1']['data'] = None
    forcebal['prs_e1']['unit'] = "keV/m^3"
    forcebal['prs_e1']['info'] = "Pressure of e1"

    forcebal['den_D1']         = {}
    forcebal['den_D1']['data'] = None
    forcebal['den_D1']['unit'] = "/m^3"
    forcebal['den_D1']['info'] = "Density of D1"

    forcebal['T_D1']         = {}
    forcebal['T_D1']['data'] = None
    forcebal['T_D1']['unit'] = "keV"
    forcebal['T_D1']['info'] = "Temperature of D1"

    forcebal['prs_D1']         = {}
    forcebal['prs_D1']['data'] = None
    forcebal['prs_D1']['unit'] = "keV/m^3"
    forcebal['prs_D1']['info'] = "Pressure of D1"

    forcebal['den_C6']         = {}
    forcebal['den_C6']['data'] = None
    forcebal['den_C6']['unit'] = "/m^3"
    forcebal['den_C6']['info'] = "Density of C6"

    forcebal['T_C6']         = {}
    forcebal['T_C6']['data'] = None
    forcebal['T_C6']['unit'] = "keV"
    forcebal['T_C6']['info'] = "Temperature of C6"

    forcebal['prs_C6']         = {}
    forcebal['prs_C6']['data'] = None
    forcebal['prs_C6']['unit'] = "keV/m^3"
    forcebal['prs_C6']['info'] = "Pressure of C6"

    forcebal['Z_eff_ex']         = {}
    forcebal['Z_eff_ex']['data'] = None
    forcebal['Z_eff_ex']['unit'] = None
    forcebal['Z_eff_ex']['info'] = "Effective Charge, Zeff - Experimental Data"

    forcebal['Z_eff']         = {}
    forcebal['Z_eff']['data'] = None
    forcebal['Z_eff']['unit'] = None
    forcebal['Z_eff']['info'] = "Effective Charge, Zeff"

    forcebal['Sqz_e1']         = {}
    forcebal['Sqz_e1']['data'] = None
    forcebal['Sqz_e1']['unit'] = None
    forcebal['Sqz_e1']['info'] = "Orbit Squeezing Factor for e1"

    forcebal['Gam_e1']         = {}
    forcebal['Gam_e1']['data'] = None
    forcebal['Gam_e1']['unit'] = "#/m^2/s"
    forcebal['Gam_e1']['info'] = "Particle Flux of e1"

    forcebal['D_eff_e1']         = {}
    forcebal['D_eff_e1']['data'] = None
    forcebal['D_eff_e1']['unit'] = "m^2/s"
    forcebal['D_eff_e1']['info'] = "Effective Particle Diffusivity for e1"

    forcebal['D_eff_g_e1']         = {}
    forcebal['D_eff_g_e1']['data'] = None
    forcebal['D_eff_g_e1']['unit'] = "m^2/s"
    forcebal['D_eff_g_e1']['info'] = "Effective Particle Diffusivity/g for e1"

    forcebal['D_n_e1']         = {}
    forcebal['D_n_e1']['data'] = None
    forcebal['D_n_e1']['unit'] = "m^2/s"
    forcebal['D_n_e1']['info'] = "Diagonal Particle Diffusivity for e1"

    forcebal['v_n_e1']         = {}
    forcebal['v_n_e1']['data'] = None
    forcebal['v_n_e1']['unit'] = "m/s"
    forcebal['v_n_e1']['info'] = "Particle Pinch from grad(T,n) for e1"

    forcebal['v_wp_e1']         = {}
    forcebal['v_wp_e1']['data'] = None
    forcebal['v_wp_e1']['unit'] = "m/s"
    forcebal['v_wp_e1']['info'] = "Ware Pinch for e1"

    forcebal['v_nwp_e1']         = {}
    forcebal['v_nwp_e1']['data'] = None
    forcebal['v_nwp_e1']['unit'] = "m/s"
    forcebal['v_nwp_e1']['info'] = "Total Pinch for e1"

    forcebal['D_n_g_e1']         = {}
    forcebal['D_n_g_e1']['data'] = None
    forcebal['D_n_g_e1']['unit'] = "m^2/s"
    forcebal['D_n_g_e1']['info'] = "Diagonal Particle Diffusivity/g for e1"

    forcebal['v_n_g_e1']         = {}
    forcebal['v_n_g_e1']['data'] = None
    forcebal['v_n_g_e1']['unit'] = "m/s"
    forcebal['v_n_g_e1']['info'] = "Particle Pinch from grad(T,n)/g for e1"

    forcebal['v_wp_g_e1']         = {}
    forcebal['v_wp_g_e1']['data'] = None
    forcebal['v_wp_g_e1']['unit'] = "m/s"
    forcebal['v_wp_g_e1']['info'] = "Ware Pinch/g for e1"

    forcebal['v_nwp_g_e1']         = {}
    forcebal['v_nwp_g_e1']['data'] = None
    forcebal['v_nwp_g_e1']['unit'] = "m/s"
    forcebal['v_nwp_g_e1']['info'] = "Total Pinch/g for e1"

    forcebal['g_Te_e1']         = {}
    forcebal['g_Te_e1']['data'] = None
    forcebal['g_Te_e1']['unit'] = None
    forcebal['g_Te_e1']['info'] = "Te Exponent in the Density Profile of e1"

    forcebal['g_Ti_e1']         = {}
    forcebal['g_Ti_e1']['data'] = None
    forcebal['g_Ti_e1']['unit'] = None
    forcebal['g_Ti_e1']['info'] = "Te Exponent in the Density Profile of e1"

    forcebal['g_ne1_e1']         = {}
    forcebal['g_ne1_e1']['data'] = None
    forcebal['g_ne1_e1']['unit'] = None
    forcebal['g_ne1_e1']['info'] = "e1 Exponent in the Density Profile of e1"

    forcebal['g_nD1_e1']         = {}
    forcebal['g_nD1_e1']['data'] = None
    forcebal['g_nD1_e1']['unit'] = None
    forcebal['g_nD1_e1']['info'] = "D1 Exponent in the Density Profile of e1"

    forcebal['g_nC6_e1']         = {}
    forcebal['g_nC6_e1']['data'] = None
    forcebal['g_nC6_e1']['unit'] = None
    forcebal['g_nC6_e1']['info'] = "C6 Exponent in the Density Profile of e1"

    forcebal['Sqz_D1']         = {}
    forcebal['Sqz_D1']['data'] = None
    forcebal['Sqz_D1']['unit'] = None
    forcebal['Sqz_D1']['info'] = "Orbit Squeezing Factor for D1"

    forcebal['Gam_D1']         = {}
    forcebal['Gam_D1']['data'] = None
    forcebal['Gam_D1']['unit'] = "#/m^2/s"
    forcebal['Gam_D1']['info'] = "Particle Flux of D1"

    forcebal['D_eff_D1']         = {}
    forcebal['D_eff_D1']['data'] = None
    forcebal['D_eff_D1']['unit'] = "m^2/s"
    forcebal['D_eff_D1']['info'] = "Effective Particle Diffusivity for D1"

    forcebal['D_eff_g_D1']         = {}
    forcebal['D_eff_g_D1']['data'] = None
    forcebal['D_eff_g_D1']['unit'] = "m^2/s"
    forcebal['D_eff_g_D1']['info'] = "Effective Particle Diffusivity/g for D1"

    forcebal['D_n_D1']         = {}
    forcebal['D_n_D1']['data'] = None
    forcebal['D_n_D1']['unit'] = "m^2/s"
    forcebal['D_n_D1']['info'] = "Diagonal Particle Diffusivity for D1"

    forcebal['v_n_D1']         = {}
    forcebal['v_n_D1']['data'] = None
    forcebal['v_n_D1']['unit'] = "m/s"
    forcebal['v_n_D1']['info'] = "Particle Pinch from grad(T,n) for D1"

    forcebal['v_wp_D1']         = {}
    forcebal['v_wp_D1']['data'] = None
    forcebal['v_wp_D1']['unit'] = "m/s"
    forcebal['v_wp_D1']['info'] = "Ware Pinch for D1"

    forcebal['v_nwp_D1']         = {}
    forcebal['v_nwp_D1']['data'] = None
    forcebal['v_nwp_D1']['unit'] = "m/s"
    forcebal['v_nwp_D1']['info'] = "Total Pinch for D1"

    forcebal['D_n_g_D1']         = {}
    forcebal['D_n_g_D1']['data'] = None
    forcebal['D_n_g_D1']['unit'] = "m^2/s"
    forcebal['D_n_g_D1']['info'] = "Diagonal Particle Diffusivity/g for D1"

    forcebal['v_n_g_D1']         = {}
    forcebal['v_n_g_D1']['data'] = None
    forcebal['v_n_g_D1']['unit'] = "m/s"
    forcebal['v_n_g_D1']['info'] = "Particle Pinch from grad(T,n)/g for D1"

    forcebal['v_wp_g_D1']         = {}
    forcebal['v_wp_g_D1']['data'] = None
    forcebal['v_wp_g_D1']['unit'] = "m/s"
    forcebal['v_wp_g_D1']['info'] = "Ware Pinch/g for D1"

    forcebal['v_nwp_g_D1']         = {}
    forcebal['v_nwp_g_D1']['data'] = None
    forcebal['v_nwp_g_D1']['unit'] = "m/s"
    forcebal['v_nwp_g_D1']['info'] = "Total Pinch/g for D1"

    forcebal['g_Te_D1']         = {}
    forcebal['g_Te_D1']['data'] = None
    forcebal['g_Te_D1']['unit'] = None
    forcebal['g_Te_D1']['info'] = "Te Exponent in the Density Profile of D1"

    forcebal['g_Ti_D1']         = {}
    forcebal['g_Ti_D1']['data'] = None
    forcebal['g_Ti_D1']['unit'] = None
    forcebal['g_Ti_D1']['info'] = "Te Exponent in the Density Profile of D1"

    forcebal['g_ne1_D1']         = {}
    forcebal['g_ne1_D1']['data'] = None
    forcebal['g_ne1_D1']['unit'] = None
    forcebal['g_ne1_D1']['info'] = "e1 Exponent in the Density Profile of D1"

    forcebal['g_nD1_D1']         = {}
    forcebal['g_nD1_D1']['data'] = None
    forcebal['g_nD1_D1']['unit'] = None
    forcebal['g_nD1_D1']['info'] = "D1 Exponent in the Density Profile of D1"

    forcebal['g_nC6_D1']         = {}
    forcebal['g_nC6_D1']['data'] = None
    forcebal['g_nC6_D1']['unit'] = None
    forcebal['g_nC6_D1']['info'] = "C6 Exponent in the Density Profile of D1"

    forcebal['Sqz_C6']         = {}
    forcebal['Sqz_C6']['data'] = None
    forcebal['Sqz_C6']['unit'] = None
    forcebal['Sqz_C6']['info'] = "Orbit Squeezing Factor for C6"

    forcebal['Gam_C6']         = {}
    forcebal['Gam_C6']['data'] = None
    forcebal['Gam_C6']['unit'] = "#/m^2/s"
    forcebal['Gam_C6']['info'] = "Particle Flux of C6"

    forcebal['D_eff_C6']         = {}
    forcebal['D_eff_C6']['data'] = None
    forcebal['D_eff_C6']['unit'] = "m^2/s"
    forcebal['D_eff_C6']['info'] = "Effective Particle Diffusivity for C6"

    forcebal['D_eff_g_C6']         = {}
    forcebal['D_eff_g_C6']['data'] = None
    forcebal['D_eff_g_C6']['unit'] = "m^2/s"
    forcebal['D_eff_g_C6']['info'] = "Effective Particle Diffusivity/g for C6"

    forcebal['D_n_C6']         = {}
    forcebal['D_n_C6']['data'] = None
    forcebal['D_n_C6']['unit'] = "m^2/s"
    forcebal['D_n_C6']['info'] = "Diagonal Particle Diffusivity for C6"

    forcebal['v_n_C6']         = {}
    forcebal['v_n_C6']['data'] = None
    forcebal['v_n_C6']['unit'] = "m/s"
    forcebal['v_n_C6']['info'] = "Particle Pinch from grad(T,n) for C6"

    forcebal['v_wp_C6']         = {}
    forcebal['v_wp_C6']['data'] = None
    forcebal['v_wp_C6']['unit'] = "m/s"
    forcebal['v_wp_C6']['info'] = "Ware Pinch for C6"

    forcebal['v_nwp_C6']         = {}
    forcebal['v_nwp_C6']['data'] = None
    forcebal['v_nwp_C6']['unit'] = "m/s"
    forcebal['v_nwp_C6']['info'] = "Total Pinch for C6"

    forcebal['D_n_g_C6']         = {}
    forcebal['D_n_g_C6']['data'] = None
    forcebal['D_n_g_C6']['unit'] = "m^2/s"
    forcebal['D_n_g_C6']['info'] = "Diagonal Particle Diffusivity/g for C6"

    forcebal['v_n_g_C6']         = {}
    forcebal['v_n_g_C6']['data'] = None
    forcebal['v_n_g_C6']['unit'] = "m/s"
    forcebal['v_n_g_C6']['info'] = "Particle Pinch from grad(T,n)/g for C6"

    forcebal['v_wp_g_C6']         = {}
    forcebal['v_wp_g_C6']['data'] = None
    forcebal['v_wp_g_C6']['unit'] = "m/s"
    forcebal['v_wp_g_C6']['info'] = "Ware Pinch/g for C6"

    forcebal['v_nwp_g_C6']         = {}
    forcebal['v_nwp_g_C6']['data'] = None
    forcebal['v_nwp_g_C6']['unit'] = "m/s"
    forcebal['v_nwp_g_C6']['info'] = "Total Pinch/g for C6"

    forcebal['g_Te_C6']         = {}
    forcebal['g_Te_C6']['data'] = None
    forcebal['g_Te_C6']['unit'] = None
    forcebal['g_Te_C6']['info'] = "Te Exponent in the Density Profile of C6"

    forcebal['g_Ti_C6']         = {}
    forcebal['g_Ti_C6']['data'] = None
    forcebal['g_Ti_C6']['unit'] = None
    forcebal['g_Ti_C6']['info'] = "Te Exponent in the Density Profile of C6"

    forcebal['g_ne1_C6']         = {}
    forcebal['g_ne1_C6']['data'] = None
    forcebal['g_ne1_C6']['unit'] = None
    forcebal['g_ne1_C6']['info'] = "e1 Exponent in the Density Profile of C6"

    forcebal['g_nD1_C6']         = {}
    forcebal['g_nD1_C6']['data'] = None
    forcebal['g_nD1_C6']['unit'] = None
    forcebal['g_nD1_C6']['info'] = "D1 Exponent in the Density Profile of C6"

    forcebal['g_nC6_C6']         = {}
    forcebal['g_nC6_C6']['data'] = None
    forcebal['g_nC6_C6']['unit'] = None
    forcebal['g_nC6_C6']['info'] = "C6 Exponent in the Density Profile of C6"
    
    forcebal['q_con_e']         = {}
    forcebal['q_con_e']['data'] = None
    forcebal['q_con_e']['unit'] = "W/m^2"
    forcebal['q_con_e']['info'] = "Electron Thermal Conduction Flux"

    forcebal['q_con_i']         = {}
    forcebal['q_con_i']['data'] = None
    forcebal['q_con_i']['unit'] = "W/m^2"
    forcebal['q_con_i']['info'] = "Ion Thermal Conduction Flux"

    forcebal['chi_eff_e']         = {}
    forcebal['chi_eff_e']['data'] = None
    forcebal['chi_eff_e']['unit'] = "m^2/s"
    forcebal['chi_eff_e']['info'] = "Effective Electron Thermal Conductivity"

    forcebal['chi_eff_i']         = {}
    forcebal['chi_eff_i']['data'] = None
    forcebal['chi_eff_i']['unit'] = "m^2/s"
    forcebal['chi_eff_i']['info'] = "Effective Ion Thermal Conductivity"

    forcebal['chi_eff_g_e']         = {}
    forcebal['chi_eff_g_e']['data'] = None
    forcebal['chi_eff_g_e']['unit'] = "m^2/s"
    forcebal['chi_eff_g_e']['info'] = "Effective Electron Thermal Conductivity/g"

    forcebal['chi_eff_g_i']         = {}
    forcebal['chi_eff_g_i']['data'] = None
    forcebal['chi_eff_g_i']['unit'] = "m^2/s"
    forcebal['chi_eff_g_i']['info'] = "Effective Ion Thermal Conductivity/g"

    forcebal['eta_par_r']         = {}
    forcebal['eta_par_r']['data'] = None
    forcebal['eta_par_r']['unit'] = "Ohm.m"
    forcebal['eta_par_r']['info'] = "Parallel Electrical Resistivity"

    forcebal['J']         = {}
    forcebal['J']['data'] = None
    forcebal['J']['unit'] = "A/m^2"
    forcebal['J']['info'] = "Total Parallel Current Density, <J.B>/Bt0"

    forcebal['I']         = {}
    forcebal['I']['data'] = None
    forcebal['I']['unit'] = "A"
    forcebal['I']['info'] = "Enclosed Total Toroidal Current"

    forcebal['J_bs']         = {}
    forcebal['J_bs']['data'] = None
    forcebal['J_bs']['unit'] = "A/m^2"
    forcebal['J_bs']['info'] = "Bootstrap Parallel Current Density, <Jbs.B>/Bt0"

    forcebal['I_bs']         = {}
    forcebal['I_bs']['data'] = None
    forcebal['I_bs']['unit'] = "A"
    forcebal['I_bs']['info'] = "Enclosed Bootstrap Toroidal Current"

    forcebal['E_par']         = {}
    forcebal['E_par']['data'] = None
    forcebal['E_par']['unit'] = "V/m"
    forcebal['E_par']['info'] = "Parallel Electric Field, <E.B/Bt0>"

    forcebal['E_rad_tot']         = {}
    forcebal['E_rad_tot']['data'] = None
    forcebal['E_rad_tot']['unit'] = "V/m"
    forcebal['E_rad_tot']['info'] = "Radial Electric Field, Erho=-d(Phi)/drho/a1, from C6 Force Balance"

    forcebal['E_rad_t']         = {}
    forcebal['E_rad_t']['data'] = None
    forcebal['E_rad_t']['unit'] = "V/m"
    forcebal['E_rad_t']['info'] = "Vtor Contribution to Erho from C6 Force Balance"

    forcebal['E_rad_p']         = {}
    forcebal['E_rad_p']['data'] = None
    forcebal['E_rad_p']['unit'] = "V/m"
    forcebal['E_rad_p']['info'] = "Vpol Contribution to Erho from C6 Force Balance"

    forcebal['E_rad_prs']         = {}
    forcebal['E_rad_prs']['data'] = None
    forcebal['E_rad_prs']['unit'] = "V/m"
    forcebal['E_rad_prs']['info'] = "Diamagnetic Contribution to Erho from C6 Force Balance"

    forcebal['E_rad_tot_o']         = {}
    forcebal['E_rad_tot_o']['data'] = None
    forcebal['E_rad_tot_o']['unit'] = "V/m"
    forcebal['E_rad_tot_o']['info'] = "Er on the Outside Midplane for C6"

    forcebal['E_rad_t_o']         = {}
    forcebal['E_rad_t_o']['data'] = None
    forcebal['E_rad_t_o']['unit'] = "V/m"
    forcebal['E_rad_t_o']['info'] = "Vtor contribution to Er on the outside midplane for C6"

    forcebal['E_rad_p_o']         = {}
    forcebal['E_rad_p_o']['data'] = None
    forcebal['E_rad_p_o']['unit'] = "V/m"
    forcebal['E_rad_p_o']['info'] = "Vpol contribution to Er on the outside midplane for C6"

    forcebal['E_rad_prs_o']         = {}
    forcebal['E_rad_prs_o']['data'] = None
    forcebal['E_rad_prs_o']['unit'] = "V/m"
    forcebal['E_rad_prs_o']['info'] = "Diamagnetic contribution to Er on the outside midplane for C6"

    forcebal['omega_exb_o']         = {}
    forcebal['omega_exb_o']['data'] = None
    forcebal['omega_exb_o']['unit'] = "/s"
    forcebal['omega_exb_o']['info'] = "Hahm-Burell ExB shear damping for C6"

    forcebal['E_rad_m_tot']         = {}
    forcebal['E_rad_m_tot']['data'] = None
    forcebal['E_rad_m_tot']['unit'] = "V/m"
    forcebal['E_rad_m_tot']['info'] = "Mass-density weighted E_rho"

    forcebal['E_rad_m_t']         = {}
    forcebal['E_rad_m_t']['data'] = None
    forcebal['E_rad_m_t']['unit'] = "V/m"
    forcebal['E_rad_m_t']['info'] = "Vtor contribution to mass-density weighted E_rho"

    forcebal['E_rad_m_p']         = {}
    forcebal['E_rad_m_p']['data'] = None
    forcebal['E_rad_m_p']['unit'] = "V/m"
    forcebal['E_rad_m_p']['info'] = "Vpol contribution to mass-density weighted E_rho"

    forcebal['E_rad_m_prs']         = {}
    forcebal['E_rad_m_prs']['data'] = None
    forcebal['E_rad_m_prs']['unit'] = "V/m"
    forcebal['E_rad_m_prs']['info'] = "Diamagnetic contribution to mass-density weighted E_rho"

    forcebal['E_rad_D1_t']         = {}
    forcebal['E_rad_D1_t']['data'] = None
    forcebal['E_rad_D1_t']['unit'] = "V/m"
    forcebal['E_rad_D1_t']['info'] = "Vtor contribution to Erho in the D1 force balance"

    forcebal['E_rad_D1_p']         = {}
    forcebal['E_rad_D1_p']['data'] = None
    forcebal['E_rad_D1_p']['unit'] = "V/m"
    forcebal['E_rad_D1_p']['info'] = "Vpol contribution to Erho in the D1 force balance"

    forcebal['E_rad_D1_prs']         = {}
    forcebal['E_rad_D1_prs']['data'] = None
    forcebal['E_rad_D1_prs']['unit'] = "V/m"
    forcebal['E_rad_D1_prs']['info'] = "Diamagnetic contribution to Erho in the D1 force balance"

    forcebal['E_rad_C6_t']         = {}
    forcebal['E_rad_C6_t']['data'] = None
    forcebal['E_rad_C6_t']['unit'] = "V/m"
    forcebal['E_rad_C6_t']['info'] = "Vtor contribution to Erho in the C6 force balance"

    forcebal['E_rad_C6_p']         = {}
    forcebal['E_rad_C6_p']['data'] = None
    forcebal['E_rad_C6_p']['unit'] = "V/m"
    forcebal['E_rad_C6_p']['info'] = "Vpol contribution to Erho in the C6 force balance"

    forcebal['E_rad_C6_prs']         = {}
    forcebal['E_rad_C6_prs']['data'] = None
    forcebal['E_rad_C6_prs']['unit'] = "V/m"
    forcebal['E_rad_C6_prs']['info'] = "Diamagnetic contribution to Erho in the C6 force balance"

    forcebal['v_t_o_e1']         = {}
    forcebal['v_t_o_e1']['data'] = None
    forcebal['v_t_o_e1']['unit'] = "m/s"
    forcebal['v_t_o_e1']['info'] = "Vtor on the outside midplane for e1"

    forcebal['M_t_o_e1']         = {}
    forcebal['M_t_o_e1']['data'] = None
    forcebal['M_t_o_e1']['unit'] = None
    forcebal['M_t_o_e1']['info'] = "Toroidal Mach Number on the Outside Midplane for e1"

    forcebal['v_p_o_e1']         = {}
    forcebal['v_p_o_e1']['data'] = None
    forcebal['v_p_o_e1']['unit'] = "m/s"
    forcebal['v_p_o_e1']['info'] = "Vpol on the outside midplane for e1"

    forcebal['M_p_o_e1']         = {}
    forcebal['M_p_o_e1']['data'] = None
    forcebal['M_p_o_e1']['unit'] = None
    forcebal['M_p_o_e1']['info'] = "Poloidal Mach Number on the Outside Midplane for e1"

    forcebal['v_par_o_e1']         = {}
    forcebal['v_par_o_e1']['data'] = None
    forcebal['v_par_o_e1']['unit'] = "m/s"
    forcebal['v_par_o_e1']['info'] = "Vpar on the outside midplane for e1"

    forcebal['v_per_o_e1']         = {}
    forcebal['v_per_o_e1']['data'] = None
    forcebal['v_per_o_e1']['unit'] = "m/s"
    forcebal['v_per_o_e1']['info'] = "Vperp on the outside midplane for e1"

    forcebal['v_t_o_D1']         = {}
    forcebal['v_t_o_D1']['data'] = None
    forcebal['v_t_o_D1']['unit'] = "m/s"
    forcebal['v_t_o_D1']['info'] = "Vtor on the outside midplane for D1"

    forcebal['M_t_o_D1']         = {}
    forcebal['M_t_o_D1']['data'] = None
    forcebal['M_t_o_D1']['unit'] = None
    forcebal['M_t_o_D1']['info'] = "Toroidal Mach Number on the outside midplane for D1"

    forcebal['v_p_o_D1']         = {}
    forcebal['v_p_o_D1']['data'] = None
    forcebal['v_p_o_D1']['unit'] = "m/s"
    forcebal['v_p_o_D1']['info'] = "Vpol on the outside midplane for D1"

    forcebal['M_p_o_D1']         = {}
    forcebal['M_p_o_D1']['data'] = None
    forcebal['M_p_o_D1']['unit'] = None
    forcebal['M_p_o_D1']['info'] = "Poloidal Mach Number on the outside midplane for D1"

    forcebal['v_par_o_D1']         = {}
    forcebal['v_par_o_D1']['data'] = None
    forcebal['v_par_o_D1']['unit'] = "m/s"
    forcebal['v_par_o_D1']['info'] = "Vpar on the outside midplane for D1"

    forcebal['v_per_o_D1']         = {}
    forcebal['v_per_o_D1']['data'] = None
    forcebal['v_per_o_D1']['unit'] = "m/s"
    forcebal['v_per_o_D1']['info'] = "Vperp on the outside midplane for D1"

    forcebal['v_t_o_C6']         = {}
    forcebal['v_t_o_C6']['data'] = None
    forcebal['v_t_o_C6']['unit'] = "m/s"
    forcebal['v_t_o_C6']['info'] = "Vtor on the outside midplane for C6"

    forcebal['M_t_o_C6']         = {}
    forcebal['M_t_o_C6']['data'] = None
    forcebal['M_t_o_C6']['unit'] = None
    forcebal['M_t_o_C6']['info'] = "Toroidal Mach Number on the outside midplane for C6"

    forcebal['v_p_o_C6']         = {}
    forcebal['v_p_o_C6']['data'] = None
    forcebal['v_p_o_C6']['unit'] = "m/s"
    forcebal['v_p_o_C6']['info'] = "Vpol on the outside midplane for C6"

    forcebal['M_p_o_C6']         = {}
    forcebal['M_p_o_C6']['data'] = None
    forcebal['M_p_o_C6']['unit'] = None
    forcebal['M_p_o_C6']['info'] = "Poloidal Mach Number on the outside midplane for C6"

    forcebal['v_par_o_C6']         = {}
    forcebal['v_par_o_C6']['data'] = None
    forcebal['v_par_o_C6']['unit'] = "m/s"
    forcebal['v_par_o_C6']['info'] = "Vpar on the outside midplane for C6"

    forcebal['v_per_o_C6']         = {}
    forcebal['v_per_o_C6']['data'] = None
    forcebal['v_per_o_C6']['unit'] = "m/s"
    forcebal['v_per_o_C6']['info'] = "Vperp on the outside midplane for C6"

    forcebal['v_t_o_C6_ex']         = {}
    forcebal['v_t_o_C6_ex']['data'] = None
    forcebal['v_t_o_C6_ex']['unit'] = "m/s"
    forcebal['v_t_o_C6_ex']['info'] = "Experimental Vtor on the outside midplane for C6"

    forcebal['M_t_o_C6_ex']         = {}
    forcebal['M_t_o_C6_ex']['data'] = None
    forcebal['M_t_o_C6_ex']['unit'] = None
    forcebal['M_t_o_C6_ex']['info'] = "Experimental Toroidal Mach Number on the outside midplane for C6"

    forcebal['v_p_o_C6_ex']         = {}
    forcebal['v_p_o_C6_ex']['data'] = None
    forcebal['v_p_o_C6_ex']['unit'] = "m/s"
    forcebal['v_p_o_C6_ex']['info'] = "Experimental Vpol on the outside midplane for C6"

    forcebal['k_p_o_C6_ex']         = {}
    forcebal['k_p_o_C6_ex']['data'] = None
    forcebal['k_p_o_C6_ex']['unit'] = "m/s"
    forcebal['k_p_o_C6_ex']['info'] = "Experimental k flux surface average for C6"

    forcebal['Omega_e1']         = {}
    forcebal['Omega_e1']['data'] = None
    forcebal['Omega_e1']['unit'] = "rad/s"
    forcebal['Omega_e1']['info'] = "Toroidal rotation frequency of e1"

    forcebal['u_par_e1']         = {}
    forcebal['u_par_e1']['data'] = None
    forcebal['u_par_e1']['unit'] = "T.m/s"
    forcebal['u_par_e1']['info'] = "Parallel velocity parameter <v.B> for e1"

    forcebal['u_p_e1']         = {}
    forcebal['u_p_e1']['data'] = None
    forcebal['u_p_e1']['unit'] = "m/T/s"
    forcebal['u_p_e1']['info'] = "Poloidal velocity parameter <v.theta>/<B.theta> for e1"

    forcebal['Omega_D1']         = {}
    forcebal['Omega_D1']['data'] = None
    forcebal['Omega_D1']['unit'] = "rad/s"
    forcebal['Omega_D1']['info'] = "Toroidal rotation frequency of D1"

    forcebal['u_par_D1']         = {}
    forcebal['u_par_D1']['data'] = None
    forcebal['u_par_D1']['unit'] = "T.m/s"
    forcebal['u_par_D1']['info'] = "Parallel velocity parameter <v.B> for D1"

    forcebal['u_p_D1']         = {}
    forcebal['u_p_D1']['data'] = None
    forcebal['u_p_D1']['unit'] = "m/T/s"
    forcebal['u_p_D1']['info'] = "Poloidal velocity parameter <v.theta>/<B.theta> for D1"

    forcebal['Omega_C6']         = {}
    forcebal['Omega_C6']['data'] = None
    forcebal['Omega_C6']['unit'] = "rad/s"
    forcebal['Omega_C6']['info'] = "Toroidal rotation frequency of C6"

    forcebal['u_par_C6']         = {}
    forcebal['u_par_C6']['data'] = None
    forcebal['u_par_C6']['unit'] = "T.m/s"
    forcebal['u_par_C6']['info'] = "Parallel velocity parameter <v.B> for C6"

    forcebal['u_p_C6']         = {}
    forcebal['u_p_C6']['data'] = None
    forcebal['u_p_C6']['unit'] = "m/T/s"
    forcebal['u_p_C6']['info'] = "Poloidal velocity parameter <v.theta>/<B.theta> for C6"

    forcebal['Omega_C6_ex']         = {}
    forcebal['Omega_C6_ex']['data'] = None
    forcebal['Omega_C6_ex']['unit'] = "rad/s"
    forcebal['Omega_C6_ex']['info'] = "Experimental Toroidal rotation frequency of C6"

    forcebal['dummy']         = {}
    forcebal['dummy']['data'] = None
    forcebal['dummy']['unit'] = None
    forcebal['dummy']['info'] = None

    return forcebal

def get_forcebal_extend_vars():
    #>============================================
    #>Notes on meaning of subscripts:
    #>============================================
    #_r             Radial grid (all variables presently on same grid)
    #_p             Poloidal
    #_t             Toroidal
    #_in            INside of torus in axial plane
    #_out           OUTside of torus in axial plane
    #_par           PARallel to B
    #_e             Electron value
    #_i             Ion (main) value
    #_im            IMpurity value
    #_con           CONduction
    #_bs            BootStrap
    #_wp            Ware Pinch
    #_g             <Grad(r_t)**2> normalization
    #_n             deNsity term

    forcebal = {}

    #>============================================
    #>============================================
    #>Radial grids:
    #>============================================

    forcebal['rho_p_r']         = {}
    forcebal['rho_p_r']['data'] = None
    forcebal['rho_p_r']['unit'] = None
    forcebal['rho_p_r']['info'] = "Normalized radial grid prop to poloidal flux"

    forcebal['rho_t_r']         = {}
    forcebal['rho_t_r']['data'] = None
    forcebal['rho_t_r']['unit'] = None
    forcebal['rho_t_r']['info'] = "Normalized radial grid prop to sqrt(toroidal flux)"

    forcebal['r_p_r']         = {}
    forcebal['r_p_r']['data'] = None
    forcebal['r_p_r']['unit'] = "m"
    forcebal['r_p_r']['info'] = "rho_p_r*a0"

    forcebal['r_t_r']         = {}
    forcebal['r_t_r']['data'] = None
    forcebal['r_t_r']['unit'] = "m"
    forcebal['r_t_r']['info'] = "rho_t_r*a0"

    forcebal['R_in_r']         = {}
    forcebal['R_in_r']['data'] = None
    forcebal['R_in_r']['unit'] = "m"
    forcebal['R_in_r']['info'] = "R of intersection of surface with axis plane on inside"

    forcebal['R_o_r']         = {}
    forcebal['R_o_r']['data'] = None
    forcebal['R_o_r']['unit'] = "m"
    forcebal['R_o_r']['info'] = "R of intersection of surface with axis plane on outside"

    #>============================================
    #>Magnetic fields and fluxes
    #>============================================

    forcebal['B_p_o_r']         = {}
    forcebal['B_p_o_r']['data'] = None
    forcebal['B_p_o_r']['unit'] = "T"
    forcebal['B_p_o_r']['info'] = "B_p at R_o_r"

    forcebal['B_t_o_r']         = {}
    forcebal['B_t_o_r']['data'] = None
    forcebal['B_t_o_r']['unit'] = "T"
    forcebal['B_t_o_r']['info'] = "B_t at R_o_r"

    forcebal['Phi_t_r']         = {}
    forcebal['Phi_t_r']['data'] = None
    forcebal['Phi_t_r']['unit'] = "Wb"
    forcebal['Phi_t_r']['info'] = "Toroidal flux"

    forcebal['Psi_r']         = {}
    forcebal['Psi_r']['data'] = None
    forcebal['Psi_r']['unit'] = "Wb/rad"
    forcebal['Psi_r']['info'] = "Poloidal flux/2*pi"

    forcebal['q_r']         = {}
    forcebal['q_r']['data'] = None
    forcebal['q_r']['unit'] = None
    forcebal['q_r']['info'] = "Safety Factor"


    #>============================================
    #>Flux functions and metrics
    #>============================================

    forcebal['F_r']         = {}
    forcebal['F_r']['data'] = None
    forcebal['F_r']['unit'] = "m*T"
    forcebal['F_r']['info'] = "R*B_t"

    forcebal['f_trap_r']         = {}
    forcebal['f_trap_r']['data'] = None
    forcebal['f_trap_r']['unit'] = None
    forcebal['f_trap_r']['info'] = "Trapped Fraction"

    forcebal['Sqz_i_r']         = {}
    forcebal['Sqz_i_r']['data'] = None
    forcebal['Sqz_i_r']['unit'] = None
    forcebal['Sqz_i_r']['info'] = "Orbit squeezing factor for i"

    forcebal['Sqz_im_r']         = {}
    forcebal['Sqz_im_r']['data'] = None
    forcebal['Sqz_im_r']['unit'] = None
    forcebal['Sqz_im_r']['info'] = "Orbit squeezing factor for im"

    forcebal['Sqz_D_r']         = {}
    forcebal['Sqz_D_r']['data'] = None
    forcebal['Sqz_D_r']['unit'] = None
    forcebal['Sqz_D_r']['info'] = "Orbit squeezing factor for i"

    forcebal['Sqz_N_r']         = {}
    forcebal['Sqz_N_r']['data'] = None
    forcebal['Sqz_N_r']['unit'] = None
    forcebal['Sqz_N_r']['info'] = "Orbit squeezing factor for i"

    forcebal['Sqz_C_r']         = {}
    forcebal['Sqz_C_r']['data'] = None
    forcebal['Sqz_C_r']['unit'] = None
    forcebal['Sqz_C_r']['info'] = "Orbit squeezing factor for i"

    forcebal['grad(r_t)**2_r']         = {}
    forcebal['grad(r_t)**2_r']['data'] = None
    forcebal['grad(r_t)**2_r']['unit'] = None
    forcebal['grad(r_t)**2_r']['info'] = "<grad(rho_t)**2>*a_0**2"

    forcebal['dV/d(r_t)_r']         = {}
    forcebal['dV/d(r_t)_r']['data'] = None
    forcebal['dV/d(r_t)_r']['unit'] = None
    forcebal['dV/d(r_t)_r']['info'] = ""

    forcebal['Vol_r']         = {}
    forcebal['Vol_r']['data'] = None
    forcebal['Vol_r']['unit'] = None
    forcebal['Vol_r']['info'] = ""


    #>============================================
    #>Plasma profiles
    #>============================================

    forcebal['den_e_r']         = {}
    forcebal['den_e_r']['data'] = None
    forcebal['den_e_r']['unit'] = "/m^3"
    forcebal['den_e_r']['info'] = "ne"

    forcebal['den_i_r']         = {}
    forcebal['den_i_r']['data'] = None
    forcebal['den_i_r']['unit'] = "/m^3"
    forcebal['den_i_r']['info'] = "Main ion density - deduced from n_e and Z_eff"

    forcebal['den_im_r']         = {}
    forcebal['den_im_r']['data'] = None
    forcebal['den_im_r']['unit'] = "/m^3"
    forcebal['den_im_r']['info'] = "Impurity density - deduced from n_e and Z_eff"

    forcebal['den_D_r']         = {}
    forcebal['den_D_r']['data'] = None
    forcebal['den_D_r']['unit'] = "/m^3"
    forcebal['den_D_r']['info'] = ""

    forcebal['den_N_r']         = {}
    forcebal['den_N_r']['data'] = None
    forcebal['den_N_r']['unit'] = "/m^3"
    forcebal['den_N_r']['info'] = ""

    forcebal['den_C_r']         = {}
    forcebal['den_C_r']['data'] = None
    forcebal['den_C_r']['unit'] = "/m^3"
    forcebal['den_C_r']['info'] = ""

    forcebal['Omega_e_r']         = {}
    forcebal['Omega_e_r']['data'] = None
    forcebal['Omega_e_r']['unit'] = None
    forcebal['Omega_e_r']['info'] = "Toroidal rotation frequency of impurity"

    forcebal['Omega_im_r']         = {}
    forcebal['Omega_im_r']['data'] = None
    forcebal['Omega_im_r']['unit'] = None
    forcebal['Omega_im_r']['info'] = "Toroidal rotation frequency of impurity"

    forcebal['Omega_D_r']         = {}
    forcebal['Omega_D_r']['data'] = None
    forcebal['Omega_D_r']['unit'] = None
    forcebal['Omega_D_r']['info'] = "Toroidal rotation frequency of impurity"

    forcebal['Omega_N_r']         = {}
    forcebal['Omega_N_r']['data'] = None
    forcebal['Omega_N_r']['unit'] = None
    forcebal['Omega_N_r']['info'] = "Toroidal rotation frequency of impurity"

    forcebal['Omega_C_r']         = {}
    forcebal['Omega_C_r']['data'] = None
    forcebal['Omega_C_r']['unit'] = None
    forcebal['Omega_C_r']['info'] = "Toroidal rotation frequency of impurity"

    forcebal['Omega_C_ex_r']         = {}
    forcebal['Omega_C_ex_r']['data'] = None
    forcebal['Omega_C_ex_r']['unit'] = None
    forcebal['Omega_C_ex_r']['info'] = "Toroidal rotation frequency of impurity"

    forcebal['T_e_r']         = {}
    forcebal['T_e_r']['data'] = None
    forcebal['T_e_r']['unit'] = "keV"
    forcebal['T_e_r']['info'] = "Electron Temperature"

    forcebal['T_i_r']         = {}
    forcebal['T_i_r']['data'] = None
    forcebal['T_i_r']['unit'] = "keV"
    forcebal['T_i_r']['info'] = "Ions Temperature"

    forcebal['Z_eff_r']         = {}
    forcebal['Z_eff_r']['data'] = None
    forcebal['Z_eff_r']['unit'] = None
    forcebal['Z_eff_r']['info'] = "Effective Atomic Number"

    forcebal['Z_eff_ex_r']         = {}
    forcebal['Z_eff_ex_r']['data'] = None
    forcebal['Z_eff_ex_r']['unit'] = None
    forcebal['Z_eff_ex_r']['info'] = "Effective Atomic Number"


    #>============================================
    #>Particle fluxes, diffusivities and pinches
    #>============================================

    forcebal['Gam_e_r']         = {}
    forcebal['Gam_e_r']['data'] = None
    forcebal['Gam_e_r']['unit'] = "/m^2/s"
    forcebal['Gam_e_r']['info'] = "Net nc radial e flux=<Gam.grad(r_t)>"

    forcebal['Gam_i_r']         = {}
    forcebal['Gam_i_r']['data'] = None
    forcebal['Gam_i_r']['unit'] = "/m^2/s"
    forcebal['Gam_i_r']['info'] = "Net nc radial i flux=<Gam.grad(r_t)>"

    forcebal['Gam_im_r']         = {}
    forcebal['Gam_im_r']['data'] = None
    forcebal['Gam_im_r']['unit'] = "/m^2/s"
    forcebal['Gam_im_r']['info'] = "Net nc radial im flux=<Gam.grad(r_t)>"

    forcebal['Gam_D_r']         = {}
    forcebal['Gam_D_r']['data'] = None
    forcebal['Gam_D_r']['unit'] = "/m^2/s"
    forcebal['Gam_D_r']['info'] = "Net nc radial i flux=<Gam.grad(r_t)>"

    forcebal['Gam_N_r']         = {}
    forcebal['Gam_N_r']['data'] = None
    forcebal['Gam_N_r']['unit'] = "/m^2/s"
    forcebal['Gam_N_r']['info'] = "Net nc radial i flux=<Gam.grad(r_t)>"

    forcebal['Gam_C_r']         = {}
    forcebal['Gam_C_r']['data'] = None
    forcebal['Gam_C_r']['unit'] = "/m^2/s"
    forcebal['Gam_C_r']['info'] = "Net nc radial i flux=<Gam.grad(r_t)>"

    forcebal['D_eff_e_r']         = {}
    forcebal['D_eff_e_r']['data'] = None
    forcebal['D_eff_e_r']['unit'] = "m^2/s"
    forcebal['D_eff_e_r']['info'] = "Eff e diffusivity=-Gam_e/(dn/(dr_t))"

    forcebal['D_eff_i_r']         = {}
    forcebal['D_eff_i_r']['data'] = None
    forcebal['D_eff_i_r']['unit'] = "m^2/s"
    forcebal['D_eff_i_r']['info'] = "Eff i diffusivity=-Gam_i/(dn/(dr_t))"

    forcebal['D_eff_im_r']         = {}
    forcebal['D_eff_im_r']['data'] = None
    forcebal['D_eff_im_r']['unit'] = "m^2/s"
    forcebal['D_eff_im_r']['info'] = "Eff im diffusivity=-Gam_im/(dn/(dr_t))"

    forcebal['D_eff_D_r']         = {}
    forcebal['D_eff_D_r']['data'] = None
    forcebal['D_eff_D_r']['unit'] = "m^2/s"
    forcebal['D_eff_D_r']['info'] = ""

    forcebal['D_eff_N_r']         = {}
    forcebal['D_eff_N_r']['data'] = None
    forcebal['D_eff_N_r']['unit'] = "m^2/s"
    forcebal['D_eff_N_r']['info'] = ""

    forcebal['D_eff_C_r']         = {}
    forcebal['D_eff_C_r']['data'] = None
    forcebal['D_eff_C_r']['unit'] = "m^2/s"
    forcebal['D_eff_C_r']['info'] = ""

    forcebal['D_eff_g_e_r']         = {}
    forcebal['D_eff_g_e_r']['data'] = None
    forcebal['D_eff_g_e_r']['unit'] = "m^2/s"
    forcebal['D_eff_g_e_r']['info'] = "Eff e diffusivity=D_eff_e_r/<grad(r_t)^2>"

    forcebal['D_eff_g_i_r']         = {}
    forcebal['D_eff_g_i_r']['data'] = None
    forcebal['D_eff_g_i_r']['unit'] = "m^2/s"
    forcebal['D_eff_g_i_r']['info'] = "Eff i diffusivity=D_eff_i_r/<grad(r_t)^2>"

    forcebal['D_eff_g_im_r']         = {}
    forcebal['D_eff_g_im_r']['data'] = None
    forcebal['D_eff_g_im_r']['unit'] = "m^2/s"
    forcebal['D_eff_g_im_r']['info'] = "Eff im diffusivity=D_eff_im_r/<grad(r_t)^2>"

    forcebal['D_eff_g_D_r']         = {}
    forcebal['D_eff_g_D_r']['data'] = None
    forcebal['D_eff_g_D_r']['unit'] = "m^2/s"
    forcebal['D_eff_g_D_r']['info'] = ""

    forcebal['D_eff_g_N_r']         = {}
    forcebal['D_eff_g_N_r']['data'] = None
    forcebal['D_eff_g_N_r']['unit'] = "m^2/s"
    forcebal['D_eff_g_N_r']['info'] = ""

    forcebal['D_eff_g_C_r']         = {}
    forcebal['D_eff_g_C_r']['data'] = None
    forcebal['D_eff_g_C_r']['unit'] = "m^2/s"
    forcebal['D_eff_g_C_r']['info'] = ""

    forcebal['D_n_e_r']         = {}
    forcebal['D_n_e_r']['data'] = None
    forcebal['D_n_e_r']['unit'] = "m^2/s"
    forcebal['D_n_e_r']['info'] = "e diagonal diffusivity"

    forcebal['D_n_i_r']         = {}
    forcebal['D_n_i_r']['data'] = None
    forcebal['D_n_i_r']['unit'] = "m^2/s"
    forcebal['D_n_i_r']['info'] = "i diagonal diffusivity"

    forcebal['D_n_im_r']         = {}
    forcebal['D_n_im_r']['data'] = None
    forcebal['D_n_im_r']['unit'] = "m^2/s"
    forcebal['D_n_im_r']['info'] = "im diagonal diffusivity"

    forcebal['D_n_D_r']         = {}
    forcebal['D_n_D_r']['data'] = None
    forcebal['D_n_D_r']['unit'] = "m^2/s"
    forcebal['D_n_D_r']['info'] = ""

    forcebal['D_n_N_r']         = {}
    forcebal['D_n_N_r']['data'] = None
    forcebal['D_n_N_r']['unit'] = "m^2/s"
    forcebal['D_n_N_r']['info'] = ""

    forcebal['D_n_C_r']         = {}
    forcebal['D_n_C_r']['data'] = None
    forcebal['D_n_C_r']['unit'] = "m^2/s"
    forcebal['D_n_C_r']['info'] = ""

    forcebal['v_n_e_r']         = {}
    forcebal['v_n_e_r']['data'] = None
    forcebal['v_n_e_r']['unit'] = "m/s"
    forcebal['v_n_e_r']['info'] = "e nondiffusive (non-WP) velocity"

    forcebal['v_n_i_r']         = {}
    forcebal['v_n_i_r']['data'] = None
    forcebal['v_n_i_r']['unit'] = "m/s"
    forcebal['v_n_i_r']['info'] = "i nondiffusive (non-WP) velocity"

    forcebal['v_n_im_r']         = {}
    forcebal['v_n_im_r']['data'] = None
    forcebal['v_n_im_r']['unit'] = "m/s"
    forcebal['v_n_im_r']['info'] = "im nondiffusive (non-WP) velocity"

    forcebal['v_n_D_r']         = {}
    forcebal['v_n_D_r']['data'] = None
    forcebal['v_n_D_r']['unit'] = "m/s"
    forcebal['v_n_D_r']['info'] = ""

    forcebal['v_n_N_r']         = {}
    forcebal['v_n_N_r']['data'] = None
    forcebal['v_n_N_r']['unit'] = "m/s"
    forcebal['v_n_N_r']['info'] = ""

    forcebal['v_n_C_r']         = {}
    forcebal['v_n_C_r']['data'] = None
    forcebal['v_n_C_r']['unit'] = "m/s"
    forcebal['v_n_C_r']['info'] = ""

    forcebal['v_wp_e_r']         = {}
    forcebal['v_wp_e_r']['data'] = None
    forcebal['v_wp_e_r']['unit'] = "m/s"
    forcebal['v_wp_e_r']['info'] = "e WP velocity"

    forcebal['v_wp_i_r']         = {}
    forcebal['v_wp_i_r']['data'] = None
    forcebal['v_wp_i_r']['unit'] = "m/s"
    forcebal['v_wp_i_r']['info'] = "i WP velocity"

    forcebal['v_wp_im_r']         = {}
    forcebal['v_wp_im_r']['data'] = None
    forcebal['v_wp_im_r']['unit'] = "m/s"
    forcebal['v_wp_im_r']['info'] = "im WP velocity"

    forcebal['v_wp_D_r']         = {}
    forcebal['v_wp_D_r']['data'] = None
    forcebal['v_wp_D_r']['unit'] = "m/s"
    forcebal['v_wp_D_r']['info'] = ""

    forcebal['v_wp_N_r']         = {}
    forcebal['v_wp_N_r']['data'] = None
    forcebal['v_wp_N_r']['unit'] = "m/s"
    forcebal['v_wp_N_r']['info'] = ""

    forcebal['v_wp_C_r']         = {}
    forcebal['v_wp_C_r']['data'] = None
    forcebal['v_wp_C_r']['unit'] = "m/s"
    forcebal['v_wp_C_r']['info'] = ""

    forcebal['D_n_g_e_r']         = {}
    forcebal['D_n_g_e_r']['data'] = None
    forcebal['D_n_g_e_r']['unit'] = "m^2/s"
    forcebal['D_n_g_e_r']['info'] = "D_n_e_r/<grad(r_t)^2>"

    forcebal['D_n_g_i_r']         = {}
    forcebal['D_n_g_i_r']['data'] = None
    forcebal['D_n_g_i_r']['unit'] = "m^2/s"
    forcebal['D_n_g_i_r']['info'] = "D_n_i_r/<grad(r_t)^2>"

    forcebal['D_n_g_im_r']         = {}
    forcebal['D_n_g_im_r']['data'] = None
    forcebal['D_n_g_im_r']['unit'] = "m^2/s"
    forcebal['D_n_g_im_r']['info'] = "D_n_im_r/<grad(r_t)^2>"

    forcebal['D_n_g_D_r']         = {}
    forcebal['D_n_g_D_r']['data'] = None
    forcebal['D_n_g_D_r']['unit'] = "m^2/s"
    forcebal['D_n_g_D_r']['info'] = ""

    forcebal['D_n_g_N_r']         = {}
    forcebal['D_n_g_N_r']['data'] = None
    forcebal['D_n_g_N_r']['unit'] = "m^2/s"
    forcebal['D_n_g_N_r']['info'] = ""

    forcebal['D_n_g_C_r']         = {}
    forcebal['D_n_g_C_r']['data'] = None
    forcebal['D_n_g_C_r']['unit'] = "m^2/s"
    forcebal['D_n_g_C_r']['info'] = ""

    forcebal['v_n_g_e_r']         = {}
    forcebal['v_n_g_e_r']['data'] = None
    forcebal['v_n_g_e_r']['unit'] = "m/s"
    forcebal['v_n_g_e_r']['info'] = "v_n_e_r/(<grad(r_t)^2>)**0.5"

    forcebal['v_n_g_i_r']         = {}
    forcebal['v_n_g_i_r']['data'] = None
    forcebal['v_n_g_i_r']['unit'] = "m/s"
    forcebal['v_n_g_i_r']['info'] = "v_n_i_r/(<grad(r_t)^2>)**0.5"

    forcebal['v_n_g_im_r']         = {}
    forcebal['v_n_g_im_r']['data'] = None
    forcebal['v_n_g_im_r']['unit'] = "m/s"
    forcebal['v_n_g_im_r']['info'] = "v_n_im_r/(<grad(r_t)^2>)**0.5"

    forcebal['v_n_g_D_r']         = {}
    forcebal['v_n_g_D_r']['data'] = None
    forcebal['v_n_g_D_r']['unit'] = "m/s"
    forcebal['v_n_g_D_r']['info'] = ""

    forcebal['v_n_g_N_r']         = {}
    forcebal['v_n_g_N_r']['data'] = None
    forcebal['v_n_g_N_r']['unit'] = "m/s"
    forcebal['v_n_g_N_r']['info'] = ""

    forcebal['v_n_g_C_r']         = {}
    forcebal['v_n_g_C_r']['data'] = None
    forcebal['v_n_g_C_r']['unit'] = "m/s"
    forcebal['v_n_g_C_r']['info'] = ""

    forcebal['v_wp_g_e_r']         = {}
    forcebal['v_wp_g_e_r']['data'] = None
    forcebal['v_wp_g_e_r']['unit'] = "m/s"
    forcebal['v_wp_g_e_r']['info'] = "v_wp_e_r/(<grad(r_t)^2>)**0.5"

    forcebal['v_wp_g_i_r']         = {}
    forcebal['v_wp_g_i_r']['data'] = None
    forcebal['v_wp_g_i_r']['unit'] = "m/s"
    forcebal['v_wp_g_i_r']['info'] = "v_wp_e_r/(<grad(r_t)^2>)**0.5"

    forcebal['v_wp_g_im_r']         = {}
    forcebal['v_wp_g_im_r']['data'] = None
    forcebal['v_wp_g_im_r']['unit'] = "m/s"
    forcebal['v_wp_g_im_r']['info'] = "v_wp_e_r/(<grad(r_t)^2>)**0.5"

    forcebal['v_wp_g_D_r']         = {}
    forcebal['v_wp_g_D_r']['data'] = None
    forcebal['v_wp_g_D_r']['unit'] = "m/s"
    forcebal['v_wp_g_D_r']['info'] = ""

    forcebal['v_wp_g_N_r']         = {}
    forcebal['v_wp_g_N_r']['data'] = None
    forcebal['v_wp_g_N_r']['unit'] = "m/s"
    forcebal['v_wp_g_N_r']['info'] = ""

    forcebal['v_wp_g_C_r']         = {}
    forcebal['v_wp_g_C_r']['data'] = None
    forcebal['v_wp_g_C_r']['unit'] = "m/s"
    forcebal['v_wp_g_C_r']['info'] = ""


    #>============================================
    #>Conduction heat fluxes and conductivities
    #>============================================

    forcebal['q_con_e_r']         = {}
    forcebal['q_con_e_r']['data'] = None
    forcebal['q_con_e_r']['unit'] = "W/m^2"
    forcebal['q_con_e_r']['info'] = "Net nc e cond flux=<qcon.grad(rt)>"

    forcebal['q_con_i_r']         = {}
    forcebal['q_con_i_r']['data'] = None
    forcebal['q_con_i_r']['unit'] = "W/m^2"
    forcebal['q_con_i_r']['info'] = "Net nc i+im cond flux=<qcon.grad(rt)>"

    forcebal['chi_eff_e_r']         = {}
    forcebal['chi_eff_e_r']['data'] = None
    forcebal['chi_eff_e_r']['unit'] = "m^2/s"
    forcebal['chi_eff_e_r']['info'] = "Eff e cond=-qcone/(ne*dkTe/(drt))"

    forcebal['chi_eff_i_r']         = {}
    forcebal['chi_eff_i_r']['data'] = None
    forcebal['chi_eff_i_r']['unit'] = "m^2/s"
    forcebal['chi_eff_i_r']['info'] = "Eff i+im cond=-qconi/((ni+nim)*dkTi/(drt))"

    forcebal['chi_eff_g_e_r']         = {}
    forcebal['chi_eff_g_e_r']['data'] = None
    forcebal['chi_eff_g_e_r']['unit'] = "m^2/s"
    forcebal['chi_eff_g_e_r']['info'] = "Eff e cond=chieffe_r/<grad(rt)^2>"

    forcebal['chi_eff_g_i_r']         = {}
    forcebal['chi_eff_g_i_r']['data'] = None
    forcebal['chi_eff_g_i_r']['unit'] = "m^2/s"
    forcebal['chi_eff_g_i_r']['info'] = "Eff i+im cond=chieffi_r/<grad(rt)^2>"


    #>============================================
    #>Resistivity and currents
    #>============================================

    forcebal['eta_par_r']         = {}
    forcebal['eta_par_r']['data'] = None
    forcebal['eta_par_r']['unit'] = "Ohm*m"
    forcebal['eta_par_r']['info'] = "Parallel electrical resistivity"

    forcebal['I_r']         = {}
    forcebal['I_r']['data'] = None
    forcebal['I_r']['unit'] = "A/m^2"
    forcebal['I_r']['info'] = ""

    forcebal['J_r']         = {}
    forcebal['J_r']['data'] = None
    forcebal['J_r']['unit'] = "A/m^2"
    forcebal['J_r']['info'] = ""

    forcebal['I_bs_r']         = {}
    forcebal['I_bs_r']['data'] = None
    forcebal['I_bs_r']['unit'] = "A"
    forcebal['I_bs_r']['info'] = "Bootstrap current"

    forcebal['J_bs_r']         = {}
    forcebal['J_bs_r']['data'] = None
    forcebal['J_bs_r']['unit'] = "A/m^2"
    forcebal['J_bs_r']['info'] = "Bootstrap current density=<Jbc.B>/Bt0"


    #>============================================
    #>Radial electric field
    #>============================================

    forcebal['E_par_r']         = {}
    forcebal['E_par_r']['data'] = None
    forcebal['E_par_r']['unit'] = "V/m"
    forcebal['E_par_r']['info'] = ""

    forcebal['E_r_tot_r']         = {}
    forcebal['E_r_tot_r']['data'] = None
    forcebal['E_r_tot_r']['unit'] = "V/m"
    forcebal['E_r_tot_r']['info'] = "Total radial electric field=-dPhi/(dr_t)"

    forcebal['E_r_tor_r']         = {}
    forcebal['E_r_tor_r']['data'] = None
    forcebal['E_r_tor_r']['unit'] = "V/m"
    forcebal['E_r_tor_r']['info'] = "E_r component from im toroidal rotation"

    forcebal['E_r_pol_r']         = {}
    forcebal['E_r_pol_r']['data'] = None
    forcebal['E_r_pol_r']['unit'] = "V/m"
    forcebal['E_r_pol_r']['info'] = "E_r component from im poloidal rotation"

    forcebal['E_r_prs_r']         = {}
    forcebal['E_r_prs_r']['data'] = None
    forcebal['E_r_prs_r']['unit'] = "V/m"
    forcebal['E_r_prs_r']['info'] = "E_r component from im pressure gradient"

    forcebal['E_rad_tot_r']         = {}
    forcebal['E_rad_tot_r']['data'] = None
    forcebal['E_rad_tot_r']['unit'] = "V/m"
    forcebal['E_rad_tot_r']['info'] = ""

    forcebal['E_rad_t_r']         = {}
    forcebal['E_rad_t_r']['data'] = None
    forcebal['E_rad_t_r']['unit'] = "V/m"
    forcebal['E_rad_t_r']['info'] = ""

    forcebal['E_rad_p_r']         = {}
    forcebal['E_rad_p_r']['data'] = None
    forcebal['E_rad_p_r']['unit'] = "V/m"
    forcebal['E_rad_p_r']['info'] = ""

    forcebal['E_rad_prs_r']         = {}
    forcebal['E_rad_prs_r']['data'] = None
    forcebal['E_rad_prs_r']['unit'] = "V/m"
    forcebal['E_rad_prs_r']['info'] = ""

    forcebal['E_rad_ex_r']         = {}
    forcebal['E_rad_ex_r']['data'] = None
    forcebal['E_rad_ex_r']['unit'] = "V/m"
    forcebal['E_rad_ex_r']['info'] = ""

    forcebal['E_rad_ex_o_r']         = {}
    forcebal['E_rad_ex_o_r']['data'] = None
    forcebal['E_rad_ex_o_r']['unit'] = "V/m"
    forcebal['E_rad_ex_o_r']['info'] = ""

    forcebal['E_rad_t_o_r']         = {}
    forcebal['E_rad_t_o_r']['data'] = None
    forcebal['E_rad_t_o_r']['unit'] = "V/m"
    forcebal['E_rad_t_o_r']['info'] = ""

    forcebal['E_rad_p_o_r']         = {}
    forcebal['E_rad_p_o_r']['data'] = None
    forcebal['E_rad_p_o_r']['unit'] = "V/m"
    forcebal['E_rad_p_o_r']['info'] = ""

    forcebal['E_rad_prs_o_r']         = {}
    forcebal['E_rad_prs_o_r']['data'] = None
    forcebal['E_rad_prs_o_r']['unit'] = "V/m"
    forcebal['E_rad_prs_o_r']['info'] = ""

    forcebal['E_rad_tot_o_r']         = {}
    forcebal['E_rad_tot_o_r']['data'] = None
    forcebal['E_rad_tot_o_r']['unit'] = "V/m"
    forcebal['E_rad_tot_o_r']['info'] = ""

    forcebal['prs_e_r']         = {}
    forcebal['prs_e_r']['data'] = None
    forcebal['prs_e_r']['unit'] = ""
    forcebal['prs_e_r']['info'] = ""

    forcebal['prs_D_r']         = {}
    forcebal['prs_D_r']['data'] = None
    forcebal['prs_D_r']['unit'] = ""
    forcebal['prs_D_r']['info'] = ""

    forcebal['prs_N_r']         = {}
    forcebal['prs_N_r']['data'] = None
    forcebal['prs_N_r']['unit'] = ""
    forcebal['prs_N_r']['info'] = ""

    forcebal['prs_C_r']         = {}
    forcebal['prs_C_r']['data'] = None
    forcebal['prs_C_r']['unit'] = ""
    forcebal['prs_C_r']['info'] = ""

    forcebal['omega_exb_o_r']         = {}
    forcebal['omega_exb_o_r']['data'] = None
    forcebal['omega_exb_o_r']['unit'] = "V/m"
    forcebal['omega_exb_o_r']['info'] = ""

    #>============================================
    #>Particle flow velocities
    #>============================================

    forcebal['v_t_o_e_r']         = {}
    forcebal['v_t_o_e_r']['data'] = None
    forcebal['v_t_o_e_r']['unit'] = "m/s"
    forcebal['v_t_o_e_r']['info'] = "Toroidal electron velocity"

    forcebal['v_t_o_i_r']         = {}
    forcebal['v_t_o_i_r']['data'] = None
    forcebal['v_t_o_i_r']['unit'] = "m/s"
    forcebal['v_t_o_i_r']['info'] = "Toroidal ion velocity"

    forcebal['v_t_o_im_r']         = {}
    forcebal['v_t_o_im_r']['data'] = None
    forcebal['v_t_o_im_r']['unit'] = "m/s"
    forcebal['v_t_o_im_r']['info'] = "Toroidal impurity velocity"

    forcebal['v_t_o_D_r']         = {}
    forcebal['v_t_o_D_r']['data'] = None
    forcebal['v_t_o_D_r']['unit'] = "m/s"
    forcebal['v_t_o_D_r']['info'] = ""

    forcebal['v_t_o_N_r']         = {}
    forcebal['v_t_o_N_r']['data'] = None
    forcebal['v_t_o_N_r']['unit'] = "m/s"
    forcebal['v_t_o_N_r']['info'] = ""

    forcebal['v_t_o_C_r']         = {}
    forcebal['v_t_o_C_r']['data'] = None
    forcebal['v_t_o_C_r']['unit'] = "m/s"
    forcebal['v_t_o_C_r']['info'] = ""

    forcebal['v_t_o_C_ex_r']         = {}
    forcebal['v_t_o_C_ex_r']['data'] = None
    forcebal['v_t_o_C_ex_r']['unit'] = "m/s"
    forcebal['v_t_o_C_ex_r']['info'] = ""

    forcebal['v_t_o_ex_r']         = {}
    forcebal['v_t_o_ex_r']['data'] = None
    forcebal['v_t_o_ex_r']['unit'] = "m/s"
    forcebal['v_t_o_ex_r']['info'] = "Toroidal impurity velocity experimental"

    forcebal['v_p_o_e_r']         = {}
    forcebal['v_p_o_e_r']['data'] = None
    forcebal['v_p_o_e_r']['unit'] = "m/s"
    forcebal['v_p_o_e_r']['info'] = "Poloidal electron velocity"

    forcebal['v_p_o_i_r']         = {}
    forcebal['v_p_o_i_r']['data'] = None
    forcebal['v_p_o_i_r']['unit'] = "m/s"
    forcebal['v_p_o_i_r']['info'] = "Poloidal ion velocity"

    forcebal['v_p_o_im_r']         = {}
    forcebal['v_p_o_im_r']['data'] = None
    forcebal['v_p_o_im_r']['unit'] = "m/s"
    forcebal['v_p_o_im_r']['info'] = "Poloidal impurity velocity"

    forcebal['v_p_o_D_r']         = {}
    forcebal['v_p_o_D_r']['data'] = None
    forcebal['v_p_o_D_r']['unit'] = "m/s"
    forcebal['v_p_o_D_r']['info'] = ""

    forcebal['v_p_o_N_r']         = {}
    forcebal['v_p_o_N_r']['data'] = None
    forcebal['v_p_o_N_r']['unit'] = "m/s"
    forcebal['v_p_o_N_r']['info'] = ""

    forcebal['v_p_o_C_r']         = {}
    forcebal['v_p_o_C_r']['data'] = None
    forcebal['v_p_o_C_r']['unit'] = "m/s"
    forcebal['v_p_o_C_r']['info'] = ""

    forcebal['v_p_o_C_ex_r']         = {}
    forcebal['v_p_o_C_ex_r']['data'] = None
    forcebal['v_p_o_C_ex_r']['unit'] = "m/s"
    forcebal['v_p_o_C_ex_r']['info'] = ""

    forcebal['v_p_o_ex_r']         = {}
    forcebal['v_p_o_ex_r']['data'] = None
    forcebal['v_p_o_ex_r']['unit'] = "m/s"
    forcebal['v_p_o_ex_r']['info'] = "Poloidal impurity velocity experimental"


    forcebal['M_t_o_e_r']         = {}
    forcebal['M_t_o_e_r']['data'] = None
    forcebal['M_t_o_e_r']['unit'] = "m/s"
    forcebal['M_t_o_e_r']['info'] = "Toroidal electron velocity"

    forcebal['M_t_o_i_r']         = {}
    forcebal['M_t_o_i_r']['data'] = None
    forcebal['M_t_o_i_r']['unit'] = "m/s"
    forcebal['M_t_o_i_r']['info'] = "Toroidal ion velocity"

    forcebal['M_t_o_im_r']         = {}
    forcebal['M_t_o_im_r']['data'] = None
    forcebal['M_t_o_im_r']['unit'] = "m/s"
    forcebal['M_t_o_im_r']['info'] = "Toroidal impurity velocity"

    forcebal['M_t_o_D_r']         = {}
    forcebal['M_t_o_D_r']['data'] = None
    forcebal['M_t_o_D_r']['unit'] = "m/s"
    forcebal['M_t_o_D_r']['info'] = ""

    forcebal['M_t_o_N_r']         = {}
    forcebal['M_t_o_N_r']['data'] = None
    forcebal['M_t_o_N_r']['unit'] = "m/s"
    forcebal['M_t_o_N_r']['info'] = ""

    forcebal['M_t_o_C_r']         = {}
    forcebal['M_t_o_C_r']['data'] = None
    forcebal['M_t_o_C_r']['unit'] = "m/s"
    forcebal['M_t_o_C_r']['info'] = ""

    forcebal['M_t_o_C_ex_r']         = {}
    forcebal['M_t_o_C_ex_r']['data'] = None
    forcebal['M_t_o_C_ex_r']['unit'] = "m/s"
    forcebal['M_t_o_C_ex_r']['info'] = ""

    forcebal['M_t_o_ex_r']         = {}
    forcebal['M_t_o_ex_r']['data'] = None
    forcebal['M_t_o_ex_r']['unit'] = "m/s"
    forcebal['M_t_o_ex_r']['info'] = "Toroidal impurity velocity experimental"

    forcebal['M_p_o_e_r']         = {}
    forcebal['M_p_o_e_r']['data'] = None
    forcebal['M_p_o_e_r']['unit'] = "m/s"
    forcebal['M_p_o_e_r']['info'] = "Poloidal electron velocity"

    forcebal['M_p_o_i_r']         = {}
    forcebal['M_p_o_i_r']['data'] = None
    forcebal['M_p_o_i_r']['unit'] = "m/s"
    forcebal['M_p_o_i_r']['info'] = "Poloidal ion velocity"

    forcebal['M_p_o_im_r']         = {}
    forcebal['M_p_o_im_r']['data'] = None
    forcebal['M_p_o_im_r']['unit'] = "m/s"
    forcebal['M_p_o_im_r']['info'] = "Poloidal impurity velocity"

    forcebal['M_p_o_D_r']         = {}
    forcebal['M_p_o_D_r']['data'] = None
    forcebal['M_p_o_D_r']['unit'] = "m/s"
    forcebal['M_p_o_D_r']['info'] = ""

    forcebal['M_p_o_N_r']         = {}
    forcebal['M_p_o_N_r']['data'] = None
    forcebal['M_p_o_N_r']['unit'] = "m/s"
    forcebal['M_p_o_N_r']['info'] = ""

    forcebal['M_p_o_C_r']         = {}
    forcebal['M_p_o_C_r']['data'] = None
    forcebal['M_p_o_C_r']['unit'] = "m/s"
    forcebal['M_p_o_C_r']['info'] = ""

    forcebal['M_p_o_C_ex_r']         = {}
    forcebal['M_p_o_C_ex_r']['data'] = None
    forcebal['M_p_o_C_ex_r']['unit'] = "m/s"
    forcebal['M_p_o_C_ex_r']['info'] = ""

    forcebal['M_p_o_ex_r']         = {}
    forcebal['M_p_o_ex_r']['data'] = None
    forcebal['M_p_o_ex_r']['unit'] = "m/s"
    forcebal['M_p_o_ex_r']['info'] = "Poloidal impurity velocity experimental"


    #>============================================
    #>Flux surface particle flow velocities
    #>============================================

    forcebal['u_par_e_r']         = {}
    forcebal['u_par_e_r']['data'] = None
    forcebal['u_par_e_r']['unit'] = "T*m/s"
    forcebal['u_par_e_r']['info'] = "<u.B> electron"

    forcebal['u_par_i_r']         = {}
    forcebal['u_par_i_r']['data'] = None
    forcebal['u_par_i_r']['unit'] = "T*m/s"
    forcebal['u_par_i_r']['info'] = "<u.B> ion"

    forcebal['u_par_im_r']         = {}
    forcebal['u_par_im_r']['data'] = None
    forcebal['u_par_im_r']['unit'] = "T*m/s"
    forcebal['u_par_im_r']['info'] = "<u.B> impurity"

    forcebal['u_par_D_r']         = {}
    forcebal['u_par_D_r']['data'] = None
    forcebal['u_par_D_r']['unit'] = "T*m/s"
    forcebal['u_par_D_r']['info'] = ""

    forcebal['u_par_N_r']         = {}
    forcebal['u_par_N_r']['data'] = None
    forcebal['u_par_N_r']['unit'] = "T*m/s"
    forcebal['u_par_N_r']['info'] = ""

    forcebal['u_par_C_r']         = {}
    forcebal['u_par_C_r']['data'] = None
    forcebal['u_par_C_r']['unit'] = "T*m/s"
    forcebal['u_par_C_r']['info'] = ""

    forcebal['u_p_e_r']         = {}
    forcebal['u_p_e_r']['data'] = None
    forcebal['u_p_e_r']['unit'] = "m/s/T"
    forcebal['u_p_e_r']['info'] = "<u.grad(theta)>/<B.grad(theta)> electron"

    forcebal['u_p_i_r']         = {}
    forcebal['u_p_i_r']['data'] = None
    forcebal['u_p_i_r']['unit'] = "m/s/T"
    forcebal['u_p_i_r']['info'] = "<u.grad(theta)>/<B.grad(theta)> ion"

    forcebal['u_p_im_r']         = {}
    forcebal['u_p_im_r']['data'] = None
    forcebal['u_p_im_r']['unit'] = "m/s/T"
    forcebal['u_p_im_r']['info'] = "<u.grad(theta)>/<B.grad(theta)> impurity"

    forcebal['u_p_D_r']         = {}
    forcebal['u_p_D_r']['data'] = None
    forcebal['u_p_D_r']['unit'] = "m/s/T"
    forcebal['u_p_D_r']['info'] = ""

    forcebal['u_p_N_r']         = {}
    forcebal['u_p_N_r']['data'] = None
    forcebal['u_p_N_r']['unit'] = "m/s/T"
    forcebal['u_p_N_r']['info'] = ""

    forcebal['u_p_C_r']         = {}
    forcebal['u_p_C_r']['data'] = None
    forcebal['u_p_C_r']['unit'] = "m/s/T"
    forcebal['u_p_C_r']['info'] = ""

    forcebal['u_p_C_ex_r']         = {}
    forcebal['u_p_C_ex_r']['data'] = None
    forcebal['u_p_C_ex_r']['unit'] = "m/s/T"
    forcebal['u_p_C_ex_r']['info'] = ""

    forcebal['dummy']         = {}
    forcebal['dummy']['data'] = None
    forcebal['dummy']['unit'] = ""
    forcebal['dummy']['info'] = ""

    return forcebal

def read_1d_file(fname):
    fid = open(fname,"r")
    lines = fid.readlines()
    nvars = int(lines[0].split()[0])
    ndata = int(lines[0].split()[1])

    nlines_vars = int((nvars+4)/5)
    nlines_data = int((ndata+4)/5)

    outbound_range = [i for i in range(ndata,5*nlines_data)]

    vars_line_bgn = 1
    unit_line_bgn = nlines_vars + 1

    forcebal = get_forcebal_vars()
    for ilines_vars in range(nlines_vars):
        line_vars = lines[vars_line_bgn + ilines_vars].split()
        line_unit = lines[unit_line_bgn + ilines_vars].split()

        nline_vars = len(line_vars)
        nline_unit = len(line_unit)

        data_line_bgn = 1 + 2 * nlines_vars + ilines_vars * nline_vars * nlines_data

        for iline_vars in range(nline_vars):
            ivars = line_vars[iline_vars]
            iunit = line_unit[iline_vars]
            forcebal[ivars]['unit'] = iunit
            forcebal[ivars]['data'] = npy.zeros(5*nlines_data)
            for ilines_data in range(nlines_data):
                idata = lines[data_line_bgn + iline_vars * nlines_data + ilines_data].split()
                forcebal[ivars]['data'][5*ilines_data+0] = idata[0]
                forcebal[ivars]['data'][5*ilines_data+1] = idata[1]
                forcebal[ivars]['data'][5*ilines_data+2] = idata[2]
                forcebal[ivars]['data'][5*ilines_data+3] = idata[3]
                forcebal[ivars]['data'][5*ilines_data+4] = idata[4]
            forcebal[ivars]['data'] = npy.delete(forcebal[ivars]['data'],outbound_range)

    return forcebal


def read_sum_file(fname):
    forcebal = get_forcebal_vars()

    fid = open(fname,"r")
    lines = fid.readlines()
    nlines = len(lines)

    forcebal['l_banana']    =     float(lines[37].split('=')[1].strip())
    forcebal['l_pfirsch']   =     float(lines[39].split('=')[1].strip())
    forcebal['l_classical'] =     float(lines[41].split('=')[1].strip())
    forcebal['l_potato']    =     float(lines[43].split('=')[1].strip())
    forcebal['l_squeeze']   =     float(lines[45].split('=')[1].strip())
    forcebal['diiid']       =           lines[47].split('=')[1].strip()
    forcebal['id_shot']     = int(float(lines[49].split('=')[1].strip()))
    forcebal['time']        =     float(lines[51].split('=')[1].strip()) 
    forcebal['nr_r']        = int(float(lines[53].split('=')[1].strip()))
    forcebal['a0']          =     float(lines[55].split('=')[1].strip())
    forcebal['a1']          =     float(lines[57].split('=')[1].strip())
    forcebal['R0']          =     float(lines[59].split('=')[1].strip())
    forcebal['Bt0']         =     float(lines[61].split('=')[1].strip())
    forcebal['elongation']  =     float(lines[63].split('=')[1].strip())

    cline = 209
    while cline < nlines:
        cline += 1
        var_name_list = lines[cline].split();     cline += 1
        var_unit_list = lines[cline].split();     cline += 1
        forcebal[var_name_list[1]]['data'] = npy.zeros(forcebal['nr_r'])
        forcebal[var_name_list[2]]['data'] = npy.zeros(forcebal['nr_r'])
        forcebal[var_name_list[3]]['data'] = npy.zeros(forcebal['nr_r'])
        forcebal[var_name_list[4]]['data'] = npy.zeros(forcebal['nr_r'])
        forcebal[var_name_list[5]]['data'] = npy.zeros(forcebal['nr_r'])
        for i in range(forcebal['nr_r']):
            var_data_list = lines[cline].split(); cline += 1
            forcebal[var_name_list[1]]['data'][i] = float(var_data_list[1])
            forcebal[var_name_list[2]]['data'][i] = float(var_data_list[2])
            forcebal[var_name_list[3]]['data'][i] = float(var_data_list[3])
            forcebal[var_name_list[4]]['data'][i] = float(var_data_list[4])
            forcebal[var_name_list[5]]['data'][i] = float(var_data_list[5])
    
    return forcebal


def read_cdf_file(fname):
    forcebal = get_forcebal_vars()

    ncfid = ncdf.Dataset(fname)
    ncfvars = ncfid.variables.keys()
    for ncfvar in ncfvars:
        forcebal[ncfvar]['data'] = ncfid.variables[ncfvar][:]

    return forcebal


def write_input_files(statefname='',statedata={},shot_id=100001,time_id=10001):
    if statefname:
       try:
           statedata = cheasefile.get_plasmastate(instatefpath=statefname)
       except IOError:
           statedata = cheasefile.get_plasmastate(statefpath=statefname)

    if type(shot_id) == str: shot_id = int(float(shot_id))
    if type(time_id) == str: time_id = int(float(time_id))
    
    fname = "te%06d.%05d" % (shot_id,time_id)
    records = npy.column_stack((statedata['rho'],statedata['Te']*1.0e-3))
    npy.savetxt(fname,records, fmt='%7.7E', delimiter='\t')

    fname = "ti%06d.%05d" % (shot_id,time_id)
    records = npy.column_stack((statedata['rho'],statedata['Ti']*1.0e-3))
    npy.savetxt(fname,records, fmt='%7.7E', delimiter='\t')

    fname = "ne%06d.%05d" % (shot_id,time_id)
    records = npy.column_stack((statedata['rho'],statedata['ne']*1.0e-19))
    npy.savetxt(fname,records, fmt='%7.7E', delimiter='\t')

    fname = "nc%06d.%05d" % (shot_id,time_id)
    records = npy.column_stack((statedata['rho'],statedata['nz']*1.0e-19))
    npy.savetxt(fname,records, fmt='%7.7E', delimiter='\t')

    fname = "zf%06d.%05d" % (shot_id,time_id)
    records = npy.column_stack((statedata['rho'],statedata['zeff']))
    npy.savetxt(fname,records, fmt='%7.7E', delimiter='\t')

    fname = "vt%06d.%05d" % (shot_id,time_id)
    records = npy.column_stack((statedata['rho'],statedata['omega']))
    npy.savetxt(fname,records, fmt='%7.7E', delimiter='\t')

    return 1

if __name__ == "__main__":
   #fname = "1d_92664_01000_a05.fbl"
   #read_1d_file(fname)
    fname = "sum_84713_02080_a02.dat"
    read_sum_file(fname)
   #fname = "cdf_84713_02080_a02.nc"
   #read_cdf_file(fname)
