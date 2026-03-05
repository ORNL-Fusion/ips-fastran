import imas
import numpy as np
import os
import datetime
import netCDF4
import imas

from fastran.equilibrium.efit_eqdsk import readg
from fastran.plasmastate.plasmastate import plasmastate
from fastran.util.zinterp import zinterp

def to_imas(f_eqdsk, f_state, shot, itime, run, provider, user):
    time = itime / 1.e3

    ps = plasmastate('ips', 1)
    ps.read(f_state)
    spec = ps.get_species()
    
    R_axis = ps['R_axis']
    B_axis_vac = ps['B_axis_vac']
    print(f'R_axis = {R_axis}, B_axis_vac = {B_axis_vac}')
    r0 = 1.69555
    b0 = R_axis * B_axis_vac / r0 
    print(f'R0 = {r0}, B_axis_vac = {b0}')
    ip = ps['curt'][-1]

    geq = readg(f_eqdsk)

    ps_xe  = 1.6022e-19
    ps_mp  = 1.6726e-27
    
    # ---------------------------------------------------------------------
    # SUMMARY
    
    summary = imas.summary()
    summary.ids_properties.homogeneous_time = imas.imasdef.IDS_TIME_MODE_HOMOGENEOUS # 1
    summary.ids_properties.provider = provider
    summary.ids_properties.creation_date = datetime.datetime.now().strftime("%y-%m-%d")
    
    summary.time.resize(1)
    summary.time[0] = float(time)
    
    # ---------------------------------------------------------------------
    # CORE PROFILES
    
    core_profiles = imas.core_profiles()
    core_profiles.ids_properties.homogeneous_time = imas.imasdef.IDS_TIME_MODE_HOMOGENEOUS # 1 
    core_profiles.time.resize(1)
    core_profiles.time[0] = float(time)
    
    core_profiles.profiles_1d.resize(1)
    
    # ---------------------------------------------------------------------
    # Grid
    rho = ps['rho'][:]
    vol = ps['vol'][:]
    nrho = len(rho)
    
    core_profiles.profiles_1d[0].grid.rho_tor_norm = rho 
    ##core_profiles.profiles_1d[0].grid.rho_tor = np.linspace(0, 1, nrho) * prof['arho']
    ##core_profiles.profiles_1d[0].grid.psi = prof['polflux'] * 2. * np.pi 
    core_profiles.profiles_1d[0].grid.volume = vol
    
    # ---------------------------------------------------------------------
    # Electron
    ne = ps['ns'][0, :]
    ne = ps.cell2node_bdry(ne)
    te = ps['Ts'][0, :] * 1.e3
    te = ps.cell2node_bdry(te)
    pe = 1.602e-19 * ne * te
    
    core_profiles.profiles_1d[0].electrons.density =  ne
    core_profiles.profiles_1d[0].electrons.temperature = te 
    core_profiles.profiles_1d[0].electrons.pressure = pe
    core_profiles.profiles_1d[0].electrons.pressure_thermal = pe
    
    zeff = ps['Zeff'][:]
    omega = ps['omegat'][:]
    zeff = ps.cell2node_bdry(zeff)
    omega = ps.cell2node_bdry(omega)
    
    # ---------------------------------------------------------------------
    # ION - thermal
    nion = 2 
    nb = 1
    map_th = {0:0, 1:2}
    map_beam = {0:3}
    
    core_profiles.profiles_1d[0].ion.resize(nion + nb)
    imas_spec_id = {'D':0, 'C':1, 'D_beam':2}
    
    pi_tot = np.zeros(nrho)
    for k in range(nion):
        ni = ps['ns'][k + 1, :]
        ni = ps.cell2node_bdry(ni)
        ti = ps['Ti'][:] * 1.e3
        ti = ps.cell2node_bdry(ti)
        pi = 1.602e-19 * ni * ti
        
        label = ps['S_name'][k + 1].strip()
    
        Z = round( ps['qatom_S'][k + 1] / ps_xe )
        A = round( ps['m_S'][k + 1] / ps_mp )
    
        print(k, Z, A)
        print('>>> label, Z, A, imas_spec_id:', label, Z, A, imas_spec_id[label])
    
        k_imas = imas_spec_id[label]
        
        core_profiles.profiles_1d[0].ion[k_imas].label = label 
        core_profiles.profiles_1d[0].ion[k_imas].element.resize(1)
        core_profiles.profiles_1d[0].ion[k_imas].element[0].a = float(A)
        core_profiles.profiles_1d[0].ion[k_imas].element[0].z_n = float(Z)
    
        core_profiles.profiles_1d[0].ion[k_imas].density = ni
        core_profiles.profiles_1d[0].ion[k_imas].density_thermal = ni
        core_profiles.profiles_1d[0].ion[k_imas].temperature = ti
        core_profiles.profiles_1d[0].ion[k_imas].pressure_thermal = pi
    
        pi_tot += pi
    
    for k in range(nb):
        nbeam = ps.dump_profile(rho, 'rho_nbi', 'nbeami', k=0)
        pbeam_pll = 1.602e-16 * nbeam * ps.dump_profile(rho, 'rho_nbi', 'epll_beami', k=0)
        pbeam_perp = 1.602e-16 * nbeam * ps.dump_profile(rho, 'rho_nbi', 'eperp_beami', k=0)
    
        label = ps['SNBI_name'][k].strip()
        Z = round( ps['qatom_SNBI'][k] / ps_xe )
        A = round( ps['m_SNBI'][k] / ps_mp )
    
        print('>>> label, Z, A, imas_spec_id:', label, Z, A, imas_spec_id[label]) #.split('_')[0]])
    
        k_imas = imas_spec_id[label] # .split('_')[0]
    
        core_profiles.profiles_1d[0].ion[k_imas].label = label 
        core_profiles.profiles_1d[0].ion[k_imas].element.resize(1)
        core_profiles.profiles_1d[0].ion[k_imas].element[0].a = float(A)
        core_profiles.profiles_1d[0].ion[k_imas].element[0].z_n = float(Z)
    
        core_profiles.profiles_1d[0].ion[k_imas].density = nbeam
        core_profiles.profiles_1d[0].ion[k_imas].density_fast = nbeam
        core_profiles.profiles_1d[0].ion[k_imas].pressure_fast_parallel  = pbeam_pll
        core_profiles.profiles_1d[0].ion[k_imas].pressure_fast_perpendicular = pbeam_perp
    
    p_th = pi_tot + pe
    core_profiles.profiles_1d[0].pressure_ion_total = pi_tot
    core_profiles.profiles_1d[0].pressure_parallel = pbeam_pll
    core_profiles.profiles_1d[0].pressure_perpendicular = pbeam_perp
    core_profiles.profiles_1d[0].pressure_thermal = p_th
    
    # ---------------------------------------------------------------------
    # Source 
    core_sources = imas.core_sources()
    core_sources.ids_properties.homogeneous_time = imas.imasdef.IDS_TIME_MODE_HOMOGENEOUS # 1 
    core_sources.time.resize(1)
    core_sources.time[0] = float(time)
    core_sources.source.resize(3)
    
    id_total = 1
    id_nb = 2
    id_ec = 3
    id_lh = 4
    id_ic = 5
    id_fusion = 6
    id_ohm = 7 
    id_br = 8
    id_synchrotron = 9
    id_line = 10
    
    pe_nb = ps.dump_vol_profile(rho, 'rho_nbi', 'pbe')
    pi_nb = ps.dump_vol_profile(rho, 'rho_nbi', 'pbi') + ps.dump_vol_profile(rho, 'rho_nbi', 'pbth')
    pth_nb = ps.dump_vol_profile(rho, 'rho_nbi', 'pbth')
    se_nb = ps.dump_vol_profile(rho, 'rho_nbi', 'sbedep') + ps.dump_vol_profile(rho, 'rho_nbi', 'sbehalo')
    pe_ec = ps.dump_vol_profile(rho, 'rho_ecrf', 'peech') 
    pe_ohm = ps.dump_vol_profile(rho, 'rho_ecrf', 'pohme') 
    
    core_sources.source[0].identifier.name = 'NB'
    core_sources.source[0].identifier.index = id_nb
    core_sources.source[0].profiles_1d.resize(1)
    core_sources.source[0].profiles_1d[0].grid.rho_tor_norm    = rho
    core_sources.source[0].profiles_1d[0].grid.volume          = vol
    core_sources.source[0].profiles_1d[0].electrons.energy     = pe_nb     
    core_sources.source[0].profiles_1d[0].total_ion_energy     = pi_nb    
    core_sources.source[0].profiles_1d[0].electrons.particles  = se_nb   
    
    core_sources.source[1].identifier.name = 'EC'
    core_sources.source[1].identifier.index = id_ec
    core_sources.source[1].profiles_1d.resize(1)
    core_sources.source[1].profiles_1d[0].grid.rho_tor_norm    = rho
    core_sources.source[1].profiles_1d[0].grid.volume          = vol
    core_sources.source[1].profiles_1d[0].electrons.energy     = pe_ec    
    
    core_sources.source[2].identifier.name = 'ohmic'
    core_sources.source[2].identifier.index = id_ohm
    core_sources.source[2].profiles_1d.resize(1)
    core_sources.source[2].profiles_1d[0].grid.rho_tor_norm    = rho
    core_sources.source[2].profiles_1d[0].grid.volume          = vol
    core_sources.source[2].profiles_1d[0].electrons.energy     = pe_ohm    
    
    # ---------------------------------------------------------------------
    # Equilibrium
    equilibrium = imas.equilibrium()
    equilibrium.ids_properties.homogeneous_time = imas.imasdef.IDS_TIME_MODE_HOMOGENEOUS # 1 
    equilibrium.time.resize(1)
    equilibrium.time[0] = float(time)
    equilibrium.time_slice.resize(1)
    
    equilibrium.vacuum_toroidal_field.r0 = r0
    equilibrium.vacuum_toroidal_field.b0 = np.array([b0])
    
    psipol = ps['psipol'][:] * 2. * np.pi # equi-drho grid
    psi = ps['psipol'][:] / ps['psipol'][-1]  
    psi_axis = psipol[0]
    psi_bdry = psipol[-1]
    
    r_axis = ps['R_axis']
    z_axis = ps['Z_axis']
    print(f'r_axis = {r_axis}, z_axis = {z_axis}')
    
    g_eq = ps['g_eq'][:]
    P_eq = ps['P_eq'][:]
    q_eq = ps['q_eq'][:]
    
    pprim = zinterp(psipol, P_eq, s=0)(psipol, der=1)
    ffprim = zinterp(psipol, g_eq, s=0)(psipol, der=1) * g_eq
    
    rgrid = ps['R_grid']
    zgrid = ps['Z_grid']
    psirz = ps['PsiRZ'][:, :].transpose() * 2. * np.pi
    
    rbdry = ps['R_geo'][:, -1]
    zbdry = ps['Z_geo'][:, -1]
    
    rlim = ps['rlim'][:]
    zlim = ps['zlim'][:]
    
    equilibrium.time_slice[0].global_quantities.ip = ip 
    
    equilibrium.time_slice[0].global_quantities.psi_axis     = psi_axis
    equilibrium.time_slice[0].global_quantities.psi_boundary = psi_bdry
    
    equilibrium.time_slice[0].boundary.outline.r               = rbdry   
    equilibrium.time_slice[0].boundary.outline.z               = zbdry    
    equilibrium.time_slice[0].profiles_1d.rho_tor_norm         = rho   
    equilibrium.time_slice[0].profiles_1d.psi                  = psipol   
    equilibrium.time_slice[0].profiles_1d.f                    = g_eq     
    equilibrium.time_slice[0].profiles_1d.f_df_dpsi            = ffprim   
    equilibrium.time_slice[0].profiles_1d.dpressure_dpsi       = pprim        
    equilibrium.time_slice[0].profiles_1d.pressure             = P_eq     
    equilibrium.time_slice[0].profiles_1d.q                    = q_eq     
    index = 1 # rectangular
    equilibrium.time_slice[0].profiles_2d.resize(1)
    equilibrium.time_slice[0].profiles_2d[0].grid_type.index   = index 
    equilibrium.time_slice[0].profiles_2d[0].grid.dim1         = rgrid 
    equilibrium.time_slice[0].profiles_2d[0].grid.dim2         = zgrid 
    equilibrium.time_slice[0].profiles_2d[0].psi               = psirz   
    
    # ---------------------------------------------------------------------
    # Save to h5 
    
    shot, run, user, database = shot, run, os.getenv('USER'), 'TEST'
    output = imas.DBEntry(imas.imasdef.HDF5_BACKEND, database, shot, run, user)
    output.create()
    
    output.put(summary)
    output.put(core_profiles)
    output.put(core_sources)
    output.put(equilibrium)
    
    output.close()
