import os
import numpy as np
import netCDF4
from fastran.equilibrium.efit_eqdsk import readg
from fastran.plasmastate.plasmastate import plasmastate
from fastran.util.zinterp import zinterp

def write_inputfiles(f_state, f_eqdsk, nexp=201):
    print('tglfep: write_inputfiles')

# torfluxa | Wb/radian      --> rhob = sqrt(torflux / pi)
# rcentr | m                --> r0
# bcentr | T                --> b0
# current | MA              --> ip
# rho | -                   --> rho
# rmin | m                  --> rminor
# polflux | Wb/radian       --> psi
# q | -                     --> qmhd
# rmaj | m                  --> rmajor
# zmag | m                  --> zaxis 
# kappa | -                 --> kappa
# delta | -                 --> delta
# zeta | -                  --> zeta 
# ne | 10^19/m^3            --> ne
# ni | 10^19/m^3            --> ps['ns'][1, :], ps['ns'][2, :], ...       
# te | keV                  --> te
# ti | keV                  --> ps['Ts'][1, :], ps['ns'][2, :], ...  
# ptot | Pa                 --> pmhd
# z_eff | -                 --> zeff

    # -- read geqdsk and plasma state file
    geq = readg(f_eqdsk)
    r0 = geq['rzero']
    b0 = np.abs(geq['bcentr'])
    ip = geq['cpasma']

    ps = plasmastate('ips', 1)
    ps.read(f_state)
    spec = ps.get_species()

    # --  profiles
    nrho = len(ps['rho'])
    rho = ps['rho'][:]
    ne = ps['ns'][0, :]*1.e-19
    te = ps['Ts'][0, :]
    ti = ps['Ti'][:]
    zeff = ps['Zeff'][:]
    omega = ps['omegat'][:]
    ne = ps.cell2node_bdry(ne)
    te = ps.cell2node_bdry(te)
    ti = ps.cell2node_bdry(ti)
    zeff = ps.cell2node_bdry(zeff)
    omega = ps.cell2node_bdry(omega)

    omega = ps['p_EQ'][:]

    density_beam = ps.dump_profile(rho, 'rho_nbi', 'nbeami', k=0)*1.e-19
    wbeam = ps.dump_profile(rho, 'rho_nbi', 'eperp_beami', k=0) \
        + ps.dump_profile(rho, 'rho_nbi', 'epll_beami', k=0)
    wbeam = density_beam*wbeam*1.602e-3 # MJ/m**3

    # -- metrics
    psi = ps['psipol'][:]/ps['psipol'][-1]  # equi-drho grid
    rho = np.sqrt(ps['phit'][:]/ps['phit'][-1])
    rhob = (ps['phit'][-1]/np.pi/b0)**0.5

    rhopsi = zinterp(psi, rho)
    ipol = ps['g_eq'][:]/(r0*b0)
    ipol = np.abs(ipol)
    volp = 4.0*np.pi*np.pi*rho*rhob*r0/ipol/(r0*r0*ps['gr2i'][:])

    g11 = volp*ps['grho2'][:]*rhob**2
    g22 = r0*volp/(4.0*np.pi*np.pi)*ps['grho2r2i'][:]*rhob**2
    g33 = r0*r0*ps['gr2i'][:]
    gradrho = ps['grho1'][:]*rhob
    area = ps['surf'][:]
    rmajor = ps['Rmajor_mean'][:]
    rminor = ps['rMinor_mean'][:]
    shift = rmajor-r0
    kappa = ps['elong'][:]
    delta = ps['triang'][:]
    pmhd = ps['P_eq'][:]
    qmhd = ps['q_eq'][:]
    
    zaxis = ps['Z_axis']
    # zeta = 0.5 * (ps['squareLO'][:] + squareUO[:]) 

    # write input.gacode
    nions = len ps['ns']-1

    f = open('input.gacode', 'w')
      f.write('# transitory statefile for TGLF-EP.'+'\n')
      f.write('# nexp'+'\n')
      f.write('{nrho}'+'\n')
      f.write('# nion'+'\n')
      f.write('{nions}'+'\n')
      f.write('# name'+'\n')
      f.write('D D C'+'\n')
      f.write('# type'+'\n')
      f.write('[therm] [fast] [therm]'+'\n')
      f.write('# masse'+'\n')
      f.write('5.4488739E-04'+'\n')
      f.write('# mass'+'\n')
      f.write('2.0000000E+00 2.0000000E+00 12.0000000E+00'+'\n')
      f.write('# ze'+'\n')
      f.write('# z'+'\n')
      f.write('1.0000000E+00 1.0000000E+00 6.0000000E+00'+'\n')
      f.write('# torfluxa | Wb/radian')
      f.write('{rhob}'+'\n')
      f.write('# rcentr | m',+'\n')
      f.write('{r0}'+'\n')
      f.write('# bcentr | T'+'\n')
      f.write('{b0}'+'\n')
      f.write('# current'+'\n')
      f.write('{ip}'+'\n')
      f.write('# rho | -'+'\n')
      for ir in range(nrho):
          f.write('{:>3} {:>14.7E}'+'\n'.format(ir,rho[ir]))
      f.write('# rmin | m'+'\n')
      for ir in range(nrho):
          f.write('{:>3} {:>14.7E}'+'\n'.format(ir,rminor[ir]))
      f.write('# polflux | Wb/radian'+'\n')
      for ir in range(nrho):
          f.write('{:>3} {:>14.7E}'+'\n'.format(ir,psi[ir]))
      f.write('# q | -'+'\n')
      for ir in range(nrho):
          f.write('{:>3} {:>14.7E}'+'\n'.format(ir,qmhd[ir]))
      f.write('# rmaj | m'+'\n')
      for ir in range(nrho):
          f.write('{:>3} {:>14.7E}'+'\n'.format(ir,rmajor[ir]))
      f.write('# zmag | m'+'\n')
      for ir in range(nrho):
          f.write('{:>3} {:>14.7E}'+'\n'.format(ir,zaxis[ir]))
      f.write('# kappa | -'+'\n')
      for ir in range(nrho):
          f.write('{:>3} {:>14.7E}'+'\n'.format(ir,kappa[ir]))
      f.write('# delta | -'+'\n')
      for ir in range(nrho):
          f.write('{:>3} {:>14.7E}'+'\n'.format(ir,delta[ir]))
      f.write('# zeta | -'+'\n')
      for ir in range(nrho):
          f.write('{:>3} {:>14.7E}'+'\n'.format(ir,zeta[ir]))
      f.write('# ne | 10^19/m^3'+'\n')
      for ir in range(nrho):
          f.write('{:>3} {:>14.7E}'+'\n'.format(ir,ne[ir]))
      f.write('# ni | 10^19/m^3'+'\n')
      for ir in range(nrho):
          f.write('{:>3}'.format(ir,ps[ir]))
          for jions in range(nions):
              f.write(' {:>14.7E}'.format(ps['ns'][jions+1,ir]))
          f.write('\n')
      f.write('# te | keV'+'\n')
      for ir in range(nrho):
          f.write('{:>3} {:>14.7E}'+'\n'.format(ir,te[ir]))
      f.write('# ti | keV'+'\n')
      for ir in range(nrho):
          f.write('{:>3}'.format(ir,ps[ir])))                          
          for jions in range(nions):
              f.write(' {:>14.7E}'.format(ps['Ts'][jions+1,ir]))
          f.write('\n')    
      f.write('# ptot | Pa'+'\n')
      for ir in range(nrho):
          f.write('{:>3} {:>14.7E}'+'\n'.format(ir,pmhd[ir]))
      f.write('# z_eff | -'+'\n')
      for ir in range(nrho):
          f.write('{:>3} {:>14.7E}'+'\n'.format(ir,zeff[ir]))

    f.close()

def update_state():
    print('tglfep update_state')
    # read tglfep output
    # update state file
   

