"""
 -----------------------------------------------------------------------
 freegs io 
 -----------------------------------------------------------------------
"""
import numpy as np
from Namelist import Namelist
import freegs
import freegs.geqdsk
import freegs.machine
from freegs.shaped_coil import ShapedCoil
import json
from fastran.util.zinterp import zinterp

def read_coil_data(f_coil_data):
    with open(f_coil_data, "r") as f:
         coils_input = json.load(f)
    coils = []
    for key in coils_input:
        coils.append([key, ShapedCoil(coils_input[key])])
        print(coils_input[key])
    return coils

def call_freegs(f_instate, f_inefit, init=False, boundary='all'):
    instate = Namelist(f_instate)

    inefit = Namelist(f_inefit)
    r0 = inefit["inefit"]["r0"][0]
    b0 = inefit["inefit"]["b0"][0]
    ip = inefit["inefit"]["ip"][0]
    print('r0, b0, ip = ', r0, b0, ip)
    rbdry = inefit["inefit"]["rbdry"]
    zbdry = inefit["inefit"]["zbdry"]
    rwall = inefit["inefit"]["rlim"]
    zwall = inefit["inefit"]["zlim"]
    press = inefit["inefit"]["press"]

    k_x1 = np.argmin(zbdry)
    k_x2 = np.argmax(zbdry)
    k_in = np.argmin(rbdry)
    k_out = np.argmax(rbdry)
    
    Rx_1, Zx_1 = rbdry[k_x1], zbdry[k_x1]
    Rx_1p, Zx_1p = rbdry[k_x1+1], zbdry[k_x1+1]
    Rx_1m, Zx_1m = rbdry[k_x1-1], zbdry[k_x1-1]

    Rx_2, Zx_2 = rbdry[k_x2], zbdry[k_x2]
    Rx_2p, Zx_2p = rbdry[k_x2+1], zbdry[k_x2+1]
    Rx_2m, Zx_2m = rbdry[k_x2-1], zbdry[k_x2-1]

    R_in, Z_in = rbdry[k_in], zbdry[k_in]
    R_out, Z_out = rbdry[k_out], zbdry[k_out]

    print('Rx_1 =', Rx_1, Rx_1p, Rx_1m)
    print('Rx_2 =', Rx_2, Rx_2p, Rx_2m)
    print('R_in, R_out =', R_in, R_out)

    if boundary == 'all':
        r_isoflux = rbdry
        z_isoflux = zbdry
    else:
        r_isoflux = [R_in, Rx_1m, Rx_1, Rx_1p, R_out, Rx_2m, Rx_2, Rx_2p, R_in]
        z_isoflux = [Z_in, Zx_1m, Zx_1, Zx_1p, Z_out, Zx_2m, Zx_2, Zx_2p, Z_in]

    xpoints = [(Rx_1, Zx_1), (Rx_2, Zx_2)]
    
    n_isoflux = len(r_isoflux)
    isoflux = []
    for k in range(n_isoflux-1):
        p = [(r_isoflux[k], z_isoflux[k], r_isoflux[k+1], z_isoflux[k+1],)]
        isoflux += p

    tokamak = freegs.machine.DIIID_Tokamak(rwall=rwall, zwall=zwall)
#   f_coil_data = 'd3d_coils'
#   coils_shape = read_coil_data(f_coil_data)
#   tokamak = freegs.machine.Machine(coils_shape, freegs.machine.Wall(rwall, zwall))

    nx = 129
    ny = 129

    eq = freegs.Equilibrium(tokamak=tokamak,
                    Rmin=0.84, Rmax=2.54, # Radial domain
                    Zmin=-1.6, Zmax=1.6, # Height range
                    nx=nx, ny=ny, # Number of grid points
                    boundary=freegs.boundary.freeBoundaryHagenow) # Boundary condition

    if init:
        profiles = freegs.jtor.ConstrainPaxisIp(press[0], # Plasma pressure on axis [Pascals]
                                                         ip, # Plasma current [Amps]
                                                         r0*b0) # Vacuum f=R*Bt
    else: 
        infreegs = Namelist("infreegs")
        pprime_ext = infreegs['profile_ext']['pprime_ext']
        ffprime_ext = infreegs['profile_ext']['ffprim_ext']
        p_ext = infreegs['profile_ext']['p_ext']
        f_ext = infreegs['profile_ext']['f_ext']
        
        npsi = infreegs['profile_ext']['npsi_ext'][0]
        psin = infreegs['profile_ext']['psin_ext']
        
        pprime = zinterp(psin, pprime_ext)
        ffprime = zinterp(psin, ffprime_ext)
        pressure =  zinterp(psin, p_ext)
        fpol =  zinterp(psin, f_ext)
        
        def pprime_func(psin):
            return pprime(psin)
        
        def ffprime_func(psin):
           return ffprime(psin)
        
        def p_func(psin):
           return pressure(psin)
        
        def f_func(psin):
           return fpol(psin)

        profiles = freegs.jtor.ProfilesPprimeFfprime(pprime_func, 
                                                              ffprime_func, 
                                                              r0*b0, 
                                                              p_func=p_func, 
                                                              f_func=f_func)
         
#   constrain = freegs.control.constrain(xpoints=xpoints, isoflux=isoflux)
    constrain = freegs.control.constrain(isoflux=isoflux)

    freegs.solve(eq, # The equilibrium to adjust
                 profiles, # The toroidal current profile function
                 constrain) # Constraint function to set coil currents
    
    # eq now contains the solution
    
    print("Done!")
    print("Plasma current: %e Amps" % (eq.plasmaCurrent()))
    print("Plasma pressure on axis: %e Pascals" % (eq.pressure(0.0)))
    print("Poloidal beta: %e" % (eq.poloidalBeta()))

    with open("lsn.geqdsk", "w") as f:
        freegs.geqdsk.write(eq, f, R0=r0)

