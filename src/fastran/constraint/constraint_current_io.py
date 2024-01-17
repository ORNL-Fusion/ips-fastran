"""
 -----------------------------------------------------------------------
 constraint component
 -----------------------------------------------------------------------
"""

import sys
import os
import shutil
from numpy import *
from fastran.util.model_profile import profile_hat, profile_spline
from fastran.plasmastate.plasmastate import plasmastate
from fastran.equilibrium.efit_eqdsk import readg

def hatJ(
    f_state,
    f_eqdsk,
    rho_jbdry=0.8,
    q0=1.1,
    wid_jhat=0.2,
    rho_jhat_min=0.1,
    rho_jhat_max=0.6,
    niter_max=100,
    jboot=1.0):

    ps = plasmastate('ips', 1)
    ps.read(f_state)

    geq = readg(f_eqdsk)
    r0  = geq["rzero" ]
    b0  = abs(geq["bcentr"])
    ip  = geq['cpasma']

    nrho = len(ps["rho"])
    rho  = ps["rho"][:]
    vol  = ps["vol" ][:]
    area = ps["area"][:]
    ipol = ps["g_eq"][:]/(r0*b0)

    j_bs = 1.e-6 * ps.dump_j_parallel(rho,"rho","curr_bootstrap", r0, b0)
    j_bs *= jboot

    print(f'jboot = {jboot}')

    jbdry = j_bs[where( rho>=rho_jbdry)[0][0]]

    rhob = (ps["phit"][-1]/pi/b0)**0.5
    ipol= ps["g_eq"][:]/(r0*b0)
    volp = 4.0*pi*pi*rho*rhob*r0/ipol/(r0*r0*ps["gr2i"][:])
    g22  = r0*volp/(4.0*pi*pi)*ps["grho2r2i"][:]*rhob**2

    # drho = rhob/(nrho-1.)
    # volp0 = 0.5*(volp[0]+volp[1])
    # ipol0 = 0.5*(ipol[0]+ipol[1])
    # jpar0 = 2.0*pi*b0*rho[1]*rhob/q0
    # jpar0 /= 0.2/g22[1]
    # jpar0 /= volp0*drho/ipol0**2

    j_tot = zeros(nrho)
    curt = zeros(nrho)

    for iter_scale in range(niter_max):
        rho_jhat = 0.5 * (rho_jhat_min + rho_jhat_max)

        jpeak_min, jpeak_max = 0.1, 10.
        for iter in range(niter_max):
            jpeak = 0.5 * (jpeak_min + jpeak_max)

            j_core = (jpeak - jbdry) * profile_hat(nrho, rho_jhat, wid_jhat)(rho) + jbdry

            for i in range(nrho):
                if rho[i] < rho_jbdry: j_tot[i] = j_core[i]
                else: j_tot[i] = j_bs[i]

            curt = zeros(nrho)
            for i in range(nrho - 1):
                dV = vol[i + 1] - vol[i]
                jparm = 0.5 * (j_tot[i] + j_tot[i + 1])
                ipolm = 0.5 * (ipol[i] + ipol[i + 1])
                curt[i + 1] = (curt[i] / ipol[i] + jparm * dV / (2. * pi * r0 * ipolm**2)) * ipol[i + 1]

            if abs(curt[-1]-ip*1.e-6) < 0.001: break

            if curt[-1] > ip * 1.e-6:
               jpeak_max = jpeak
            else:
               jpeak_min = jpeak

        drho = rhob / (nrho - 1.)
        jpar0 = 0.5 * (j_tot[0] + j_tot[1])
        volp0 = 0.5 * (volp[0] + volp[1])
        ipol0 = 0.5 * (ipol[0] + ipol[1])
        dpsidrho = 0.2 / g22[1] * (jpar0 * volp0 * drho / ipol0**2)
        q0_cal = 2. * pi * b0* rho[1] * rhob / dpsidrho

        if q0_cal > q0:
            rho_jhat_max = rho_jhat
        else:
            rho_jhat_min = rho_jhat

        print('j0, rho_jhat, ip_cal =',jpar0, rho_jhat, curt[-1])

    ps.load_j_parallel(rho, 1.0e6*j_tot, "rho_eq", "curt", r0, b0, tot=True)

    ps.store(f_state)

def parabolicJ(f_state, f_eqdsk, rho_jbdry=0.8, q0=1.1, alpha=3.0, niter_max=100, jboot=1.0):
    ps = plasmastate('ips', 1)
    ps.read(f_state)

    geq = readg(f_eqdsk)
    r0  = geq["rzero" ]
    b0  = abs(geq["bcentr"])
    ip  = geq['cpasma']

    nrho = len(ps["rho"])
    rho  = ps["rho"][:]
    vol  = ps["vol" ][:]
    area = ps["area"][:]
    ipol = ps["g_eq"][:]/(r0*b0)

    j_bs = 1.e-6 * ps.dump_j_parallel(rho,"rho","curr_bootstrap", r0, b0)
    j_bs *= jboot

    print(f'jboot = {jboot}')

    jbdry = j_bs[where( rho>=rho_jbdry)[0][0]]

    rhob = (ps["phit"][-1]/pi/b0)**0.5
    ipol= ps["g_eq"][:]/(r0*b0)
    volp = 4.0*pi*pi*rho*rhob*r0/ipol/(r0*r0*ps["gr2i"][:])
    g22  = r0*volp/(4.0*pi*pi)*ps["grho2r2i"][:]*rhob**2

    drho = rhob/(nrho-1.)
    volp0 = 0.5*(volp[0]+volp[1])
    ipol0 = 0.5*(ipol[0]+ipol[1])
    jpar0 = 2.0*pi*b0*rho[1]*rhob/q0
    jpar0 /= 0.2/g22[1]
    jpar0 /= volp0*drho/ipol0**2

    beta_min = 0.1
    beta_max = 5.0
    j_tot = zeros(nrho)
    curt = zeros(nrho)

    for iter_scale in range(niter_max):
        beta = 0.5*(beta_min+beta_max)
        j_core = jpar0*(1.0-(rho/rho_jbdry)**alpha)**beta + jbdry

        for i in range(nrho):
            if rho[i]<rho_jbdry: j_tot[i] = j_core[i]
            else: j_tot[i] = j_bs[i]

        curt = zeros(nrho)
        for i in range(nrho-1):
            dV = vol[i+1]-vol[i]
            jparm = 0.5*(j_tot[i]+j_tot[i+1])
            ipolm = 0.5*(ipol[i]+ipol[i+1])
            curt[i+1] = (curt[i]/ipol[i]+jparm*dV/(2.0*pi*r0*ipolm**2))*ipol[i+1]

        if abs(curt[-1]-ip*1.e-6) < 0.001: break

        if curt[-1] > ip*1.e-6:
           beta_min = beta
        else:
           beta_max = beta

        print('j0, beta, ip_cal =',jpar0, beta, curt[-1])

    ps.load_j_parallel(rho, 1.0e6*j_tot, "rho_eq", "curt", r0, b0, tot=True)

    ps.store(f_state)

def gauss_asym(p, x):
    rvec = zeros(len(x))
    for k in range(len(x)):
        if x[k]>=p[1]: rvec[k] = p[0]*exp(-(x[k]-p[1])**2/(2*p[2]**2))
        if x[k]< p[1]: rvec[k] = p[0]*exp(-(x[k]-p[1])**2/(2*p[3]**2))
    return rvec

def broadJ(f_state, f_eqdsk, rho_j, j_in, j_out, rho_jbdry, niter_max=100):

    gauss = lambda p, x: exp(-(x-p[0])**2/(2*p[1]**2))

    ps = plasmastate('ips',1)
    ps.read(f_state)

    geq = readg(f_eqdsk)
    r0  = geq["rzero" ]
    b0  = abs(geq["bcentr"])
    ip  = geq['cpasma']

    nrho = len(ps["rho"])
    rho  = ps["rho"][:]
    vol  = ps["vol" ][:]
    area = ps["area"][:]
    ipol = ps["g_eq"][:]/(r0*b0)

    j_bs = 1.e-6*ps.dump_j_parallel(rho,"rho","curr_bootstrap", r0, b0)
    j_nb = 1.e-6*ps.dump_j_parallel(rho,"rho_nbi","curbeam", r0, b0)

    jaxis = 0.5
    jpeak = 1.0
    jaxis1 = 0.0
    jbdry1 = -1.0
    rho_jpeak = 0.7
    # rho_jbdry = 0.85

    jbdry = j_bs[where( rho>=rho_jbdry)[0][0]]

    jpeak_min = 0.1
    jpeak_max = 10.0

    j_tot = zeros(nrho)
    curt = zeros(nrho)

    for iter_scale in range(niter_max):

        x = 0.5*(jpeak_min+jpeak_max)

        j_tot = zeros(nrho)
        for i in range(nrho):
            if rho[i]<rho_jbdry:
                j_tot[i] = gauss_asym([x,rho_j,j_out,j_in],[rho[i]])[0] - gauss_asym([x,rho_j,j_out,j_in],[rho_jbdry])[0] + jbdry
            else:
                j_tot[i] = j_bs[i]

        curt = zeros(nrho)
        for i in range(nrho-1):
            dV = vol[i+1]-vol[i]
            jparm = 0.5*(j_tot[i]+j_tot[i+1])
            ipolm = 0.5*(ipol[i]+ipol[i+1])
            curt[i+1] = (curt[i]/ipol[i]+jparm*dV/(2.0*pi*r0*ipolm**2))*ipol[i+1]

        if abs(curt[-1]-ip*1.e-6) < 0.001: break

        if curt[-1] > ip*1.e-6:
           jpeak_max = x
        else:
           jpeak_min = x

    print('jpeak =',x)
    print('ip_cal = ',curt[-1])
    ps.load_j_parallel(rho, 1.0e6*j_tot, "rho_eq", "curt", r0, b0, tot=True)

    ps.store(f_state)

def broadJ_spline(f_state, f_eqdsk, rho_j, jaxis, rho_jbdry, niter_max=100):

    gauss = lambda p, x: exp(-(x-p[0])**2/(2*p[1]**2))

    ps = plasmastate('ips',1)
    ps.read(f_state)

    geq = readg(f_eqdsk)
    r0  = geq["rzero" ]
    b0  = abs(geq["bcentr"])
    ip  = geq['cpasma']

    nrho = len(ps["rho"])
    rho  = ps["rho"][:]
    vol  = ps["vol" ][:]
    area = ps["area"][:]
    ipol = ps["g_eq"][:]/(r0*b0)

    j_bs = 1.e-6*ps.dump_j_parallel(rho, "rho", "curr_bootstrap", r0, b0)
    j_nb = 1.e-6*ps.dump_j_parallel(rho, "rho_nbi", "curbeam", r0, b0)

    # jaxis = 0.5
    jaxis1 = 0.0
    jbdry1 = -1.0
    rho_jpeak = rho_j

    jbdry = j_bs[where( rho>=rho_jbdry)[0][0]]

    jpeak_min = 0.1
    jpeak_max = 10.0

    j_tot = zeros(nrho)
    curt = zeros(nrho)

    for iter_scale in range(niter_max):

        x = 0.5*(jpeak_min+jpeak_max)

        jaxis1 = x*(1.0 - jaxis)/rho_jpeak

        # x0 = [0.0     , 0.5*rho_jpeak, rho_jpeak, rho_jbdry]
        # #y0 = [jaxis*x , 0.6*x        , x        , jbdry    ]
        # y0 = [jaxis*x , 0.5*(1.+jaxis)*x    , x        , jbdry    ]
        # yp = [jaxis1  , 2.0*jaxis1   , 0.0      , jbdry1   ]
        # spl = profile_spline(x=x0,y=y0,yp=yp)

        x0 = [0.0     , rho_jpeak, rho_jbdry]
        y0 = [jaxis*x ,      x        , jbdry    ]
        yp = [jaxis1  ,    0.0      , jbdry1   ]
        spl = profile_spline(x=x0,y=y0,yp=yp)

        j_tot = zeros(nrho)
        for i in range(nrho):
            if rho[i]<rho_jbdry: j_tot[i] = spl(rho[i])
            else: j_tot[i] = j_bs[i]

        curt = zeros(nrho)
        for i in range(nrho-1):
            dV = vol[i+1]-vol[i]
            jparm = 0.5*(j_tot[i]+j_tot[i+1])
            ipolm = 0.5*(ipol[i]+ipol[i+1])
            curt[i+1] = (curt[i]/ipol[i]+jparm*dV/(2.0*pi*r0*ipolm**2))*ipol[i+1]

        if abs(curt[-1]-ip*1.e-6) < 0.001: break

        if curt[-1] > ip*1.e-6:
           jpeak_max = x
        else:
           jpeak_min = x

    print('jpeak =', x)
    print('ip_cal = ',curt[-1])
    ps.load_j_parallel(rho, 1.0e6*j_tot, "rho_eq", "curt", r0, b0, tot=True)

    ps.store(f_state)
