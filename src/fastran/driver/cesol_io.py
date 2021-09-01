"""
 -----------------------------------------------------------------------
 cesol driver io
 -----------------------------------------------------------------------
"""

from numpy import *
import netCDF4

from Namelist import Namelist
from plasmastate import plasmastate
import eped_state
import eped_io
from efit_eqdsk import readg
from zinterp import zinterp

def eped_to_plasmastate(f_state, f_eped, f_instate):
    print('EPED TO PLASMA STATE')

    #-- plasma state
    ps = plasmastate('ips',1)
    ps.read(f_state)

    rho   = ps["rho"][:]
    nrho  = len(rho)
    ne    = ps["ns"][0,:]*1.0e-19
    te    = ps["Ts"][0,:]
    ti    = ps["Ti"][:]
    zeff  = ps["Zeff"][:]
    omega = ps["omegat"][:]
    ne    = ps.cell2node_bdry(ne)
    te    = ps.cell2node_bdry(te)
    ti    = ps.cell2node_bdry(ti)
    zeff  = ps.cell2node_bdry(zeff)
    omega = ps.cell2node_bdry(omega)

    #-- instate
    instate = Namelist(f_instate)
    tesep = instate["flux"]["tesep"][0]
    tisep = instate["flux"]["tisep"][0]

    #-- eped
    eped_data = eped_io.eped_io()
    eped_data.read(".",f_eped)

    print('check:', eped_data.neped, eped_data.nesep)

    eped_data.make_profile(
        nrho = nrho,
        ne0 = ne[0],
        te0 = te[0], tesep = tesep*1.e-3,
        ti0 = ti[0], tisep = tisep*1.e-3,
        alpha = 1.5, beta = 1.5
    )

    ne_update, te_update, ti_update = eped_data.patch_pedestal(rho, ne, te, ti)

    ne_update = eped_data.ne  #----------
    ps["ns"][0] = 1.0e19*ps.node2cell(ne_update)
    ps["Ts"][0] = ps.node2cell(te_update)
    nspec_th = len(ps["Ts"])-1
    print('nspec_th =',nspec_th)
    for k in range(nspec_th):
        ps["Ts"][k+1] = ps.node2cell(ti_update)
    ps["Ti"][:] = ps.node2cell(ti_update)

    ps.update_particle_balance()

    ps.store(f_state)

    check = Namelist()
    check["data"]["ne"] = ne
    check["data"]["ne_update"] = ne_update
    check["data"]["te"] = te
    check["data"]["te_update"] = te_update
    check["data"]["ti"] = te
    check["data"]["ti_update"] = ti_update
    check["data"]["te_ps"] = ps["Ts"][0]
    check["data"]["ti_ps"] = ps["Ts"][1]
    check.write("check.dat")

def plasmastate_to_eped(f_state, f_eqdsk, f_eped, f_instate):
    #-- plasmastate
    ps = plasmastate('ips',1)
    ps.read(f_state)

    rho   = ps["rho"][:]
    nrho  = len(rho)
    ne    = ps["ns"][0,:]*1.0e-19
    te    = ps["Ts"][0,:]
    ti    = ps["Ti"][:]
    zeff  = ps["Zeff"][:]
    ne    = ps.cell2node_bdry(ne)
    te    = ps.cell2node_bdry(te)
    ti    = ps.cell2node_bdry(ti)
    zeff  = ps.cell2node_bdry(zeff)

    #-- geqdsk
    geq = readg(f_eqdsk)
    r0  = geq["rzero" ]
    b0  = abs(geq["bcentr"])
    ip  = geq['cpasma']*1.0e-6
    rb  = geq['rbdry']
    zb  = geq['zbdry']

    #-- instate
    instate = Namelist(f_instate)
    betan = instate["flux"]["betan"][0]
    neped = instate["flux"]["neped"][0]
    nesep = instate["flux"]["nesep"][0]

    #--
    # eped_data = eped_io.eped_io()
    # eped_data.read(".",f_eped)
    # zeffped = zinterp(rho,zeff)(1.0-eped_data.wped)
    # print 'zeff: ', eped_data.zeffped, zeffped

    #-- eped
    eped = eped_state.eped_state()
    eped.load(f_eped)

    R = 0.5*( max(rb) + min(rb) )
    Z = 0.5*( max(zb) + min(zb) )
    a = 0.5*( max(rb) - min(rb) )
    kappa = 0.5*( ( max(zb) - min(zb) )/ a )
    delta_u = ( R - rb[argmax(zb)] )/a
    delta_l = ( R - rb[argmin(zb)] )/a
    delta = 0.5*( delta_u + delta_l )

    bt = r0*b0/R
    print(R, a, kappa, delta, bt, ip)

    eped["ip"      ][0] = ip
    eped["bt"      ][0] = bt
    eped["r"       ][0] = R
    eped["a"       ][0] = a
    eped["kappa"   ][0] = kappa
    eped["delta"   ][0] = delta
    eped["zeta"    ][0] = 0.0
    eped["neped"   ][0] = neped
    eped["nesep"   ][0] = nesep
    eped["betan"   ][0] = betan
    # eped["zeffped" ][0] = zeffped

    eped.close()

def fastran_to_instate(f_fastran, f_instate):
    print('f_fastran =',f_fastran)

    fastran = netCDF4.Dataset(f_fastran,'r',format='NETCDF4')
    betan = fastran.variables["betan"][:][-1]
    pnbe = fastran.variables["pnbe"][:][-1]
    pnbi = fastran.variables["pnbi"][:][-1]
    prfe = fastran.variables["prfe"][:][-1]
    prfi = fastran.variables["prfi"][:][-1]
    sn = fastran.variables["sn"][:][-1]
    pei  = 0.0 #fastran.variables["pei"][:][-1]

    b0 = fastran.variables["b0"][-1]
    a0 = fastran.variables["a0"][-1]
    ip = fastran.variables["ip"][-1]
    rho = fastran.variables["rho"][:]
    rhob = fastran.variables["rhob"][:,:][-1]
    volp = fastran.variables["volp"][:,:][-1]
    ne = fastran.variables["ne"][:,:][-1]

    drho = rhob[1]-rhob[0]
    nrho = len(rho)
    vol = zeros(nrho)
    for k in range(1,nrho):
        vol[k] = vol[k-1]+drho*(volp[k-1]+volp[k])*0.5
    def integrate(f):
        fsum = sum([(vol[i+1]-vol[i])*f[i] for i in range(nrho-1)])
        return fsum
    nesum = integrate(ne)

    tauth = fastran.variables["tauth"][-1]
    if isnan(tauth): tauth = 0.0

    instate = Namelist(f_instate)

    instate["flux"]["betan"] = [betan]
    instate["flux"]["pnbe"] = [pnbe]
    instate["flux"]["pnbi"] = [pnbi]
    instate["flux"]["prfe"] = [prfe]
    instate["flux"]["prfi"] = [prfi]
    instate["flux"]["pei" ] = [pei ]
    instate["flux"]["sn" ] = [sn*1.0e19 ]

    taup = instate["flux"]["taup"][0] * tauth

    stot =  instate["flux"]["puff"][0] + instate["flux"]["particle_core_src"][0] + instate["flux"]["sn"][0]
#   stot /= 2.0*pi

    nesum_predict = stot*1.0e-19*taup

    instate["flux"]["nesum"] = [ nesum ]
    instate["flux"]["nesum_predict"] = [ nesum_predict ]
    instate["flux"]["tauth"] = [ tauth  ]

    instate.write(f_instate)
