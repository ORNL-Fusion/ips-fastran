#!/usr/bin/env python

import sys,os,glob,shutil,cPickle
from numpy import *
import netCDF4

from Namelist import Namelist
from zinterp import zinterp
import scipy.optimize

from toq_io import read_peddata_profile
from zmodelprof import profile_pedestal
from zprofile import zprofile
#from plasmastate import plasmastate

def get_ni(ne,zeff,zmain=1.,zimp=6.,nbfast=0.,nalpha=0.,nhe=0.):

    tmp = zmain*zimp*(zimp-zmain)
    ni = (zimp**2*(ne-nbfast-2.0*(nalpha+nhe))-zimp*(zeff*ne-nbfast-4.0*(nalpha+nhe)))/tmp
    nz = (ne*(zeff-1.0)-2.0*(nalpha+nhe))/tmp
    return ni, nz

class eped_io():

    def __init__(self):

        pass

    def read(self,rdir,fn_eped_state,fn_peddata='',dskbal=False):

        #-- from eped.nc

        f_ps = os.path.join(rdir,fn_eped_state)
        nc = netCDF4.Dataset(f_ps,'r',format='NETCDF4')

        gamma_pb = nc.variables["gamma_PB"][:,:]
        gamma = array([ max(max(g),0.001) for g in gamma_pb ])

        k_start = 0
        k_EPED = k_start+where(gamma[k_start:] > 1.0)[0][0]

        #k_EPED = nc.variables["k_EPED"][0]
        wped = nc.variables["eq_wped_rho"][k_EPED]
        tped = nc.variables["eq_tped"][k_EPED]
        pped = nc.variables["eq_pped"][k_EPED]*1.0e3
        neped = nc.variables["neped"][0]
        nesep = nc.variables["nesep"][0]
        zeffped = nc.variables["zeffped"][0]
        ptop = nc.variables["eq_ptop"][k_EPED]*1.0e3
        ttop = nc.variables["eq_ttop"][k_EPED]

        betanped = nc.variables["eq_betanped"][k_EPED]
        betantop = nc.variables["eq_bedatntop"][k_EPED]
        nustar = nc.variables["eq_nustar"][k_EPED]

        #-- from peddata

        if fn_peddata:

            peddata = read_peddata_profile(os.path.join(rdir,fn_peddata))

            nrho_eped = len(peddata["rho"])

            rho_eped = [ peddata["rho"][k] for k in range(nrho_eped-1,-1,-1) ]
            ne_eped = [ peddata["ne"][k]  for k in range(nrho_eped-1,-1,-1) ]
            te_eped = [ peddata["Te"][k]  for k in range(nrho_eped-1,-1,-1) ]
            ptot_eped = [ peddata["ptot"][k]  for k in range(nrho_eped-1,-1,-1) ]

            ne_eped_spl = zinterp(rho_eped,ne_eped)
            te_eped_spl = zinterp(rho_eped,te_eped)
            ptot_eped_spl = zinterp(rho_eped,ptot_eped)

            rho_top = 1.0-1.5*wped
            ne_top = ne_eped_spl(rho_top)
            ptot_top = ptot_eped_spl(rho_top)

            ni_top, nz_top = get_ni(ne_top,zeff=zeffped)
            te_top = 0.5*ptot_top/(1.602*ne_top)
            ti_top = 0.5*ptot_top/(1.602*(ni_top+nz_top))

            ni_ped, nz_ped = get_ni(neped,zeff=zeffped)
            te_ped = 0.5*pped/(1.602*neped)
            ti_ped = 0.5*pped/(1.602*(ni_ped+nz_ped))

        else:

            rho_top = 1.0-1.5*wped

            ptot_top = ptop
            te_top = ttop
            ne_top = ptop/(2*1.602*ttop)

            ni_top, nz_top = get_ni(ne_top,zeff=zeffped)
            te_top = 0.5*ptot_top/(1.602*ne_top)
            ti_top = 0.5*ptot_top/(1.602*(ni_top+nz_top))

            ni_ped, nz_ped = get_ni(neped,zeff=zeffped)
            te_ped = 0.5*pped/(1.602*neped)
            ti_ped = 0.5*pped/(1.602*(ni_ped+nz_ped))

            rho_eped = []
            ne_eped = []
            te_eped = []

        print "TOP, PED"
        print 'rho  = %6.3f %6.3f'%(rho_top,1.0-wped)
        print 'ptot = %6.3f %6.3f'%( ptot_top, pped )
        print 'ne   = %6.3f %6.3f'%( ne_top, neped )
        print 'ni   = %6.3f %6.3f'%( ni_top, ni_ped )
        print 'nz   = %6.3f %6.3f'%( nz_top, nz_ped )
        print 'te   = %6.3f %6.3f'%( te_top, te_ped )
        print 'ti   = %6.3f %6.3f'%( ti_top, ti_ped )

        self.wped     = wped
        self.neped    = neped
        self.nesep    = nesep
        self.netop    = ne_top
        self.niped    = ni_ped
        self.nitop    = ni_top
        self.teped    = te_ped
        self.tetop    = te_top
        self.tiped    = ti_ped
        self.titop    = ti_top
        self.ptotped  = pped
        self.ptottop  = ptot_top
        self.zeffped  = zeffped
        self.betanped = betanped
        self.betantop = betantop
        self.nustar = nustar

        self.ip       = nc.variables["ip"     ][0]
        self.bt       = nc.variables["bt"     ][0]
        self.r        = nc.variables["r"      ][0]
        self.a        = nc.variables["a"      ][0]
        self.kappa    = nc.variables["kappa"  ][0]
        self.delta    = nc.variables["delta"  ][0]
        self.zeta     = nc.variables["zeta"   ][0]
        self.neped    = nc.variables["neped"  ][0]
        self.betan    = nc.variables["betan"  ][0]
        self.zeffped  = nc.variables["zeffped"][0]
        self.m        = nc.variables["m"      ][0]
        self.z        = nc.variables["z"      ][0]
        self.mi       = nc.variables["mi"     ][0]
        self.zi       = nc.variables["zi"     ][0]
        self.rho_eped = rho_eped
        self.ne_eped  = ne_eped
        self.te_eped  = te_eped

        if dskbal:
            self.q95 = nc.variables["eq_q95"][k_EPED]
            self.li = nc.variables["eq_li"][k_EPED]
            self.betap = nc.variables["eq_betap"][k_EPED]

    def dump_param(self):

        return {
            "wped"    : self.wped    ,
            "neped"   : self.neped   ,
            "netop"   : self.netop   ,
            "niped"   : self.niped   ,
            "nitop"   : self.nitop   ,
            "teped"   : self.teped   ,
            "tetop"   : self.tetop   ,
            "tiped"   : self.tiped   ,
            "titop"   : self.titop   ,
            "ptotped" : self.ptotped ,
            "ptottop" : self.ptottop ,
            "zeffped" : self.zeffped ,
            "ip"      : self.ip      ,
            "bt"      : self.bt      ,
            "r"       : self.r       ,
            "a"       : self.a       ,
            "kappa"   : self.kappa   ,
            "delta"   : self.delta   ,
            "zeta"    : self.zeta    ,
            "neped"   : self.neped   ,
            "betan"   : self.betan   ,
            "zeffped" : self.zeffped ,
            "m"       : self.m       ,
            "z"       : self.z       ,
            "mi"      : self.mi      ,
            "zi"      : self.zi      ,
            "rho_eped": self.rho_eped,
            "ne_eped" : self.ne_eped ,
            "te_eped" : self.te_eped ,
        }

    def make_profile(self, nrho = 101, ne0 = -1, te0 = -1, tesep = 0.075, ti0 = -1, tisep = 0.075, alpha = 1.5, beta = 1.5):

        wped  = self.wped
        neped = self.neped
        netop = self.netop
        teped = self.teped
        tetop = self.tetop
        tiped = self.tiped
        titop = self.titop
        nesep = self.nesep

        if nesep < 0: nesep = 0.25*neped

        rho = arange(nrho)/(nrho-1.)

        if ne0 < 0:
            ne = profile_pedestal(nrho,1.0-0.5*wped,wped,neped,nesep,ne0,alpha,beta,ytop=netop)(rho)
        else:
            ne = profile_pedestal(nrho,1.0-0.5*wped,wped,neped,nesep,ne0,alpha,beta)(rho)
            ne_spl = zinterp(rho,ne)
            self.netop = ne_spl(1.0-1.5*wped)

        self.ne_model = ne

        if te0 < 0:
            te = profile_pedestal(nrho,1.0-0.5*wped,wped,teped,tesep,te0,alpha,beta,ytop=tetop)(rho)
        else:
            te = profile_pedestal(nrho,1.0-0.5*wped,wped,teped,tesep,te0,alpha,beta)(rho)
            te_spl = zinterp(rho,te)
            self.tetop = te_spl(1.0-1.5*wped)

        if ti0 < 0:
            ti = profile_pedestal(nrho,1.0-0.5*wped,wped,tiped,tisep,ti0,alpha,beta,ytop=titop)(rho)
        else:
            ti = profile_pedestal(nrho,1.0-0.5*wped,wped,tiped,tisep,ti0,alpha,beta)(rho)
            ti_spl = zinterp(rho,ti)
            self.titop = ti_spl(1.0-1.5*wped)

        self.ne = ne
        self.te = te
        self.ti = ti

    def patch_pedestal(self, rho, nein, tein, tiin):

        nrho = len(rho)

        ne = zeros(nrho)
        te = zeros(nrho)
        ti = zeros(nrho)

        xtop = 1.0 - 1.5*self.wped

        nein_spl = zinterp(rho,nein)
        tein_spl = zinterp(rho,tein)
        tiin_spl = zinterp(rho,tiin)

        nein_top = nein_spl(xtop)
        tein_top = tein_spl(xtop)
        tiin_top = tiin_spl(xtop)

        print self.netop - nein_top


        for i in range(nrho):
            if rho[i] < xtop:
                ne[i] = nein[i] + self.netop - nein_top
                te[i] = tein[i] + self.tetop - tein_top
                ti[i] = tiin[i] + self.titop - tiin_top
            else:
                ne[i] = self.ne[i]
                te[i] = self.te[i]
                ti[i] = self.ti[i]

        return ne, te, ti

    def write(self,inp):

        rb = array(inp["rb"])
        zb = array(inp["zb"])

        R = 0.5*( max(rb) + min(rb) )
        Z = 0.5*( max(zb) + min(zb) )
        a = 0.5*( max(rb) - min(rb) )
        kappa = 0.5*( ( max(zb) - min(zb) )/ a )
        delta_u = ( R - rb[argmax(zb)] )/a
        delta_l = ( R - rb[argmin(zb)] )/a
        delta = 0.5*( delta_u + delta_l )

        bt = inp["r0"]*inp["b0"]/R
        betan = inp["betan"]
        ip = inp["ip"]*1.0e-6

        print R, a, kappa, delta, bt, ip

        rho = array(inp["rho"])
        ne = array(inp["ne"])
        te = array(inp["te"])
        zeff = array(inp["zeff"])

        te_fit = zprofile(rho,te,"mtanh_parab")
        te_fit.fit()

        print "xmid = ", te_fit.xmid
        print "xwid = ", te_fit.xwid
        xped = te_fit.xmid-0.5*te_fit.xwid

        ne_spl = zinterp(rho,ne)
        zeff_spl = zinterp(rho,zeff)

        neped = ne_spl(xped)
        zeffped = zeff_spl(xped)

        print 'neped = ',neped
        print 'zeffped = ',zeffped

        eped_input = Namelist()

        eped_input["eped_input"]["num_scan"] = [1]
        eped_input["eped_input"]["shot"    ] = [0]
        eped_input["eped_input"]["timeid"  ] = [0]
        eped_input["eped_input"]["runid"   ] = [1]
        eped_input["eped_input"]["ip"      ] = [ip]
        eped_input["eped_input"]["bt"      ] = [bt]
        eped_input["eped_input"]["r"       ] = [R]
        eped_input["eped_input"]["a"       ] = [a]
        eped_input["eped_input"]["kappa"   ] = [kappa]
        eped_input["eped_input"]["delta"   ] = [delta]
        eped_input["eped_input"]["zeta"    ] = [0.0]
        eped_input["eped_input"]["neped"   ] = [neped]
        eped_input["eped_input"]["betan"   ] = [betan]
        eped_input["eped_input"]["zeffped" ] = [zeffped]
        eped_input["eped_input"]["m"       ] = [2]
        eped_input["eped_input"]["z"       ] = [1]
        eped_input["eped_input"]["mi"      ] = [12]
        eped_input["eped_input"]["zi"      ] = [6]
        eped_input["eped_input"]["teped"   ] = [-1.0]
        eped_input["eped_input"]["tewid"   ] = [-1.0]
        eped_input["eped_input"]["ptotped" ] = [-1.0]
        eped_input["eped_input"]["ptotwid" ] = [-1.0]

        eped_input.write("eped.input")

    """
    def write_from_state(self,f_state,f_geqdsk,f_inbc,f_eped_state):

        # entry

        print "\n"+72*"="
        print '= make eped input '

        ps = plasmastate('ips',1)
        ps.read(f_state)

        geq = readg(f_eqdsk)
        r0  = geq["rzero" ]
        b0  = abs(geq["bcentr"])
        ip  = geq['cpasma']

        # calculate betan

        vol  = ps["vol" ][:]
        area = ps["area"][:]
        ipol = ps["g_eq"][:]/(r0*b0)
        rho  = ps["rho"][:]
        nrho = len(rho)
        ne   = ps["ns"][0,:]*1.0e-19
        ni   = ps["ni"][:]*1.0e-19
        te   = ps["Ts"][0,:]
        ti   = ps["Ti"][:]
        a0   = ps["rMinor_mean"][-1]

        wth  = 1.5*1.602e3*(ne*te+ni*ti)
        wth_sum = sum([(vol[i+1]-vol[i])*wth[i] for i in range(nrho-1)])

        density_beam = ps.dump_profile(rho,"rho_nbi","nbeami",k=0)*1.e-19
        wbeam = ps.dump_profile(rho,"rho_nbi","eperp_beami",k=0) \
              + ps.dump_profile(rho,"rho_nbi","epll_beami" ,k=0)
        wbeam = density_beam*wbeam*1.602e3
        wbeam_sum = sum([(vol[i+1]-vol[i])*0.5*(wbeam[i+1]-wbeam[i]) for i in range(nrho-1)])

        betan_th = wth_sum/vol[-1]
        betan_th*= 2.0*4.0*pi*1.0e-7/b0**2
        betan_th/= 1.5
        betan_th/= fabs(ip*1.0e-6)/(a0*abs(b0))
        betan_th*= 1.0e2

        betan_beam = wbeam_sum/vol[-1]
        betan_beam*= 2.0*4.0*pi*1.0e-7/b0**2
        betan_beam/= 1.5
        betan_beam/= fabs(ip*1.0e-6)/(a0*abs(b0))
        betan_beam*= 1.0e2

        print 'betan_th   =', betan_th
        print 'betan_beam =', betan_beam

        # update eped_state

        eped = eped_state()
        eped.load(f_eped_state)

        rmajor  = ps["Rmajor_mean"][:]
        rminor  = ps["rMinor_mean"][:]
        kappa   = ps["elong"][:]
        delta   = ps["triang"][:]

        eped["ip"      ][0] = ip*1.0e-6
        eped["bt"      ][0] = b0
        eped["r"       ][0] = rmajor[-1]
        eped["a"       ][0] = rminor[-1]
        eped["kappa"   ][0] = kappa[-1]
        eped["delta"   ][0] = delta[-1]
        eped["zeta"    ][0] = 0.0
       #eped["neped"   ][0] = 6.0
        eped["betan"   ][0] = betan_th+betan_beam
       #eped["zeffped" ][0] = 1.6

        eped.close()
    """

if __name__ == "__main__":

    pass

    #from zplotbase import *

    #NEMAX, TEMAX, TIMAX = 5., 6., 6.

    #rdir = "eped01"
    #fn_eped_state = "eped_state.nc"
    #fn_peddata = "peddata_26"
    #print fn_eped_state, fn_peddata
    #eped_data = eped_io()
    #eped_data.read(rdir,fn_eped_state,fn_peddata)
    #eped_data.make_profile()
    #wped = eped_data.wped

    #instate = Namelist("instate")

    #nein = array(instate["instate"]["ne"])
    #tein = array(instate["instate"]["te"])
    #tiin = array(instate["instate"]["ti"])

    #ne,te,ti = eped_data.patch_pedestal(nein,tein,tiin)
