#!/usr/bin/env python

"""
 -----------------------------------------------------------------------
 plasma state utility
 JM
 -----------------------------------------------------------------------
"""

import os,subprocess,glob
from numpy import *
import netCDF4
from Namelist import Namelist
from zinterp import zinterp

class plasma_state_file():

    def __init__(self,f_ps,mode='r+',r0=0.0,b0=0.0,ip=0.0):

        self.data = netCDF4.Dataset(f_ps,mode,format='NETCDF4')

        self.nrho = len(self.data.dimensions["dim_nrho"])
        self.nrho_eq = len(self.data.dimensions["dim_nrho_eq"])
        self.nrho_eq_geo = len(self.data.dimensions["dim_nrho_eq_geo"])

        if self.nrho_eq != self.nrho:
            raise Exception('nrho_eq != nrho')
        if self.nrho_eq_geo != self.nrho:
            raise Exception('nrho_eq_geo != nrho')

        self.r0 = r0
        self.b0 = b0
        self.ip = ip

        self.map = { 
            "bs":["curr_bootstrap","rho"],
            "oh":["curr_ohmic","rho"],
            "nb":["curbeam","rho_nbi"],
            "ec":["curech","rho_ecrf"],
            "ic":["curich","rho_icrf"] }

        try:
            FASTRAN_ROOT = os.environ["FASTRAN_ROOT"]
            DIR_BIN = os.path.join(FASTRAN_ROOT,"bin")
            self.pstool_bin = os.path.join(DIR_BIN,"pstool")
        except:
            self.pstool_bin = 'pstool'

        self.f_ps = f_ps
        self.mode = mode

        self.geqdsk = ''
    
    def info(self):

        pass

    def close(self):

        self.data.close()
        #log=open("pstool.log","w")
        #retcode = subprocess.call([self.pstool_bin, "rehash",self.f_ps],
        #    stdout=log,stderr=log)
        #log.close()

    def __getitem__(self, key):

        if key in self.data.variables: return self.data.variables[key]
        else: return []

    def put_sum_profile(self,prof,key,mode,sum=0.0):

        zone = self[mode]
        self[key][0] = 0.0 
        for i in range(1,self.nrho):
            self[key][i] = self[key][i-1]+prof[i]*(zone[i]-zone[i-1])
        if sum>0.0:
            self[key][:] *= sum/self[key][-1]

    def load_profile(self,rho,prof,key,mode,k=-1,sum=0.0):

        zone_spl = zinterp(self["rho"][:],self[mode][:])
        prof_spl = zinterp(rho,prof)

        rho_key = self["rho"+self[key].dimensions[-1].split("rho")[-1]][:]

        if k<0:
            tmp = 0.0
            for i in range(len(rho_key)-1):
                self[key][i] = 0.5*(prof_spl(rho_key[i+1])+prof_spl(rho_key[i])) \
                                  *(zone_spl(rho_key[i+1])-zone_spl(rho_key[i]))
                tmp += self[key][i]
            if sum>0.0:
                self[key][:] *= sum/tmp
        else:
            tmp = 0.0
            for i in range(len(rho_key)-1):
                self[key][k][i] = 0.5*(prof_spl(rho_key[i+1])+prof_spl(rho_key[i])) \
                                  *(zone_spl(rho_key[i+1])-zone_spl(rho_key[i]))
                tmp += self[key][k][i]
            if sum>0.0:
                self[key][k][:] *= sum/tmp

        print 'integrated : ',tmp

    def dump_profile(self,rho,key,mode=None,k=-1):

        if mode in ["vol","area"]:

            zone_spl = zinterp(self["rho"][:],self[mode][:])
    
            rho_key = self["rho"+self[key].dimensions[-1].split("rho")[-1]][:]
            nrho_key = len(rho_key)
    
            if k>=0: 
                prof = self[key][k][:]
            else:
                prof = self[key][:]
    
            cell = zeros(nrho_key-1)
            for i in range(nrho_key-1):
                cell[i] = prof[i]/(zone_spl(rho_key[i+1])-zone_spl(rho_key[i]))
            node = self.cell2node(cell)
    
            rvec = zinterp(rho_key,node)(rho)

        else:

            rho_key = self["rho"+self[key].dimensions[-1].split("rho")[-1]][:]

            if k>=0: 
                prof = self[key][k][:]
            else:
                prof = self[key][:]

            rvec = zinterp(rho_key,self.cell2node(prof))(rho)

        return rvec

    def load_j_parallel(self,jpar):

        curt = zeros(self.nrho)

        r0   = self.r0
        b0   = self.b0
        ip   = self.ip
        rho  = self["rho"][:]
        ipol = self["g_eq"][:]/(r0*b0)
        vol  = self["vol"][:]
        area = self["area"][:]

        curt[0] = 0.0
        for i in range(self.nrho-1):
            dV = vol[i+1]-vol[i]
            jparm = 0.5*(jpar[i]+jpar[i+1])
            ipolm = 0.5*(ipol[i]+ipol[i+1])
            curt[i+1] = (curt[i]/ipol[i]+jparm*dV/(2.0*pi*r0*ipolm**2))*ipol[i+1]

        self["curt"][:] = curt[:]*ip/curt[-1]

    def dump_j_parallel(self):
        
        jpar = zeros(self.nrho-1)
        curt = self["curt"][:] 

        r0   = self.r0
        b0   = self.b0
        ip   = self.ip
        rho  = self["rho"][:]
        ipol = self["g_eq"][:]/(r0*b0)
        vol  = self["vol"][:]
        area = self["area"][:]

        for i in range(self.nrho-1):
            dV = vol[i+1]-vol[i]
            ipolm = 0.5*(ipol[i]+ipol[i+1])
            jpar[i] = 2.0*pi*r0*ipolm**2/dV*(curt[i+1]/ipol[i+1]-curt[i]/ipol[i])

        return self.cell2node(jpar)

    def load_j_parallel_CD(self,rhoin,jparin,key):

        xkey = self.map[key][1]
        ykey = self.map[key][0]

        rho = self[xkey][:]
        curt = zeros(len(rho))

        r0   = self.r0
        b0   = self.b0
        ipol = zinterp(self["rho"][:],self["g_eq"][:]/(r0*b0))(rho)
        vol  = zinterp(self["rho"][:],self["vol"][:])(rho)
        area = zinterp(self["rho"][:],self["area"][:])(rho)

        jpar = zinterp(rhoin,jparin)(rho)

        curt[0] = 0.0
        for i in range(len(rho)-1):
            dV = vol[i+1]-vol[i]
            jparm = 0.5*(jpar[i]+jpar[i+1])
            ipolm = 0.5*(ipol[i]+ipol[i+1])
            curt[i+1] = (curt[i]/ipol[i]+jparm*dV/(2.0*pi*r0*ipolm**2))*ipol[i+1]

        print "I%s=%5.3e"%(key,1.0e-6*curt[-1])+"MA"

        for i in range(len(rho)-1):
            self[ykey][i] = curt[i+1]-curt[i]

    def dump_j_parallel_CD(self,rhoin,key):
        
        ykey = self.map[key][0]
        xkey = self.map[key][1]

        rho = self[xkey][:]
        jpar = self[ykey][:]
        curt = zeros(len(rho))

        r0 = self.r0
        b0 = self.b0
        ipol = zinterp(self["rho"][:],self["g_eq"][:]/(r0*b0))(rho)
        vol  = zinterp(self["rho"][:],self["vol"][:])(rho)
        area = zinterp(self["rho"][:],self["area"][:])(rho)

        curt[0] = 0.0
        for i in range(len(rho)-1): 
            curt[i+1] = curt[i]+jpar[i]

        for i in range(len(rho)-1):
            dV = vol[i+1]-vol[i]
            ipolm = 0.5*(ipol[i]+ipol[i+1])
            jpar[i] = 2.0*pi*r0*ipolm**2/dV*(curt[i+1]/ipol[i+1]-curt[i]/ipol[i])

        return zinterp(rho,self.cell2node(jpar))(rhoin)

    def node2cell(self,node):
    
        nrho = len(node)
        cell = zeros(nrho-1)
        for i in range(nrho-1):
            cell[i] = 0.5*(node[i]+node[i+1])
        return cell

    def cell2node(self,cell):
    
        nrho = len(cell)+1
        node = zeros(nrho)
        node[0] = cell[0]
        for i in range(1,nrho-1):
            node[i] = 0.5*(cell[i-1]+cell[i])
        node[-1] = cell[-1]
        return node

    #def cell2node_bdry(self,cell):
    #
    #    nrho = len(cell)+1
    #    node = zeros(nrho)
    #    node[0] = cell[0]
    #    for i in range(1,nrho-1):
    #        node[i] = 0.5*(cell[i-1]+cell[i])
    #    node[-1] = 2.0*cell[-1]-node[-2] #<=====
    #    return node

    def cell2node_bdry(self,cell):
    
        nrho = len(cell)+1
        node = zeros(nrho)
        node[0] = cell[0]
        for i in range(1,nrho):
            node[i] = 2.0*cell[i-1]-node[i-1]
        return node

def instate2ps(instate,ps):

    #------------------------------------------------------------------
    # from instate

    for key in instate.keys(): instate[key] = array(instate[key])

    nrho  = instate["nrho"][0]
    rho   = instate["rho"]

    n_ion = instate["n_ion"][0]
    z_ion = instate["z_ion"]
    a_ion = instate["a_ion"]
    f_ion = instate["f_ion"]

    n_imp = instate["n_imp"][0]
    z_imp = instate["z_imp"]
    a_imp = instate["a_imp"]
    f_imp = instate["f_imp"]

    ne = instate["ne"]
    te = instate["te"]
    ti = instate["ti"]
    omega = instate["omega"]
    zeff = instate["zeff"]

    ne[0] = ne[1]
    te[0] = te[1]
    ti[0] = ti[1]
    omega[0] = omega[1]
    zeff[0] = zeff[1]

    #------------------------------------------------------------------
    # put zeros if not defined

    for key in [
        "j_oh", "j_bs", "j_nb", "j_ec", "j_ic", \
        "pe_nb", "pe_ec", "pe_ic", "pe_fus", "pe_ionization", "p_rad", \
        "pi_nb", "pi_ec", "pi_ic", "pi_fus", "pi_ionization", "pi_cx", "p_ohm", "p_ei", \
        "torque_nb", "torque_in", "se_nb", "se_ionization", "si_nb", "si_ionization", \
        "q", "psipol", \
        "density_beam", "wbeam", "density_alpha", "walpha", \
        "chie",
        "chii" ]:
        if key.upper() not in instate.keys(): instate[key] = zeros(nrho)


    density_beam = instate["density_beam"]
    density_alpha = instate["density_alpha"]
    wbeam = instate["wbeam"]

    rho = ps["rho"][:]

    #------------------------------------------------------------------
    # density

    a=0; b=0; c=0; d=0
    for k in range(n_imp):
        b = b+f_imp[k]*z_imp[k]
        d = d+f_imp[k]*z_imp[k]*z_imp[k]
    for k in range(n_ion):
        a = a+f_ion[k]*z_ion[k]
        c = c+f_ion[k]*z_ion[k]*z_ion[k]

    density_ion = {}
    for k in range(n_ion):
        density_ion[k] = zeros(nrho)

    density_imp = {}
    for k in range(n_imp):
        density_imp[k] = zeros(nrho)

    density_th = zeros(nrho)

    for i in range(nrho):

        zne_adj = ne[i]
        zzne_adj = ne[i]*zeff[i]

        # depletion due to beam ions

        zne_adj = zne_adj - 1.0*density_beam[i]
        zzne_adj = zzne_adj - 1.0**2*density_beam[i]

        # effective main ion and impurity densities

        nion = (zne_adj *d-zzne_adj*b)/(a*d-b*c)
        nimp = (zzne_adj*a-zne_adj *c)/(a*d-b*c)

        for k in range(n_ion):
            density_ion[k][i] = f_ion[k]*nion
        for k in range(n_imp):
            density_imp[k][i] = f_imp[k]*nimp

        for k in range(n_ion):
            density_th[i] = density_th[i] + density_ion[k][i]
        for k in range(n_imp):
            density_th[i] = density_th[i] + density_imp[k][i]

    ps["ns"][0] = 1.0e19*ps.node2cell(ne)

    for k in range(n_ion):
        ps["ns"][k+1] = 1.0e19*ps.node2cell(density_ion[k])

    for k in range(n_imp):
        ps["ns"][k+n_ion+1] = 1.0e19*ps.node2cell(density_imp[k])

    ps["ni"][:] = 1.0e19*ps.node2cell(density_th)

    #------------------------------------------------------------------
    # beam

    ps["nbeami"][0] = 1.0e19*ps.node2cell(density_beam)
    ps["eperp_beami"][0] = 2.0*20.0*ones(nrho-1)
    ps["epll_beami"][0] = 20.0*ones(nrho-1)

    tbeam = wbeam/1.602e-3/(density_beam+1.0e-6)
    ps["eperp_beami"][0] = ps.node2cell(2.0*tbeam/3.0)
    ps["epll_beami"][0] = ps.node2cell(tbeam/3.0)

    #------------------------------------------------------------------
    # temperature

    ps["Ts"][0] =ps.node2cell(te)

    for k in range(n_ion):
        ps["Ts"][k+1] = ps.node2cell(ti)
    for k in range(n_imp):
        ps["Ts"][k+n_ion+1] = ps.node2cell(ti)

    ps["Ti"][:] = ps.node2cell(ti)

    #------------------------------------------------------------------
    # zeff

    ps["Zeff"][:] = ps.node2cell(zeff)
    ps["Zeff_th"][:] = ps.node2cell(zeff)

    #------------------------------------------------------------------
    # rotation

    ps["omegat"][:] = ps.node2cell(omega)

    #--------------------------------------------------------------
    # current
    
    j_tot = 1.e6*instate["j_tot"]
    j_tot[0] = j_tot[1]
    ps.load_j_parallel(j_tot)
    
    for key in ["j_nb","j_ec","j_ic","j_bs","j_oh"]:
        if key.upper() not in instate.keys(): instate[key] = zeros(nrho)
    
    j_nb  = 1.e6*instate["j_nb"]
    j_ec  = 1.e6*instate["j_ec"]
    j_ic  = 1.e6*instate["j_ic"]
    j_bs  = 1.e6*instate["j_bs"]
    j_oh  = 1.e6*instate["j_oh"]
    
    ps.load_j_parallel_CD(rho,j_nb,"nb")
    ps.load_j_parallel_CD(rho,j_ec,"ec")
    ps.load_j_parallel_CD(rho,j_ic,"ic")
    ps.load_j_parallel_CD(rho,j_bs,"bs")
    ps.load_j_parallel_CD(rho,j_oh,"oh")
    
    #-------------------------------------------------------------
    # heating
    
    for key in ["pe_nb","pi_nb","pe_ec","pe_ic","pi_ic","pe_fus","pi_fus"]:
        if key.upper() not in instate.keys(): instate[key] = zeros(nrho)
    
    pe_nb  = 1.e6*instate["pe_nb" ]
    pi_nb  = 1.e6*instate["pi_nb" ]
    pe_ec  = 1.e6*instate["pe_ec" ]
    pe_ic  = 1.e6*instate["pe_ic" ]
    pi_ic  = 1.e6*instate["pi_ic" ]
    pe_fus = 1.e6*instate["pe_fus"]
    pi_fus = 1.e6*instate["pi_fus"]
    
    ps.load_profile (rho,pe_nb,"pbe","vol")
    ps.load_profile (rho,pi_nb,"pbi","vol")
    ps.load_profile (rho,pe_ec,"peech","vol")
    ps.load_profile (rho,pe_ic,"picrf_totals","vol",k=0)
    ps.load_profile (rho,pi_ic,"picth","vol")

    #------------------------------------------------------------------
    # temp

    ps["sc0"][:] = 2.0e21
    ps["n0norm"][:] = 1.0e-10 
    ps["T0sc0"][:] = 0.01 
    ps["sc0_to_sgas"][:] = 1

########################################################################
#  test

if __name__ == "__main__":

    #load_j_parallel(j_par)
    #dump_j_parallel()
    #load_j_parallel_CD(j_ec,"ec")
    #dump_j_parallel_CD("ec")
    
    instate = Namelist("instate")["instate"]
    r0 = instate["r0"][0]
    b0 = abs(instate["b0"][0])
    ip = 1.0e6*instate["ip"][0]
    j_tot = 1.0e6*array(instate["j_tot"])
    j_ec = 1.0e6*array(instate["j_tot"])
    
    ps = plasma_state_file("ps.nc",r0=r0,b0=b0,ip=ip)
    ps.info()
    
    ps.load_j_parallel(j_tot)
    
    jpar = ps.dump_j_parallel()
    
    for i in range(51):
        print i, j_tot[i]*1.0e-6,jpar[i]*1.0e-6
    
    ps.load_j_parallel_CD(j_tot,"ec")
    jec = ps.dump_j_parallel_CD("ec")

    for i in range(51):
        print i, j_tot[i]*1.0e-6,jec[i]*1.0e-6
    
    #ps["ns"][0]=1.0

    temp = ps["Ts"][0][:]
    print temp
    
    print ps["eperp_beami"][:] #(~nrho_nbi,nspec_beam)

    ps.close()

