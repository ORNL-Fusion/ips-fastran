#!/usr/bin/env python

"""
 -----------------------------------------------------------------------
 plasma state utility
 JM
 -----------------------------------------------------------------------
"""

import os
from numpy import *
import netCDF4
from Namelist import Namelist

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
            "bs":"curr_bootstrap",
            "oh":"curr_ohmic",
            "nb":"curbeam",
            "ec":"curech",
            "ic":"curich" }
    
    def info(self):
        pass

    def close(self):

        self.data.close()

    def __getitem__(self, key):

        if key in self.data.variables: return self.data.variables[key]
        else: return []

    def put_sum_profile(self,prof,key,mode,sum=0.0):

        zone = self[mode]
        self[key][0] = 0.0 
        for i in range(1,self.nrho):
            self[key][i] = self[key][i-1]+prof[i]*(zone[i]-zone[i-1])
        print self[key][:]
        if sum>0.0:
            self[key][:] *= sum/self[key][-1]
        print self[key][:]

    def load_profile(self,prof,key,mode,sum=0.0):

        zone = self[mode]
        tmp = 0.0
        for i in range(self.nrho-1):
            self[key][i] = 0.5*(prof[i+1]+prof[i])*(zone[i+1]-zone[i])
            tmp += self[key][i]
        if sum>0.0:
            self[key][:] *= sum/tmp
        print 'integrated : ',tmp

    def dump_profile(self,key,mode,k=-1):

        zone = self[mode][:]
        if k>=0: 
            prof = self[key][k][:]
        else:
            prof = self[key][:]
        rvec = zeros(self.nrho-1)
        for i in range(self.nrho-1):
            rvec[i] = prof[i]/(zone[i+1]-zone[i])
        return self.cell2node(rvec)

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

    def load_j_parallel_CD(self,jpar,mode):

        curt = zeros(self.nrho)

        r0   = self.r0
        b0   = self.b0
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

        print "I%s=%5.3e"%(mode,1.0e-6*curt[-1])+"MA"

        for i in range(self.nrho-1):
            self[self.map[mode]][i] = curt[i+1]-curt[i]

    def dump_j_parallel_CD(self,mode):
        
        curt = zeros(self.nrho)
        jpar = self[self.map[mode]][:] 

        r0 = self.r0
        b0 = self.b0
        rho  = self["rho"][:]
        ipol = self["g_eq"][:]/(r0*b0)
        vol  = self["vol"][:]
        area = self["area"][:]

        curt[0] = 0.0
        for i in range(self.nrho-1): 
            curt[i+1] = curt[i]+jpar[i]

        for i in range(self.nrho-1):
            dV = vol[i+1]-vol[i]
            ipolm = 0.5*(ipol[i]+ipol[i+1])
            jpar[i] = 2.0*pi*r0*ipolm**2/dV*(curt[i+1]/ipol[i+1]-curt[i]/ipol[i])

        return self.cell2node(jpar)

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
       #node[-1] = cell[-1]
        node[-1] = 2.0*cell[-1]-node[-2]
        return node

#    def node2cell(self,node):
#    
#        nrho = len(node)
#        cell = zeros(nrho-1)
#        cell[0] = node[1]
#        for i in range(1,nrho-1):
#            cell[i] = 0.5*(node[i]+node[i+1])
#        return cell
#
#    def cell2node(self,cell):
#    
#        nrho = len(cell)+1
#        node = zeros(nrho)
#        node[0] = cell[0]
#        for i in range(nrho-1):
#            node[i+1] = 2.0*cell[i]-node[i]
#        return node


def instate2ps(instate,ps):

    #------------------------------------------------------------------
    # from instate

    nrho  = instate["nrho" ][0]

    n_ion = instate["n_ion"][0]
    z_ion = instate["z_ion"]
    a_ion = instate["a_ion"]
    f_ion = instate["f_ion"]

    n_imp = instate["n_imp"][0]
    z_imp = instate["z_imp"]
    a_imp = instate["a_imp"]
    f_imp = instate["f_imp"]

    ne = array(instate["ne"])
    te = array(instate["te"])
    ti = array(instate["ti"])
    omega = array(instate["omega"])
    zeff = array(instate["zeff"])
    density_beam = array(instate["density_beam"])

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

    #------------------------------------------------------------------
    # current

    j_tot = 1.0e6*array(instate["j_tot"])

    ps.load_j_parallel(j_tot)

    j_nb  = 1.0e6*array(instate["j_nb"])
    j_ec  = 1.0e6*array(instate["j_ec"])
    j_ic  = 1.0e6*array(instate["j_ic"])
    j_bs  = 1.0e6*array(instate["j_bs"])
    j_oh  = 1.0e6*array(instate["j_oh"])

    ps.load_j_parallel_CD(j_nb,"nb")
    ps.load_j_parallel_CD(j_ec,"ec")
    ps.load_j_parallel_CD(j_ic,"ic")
    ps.load_j_parallel_CD(j_bs,"bs")
    ps.load_j_parallel_CD(j_oh,"oh")

    #------------------------------------------------------------------
    # heating


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

