#! /usr/bin/env python

"""
 -----------------------------------------------------------------------
 fastran plasma state file utility 
 -----------------------------------------------------------------------
"""

from plasmastate import *
from numpy import *
from Namelist import Namelist
from zinterp import zinterp

class zplasmastate (PlasmaState):

    def init_from_instate (self,fn="instate",runid='fastran',shot=0,t0=0.0,t1=0.0):

        instate = Namelist(fn)["instate"]

        nrho = instate["nrho"][0]
        nth  = 101

        n_ion = instate["n_ion"][0]
        n_imp = instate["n_imp"][0]

        z_spec = instate["z_ion"]+instate["z_imp"]
        a_spec = instate["a_ion"]+instate["a_imp"]

        nspec_th = instate["n_ion"][0]+instate["n_imp"][0]

        nspec_rfmin = instate["n_min"][0]
        z_rfmin = instate["z_min"]
        a_rfmin = instate["a_min"]

        nspec_beam = instate["n_beam"][0]
        z_beam = instate["z_beam"]
        a_beam = instate["a_beam"]

        nspec_fusion = instate["n_fusion"][0]
        z_fusion = instate["z_beam"]
        a_fusion = instate["a_beam"]

        self["nspec_th"] = nspec_th 
        self["nspec_rfmin"] = nspec_rfmin
        self["nrho"] =  nrho
        self["nrho_eq"] = nrho
        self["nth_eq"] = nth
        self["nrho_eq_geo"] = nrho
        self["nr"] = 0
        self["nz"] = 0

        self["nrho_gas"] = nrho
        self["nrho_nbi"] = nrho
        self["nrho_ecrf"] = nrho
        self["nrho_icrf"] = nrho
        self["nrho_fus"] = nrho
        self["nrho_anom"] = nrho

        self["global_label"] = runid
        self["tokamak_id"] = instate["tokamak_id"][0]
        self["runid"] = runid
        self["shot_number"] = shot

        self["t0"] = t0
        self["t1"] = t1

        self.alloc()

        self.setThermalSpecies(-1, -1, 0)

        for k in range(nspec_th):
            self.setThermalSpecies(z_spec[k], z_spec[k], a_spec[k])

        for k in range(nspec_rfmin):
            self.setRFMinoritySpecies(z_rfmin[k], z_rfmin[k], a_rfmin[k])

        for k in range(nspec_beam):
            self.setNeutralBeamSpecies(z_beam[k], z_beam[k], a_beam[k])

        for k in range(nspec_fusion):
            self.setFusionSpecies(z_fusion[k], z_fusion[k], a_fusion[k])

        self.finishSpecies()

        rho =linspace(0.0,1.0,num=nrho)
        self["rho"] = rho
        self["rho_eq"] = rho
        self["th_eq"] = linspace(0.0,2.0*pi,num=nth)
        self["rho_eq_geo"] = rho

        self["rho_gas"] = rho
        self["rho_nbi"] = rho
        self["rho_ecrf"] = rho
        self["rho_icrf"] = rho
        self["rho_fus"] = rho
        self["rho_anom"] = rho

    def load_innubeam (self, fn="innubeam"):

        innubeam = Namelist(fn)
        nbi = innubeam["nbi_config"]

        self["t0"] = 0.0
        self["t1"] = innubeam["nubeam_run"]["dt_nubeam"][0]

        nbeam = nbi["nbeam"][0]
        self["nbeam"] = nbeam
        self.alloc()

        self["nbi_src_name"] = ['NB%02d'%k for k in range(nbeam)]

        nbion_table = { (1,1):'H',(1,2):'D',(2,3):'He3',(2,4):'He4'}
        self["nbion"] = [ nbion_table[(int(nbi["xzbeama"][k]),int(nbi["abeama"][k]))] for k in range(nbeam) ]

        srtcen = 1.0e-2*array(nbi["rtcena"])
        for k in range(nbeam):
            if not nbi["nlco"]: srtcen[k] -= 1.
        self["srtcen"] = srtcen # (signed: + means ccw momentum inj.)

        self["lbsctan"] = 1.e-2*array(nbi["xlbtna"]) # dist., sce to tangency pt.
        self["zbsc"] = 1.e-2*array(nbi["xybsca"]) # height, sce above midplane
        self["phibsc"] =  nbeam*[0.0] #nbi["xbzeta"] # toroidal angle of sce

        self["lbscap"] = 1.e-2*array(nbi["xlbapa"]) # dist., sce to aperture
        self["zbap"] = 1.e-2*array(nbi["xybapa"]) # height, aperture center


        zdivcon = 2.*pi/(360.0*sqrt(2.0))
        self["nbshape"] = [ ('','rectangle','circle')[s] for s in nbi["nbshapa"] ] # shape of source

        self["b_halfwidth"] = 1.0e-2*array(nbi["bmwidra"]) # half-width
        self["b_halfheight"] = 1.0e-2*array(nbi["bmwidza"]) # half-height
        self["b_hfocal_length"] = 1.0e-2*array(nbi["foclra"]) # horiz. focal len
        self["b_vfocal_length"] = 1.0e-2*array(nbi["foclza"]) # vert. focal len
        self["b_hdivergence"] = array(nbi["divra"])/zdivcon # horiz. diverg.
        self["b_vdivergence"] = array(nbi["divza"])/zdivcon # vert.  diverg.

        self["nbap_shape"] = [ ('','rectangle','circle')[s] for s in nbi["nbapsha"] ] # shape of aperture
        self["ap_halfwidth"] = 1.e-2*array(nbi["rapedga"])
        self["ap_halfheight"] = 1.e-2*array(nbi["xzpedga"])
        self["ap_horiz_offset"] = 1.e-2*array(nbi["xrapoffa"])
        self["ap_vert_offset" ] = 1.e-2*array(nbi["xzapoffa"])

        self["nbap2_shape"] = [ ('none','rectangle','circle')[s] for s in nbi["nbapsh2"] ] # shape of 2nd aperture

        self["Lbscap2"] = 1.e-2*array(nbi["xlbapa2"])
        self["ap2_halfwidth"] = 1.e-2*array(nbi["rapedg2"])
        self["ap2_halfheight"]  = 1.e-2*array(nbi["xzpedg2"])
        self["ap2_horiz_offset"]= 1.e-2*array(nbi["xrapoff2"])
        self["ap2_vert_offset"] = 1.e-2*array(nbi["xzapoff2"])

        self["beam_type"] = nbeam*['Standard']
        
        self["power_nbi"] = nbi["pinja"]
        self["kvolt_nbi"] = 1.e-3*array(nbi["einja"])
        self["frac_full"] = nbi["ffulla"]  
        self["frac_half"] = nbi["fhalfa"]  

        self["einj_max"] = 1.25*1.0e-3*array(nbi["einja"])
        self["einj_min"] = nbeam*[20.0]

        self["dn0out"] = 0.5e18  #defalut:0.5e18

    def load_geqdsk (self,fn_geqdsk,bdy_crat=1.e-6):

        kcur_option = 0
        rho_curbrk = 0.9 

        self["EQ_Code_Info"] = 'efit'
        self.updateEquilibrium(fn_geqdsk, bdy_crat, kcur_option, rho_curbrk)
        self.deriveMhdEq('everything')

    def init_from_geqdsk (self,fn_geqdsk,nrho=101,nth=101,shot=0,time=0,bdy_crat=1.e-6):

        self["nrho_eq"] = nrho
        self["nth_eq"] = nth
        self["nrho_eq_geo"] = nrho
        self["nr"] = 0
        self["nz"] = 0

        self.alloc()

        rho =linspace(0.0,1.0,num=nrho)
        self["rho_eq"] = rho
        self["th_eq"] = linspace(0.0,2.0*pi,num=nth)
        self["rho_eq_geo"] = rho

        kcur_option = 0
        rho_curbrk = 0.9 

        self["EQ_Code_Info"] = 'efit'
        self.updateEquilibrium(fn_geqdsk, bdy_crat, kcur_option, rho_curbrk)
        self.deriveMhdEq('everything')

    def load_vol_profile (self,rho,prof,xkey,ykey,k=-1,sum=0):

        zone_spl = zinterp(self["rho_eq"],self["vol"])
        prof_spl = zinterp(rho,prof)

        rho_ps = self[xkey]
     
        dat = zeros(len(rho_ps)-1)
        vol_integral = 0.
        for i in range(len(rho_ps)-1):
            dat[i] = 0.5*(prof_spl(rho_ps[i+1])+prof_spl(rho_ps[i])) \
                        *(zone_spl(rho_ps[i+1])-zone_spl(rho_ps[i]))
            vol_integral += dat[i]
        if sum>0:
            dat *= sum/vol_integral

        if k<0:
            self[ykey] = dat
        else:
            self[ykey][k] = dat

        print 'integrated : ', vol_integral

    def dump_vol_profile (self,rho,xkey,ykey,k=-1):

        zone_spl = zinterp(self["rho_eq"],self["vol"])

        rho_ps = self[xkey]
        if k<0:
           prof_ps = self[ykey]
        else:
           prof_ps = self[ykey][k]

        cell = zeros(len(rho_ps)-1)
        for i in range(len(rho_ps)-1):
            cell[i] = prof_ps[i]/(zone_spl(rho_ps[i+1])-zone_spl(rho_ps[i]))
        node = self.cell2node(cell)

        return zinterp(rho_ps,node)(rho)

    def load_profile (self,rho,prof,xkey,ykey,k=-1):

        prof_spl = zinterp(rho,prof)

        rho_ps = self[xkey]
        prof_ps = prof_spl(rho_ps)

        if k<0:
            self[ykey] = prof_ps
        else:
            self[ykey][k] = prof_ps
       
    def dump_profile (self,rho,xkey,ykey,k=-1):

        rho_ps = self[xkey]

        if k<0:
           prof_ps = self[ykey]
        else:
           prof_ps = self[ykey][k]

        node = self.cell2node(prof_ps)

        return zinterp(rho_ps,node)(rho)

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

    def cell2node_bdry(self,cell):
    
        nrho = len(cell)+1
        node = zeros(nrho)
        node[0] = cell[0]
        for i in range(1,nrho-1):
            node[i] = 0.5*(cell[i-1]+cell[i])
        node[-1] = 2.0*cell[-1]-node[-2] #<=====
        return node

    def load_j_parallel(self,rho,jpar,xkey,ykey,r0,b0,tot=False):

        rho_ps = self[xkey][:]

        ipol = zinterp(self["rho_eq"],self["g_eq"]/(r0*b0))(rho_ps)
        vol  = zinterp(self["rho_eq"],self["vol"])(rho_ps)
        area = zinterp(self["rho_eq"],self["area"])(rho_ps)

        jpar_ps = zinterp(rho,jpar)(rho_ps)

        curt = zeros(len(rho_ps))
        curt[0] = 0.0
        for i in range(len(rho_ps)-1):
            dV = vol[i+1]-vol[i]
            jparm = 0.5*(jpar_ps[i]+jpar_ps[i+1])
            ipolm = 0.5*(ipol[i]+ipol[i+1])
            curt[i+1] = (curt[i]/ipol[i]+jparm*dV/(2.0*pi*r0*ipolm**2))*ipol[i+1]

        print "I%s=%5.3e"%(ykey,1.0e-6*curt[-1])+"MA"

        if tot:
            self[ykey] = curt
        else:
            for i in range(len(rho)-1):
                self[ykey][i] = curt[i+1]-curt[i]

    def dump_j_parallel(self,rho,xkey,ykey,r0,b0,tot=False):
        
        rho_ps = self[xkey]
        jpar_ps = self[ykey]

        ipol = zinterp(self["rho_eq"],self["g_eq"][:]/(r0*b0))(rho_ps)
        vol  = zinterp(self["rho_eq"],self["vol"][:])(rho_ps)
        area = zinterp(self["rho_eq"],self["area"][:])(rho_ps)

        if tot:
           curt = jpar_ps
        else:
           curt = zeros(len(rho_ps))
           curt[0] = 0.0
           for i in range(len(rho)-1): 
               curt[i+1] = curt[i]+jpar_ps[i]

        jpar = zeros(len(rho_ps)-1)
        for i in range(len(rho_ps)-1):
            dV = vol[i+1]-vol[i]
            ipolm = 0.5*(ipol[i]+ipol[i+1])
            jpar[i] = 2.0*pi*r0*ipolm**2/dV*(curt[i+1]/ipol[i+1]-curt[i]/ipol[i])

        return zinterp(rho,self.cell2node(jpar))(rho)

#--test

if __name__ == "__main__":

    ps = zplasmastate("test",1)
    ps.init_from_instate()
    ps.load_innubeam()
    ps.load_geqdsk("g147634.03365")
    print  ps["power_nbi"]
    ps.store("out.nc")
    
    nrho = 20
    rho = arange(nrho)/(nrho-1.0)
    pbe = rho**2
    print pbe
    ps.load_vol_profile(rho,pbe,"rho_nbi","pbe")
    pbe2 = ps.dump_vol_profile(rho,"rho_nbi","pbe")
    print pbe2
    
    tmp = ps.dump_profile(rho,"rho_nbi","pbe")
    print tmp
    
