"""
 -----------------------------------------------------------------------
plasma state, pyps backed
 -----------------------------------------------------------------------
"""

from numpy import *
from Namelist import Namelist
from zinterp import zinterp
from efit_eqdsk import readg

from PlasmaState import PlasmaState
from plasmastate_base import plasmastate_base

class plasmastate(PlasmaState, plasmastate_base):

    def init(self, f_ps, **keyargs):
        nrho = keyargs['nrho'][0]
        nth  = 101
        nlim = keyargs['nlim'][0]

        nrho_nbi = nrho

        nspec_th = keyargs['nspec_th'][0]
        z_spec = keyargs["z_ion"]
        a_spec = keyargs["a_ion"]

        nspec_beam = keyargs['nspec_beam'][0]
        if nspec_beam:
            z_beam = keyargs['z_beam']
            a_beam = keyargs['a_beam']

        nspec_fusion = keyargs['nspec_fusion'][0]
        if nspec_fusion:
            z_fusion = keyargs['z_fusion']
            a_fusion = keyargs['a_fusion']

        nspec_rfmin = keyargs['nspec_rfmin'][0]
        if nspec_rfmin:
            z_rfmin = keyargs['z_rfmin']
            a_rfmin = keyargs['a_rfmin']

        nicrf_src = keyargs['nicrf_src'][0]
        nlhrf_src = keyargs['nlhrf_src'][0]

        self["nspec_th"] = nspec_th
        self["nspec_beam"] = nspec_beam
        self["nspec_fusion"] = nspec_fusion
        self["nspec_rfmin"] = nspec_rfmin

        self["nrho"] =  nrho
        self["nrho_eq"] = nrho
        self["nth_eq"] = nth
        self["nrho_eq_geo"] = nrho
        self["nr"] = 0
        self["nz"] = 0

        self["nrho_gas"] = nrho
        self["nrho_nbi"] = nrho_nbi #nrho
        self["nrho_ecrf"] = nrho
        self["nrho_icrf"] = nrho
        self["nrho_lhrf"] = nrho
        self["nrho_fus"] = nrho
        self["nrho_anom"] = nrho

        self["nicrf_src"] = nicrf_src
        self["nlhrf_src"] = nlhrf_src

        #-------
        self["global_label"] = keyargs['global_label'][0]
        self["tokamak_id"] = ''
        self["runid"] = keyargs['runid'][0]
        self["shot_number"] = keyargs['shot_number'][0]

        self["t0"] = keyargs['time'][0]
        self["t1"] = keyargs['time'][0]

        self["num_rzlim"] = nlim

        self.alloc()

        self["icrf_src_name"] = ['IC%02d'%k for k in range(nicrf_src)]
        self["lhrf_src_name"] = ['LH%02d'%k for k in range(nlhrf_src)]

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

        rho =linspace(0.0, 1.0, num=nrho)
        rho_nbi =linspace(0.0, 1.0, num=nrho_nbi)

        self["rho"] = rho
        self["rho_eq"] = rho
        self["th_eq"] = linspace(0.0, 2.0*pi, num=nth)
        self["rho_eq_geo"] = rho

        self["rho_gas"] = rho
        self["rho_nbi"] = rho_nbi #rho
        self["rho_ecrf"] = rho
        self["rho_icrf"] = rho
        self["rho_lhrf"] = rho
        self["rho_fus"] = rho
        self["rho_anom"] = rho

    def load_innubeam(self, fn="innubeam"):
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
            if not nbi["nlco"][k]: srtcen[k] *= -1.
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

    def load_geqdsk(self, fn_geqdsk, keep_cur_profile=True, bdy_crat=1.e-6, kcur_option=1, rho_curbrk=0.9):
        if keep_cur_profile:
            geq = readg(fn_geqdsk)
            r0  = geq["rzero" ]
            b0  = abs(geq["bcentr"])
            ip  = geq['cpasma']

            j_tot = self.dump_j_parallel(self["rho"],"rho_eq","curt",r0,b0,tot=True)

        self["EQ_Code_Info"] = 'efit'
        self.updateEquilibrium(fn_geqdsk, bdy_crat, kcur_option, rho_curbrk)
        self.deriveMhdEq('everything')

        if keep_cur_profile:
            ps_geqdsk = plasmastate('ips',1)
            ps_geqdsk.init_from_geqdsk (fn_geqdsk, nrho=self["nrho"], nth=101)

            self["vol"][:] = ps_geqdsk["vol"][:]
            self["g_eq"][:] = ps_geqdsk["g_eq"][:]
            self.load_j_parallel(self["rho"], j_tot, "rho_eq", "curt", r0, b0, tot=True)
            j_tot = self.dump_j_parallel(self["rho"],"rho_eq","curt",r0,b0,tot=True)

    def init_from_geqdsk(self, fn_geqdsk, nrho=101, nth=101, shot=0, time=0, bdy_crat=1.e-6):
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

        kcur_option = 1
        rho_curbrk = 0.9

        self["EQ_Code_Info"] = 'efit'
        self.updateEquilibrium(fn_geqdsk, bdy_crat, kcur_option, rho_curbrk)
        self.deriveMhdEq('everything')

#--test
if __name__ == "__main__":
    ps = plasmastate("test",1)
    ps.init_from_instate()
    ps.load_innubeam()
    ps.load_geqdsk("g147634.03365")
    print  (ps["power_nbi"])
    ps.store("out.nc")

    nrho = 20
    rho = arange(nrho)/(nrho-1.0)
    pbe = rho**2
    print (pbe)
    ps.load_vol_profile(rho,pbe,"rho_nbi","pbe")
    pbe2 = ps.dump_vol_profile(rho,"rho_nbi","pbe")
    print (pbe2)

    tmp = ps.dump_profile(rho,"rho_nbi","pbe")
    print (tmp)
