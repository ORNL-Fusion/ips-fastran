import numpy as np
from Namelist import Namelist
from fastran.util.modelprofile import profile_pedestal, profile_parab
from fastran.util import shape_io 
from fastran.util.formula import mu0

instate_variables = [
    "ne",
    "te",
    "ti",
    "omega",
    "zeff",
    "j_oh", 
    "j_bs", 
    "j_nb", 
    "j_ec", 
    "j_ic", 
    "pe_nb", 
    "pe_ec", 
    "pe_ic", 
    "pe_fus", 
    "pe_ionization", 
    "p_rad", 
    "pi_nb", 
    "pi_ec", 
    "pi_ic", 
    "pi_fus",
    "pi_ionization", 
    "pi_cx", 
    "p_ohm", 
    "p_ei", 
    "torque_nb", 
    "torque_in", 
    "se_nb", 
    "se_ionization", 
    "si_nb", 
    "si_ionization", 
    "q", 
    "psipol", 
    "density_beam", 
    "wbeam", 
    "tbeam",
    "density_alpha", 
    "walpha", 
    "talpha",
    "chie", 
    "chii", 
    "p_eq",
    "j_tot", 
]

class Instate():
    def __init__(self, fname):
        self.data = Namelist(fname)

        #-- check error
        self.check()

        #-- convert to numpy array
        for key in self.keys():
            if key.upper() not in ['TOKAMAK_ID', 'PRESSURE_MODEL', 'CURRENT_MODEL']:
                self[key] = np.array(self[key])

    def check(self):
        self.density_model = self["density_model"][0]
        if self.density_model not in [0, 1]:
            raise Exception("density_model error")
        self.model_shape = self["model_shape"][0]
        if self.model_shape not in [0, 1, 2]:
            raise Exception("model_shape error")
        self.model_beam = 0

    def __getitem__(self, key):
        return self.data["instate"][key]

    def __setitem__(self, key, value):
        self.data["instate"][key] = value

    def __contains__(self, key):
        return True if key.upper() in self.keys() else False 

    def keys(self):
        return self.data["instate"].keys()

    def default(self, key, value):
        if key in self: return self[key]
        else: return value

    def write(self, fname):
        self.data.write(fname)

    def write_netcdf(self, fname):
        self.data.write(fname)
            
    def set_shape(self):
        if self.model_shape  == 0:
            print("Shape from instate data array")
            rb, zb, rlim, zlim = \
                self["rbdry"], \
                self["zbdry"], \
                self["rlim"], \
                self["zlim"]
        elif self.model_shape == 1:
            print("Simple Shape")
            rb, zb, rlim, zlim = \
                shape_io.set_shape(
                    R0 = self["r0"][0],
                    a0 = self["a0"][0],
                    kappa = self["kappa"][0],
                    delta = self["delta"][0],
                    nt = self["nbdry"][0])
        elif self.model_shape == 2:
            print("Luce Shape")
            R0 = self["r0"][0]
            a0 = self["a0"][0]
            eps = a0/R0
            kapu = self["kappa"][0]
            kapl = self["kappa"][0]
            delu = self["delta"][0]
            dell = self["delta"][0]
            z0 = self.default("z0", [0.])[0]
            nbdry = self["nbdry"][0]

            zetaou = 0.
            zetaiu = 0.
            zetail = 0.
            zetaol = 0.
            zoffset = 0.

            rb, zb, zref = shape_io.boundaryShape(a0, eps, kapu, kapl, delu, dell, zetaou, zetaiu, zetail, zetaol, zoffset,
                              upnull=True, lonull=True, npts=nbdry, doPlot=False)

            print(rb)
            rb = np.append(rb, rb[0])
            zb = np.append(zb, zb[0])

            zb += z0

            dlim = 0.05
            rmax = max(rb) + dlim
            rmin = min(rb) - dlim
            zmax = max(zb) + dlim
            zmin = min(zb) - dlim
            rlim = [ rmax, rmin, rmin, rmax, rmax ]
            zlim = [ zmax, zmax, zmin, zmin, zmax ]

        self["nbdry"] = [len(rb)]
        self["rbdry"] = rb
        self["zbdry"] = zb
        self["nlim"] = [len(rlim)]
        self["rlim"] = rlim
        self["zlim"] = zlim

    def from_timetrace(self, itime, extract=["ip", "bt", "rbdry", "zbdry", "ne", "te", "ti", "omega", "zeff"]):
        #-- extract profile timetrace data 
        pass

    def zeros(self):
        for key in instate_variables:
            if key.upper() not in self.data["instate"].keys():
                self[key] = np.zeros(self["nrho"][0])

    def particle_balance(self):
        nrho = self["nrho"][0]
        rho = self["rho"]

        n_ion = self["n_ion"][0]
        z_ion = self["z_ion"]
        a_ion = self["a_ion"]
        f_ion = self["f_ion"]

        n_imp = self["n_imp"][0]
        z_imp = self["z_imp"]
        a_imp = self["a_imp"]
        f_imp = self["f_imp"]

        z_beam = self["z_beam"][0]

        density_ion = np.zeros((n_ion, nrho))
        density_imp = np.zeros((n_imp, nrho))
        density_th = np.zeros(nrho)

        if self.density_model == 0:
            print('density_model = 0')

            a=0; b=0; c=0; d=0
            for k in range(n_imp):
                b = b + f_imp[k]*z_imp[k]
                d = d + f_imp[k]*z_imp[k]*z_imp[k]
            for k in range(n_ion):
                a = a + f_ion[k]*z_ion[k]
                c = c + f_ion[k]*z_ion[k]*z_ion[k]

            for i in range(nrho):
                zne_adj = self["ne"][i]
                zzne_adj = self["ne"][i]*self["zeff"][i]

                # depletion due to beam ions
                zne_adj = zne_adj - z_beam*self["density_beam"][i] - 2.0*self["density_alpha"][i]
                zzne_adj = zzne_adj - z_beam**2*self["density_beam"][i] - 2.0**2*self["density_alpha"][i]

                # effective main ion and impurity densities
                nion = (zne_adj *d-zzne_adj*b)/(a*d-b*c)
                nimp = (zzne_adj*a-zne_adj *c)/(a*d-b*c)

                for k in range(n_ion):
                    density_ion[k][i] = f_ion[k]*nion
                for k in range(n_imp):
                    density_imp[k][i] = f_imp[k]*nimp

        elif self.density_model == 1:
            print('density_model = 1: impurity = f*ne')

            for k in range(n_imp):
                density_imp[k] = self["ne"]*f_imp[k]
            nith = self["ne"] - self["density_beam"] - 2.*self["density_alpha"]
            for k in range(n_imp):
                nith = nith - z_imp[k]*density_imp[k]
            for k in range(n_ion):
                density_ion[k] = nith*f_ion[k]

            zeff = nith + z_beam**2*self["density_beam"] + 2.**2*self["density_alpha"] 
            for k in range(n_imp):
                zeff = zeff + z_imp[k]**2*density_imp[k]
            self["zeff"] = zeff/self["ne"]

        for k in range(n_ion):
            self["density_ion_{}".format(k)] = density_ion[k]
            density_th = density_th + density_ion[k]
        for k in range(n_imp):
            self["density_imp_{}".format(k)] = density_imp[k]
            density_th = density_th + density_imp[k]
        self["density_th"] = density_th

        self["ni"] = np.array([sum(tmp) for tmp in density_ion.transpose()])
        self["nz"] = np.array([sum(tmp) for tmp in density_imp.transpose()])

        """
        zne_adj = ne
        zzne_adj = ne*zeff

        zne_adj = zne_adj - 1.0*density_beam
        zzne_adj = zzne_adj - 1.0**2*density_beam

        nion = (zne_adj *d-zzne_adj*b)/(a*d-b*c)
        nimp = (zzne_adj*a-zne_adj *c)/(a*d-b*c)

        density_ion = array([f_ion[k]*nion for k in range(n_ion)])
        density_imp = array([f_imp[k]*nimp for k in range(n_imp)])

        density_th = array([sum(tmp) for tmp in density_ion.transpose()])
        density_th+= array([sum(tmp) for tmp in density_imp.transpose()])

        ni = array([sum(tmp) for tmp in density_ion.transpose()])
        nz = array([sum(tmp) for tmp in density_imp.transpose()])
        """

    def model_profile_(self):
        nrho = self["nrho"][0]
        rho = np.arange(nrho)/(nrho-1.0)

        self["rho"] = rho
        for key in ["ne", "te", "ti"]:
            self[key] = profile_pedestal(
                nrho, 
                self["{}_xmid".format(key)][0], 
                self["{}_xwid".format(key)][0],
                self["{}_ped".format(key)][0], 
                self["{}_sep".format(key)][0],
                self["{}_axis".format(key)][0], 
                self["{}_alpha".format(key)][0], 
                self["{}_beta".format(key)][0], 
                ifit=self["{}_fit".format(key)][0])(rho)

        for key in ["omega", "zeff", "density_beam", "tbeam", "density_alpha", "jpar"]:
            self[key] = profile_parab(
                self["{}_axis".format(key)][0], 
                self["{}_sep".format(key)][0],
                self["{}_alpha".format(key)][0], 
                self["{}_beta".format(key)][0])(rho)

        self["j_tot"] = self["jpar"][:]

    def model_profile(self):
        #-- generate profies from the model profile parameters
        nrho = self["nrho"][0]
        rho = np.arange(nrho)/(nrho-1.0)

        xmid = self.default("xmid", [-1.0])[0]
        xwid = self.default("xwid", [-1.0])[0]

        if xmid > 0:
            self["ne_xmid"] = [xmid]
            self["te_xmid"] = [xmid]
            self["ti_xmid"] = [xmid]
        if xwid > 0:
            self["ne_xwid"] = [xwid]
            self["te_xwid"] = [xwid]
            self["ti_xwid"] = [xwid]

        betan_ped = self.default("betan_ped", [-1.0])[0]
        if betan_ped > 0.0:
            rb = np.array(self["rbdry"])
            zb = np.array(self["zbdry"])
            a0 = 0.5*( np.max(rb) - np.min(rb) )

            b0 = abs(self["b0"][0])
            c_betan = 4.0*1.602e5*self["ne_ped"][0]*mu0/b0**2*(a0*b0/self["ip"][0])
            te_ped = betan_ped/c_betan
            ti_ped = te_ped
            print ("PEDEDSTAL BETAN = ", betan_ped, te_ped)

            self["te_ped" ] = [te_ped ]
            self["te_xmid"] = [xmid]
            self["te_xwid"] = [xwid]

            self["ti_ped" ] = [ti_ped ]
            self["ti_xmid"] = [xmid]
            self["ti_xwid"] = [xwid]

        self["rho"] = rho
        self["ne"] = profile_pedestal(nrho, self["ne_xmid"][0], self["ne_xwid"][0], self["ne_ped"][0], self["ne_sep"][0], self["ne_axis"][0], self["ne_alpha"][0], self["ne_beta"][0], ifit=self["ne_fit"][0])(rho)
        self["te"] = profile_pedestal(nrho, self["te_xmid"][0], self["te_xwid"][0], self["te_ped"][0], self["te_sep"][0], self["te_axis"][0], self["te_alpha"][0], self["te_beta"][0], ifit=self["te_fit"][0])(rho)
        self["ti"] = profile_pedestal(nrho, self["ti_xmid"][0], self["ti_xwid"][0], self["ti_ped"][0], self["ti_sep"][0], self["ti_axis"][0], self["ti_alpha"][0], self["ti_beta"][0], ifit=self["ti_fit"][0])(rho)

        self["omega"] = (self["omega_axis"][0] - self["omega_sep"][0])*(1. - rho**self["omega_alpha"][0])**self["omega_beta"][0] + self["omega_sep"][0]
        self["zeff"] = self["zeff_axis"][0]*np.ones(nrho)
        self["density_beam"] = (self["nbeam_axis"][0] - self["nbeam_sep"][0])*(1. - rho**self["nbeam_alpha"])**self["nbeam_beta"][0] + self["nbeam_sep"][0]
        self["wbeam"] = 1.5*1.602e3*self["density_beam"]*self["tbeami"][0]*1.e-6
        self["density_alpha"] = np.zeros(nrho)
        self["walpha"] = np.zeros(nrho)

        self["j_tot"] = (self["jpar_axis"][0] - self["jpar_sep"][0])*(1.0 - rho**self["jpar_alpha"][0])**self["jpar_beta"][0] + self["jpar_sep"][0]

    def to_ps(self, ps, **keyarg):
        """
        required input profiles in instate["instate"] 
            ne, te, ti, omega, zeff 
            tbeam, density_beam (model_beam = 0) or wbeam (model_beam = 1)
      
        """
        nrho = self["nrho"][0]
        r0 = self["r0"][0]
        b0 = self["b0"][0]

        #ne[0] = ne[1]
        #te[0] = te[1]
        #ti[0] = ti[1]
        #omega[0] = omega[1]
        #zeff[0] = zeff[1]

        #-- density
        ps["ns"][0,:] = 1.e19*ps.node2cell(self["ne"])
        for k in range(self["n_ion"][0]):
            ps["ns"][k + 1, :] = 1.e19*ps.node2cell(self["density_ion_{}".format(k)])
        for k in range(self["n_imp"][0]):
            ps["ns"][k + self["n_ion"][0] + 1,:] = 1.e19*ps.node2cell(self["density_imp_{}".format(k)])
        ps["ni"][:] = 1.e19*ps.node2cell(self["density_th"])
    
        #-- beam
        if self.model_beam == 0:
            ps["nbeami"][0][:] = 1.e19*ps.node2cell(self["density_beam"])
            ps["eperp_beami"][0][:] = 2.*self["tbeam"][0]*np.ones(nrho - 1)
            ps["epll_beami"][0][:] = self["tbeam"][0]*np.ones(nrho - 1)
        elif self.model_beam == 1:
            tbeam = self["wbeam"]/1.602e-3/(delf["density_beam"]+1.e-6)
            ps["eperp_beami"][0][:] = ps.node2cell(self["zbeam"][0]*self["tbeam"]/3.0)
            ps["epll_beami"][0][:] = ps.node2cell(self["tbeam"]/3.0)
       
        #-- temperature
        ps["Ts"][0,:] = ps.node2cell(self["te"])
    
        for k in range(self["n_ion"][0]):
            ps["Ts"][k + 1, :] = ps.node2cell(self["ti"])
        for k in range(self["n_imp"][0]):
            ps["Ts"][k + self["n_ion"][0] + 1, :] = ps.node2cell(self["ti"])
    
        ps["Ti"][:] = ps.node2cell(self["ti"])
    
        #-- zeff
        ps["Zeff"][:] = ps.node2cell(self["zeff"])
        ps["Zeff_th"][:] = ps.node2cell(self["zeff"])
    
        #-- rotation
        ps["omegat"][:] = ps.node2cell(self["omega"])
    
        #-- current
        self["j_tot"][0] = self["j_tot"][1]
        ps.load_j_parallel(self["rho"], 1.e6*self["j_tot"], "rho_eq", "curt", r0, b0, tot=True)

        ps.load_j_parallel(self["rho"], 1.e6*self["j_nb"], "rho_nbi", "curbeam", r0, b0)
        ps.load_j_parallel(self["rho"], 1.e6*self["j_ec"], "rho_ecrf", "curech", r0, b0)
        ps.load_j_parallel(self["rho"], 1.e6*self["j_ic"], "rho_icrf", "curich", r0, b0)
        ps.load_j_parallel(self["rho"], 1.e6*self["j_bs"], "rho", "curr_bootstrap", r0, b0)
        ps.load_j_parallel(self["rho"], 1.e6*self["j_oh"], "rho", "curr_ohmic", r0 ,b0)
    
        #-- MHD pressure
        ps["P_eq"][:] = self["p_eq"]
    
        #-- heating
        pe_fus = 1.e6*self["pe_fus"]
        pi_fus = 1.e6*self["pi_fus"]
    
        ps.load_vol_profile(self["rho"], 1.e6*self["pe_nb"], "rho_nbi", "pbe")
        ps.load_vol_profile(self["rho"], 1.e6*self["pi_nb"], "rho_nbi", "pbi")
        ps.load_vol_profile(self["rho"], 1.e6*self["pe_ec"], "rho_ecrf", "peech")
        ps.load_vol_profile(self["rho"], 1.e6*self["pe_ic"], "rho_icrf", "picrf_totals", k=0)
        ps.load_vol_profile(self["rho"], 1.e6*self["pi_ic"], "rho_icrf", "picth")
    
        #-- particle source
        ps.load_vol_profile(self["rho"], 1.e19*self["se_nb"], "rho_nbi", "sbedep")
    
        #-- torque
        ps.load_vol_profile(self["rho"], self["torque_nb"], "rho_nbi", "tqbe")
    
    def from_ps(self):
        pass
    
if __name__=="__main__":
    instate = Instate("instate0")

    instate.model_profile()
    instate.zeros()
    instate.particle_balance()
    print(instate.data)

    from fastran.plasmastate.plasmastate import plasmastate

    ps = plasmastate('ips', 1)

    ps.init(
        ps,
        global_label = ['ips'],
        runid = ['ips'],
        shot_number = [123456],
        nspec_th = [instate["n_ion"][0] + instate["n_imp"][0]] ,
        nspec_beam = [instate["n_beam"][0]],
        nspec_fusion = [instate["n_fusion"][0]],
        nspec_rfmin = [instate["n_min"][0]],
        nspec_gas = [instate["n_ion"][0]],
        z_ion = instate["z_ion"] + instate["z_imp"],
        a_ion = instate["a_ion"] + instate["a_imp"],
        z_beam = instate["z_beam"],
        a_beam = instate["a_beam"],
        z_fusion = instate["z_fusion"],
        a_fusion = instate["a_fusion"],
        z_rfmin = instate["z_min"],
        a_rfmin = instate["a_min"],
        z_gas = instate["z_ion"],
        a_gas = instate["a_ion"],
        nicrf_src = [1],
        nlhrf_src = [1],
        nrho = [101],
        time = [0.0],
        nlim = instate["nlim"]
    )

    instate.to_ps(ps)

    ps.store("ps.nc")
