"""
 -----------------------------------------------------------------------
 model h/cd component
 -----------------------------------------------------------------------
"""

from numpy import *
from scipy import optimize

from component import Component

from plasmastate import plasmastate
from Namelist import Namelist
from efit_eqdsk import readg

class hcd_model(Component):
    def __init__(self, services, config):
        Component.__init__(self, services, config)
        print('Created %s' % (self.__class__))

    def init(self, timeid=0):
        print ('hcd_model.init() called')

    def step(self, timeid=0):
        #--- entry
        print ('hcd_model.init() started')
        services = self.services

        #--- stage plasma state files
        services.stage_state()

        #--- get plasma state file names
        cur_state_file = services.get_config_param('CURRENT_STATE')
        cur_eqdsk_file = services.get_config_param('CURRENT_EQDSK')

        ps_backend = getattr(self, 'PS_BACKEND', 'PS')
        if ps_backend == 'INSTATE':
            cur_instate_file = services.get_config_param('CURRENT_INSTATE')

        #--- stage input files
        services.stage_input_files(self.INPUT_FILES)

        #--- read inhcd
        inhcd = Namelist("inhcd")
        var_list = [ var.upper() for var in inhcd["inhcd"].keys() ]
        for key in var_list:
            try:
                inhcd["inhcd"][key][0] = float(getattr(self, key))
                print(key, 'updated')
            except AttributeError:
                pass
        inhcd.write("inhcd")

        add = int(getattr(self, "ADD", "0"))
        print('add = ',add)

        if ps_backend == 'PS':
            self.update_state(cur_state_file, cur_eqdsk_file, inhcd, add)
        elif ps_backend == 'INSTATE':
            self.update_instate(cur_instate_file, cur_eqdsk_file, inhcd, add)

        #--- update plasma state files
        services.update_state()

        #--- archive output files
        services.stage_output_files(timeid, self.OUTPUT_FILES)

    def finalize(self, timeid=0):
        print ('hcd_model.finalize() called')

    def gauss_profile(self, rho, vol, xmid, xwid):
        nrho = len(rho)
        gauss = lambda p, x: p[0]*exp(-(x-p[1])**2/(2*p[2]**2))
        def func(y0):
            y = gauss([y0,xmid,xwid], rho )
            P = 0.0
            for i in range(nrho-1):
                P +=  0.5*(y[i+1]+y[i])*(vol[i+1]-vol[i])
            return P-1.0
        f0 = optimize.bisect(func, 0.0, 10.0, xtol=1.0e-6)
        profile = gauss( [f0,xmid,xwid], rho )
        return profile

    def update_state(self, cur_state_file, cur_eqdsk_file, inhcd, add):
        nsrc = inhcd["inhcd"]["nsrc"][0]
        Pe = inhcd["inhcd"]["Pe"]
        Pi = inhcd["inhcd"]["Pi"]
        xmid = inhcd["inhcd"]["xmid"]
        xwid = inhcd["inhcd"]["xwid"]

        j0_seed = inhcd["inhcd"]["j0_seed"][0]
        x0_seed = inhcd["inhcd"]["x0_seed"][0]
        drho_seed = inhcd["inhcd"]["drho_seed"][0]

        print("Pe =", Pe)
        print("Pi =", Pi)
        print("xmid =", xmid)
        print("xwid =", xwid)

        #--- read ps
        ps = plasmastate('ips', 1)
        ps.read(cur_state_file)

        geq = readg(cur_eqdsk_file)
        r0  = geq["rzero" ]
        b0  = abs(geq["bcentr"])
        ip  = geq['cpasma']

        rho = ps["rho"][:]
        vol = ps["vol"][:]
        nrho = len(rho)

        #--- heating
        def vol_integral(f):
            val = 0.0
            for i in range(nrho-1):
                f_m = 0.5*(f[i+1]+f[i])
                dvol = vol[i+1]-vol[i]
                val += f_m*dvol
            return val

        pe_sum = zeros(nrho)
        pi_sum = zeros(nrho)
        for k in range(nsrc):
            profile = self.gauss_profile(rho, vol, xmid[k], xwid[k])
            pe_sum = pe_sum + profile*Pe[k]*1.0e6
            pi_sum = pi_sum + profile*Pi[k]*1.0e6

        ps.load_vol_profile(rho,pe_sum, "rho_icrf", "picrf_totals", k=0, add=add)
        ps.load_vol_profile(rho,pi_sum, "rho_icrf", "picrf_totals", k=1, add=add)

        #--- current
        if j0_seed > 0.0:
            jec_seed = j0_seed*exp(-(rho-x0_seed)**2/(2*drho_seed**2))
        else:
            jec_seed = zeros(nrho)
        ps.load_j_parallel(rho,jec_seed*1.0e6,"rho_icrf","curich",r0,b0,add=add)

        #--- write to ps
        ps.store(cur_state_file)

    def update_instate(self, cur_instate_file, cur_eqdsk_file, inhcd, add):
        nsrc = inhcd["inhcd"]["nsrc"][0]
        Pe = inhcd["inhcd"]["Pe"]
        Pi = inhcd["inhcd"]["Pi"]
        xmid = inhcd["inhcd"]["xmid"]
        xwid = inhcd["inhcd"]["xwid"]

        j0_seed = inhcd["inhcd"]["j0_seed"][0]
        x0_seed = inhcd["inhcd"]["x0_seed"][0]
        drho_seed = inhcd["inhcd"]["drho_seed"][0]

        print("Pe =", Pe)
        print("Pi =", Pi)
        print("xmid =", xmid)
        print("xwid =", xwid)

        #--- read instate
        instate = Namelist(cur_instate_file)

        rho = array(instate["inmetric"]["rho"])
        rhob = instate["inmetric"]["rhob"][0]
        volp = instate["inmetric"]["volp"]
        nrho = len(rho)

        vol = zeros(nrho)
        for i in range(nrho-1):
            vol[i+1] = vol[i] + 0.5*(volp[i]+volp[i+1])*(rho[i+1]-rho[i])*rhob

        def vol_integral(f):
            val = 0.0
            for i in range(nrho-1):
                f_m = 0.5*(f[i+1]+f[i])
                dvol = vol[i+1]-vol[i]
                val += f_m*dvol
            return val

        #--- heating
        pe_sum = zeros(nrho)
        pi_sum = zeros(nrho)
        for k in range(nsrc):
            profile = self.gauss_profile(rho,vol,xmid[k],xwid[k])
            pe_sum = pe_sum + profile*Pe[k]
            pi_sum = pi_sum + profile*Pi[k]

        if add:
            instate["instate"]["pe_ic"] = array(instate["instate"]["pe_ic"]) + pe_sum
            instate["instate"]["pi_ic"] = array(instate["instate"]["pi_ic"]) + pi_sum

        #--- current
        if j0_seed > 0.0:
            jec_seed = j0_seed*exp(-(rho-x0_seed)**2/(2*drho_seed**2))
        else:
            jec_seed = zeros(nrho)
        instate["instate"]["j_ic"] = jec_seed

        #--- write to instate
        instate.write(cur_instate_file)
