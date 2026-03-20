"""
calculate critical_deltap with equilibrium, profile by Shaing model
"""

from ipsframework import Component
from Namelist import Namelist
from fastran.util.zinterp import zinterp
from numpy import *
from scipy import constants, interpolate, optimize, stats
from fastran.equilibrium import efit_eqdsk
from fastran.plasmastate.plasmastate import plasmastate


class critical_deltap(Component):
    def __init__(self, services, config):
        Component.__init__(self, services, config)
        print('Created %s' % (self.__class__))

    def init(self, timeid=0):
        print('critical_deltap.init() called')

    def step(self, timeid=0):
        print('critical_deltap.step() started')

        self.services.stage_state()
        cur_instate_file = self.services.get_config_param('CURRENT_INSTATE')
        cur_eqdsk_file = self.services.get_config_param('CURRENT_EQDSK')
        cur_aeqdsk_file = self.services.get_config_param('CURRENT_AEQDSK')

        instate = Namelist(cur_instate_file)
        #gfile = efit_eqdsk.readg(cur_eqdsk_file)
        ps = plasmastate.plasmastate('ips', 1)
        ps.init_from_geqdsk(cur_eqdsk_file)
        #ps.load_geqdsk(cur_eqdsk_file)

        afile = efit_eqdsk.reada(cur_aeqdsk_file)

        R0 = float(afile['r0'])
        B0 = float(abs(afile['b0']))
        betap = float(afile['betap'])

        rho = instate["instate"]["rho"]
        qr_spl=interpolate.UnivariateSpline(rho, ps["q_eq"][:])
        q2_spl=interpolate.UnivariateSpline(rho, ps["q_eq"][:] - 2.)
        pr_spl=interpolate.UnivariateSpline(rho, ps["P_eq"][:])
        ar_spl=interpolate.UnivariateSpline(rho, ps["rMinor_mean"][:])
        q2_rho = q2_spl.roots()[-1]
        rs = ar_spl(q2_rho)

        Lq = qr_spl(q2_rho)/(qr_spl.derivatives(q2_rho)[1]/ar_spl.derivatives(q2_rho)[1])
        Lp = -1.*pr_spl(q2_rho)/(pr_spl.derivatives(q2_rho)[1]/ar_spl.derivatives(q2_rho)[1])

        ne = instate["instate"]["ne"]
        te = instate["instate"]["te"]
        pe = map(lambda x, y:1.e19*x*constants.value('atomic unit of charge')*1.e3*y, ne, te)
        zeff = instate["instate"]["zeff"]

        ne_spl = interpolate.UnivariateSpline(rho,ne)
        te_spl = interpolate.UnivariateSpline(rho,te)
        pe_spl = interpolate.UnivariateSpline(rho,pe)
        zeff_spl = interpolate.UnivariateSpline(rho,zeff)

        q2 = 2.
        #zeff = 30.*nz_spl(q2_rho)/ne_spl(q2_rho) + 1.

        nustare = (6.921*1.e-18)*R0*q2*1.e19*ne_spl(q2_rho)*zeff_spl(q2_rho)*(24.-log((1.e19*ne_spl(q2_rho))**0.5/(1.e3*te_spl(q2_rho))))/((rs/R0)**1.5*(1.e3*te_spl(q2_rho))**2)

        cridelp = 1.2/nustare**0.5*(1/1.7)**0.5*(Lq/rs)*betap*(rs/R0)**0.5/Lp

        print(Lq, Lp, rs, cridelp)

        instate["stab"]["Lq"] = [Lq]
        instate["stab"]["Lp"] = [Lp]
        instate["stab"]["cridelp"] = [rs*cridelp]
        instate.write(cur_instate_file)

        self.services.update_state()
        self.services.stage_output_files(timeid, self.OUTPUT_FILES, save_plasma_state=False)

        return

    def finalize(self, timeid=0):
        return

"""
    def calculate_delpcri(self, f_eqdsk,f_aeqdsk,f_instate):

        instate = Namelist(f_instate)
        gfile = efit_eqdsk.readg(f_eqdsk)
        ps = plasmastate.plasmastate('ips',1)
        ps.init_from_geqdsk(gfile)

        afile=efit_eqdsk.reada(f_aeqdsk)

        R0=float(aeqdsk['r0'])
        B0=float(abs(aeqdsk['b0']))
        betap=float(aeqdsk['betap'])

        qr_spl=interpolate.UnivariateSpline(rho,ps["q_eq"][:])
        q2_spl=interpolate.UnivariateSpline(rho,ps["q_eq"][:]-2.)
        pr_spl=interpolate.UnivariateSpline(rho,ps["P_eq"][:])
        ar_spl=interpolate.UnivariateSpline(rho,ps["rMinor_mean"][:])
        q2_rho=q2_spl.roots()[-1]
        rs=ar_spl(q2_rho)

        Lq=qr_spl(q2_rho)/(qr_spl.derivatives(q2_rho)[1]/ar_spl.derivatives(q2_rho)[1])
        Lp=-1.*pr_spl(q2_rho)/(pr_spl.derivatives(q2_rho)[1]/ar_spl.derivatives(q2_rho)[1])

        rho = instate["instate"]["rho"]
        ne = instate["instate"]["ne"]
        te = instate["instate"]["te"]
        pe=1E19*ne*constants.value('atomic unit of charge')*1E3*te
        zeff = instate["instate"]["zeff"]

        ne_spl=interpolate.UnivariateSpline(rho,ne)
        te_spl=interpolate.UnivariateSpline(rho,te)
        pe_spl=interpolate.UnivariateSpline(rho,pe)
        zeff_spl=interpolate.UnivariateSpline(rho,zeff)

        q2=2.
        #zeff=30.0*nz_spl(q2_rho)/ne_spl(q2_rho)+1.0

        nustare=(6.921*1E-18)*R0*q2*1E19*ne_spl(q2_rho)*zeff_spl(q2_rho)*(24-numpy.log((1E19*ne_spl(q2_rho))**0.5/(1E3*te_spl(q2_rho))))/((rs/R0)**1.5*(1E3*te_spl(q2_rho))**2)

        cridelp=1.2/nustare**0.5*(1/1.7)**0.5*(Lq/rs)*betap*(rs/R0)**0.5/Lp

        print Lq, Lp, rs, cridelp

        return Lq, Lp, rs, cridelp
"""
