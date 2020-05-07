"""
 -----------------------------------------------------------------------
 constrain component
 -----------------------------------------------------------------------
"""

from component import Component
from plasmastate import plasmastate

from numpy import *
from Namelist import Namelist
from zmodelprof import profile_pedestal
from zinterp import zinterp
from formula import get_ni

class model_pedestal(Component):

    def __init__(self, services, config):
        Component.__init__(self, services, config)
        print('Created %s' % (self.__class__))

    def init(self, timeStamp=0):
        print('pedestal.init() called')

    def step(self, timeStamp=0):
        print('pedestal.step() started')

        #--- entry
        services = self.services

        #--- stage plasma state files
        services.stage_plasma_state()

        #--- get plasma state file name
        cur_state_file = services.get_config_param('CURRENT_STATE')
        cur_bc_file = services.get_config_param('CURRENT_BC')
        cur_instate_file = services.get_config_param('CURRENT_INSTATE')
        cur_eqdsk_file = services.get_config_param('CURRENT_EQDSK')

        #--- stage input files
        services.stage_input_files(self.INPUT_FILES)

        #--- apply constraint to local plasma state
        inbc = Namelist(cur_bc_file)
        r0  = inbc["inbc"]["r0"][0]
        b0  = inbc["inbc"]["b0"][0]
        ip  = inbc["inbc"]['ip'][0]

        instate = Namelist(cur_instate_file)
        neped = instate["instate"]["ne_ped"][0]
        nesep = instate["instate"]["ne_sep"][0]/neped

        inped = self.INPED
        wped_fit = loglinear(inped,'wped')
        teped_fit = loglinear(inped,'teped')

        wped = wped_fit.get(ip=ip, neped=neped, nesep=nesep)
        pped = pped_fit.get(ip=ip, neped=neped, nesep=nesep)

        zeffped = instate["instate"]["zeff_axis"][0]

        niped, nzped = get_ni(neped, zeff=zeffped)
        teped = 0.5*pped/(1.602*neped)
        tiped = 0.5*pped/(1.602*(niped+nzped))

        print('neped:', neped)
        print('niped:', niped)
        print('nzped:', nzped)
        print('teped:', teped)
        print('tiped:', tiped)

        ps = plasmastate('ips', 1)
        ps.read(cur_state_file)

        rho   = ps["rho"][:]
        nrho  = len(rho)
        ne = ps["ns"][0,:]*1.0e-19
        te = ps["Ts"][0,:]
        ti = ps["Ti"][:]
        ne = ps.cell2node_bdry(ne)
        te = ps.cell2node_bdry(te)
        ti = ps.cell2node_bdry(ti)

        alpha = 1.5
        beta = 1.5
        xtop = 1.0-1.5*wped

        ne_model = profile_pedestal(nrho, 1.0-0.5*wped, wped, neped, nesep, ne[0], alpha, beta)
        te_model = profile_pedestal(nrho, 1.0-0.5*wped, wped, teped, tesep, te[0], alpha, beta)
        ti_model = profile_pedestal(nrho, 1.0-0.5*wped, wped, tiped, tisep, ti[0], alpha, beta)

        ne_model_top = ne_model(xtop)
        te_model_top = te_model(xtop)
        ti_model_top = ti_model(xtop)

        ne_top = zinterp(rho,ne)(xtop)
        te_top = zinterp(rho,te)(xtop)
        ti_top = zinterp(rho,ti)(xtop)

        ne_update = zeros(nrho)
        te_update = zeros(nrho)
        ti_update = zeros(nrho)

        print("ne top", ne_model_top , ne_top)
        print("te top", te_model_top , te_top)
        print("ti top", ti_model_top , ti_top)

        for i in range(nrho):
            if rho[i] < xtop:
                ne_update[i] = ne[i] + ne_model_top - ne_top
                te_update[i] = te[i] + te_model_top - te_top
                ti_update[i] = ti[i] + ti_model_top - ti_top
            else:
                ne_update[i] = ne_model(rho[i])
                te_update[i] = te_model(rho[i])
                ti_update[i] = ti_model(rho[i])

        ps["ns"][0] = 1.0e19*ps.node2cell(ne_update)
        ps["Ts"][0] = ps.node2cell(te_update)
        nspec_th = len(ps["Ts"])-1
        print('nspec_th =',nspec_th)
        for k in range(nspec_th):
            ps["Ts"][k+1] = ps.node2cell(ti_update)
        ps["Ti"][:] = ps.node2cell(ti_update)

        ps.store(cur_state_file)

        #--- update plasma state files
        services.update_plasma_state()

        #--- archive output files
        services.stage_output_files(timeStamp, self.OUTPUT_FILES)

    def finalize(self, timeStamp=0.0):
        print('pedestal.finalize() called')
