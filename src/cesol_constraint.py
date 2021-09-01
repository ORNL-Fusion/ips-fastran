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

class cesol_constraint(Component):
    def __init__(self, services, config):
        Component.__init__(self, services, config)
        print('Created %s' % (self.__class__))

    def init(self, timeid=0):
        print('cesol_constraint.init() called')

    def step(self, timeid=0):
        print('cesol_constraint.step() started')

        #--- entry
        services = self.services

        #--- stage plasma state files
        services.stage_state()

        #--- get plasma state file name
        cur_state_file = services.get_config_param('CURRENT_STATE')
        cur_instate_file = services.get_config_param('CURRENT_INSTATE')
        cur_eqdsk_file = services.get_config_param('CURRENT_EQDSK')

        #--- stage input files
        services.stage_input_files(self.INPUT_FILES)

        #--- apply constraint to local plasma state
        instate = Namelist(cur_instate_file)

        neped = instate["flux"]["neped"][0]
        nesep = instate["flux"]["nesep"][0]
        tesep = instate["flux"]["tesep"][0]*1.0e-3
        tisep = instate["flux"]["tisep"][0]*1.0e-3

        wped  = instate["flux"]["wped"][0]
        ptop  = instate["flux"]["ptop"][0]
        pped  = instate["flux"]["pped"][0]

        zeffped = 1.7 #----------

        niped, nzped = get_ni(neped, zeff=zeffped)
        teped = 0.5*pped/(1.602*neped)
        tiped = 0.5*pped/(1.602*(niped+nzped))

        print('neped:', neped)
        print('niped:', niped)
        print('nzped:', nzped)
        print('teped:', teped)
        print('tiped:', tiped)

        ps = plasmastate('ips',1)
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
        xtop = 1.0 - 1.5*wped

        ne_model = profile_pedestal(nrho, 1.0-0.5*wped, wped, neped, nesep, ne[0], alpha, beta)
        te_model = profile_pedestal(nrho, 1.0-0.5*wped, wped, teped, tesep, te[0], alpha, beta)
        ti_model = profile_pedestal(nrho, 1.0-0.5*wped, wped, tiped, tisep, ti[0], alpha, beta)

        ne_model_top = ne_model(xtop)
        te_model_top = te_model(xtop)
        ti_model_top = ti_model(xtop)

        ne_top = zinterp(rho, ne)(xtop)
        te_top = zinterp(rho, te)(xtop)
        ti_top = zinterp(rho, ti)(xtop)

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
        services.update_state()

        #--- archive output files
        services.stage_output_files(timeid, self.OUTPUT_FILES)

    def finalize(self, timeid=0.0):
        print('cesol_constraint.finalize() called')
