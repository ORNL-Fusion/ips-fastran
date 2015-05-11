#! /usr/bin/env python

"""
 -----------------------------------------------------------------------
 fastran model equilibrium driver
 JM
 -----------------------------------------------------------------------
"""

import sys,os,shutil
from numpy import *

from component import Component

from Namelist import Namelist
from zplasmastate import plasma_state_file
from zmodelprof import profile_pedestal, spline_profile 
import scipy.optimize

class fastran_modeleq(Component):

    def __init__(self, services, config):

        Component.__init__(self, services, config)
        print 'Created %s' % (self.__class__)

    # ------------------------------------------------------------------
    # init function

    def init(self, timestamp=0):
        print 80*'*'
        print 'DIRVER INIT'
        return

    # ------------------------------------------------------------------
    # step function

    def step(self, timestamp=0):

        services = self.services

        #-- stage input and plasma state files

        services.stage_input_files(self.INPUT_FILES)
        services.stage_plasma_state()

        cur_state_file = services.get_config_param('CURRENT_STATE')
        cur_eqdsk_file = services.get_config_param('CURRENT_EQDSK')
        cur_bc_file = services.get_config_param('CURRENT_BC')

        #-- get Portal RUNID and save to a file

        run_id = services.get_config_param('PORTAL_RUNID')
        sym_root = services.getGlobalConfigParameter('SIM_ROOT')
        path = os.path.join(sym_root, 'PORTAL_RUNID')
        runid_file = open(path, 'a')
        runid_file.writelines(run_id + '\n')
        runid_file.close()

        #-- get list of ports

        ports = services.getGlobalConfigParameter('PORTS')
        port_names = ports['NAMES'].split()
        print 'PORTS =', port_names
        port_dict = {}
        port_id_list = []

        #-- instantiate components in port_names list, except DRIVER itself

        for port_name in port_names:
            if port_name in ["DRIVER"]: continue
            port = services.get_port(port_name) 
            if(port == None):
                print 'Error accessing %s component'%port_name
                raise
            port_dict[port_name] = port
            port_id_list.append(port)

        #-- initialize components in PORTS list for startup

        init_mode = 'init'
        for port_name in port_names:
            if port_name in ['INIT','DRIVER']: continue 
            self.component_call(services,port_name,port_dict[port_name],init_mode,0)
        
        #-- get plasma state files into driver work directory 

        services.stage_plasma_state()
 
        #-- post init processing: stage output

        services.stage_output_files(0, self.OUTPUT_FILES)

        #-- iteration loop

        nstep_eq = int(services.sim_conf["ITERATION_LOOP"]["NSTEP_EQ"])
        nstep = int(services.sim_conf["ITERATION_LOOP"]["NSTEP"])
        print "number of interation :",nstep

        for k in range(nstep_eq):

            print ''
            print 72*"="
            print '= model driver: iteration number = ', k
            print ''
            services.update_time_stamp(k)

            # call step for each component

            if 'EQ0' in port_names:
                self.component_call(services,'EQ0',port_dict['EQ0'],'step',k) 
            if 'TR0' in port_names:
                self.component_call(services,'TR0',port_dict['TR0'],'step',k)
            if 'MONITOR' in port_names:
                self.component_call(services,'MOINTOR',port_dict['MONITOR'],'step',k) 

            services.stage_plasma_state()

            services.stage_output_files(k, self.OUTPUT_FILES)

            io_update_state(cur_state_file,cur_bc_file)

            services.update_plasma_state()


        for k in range(nstep_eq,nstep_eq+nstep):

            print ''
            print 72*"="
            print '= model driver: iteration number = ', k
            print ''
            services.update_time_stamp(k)

            if 'EQ' in port_names:
                self.component_call(services,'EQ',port_dict['EQ'],'step',k)
            if 'EC' in port_names:
                self.component_call(services,'EC',port_dict['EC'],'step',k)
            if 'IC' in port_names:
                self.component_call(services,'IC',port_dict['IC'],'step',k)
            if 'NB' in port_names:
                self.component_call(services,'NB',port_dict['NB'],'step',k)
            if 'TR' in port_names:
                self.component_call(services,'TR',port_dict['TR'],'step',k)
            if 'MONITOR' in port_names:
                self.component_call(services,'MOINTOR',port_dict['MONITOR'],'step',k)

            services.stage_plasma_state()
            services.stage_output_files(k, self.OUTPUT_FILES)
            services.update_plasma_state()
      
        #-- post simulation: call finalize on each component

        print ''
        print 72*"="
        print '= modeleq driver: finalize'
        print ''
        for port_name in port_names:
            if port_name in ['INIT','DRIVER']: continue 
            self.component_call(services, port_name, port_dict[port_name], 'finalize', k)

    # ------------------------------------------------------------------
    # checkpoint function

    def checkpoint(self, timestamp=0):
        print 'fastran_modeleq.checkpoint() called'
        
    # ------------------------------------------------------------------
    # finalize function

    def finalize(self, timestamp = 0):
      # Driver finalize - nothing to be done
        print 80*'*'
        print 'DIRVER FINALIZE'

        sym_root = services.getGlobalConfigParameter('SIM_ROOT')
        outfile = os.path.join(sym_root,'RESULT')
        f = open(outfile,"w")
        f.write("%d scan_id\n"%0)
        f.close()
        pass

    # ----------------------------------------------------------------------
    # "Private" driver methods

    # Component call - wraps the exception handling for all component calls
    def component_call(self, services, port_name, comp, mode, time):
            comp_mode_string = port_name + ' ' + mode
            print (comp_mode_string)
            try:
                services.call(comp, mode, time)
            except Exception:
                services.exception(comp_mode_string + ' failed')
                raise
            else:
                print comp_mode_string + ' finished'
            
            return
    
    # preprocess

    def adjust_profile(self,timeStamp):

        services = self.services
        cur_fastran_file = services.get_config_param('CURRENT_FASTRAN')

#------------
# Local

def io_update_state(f_state,f_inbc):

    f_instate = "instate"

    #------------------------------------------------------------------ 
    # read plasma state

    inbc = Namelist(f_inbc)["inbc"]
    r0 = inbc["r0"][0]
    b0 = abs(inbc["b0"][0])
    ip = inbc['ip'][0]*1.0e6
    ps = plasma_state_file(f_state,r0=r0,b0=b0,ip=ip)

    #------------------------------------------------------------------ 
    # calculate betan

    vol  = ps["vol" ][:]
    area = ps["area"][:]
    ipol = ps["g_eq"][:]/(r0*b0)
    rho  = ps["rho"][:]
    ne   = ps["ns"][0,:]*1.0e-19
    ni   = ps["ni"][:]*1.0e-19
    te   = ps["Ts"][0,:]
    ti   = ps["Ti"][:]
    a0   = ps["rMinor_mean"][-1]

    wth  = 1.5*1.602e3*(ne*te+ni*ti)
    wth_sum = sum([(vol[i+1]-vol[i])*wth[i] for i in range(ps.nrho-1)])

    density_beam = ps.dump_profile(rho,"nbeami",k=0)*1.e-19
    wbeam = ps.dump_profile(rho,"eperp_beami",k=0) \
          + ps.dump_profile(rho,"epll_beami" ,k=0)
    wbeam = 1.602e3*density_beam*wbeam
    wbeam_sum = sum([(vol[i+1]-vol[i])*wbeam[i] for i in range(ps.nrho-1)])

    betan = wth_sum/vol[-1]
    betan*= 2.0*4.0*pi*1.0e-7/b0**2
    betan/= 1.5
    betan/= fabs(ip*1.0e-6)/(a0*abs(b0))
    betan*= 1.0e2

    betan_beam = wbeam_sum/vol[-1]
    betan_beam*= 2.0*4.0*pi*1.0e-7/b0**2
    betan_beam/= 1.5
    betan_beam/= fabs(ip*1.0e-6)/(a0*abs(b0))
    betan_beam*= 1.0e2

    print 'betan_th   =', betan
    print 'betan_beam =', betan_beam

    #------------------------------------------------------------------ 
    # scale : betan

    instate = Namelist(f_instate)["instate"]
    for key in instate.keys(): instate[key] = array(instate[key])

    nrho       = instate["nrho"       ][0]
    n_ion      = instate["n_ion"      ][0]
    z_ion      = instate["z_ion"      ]
    a_ion      = instate["a_ion"      ]
    f_ion      = instate["f_ion"      ]
    n_imp      = instate["n_imp"      ][0]
    z_imp      = instate["z_imp"      ]
    a_imp      = instate["a_imp"      ]
    f_imp      = instate["f_imp"      ]
    xmid       = instate["xmid"       ][0]
    xwid       = instate["xwid"       ][0]
    ne_axis    = instate["ne_axis"    ][0]
    ne_ped     = instate["ne_ped"     ][0]
    ne_sep     = instate["ne_sep"     ][0]
    ne_alpha   = instate["ne_alpha"   ][0]
    ne_beta    = instate["ne_beta"    ][0]
    te_axis    = instate["te_axis"    ][0]
    te_ped     = instate["te_ped"     ][0]
    te_sep     = instate["te_sep"     ][0]
    te_alpha   = instate["te_alpha"   ][0]
    te_beta    = instate["te_beta"    ][0]
    ti_axis    = instate["ti_axis"    ][0]
    ti_ped     = instate["ti_ped"     ][0]
    ti_sep     = instate["ti_sep"     ][0]
    ti_alpha   = instate["ti_alpha"   ][0]
    ti_beta    = instate["ti_beta"    ][0]
    jpar_axis  = instate["jpar_axis"  ][0]
    jpar_sep   = instate["jpar_sep"   ][0]
    jpar_alpha = instate["jpar_alpha" ][0]
    jpar_beta  = instate["jpar_beta"  ][0]
    zeff       = instate["zeff"       ][0]
    nbeam_axis = instate["nbeam_axis" ][0]
    nbeam_sep  = instate["nbeam_sep"  ][0]
    nbeam_alpha= instate["nbeam_alpha"][0]
    nbeam_beta = instate["nbeam_beta" ][0]
    tbeami     = instate["tbeami"     ][0]
    
    rho = arange(nrho)/(nrho-1.0)
    ne = profile_pedestal(nrho,xmid,xwid,ne_ped,ne_sep,ne_axis,ne_alpha,ne_beta)(rho)
    te = profile_pedestal(nrho,xmid,xwid,te_ped,te_sep,te_axis,te_alpha,te_beta)(rho)
    ti = profile_pedestal(nrho,xmid,xwid,ti_ped,ti_sep,ti_axis,ti_alpha,ti_beta)(rho)
    omega = zeros(nrho)
    zeff = zeff*ones(nrho)
    jtot = (jpar_axis-jpar_sep)*(1.0-rho**jpar_alpha)**jpar_beta+jpar_sep
    density_beam = (nbeam_axis-nbeam_sep)*(1.0-rho**nbeam_alpha)**nbeam_beta+nbeam_sep
    density_alpha = zeros(nrho)
    tbeami = tbeami*ones(nrho)

    betan_target  = instate["betan"][0]
    betan_beam_target  = instate["betan_beam"][0]

    teaxis  = te[0]*betan_target/betan  
    tiaxis  = ti[0]*betan_target/betan  
    te = profile_pedestal(nrho,xmid,xwid,te_ped,te_sep,te_axis,te_alpha,te_beta)(rho)
    ti = profile_pedestal(nrho,xmid,xwid,ti_ped,ti_sep,ti_axis,ti_alpha,ti_beta)(rho)

    if betan_beam_target > 0.0:
        density_beam *=  betan_beam0/betan_beam

    ps["Ts"][0] =ps.node2cell(te)
    for k in range(n_ion): ps["Ts"][k+1] = ps.node2cell(ti)
    for k in range(n_imp): ps["Ts"][k+n_ion+1] = ps.node2cell(ti)
    ps["Ti"][:] = ps.node2cell(ti)

    ps["nbeami"][0] = 1.0e19*ps.node2cell(density_beam)

    #------------------------------------------------------------------ 
    # scale : current

    rho_jpeak = instate["rho_jpeak"][0]    
    rho_jbdry = instate["rho_jbdry"][0]    
    jaxis     = instate["jaxis"][0]      
    jpeak     = instate["jpeak"][0]      
    jaxis1    = instate["jaxis_prime"][0]
    jbdry1    = instate["jbdry_prime"][0]

    j_bs  = 1.0e-6*ps.dump_j_parallel_CD(rho,"bs")
    jbdry = j_bs[where( rho>=rho_jbdry)[0][0]]

    def _cal_ip(x):

         x0 = [0.0,rho_jpeak,rho_jbdry]
         y0 = [jaxis/jpeak*x,x,jbdry]
         yp = [jaxis1,0.0,jbdry1]
         spl = spline_profile(x=x0,y=y0,yp=yp)

         j_tot = zeros(nrho)
         for i in range(nrho):
             if rho[i]<rho_jbdry: j_tot[i] = spl(rho[i])
             else: j_tot[i] = j_bs[i]

         curt = zeros(nrho)
         for i in range(nrho-1):
             dV = vol[i+1]-vol[i]
             jparm = 0.5*(j_tot[i]+j_tot[i+1])
             ipolm = 0.5*(ipol[i]+ipol[i+1])
             curt[i+1] = (curt[i]/ipol[i]+jparm*dV/(2.0*pi*r0*ipolm**2))*ipol[i+1]
         return curt[-1]-ip*1.e-6

    x_find = scipy.optimize.bisect(
        _cal_ip, 0.5, 10.0, args=(),
        xtol=1e-12, rtol=1.e-3, maxiter=100, full_output=False, disp=True)
    print 'j_peak =', x_find

    x0 = [0.0,rho_jpeak,rho_jbdry]
    y0 = [jaxis/jpeak*x_find,x_find,jbdry]
    yp = [jaxis1,0.0,jbdry1]
    spl = spline_profile(x=x0,y=y0,yp=yp)

    j_tot = zeros(nrho)
    for i in range(nrho):
        if rho[i]<rho_jbdry: j_tot[i] = spl(rho[i])
        else: j_tot[i] = j_bs[i]

    ps.load_j_parallel(j_tot)

    #------------------------------------------------------------------ 
    # close

    ps.close()
