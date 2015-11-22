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
from zmodelprof import profile_pedestal, spline_profile, profile_hat 
import scipy.optimize
from zinterp import zinterp

class fastran_modeleq(Component):

    def __init__(self, services, config):

        Component.__init__(self, services, config)
        print 'Created %s' % (self.__class__)

    # ------------------------------------------------------------------
    # init function

    def init(self, timestamp=0):
        print 80*'*'
        print 'DIRVER INIT'
        self.kiter = 0
        self.rho_hat_now = -1.0
        self.rho_hat_list = []
        self.qmin_list = []
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

        # run_id = services.get_config_param('PORTAL_RUNID')
        # sym_root = services.getGlobalConfigParameter('SIM_ROOT')
        # path = os.path.join(sym_root, 'PORTAL_RUNID')
        # runid_file = open(path, 'a')
        # runid_file.writelines(run_id + '\n')
        # runid_file.close()

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

        #-- load input

        f_instate = "instate"
        self.instate = Namelist(f_instate)["instate"]

        #-- iteration loop

        nstep_eq = int(services.sim_conf["ITERATION_LOOP"]["NSTEP_EQ"])
        nstep = int(services.sim_conf["ITERATION_LOOP"]["NSTEP"])
        print "number of interation :",nstep

        for k in range(nstep_eq):
    
            print ''
            print 72*"="
            print '= model driver: iteration number = ', k
            print ''
            
            time_id = k

            services.update_time_stamp(time_id)
    
            # call step for each component
    
            if 'PED' in port_names:
                self.component_call(services,'PED',port_dict['PED'],'step',time_id) 
            if 'EQ0' in port_names:
                self.component_call(services,'EQ0',port_dict['EQ0'],'step',time_id) 
            if 'TR0' in port_names:
                self.component_call(services,'TR0',port_dict['TR0'],'step',time_id)
            if 'MONITOR' in port_names:
                self.component_call(services,'MOINTOR',port_dict['MONITOR'],'step',time_id) 
    
            services.stage_plasma_state()
    
            services.stage_output_files(time_id, self.OUTPUT_FILES)
    
            if k < nstep_eq-1:
                iconv = io_update_state(self,cur_state_file,cur_bc_file,0.5)
                print 'iconv = ',iconv
            #if iconv > 0: break
        
            services.update_plasma_state()

        if 'STAB' in port_names:
            self.component_call(services,'STAB',port_dict['STAB'],'step',time_id) 

        for k in range(nstep_eq,nstep_eq+nstep):
    
            print ''
            print 72*"="
            print '= model driver: iteration number = ', k
            print ''
            
            time_id = k

            services.update_time_stamp(time_id)
    
            # call step for each component
    
            iStamp = time_id

            if 'EC' in port_names:
                self.component_call(services,'EC',port_dict['EC'],'step',iStamp) #t)
            if 'IC' in port_names:
                self.component_call(services,'IC',port_dict['IC'],'step',iStamp) #t)
            if 'NB' in port_names:
                self.component_call(services,'NB',port_dict['NB'],'step',iStamp) #t)
            if 'EQ' in port_names:
                self.component_call(services,'EQ',port_dict['EQ'],'step',iStamp) #t)
            if 'TR' in port_names:
                self.component_call(services,'TR',port_dict['TR'],'step',iStamp) #t)
            if 'MONITOR' in port_names:
                self.component_call(services,'MOINTOR',port_dict['MONITOR'],'step',iStamp) #t)
    
            services.stage_plasma_state()
            services.stage_output_files(time_id, self.OUTPUT_FILES)
            services.update_plasma_state()

        #-- post simulation: call finalize on each component

        print ''
        print 72*"="
        print '= modeleq driver: finalize'
        print ''
        for port_name in port_names:
            if port_name in ['INIT','DRIVER']: continue 
            self.component_call(services, port_name, port_dict[port_name], 'finalize', time_id)

        self.finalize(time_id)

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

        services = self.services
        sym_root = services.getGlobalConfigParameter('SIM_ROOT')
        outfile = os.path.join(sym_root,'RESULT')
        f = open(outfile,"w")
       #f.write("%d scan_id\n"%0)
        f.write("%f rho_jhat\n"%self.rho_hat_now)
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

def io_update_state(self,f_state,f_inbc,rho_hat=-1,nmax_iter=100):

    #------------------------------------------------------------------ 
    # entry

    print "\n"+72*"="
    print '= model driver: adjust profiles '

    f_instate = "instate"
    instate = Namelist(f_instate)["instate"]
    iconv = 0

    #------------------------------------------------------------------ 
    # read plasma state

    inbc = Namelist(f_inbc)["inbc"]
    r0 = inbc["r0"][0]
    b0 = abs(inbc["b0"][0])
    ip = inbc['ip'][0]*1.0e6
    ps = plasma_state_file(f_state,r0=r0,b0=b0,ip=ip)

    vol  = ps["vol" ][:]
    area = ps["area"][:]
    ipol = ps["g_eq"][:]/(r0*b0)
    rho  = ps["rho"][:]
    ne   = ps["ns"][0,:]*1.0e-19
    ni   = ps["ni"][:]*1.0e-19
    te   = ps["Ts"][0,:]
    ti   = ps["Ti"][:]
    a0   = ps["rMinor_mean"][-1]

    qmhd = ps["q_eq"][:]

    #------------------------------------------------------------------ 
    # calculate betan

    wth  = 1.5*1.602e3*(ne*te+ni*ti)
    wth_sum = sum([(vol[i+1]-vol[i])*wth[i] for i in range(ps.nrho-1)])

    density_beam = ps.dump_profile(rho,"nbeami",k=0)*1.e-19
    wbeam = ps.dump_profile(rho,"eperp_beami",k=0) \
          + ps.dump_profile(rho,"epll_beami" ,k=0)
    wbeam = 1.602e3*density_beam*wbeam
    wbeam_sum = sum([(vol[i+1]-vol[i])*wbeam[i] for i in range(ps.nrho-1)])

    betan_th = wth_sum/vol[-1]
    betan_th*= 2.0*4.0*pi*1.0e-7/b0**2
    betan_th/= 1.5
    betan_th/= fabs(ip*1.0e-6)/(a0*abs(b0))
    betan_th*= 1.0e2

    betan_beam = wbeam_sum/vol[-1]
    betan_beam*= 2.0*4.0*pi*1.0e-7/b0**2
    betan_beam/= 1.5
    betan_beam/= fabs(ip*1.0e-6)/(a0*abs(b0))
    betan_beam*= 1.0e2

    print 'betan_th   =', betan_th
    print 'betan_beam =', betan_beam
    print 'a0 = ', a0

    #------------------------------------------------------------------ 
    # scale : betan

    for key in instate.keys(): instate[key] = array(instate[key])

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

    ptot_alpha = instate["ptot_alpha" ][0]
    ptot_beta  = instate["ptot_beta"  ][0]

    #------------------------------------------------------------------ 
    # grid and ..
    
    nrho = instate["nrho"][0]
    rho = arange(nrho)/(nrho-1.0)

    n_ion = instate["n_ion"][0]
    n_imp = instate["n_imp"][0]

    xmid  = instate["xmid"][0]
    xwid  = instate["xwid"][0]

    ne = ps.cell2node(ne)
    ni = ps.cell2node(ni)
    te = ps.cell2node(te)
    ti = ps.cell2node(ti)

    #------------------------------------------------------------------ 
    # scale : thermal pressure
 
    betan_target = instate["betan_th"][0]

    if betan_target > 0.0:

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

        tscale_min = 0.0
        tscale_max = 10.0

        for iter_scale in range(nmax_iter):

            tscale = 0.5*(tscale_min+tscale_max)

            te = profile_pedestal(nrho,xmid,xwid,te_ped,te_sep,tscale*te_axis,te_alpha,te_beta)(rho)
            ti = profile_pedestal(nrho,xmid,xwid,ti_ped,ti_sep,tscale*ti_axis,ti_alpha,ti_beta)(rho)

            pth  = 1.602e3*(ne*te+ni*ti)
            pth_sum = sum([(vol[i+1]-vol[i])*(pth[i+1]+pth[i])*0.5 for i in range(ps.nrho-1)])
     
            betan_th = pth_sum/vol[-1]
            betan_th*= 2.0*4.0*pi*1.0e-7/b0**2
            betan_th/= fabs(ip*1.0e-6)/(a0*abs(b0))
            betan_th*= 1.0e2
     
            if abs(betan_th-betan_target) < 0.01: break
     
            if betan_th > betan_target:
               tscale_max = tscale
            else:
               tscale_min = tscale
        
            print "scale to match betan_th : %5.3f,%5.3f"%(tscale,betan_th)
    else:

        pth  = 1.602e3*(ne*te+ni*ti)


   #pth_ped = pth[where( rho>=xmid-xwid/2)[0][0]]
    pth_ped = zinterp(rho,pth)(xmid-xwid/2)

    ps["Ts"][0] =ps.node2cell(te)
    for k in range(n_ion): ps["Ts"][k+1] = ps.node2cell(ti)
    for k in range(n_imp): ps["Ts"][k+n_ion+1] = ps.node2cell(ti)
    ps["Ti"][:] = ps.node2cell(ti)

    print 'Temperature scale:',te[0],ti[0],betan_target/betan_th,tscale

    #------------------------------------------------------------------ 
    # scale : total pressure

    try:
        pressure_model =  instate["pressure_model"][0].strip().lower()
    except:
        pressure_model = "total"
    print "pressure_model = ",pressure_model

    betan_target = instate["betan"][0]

    if pressure_model == "kinetic":

        p_axis_min = 0.0
        p_axis_max = 1000.0e3
    
        pmhd = zeros(nrho)
    
        for iter_scale in range(nmax_iter):
    
            p_axis = 0.5*(p_axis_min+p_axis_max)
    
            xped = xmid-xwid/2
            for k,xval in enumerate(rho):
                if 1.0-xval/xped > 0.0: pmhd[k] = p_axis*(1.0-(xval/xped)**ptot_alpha)**ptot_beta + pth_ped
                else: pmhd [k] = pth[k]
    
            pmhd_sum = sum([(vol[i+1]-vol[i])*(pmhd[i+1]+pmhd[i])*0.5 for i in range(ps.nrho-1)])
            betan = pmhd_sum/vol[-1]
            betan*= 2.0*4.0*pi*1.0e-7/b0**2
            betan/= fabs(ip*1.0e-6)/(a0*abs(b0))
            betan*= 1.0e2
    
            if abs(betan-betan_target) < 0.01: break
    
            if betan > betan_target:
               p_axis_max = p_axis
            else:
               p_axis_min = p_axis
    
            print "scale to match betan_th : %5.3f,%5.3f"%(p_axis*1.0e-3,betan)

        print 'p_axis = ',p_axis
        print 'betan = ',betan
    
    elif pressure_model == "total":

        ptot_axis    = instate["ptot_axis" ][0]
        ptot_ped     = instate["ptot_ped"  ][0]
        ptot_sep     = instate["ptot_sep"  ][0]
        ptot_alpha   = instate["ptot_alpha"][0]
        ptot_beta    = instate["ptot_beta" ][0]

        p_axis_min = 0.0 
        p_axis_max = 100.0*ptot_axis
    
        pmhd = zeros(nrho)
    
        for iter_scale in range(100):
    
            p_axis = 0.5*(p_axis_min+p_axis_max)
    
            pmhd = profile_pedestal(nrho,xmid,xwid,ptot_ped,ptot_sep,p_axis,ptot_alpha,ptot_beta)(rho)

            pmhd_sum = sum([(vol[i+1]-vol[i])*(pmhd[i+1]+pmhd[i])*0.5 for i in range(ps.nrho-1)])
            betan = pmhd_sum/vol[-1]
            betan*= 2.0*4.0*pi*1.0e-7/b0**2
            betan/= fabs(ip*1.0e-6)/(a0*abs(b0))
            betan*= 1.0e2
    
            if abs(betan-betan_target) < 0.01: break
    
            if betan > betan_target:
               p_axis_max = p_axis
            else:
               p_axis_min = p_axis

            print "scale to match betan_th : %5.3e,%5.3f"%(p_axis*1.0e-3,betan)
    
        print 'p_axis = ',p_axis
        print 'betan = ',betan

    else:
        print "check pressumre_model"
        raise

    ps["P_eq"][:] = pmhd

    #------------------------------------------------------------------ 
    # scale : current

    current_model = instate["current_model"][0].strip().lower()
    print 'current_mode = ',current_model

    try:
        cboot = instate["cboot"][0]
    except:
        cboot = 1.0
    print 'cboot = ',cboot

    try:
        q0 = instate["q0"][0]
    except:
        q0 = -1.0
    print 'q0 = ', q0

    if current_model == "hat":

        rho_jbdry = instate["rho_jbdry"][0]    
        rho_jhat  = instate["rho_jhat"][0]    
        wid_jhat  = instate["wid_jhat"][0]    

        if q0>0.0:

            if self.rho_hat_now < 0.0:
                rho_jhat  = instate["rho_jhat"][0]  
            else:
                self.rho_hat_list.append(self.rho_hat_now)
                self.qmin_list.append(qmhd[0])
                if len(self.rho_hat_list) < 2:
                    if qmhd[0] > 1.0:
                        rho_jhat = self.rho_hat_now*0.9
                    else:
                        rho_jhat = self.rho_hat_now*1.1
                else:
                    rho_1 = self.rho_hat_list[-1]
                    rho_2 = self.rho_hat_list[-2]
                    qmin_1 = self.qmin_list[-1]
                    qmin_2 = self.qmin_list[-2]
                    rho_jhat = (q0-qmin_1)/(qmin_2-qmin_1)*(rho_2-rho_1)+rho_1
                    print '************',rho_1,rho_2,qmin_1,qmin_2
                    print '************',rho_jhat
            self.rho_hat_now  = rho_jhat

        j_bs  = cboot*1.0e-6*ps.dump_j_parallel_CD(rho,"bs")
        jbdry = j_bs[where( rho>=rho_jbdry)[0][0]]
    
        jpeak_min = 0.1
        jpeak_max = 10.0
    
        j_tot = zeros(nrho)
        curt = zeros(nrho)
    
        for iter_scale in range(100):
    
            jpeak = 0.5*(jpeak_min+jpeak_max)   
    
            j_hat = profile_hat(nrho,rho_jhat,wid_jhat)
            j_tmp = j_hat[where(rho>=rho_jbdry)[0][0]]
    
            for i in range(nrho):
                if rho[i]<rho_jbdry: j_tot[i] = jpeak*(j_hat(rho[i])-j_tmp)+jbdry
                else: j_tot[i] = j_bs[i]
    
            for i in range(nrho-1):
                dV = vol[i+1]-vol[i]
                jparm = 0.5*(j_tot[i]+j_tot[i+1])
                ipolm = 0.5*(ipol[i]+ipol[i+1])
                curt[i+1] = (curt[i]/ipol[i]+jparm*dV/(2.0*pi*r0*ipolm**2))*ipol[i+1]
    
            if abs(curt[-1]-ip*1.e-6) < 0.001: break
    
            if curt[-1] > ip*1.e-6:
               jpeak_max = jpeak
            else:
               jpeak_min = jpeak
    
        print 'jpeak =',jpeak
        print 'ip_cal = ',curt[-1]
    
        ps.load_j_parallel(j_tot)

    elif current_model == "broad":

        rho_jpeak = instate["rho_jpeak"][0]    
        rho_jbdry = instate["rho_jbdry"][0]    
        jaxis     = instate["jaxis"][0]      
        jpeak     = instate["jpeak"][0]      
        jaxis1    = instate["jaxis_prime"][0]
        jbdry1    = instate["jbdry_prime"][0]

        j_bs  = cboot*1.0e-6*ps.dump_j_parallel_CD(rho,"bs")
        jbdry = j_bs[where( rho>=rho_jbdry)[0][0]]
    
        jpeak_min = 0.1
        jpeak_max = 10.0
    
        j_tot = zeros(nrho)
        curt = zeros(nrho)
    
        for iter_scale in range(100):
    
            x = 0.5*(jpeak_min+jpeak_max)   

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

            if abs(curt[-1]-ip*1.e-6) < 0.001: break
    
            if curt[-1] > ip*1.e-6:
               jpeak_max = x
            else:
               jpeak_min = x

        print 'jpeak =',jpeak
        print 'ip_cal = ',curt[-1]
    
        ps.load_j_parallel(j_tot)

    else:
        print "check current_model"
        raise
        
    #------------------------------------------------------------------ 
    # close

    ps.close()

    if abs(betan-betan_target)/betan_target < 0.005: iconv = 1
    else: iconv = 0

    return iconv

