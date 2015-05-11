#! /usr/bin/env python

"""
 -----------------------------------------------------------------------
 fastran model equilibrium init component
 JM
 -----------------------------------------------------------------------
"""

import sys,os,os.path,shutil,pickle
import subprocess
from numpy import *

from component import Component

from Namelist import Namelist
import zefit

import netCDF4
import zefitutil
from zplasmastate import plasma_state_file,instate2ps
from zmodelprof import profile_pedestal
from zefit import fixbdry_kfile_init

class modeleq_init(Component):

    def __init__(self, services, config):

        Component.__init__(self, services, config)
        print 'Created %s' % (self.__class__)

    def init(self, timestamp=0.0):

        print ('modeleq_init.init() called')
        services = self.services

        input_dir_fastran = services.get_config_param('INPUT_DIR_FASTRAN')
        input_dir_sim = services.get_config_param('INPUT_DIR_SIM')

        print 'INPUT_DIR_FASTRAN = ',input_dir_fastran
        print 'INPUT_DIR_SIM = ',input_dir_sim

        if not os.path.exists(input_dir_sim):
            print 'making directory', input_dir_sim
            os.mkdir(input_dir_sim)

        try:
            input_id = int(float(self.INPUT_ID))
        except:
            input_id = 0
        print 'input_id = ',input_id


        for input_file in self.INPUT_FILES.split():
 
            src = os.path.join(input_dir_fastran,input_file)
            if input_id > 0:
                src = src+".%d"%input_id
            target = os.path.join(input_dir_sim,input_file)
            print 'src    :',src
            print 'target :',target
            shutil.copyfile (src,target)

    def step(self, timeStamp):

        #----------------------------------------------------------
        #-- entry

        print ('modeleq_init.step() called')
        services = self.services

        tokamak_id = services.get_config_param('TOKAMAK_ID')
        shot_number = services.get_config_param('SHOT_NUMBER')
        run_id = services.get_config_param('RUN_ID')

        #----------------------------------------------------------
        #-- stage input files

        services.stage_input_files(self.INPUT_FILES)

        #----------------------------------------------------------
        #-- get plasma state file names

        cur_state_file = services.get_config_param('CURRENT_STATE')
        cur_eqdsk_file = services.get_config_param('CURRENT_EQDSK')
        cur_bc_file = services.get_config_param('CURRENT_BC')

        #----------------------------------------------------------
        #-- read instate file

        instate = Namelist("instate")["instate"]
        nrho = instate["nrho"][0]
        rho = arange(nrho)/(nrho-1.0)

        #----------------------------------------------------------
        #-- get pstool excutable name

        pstool_bin = os.path.join(self.BIN_PATH, 'pstool')

        #----------------------------------------------------------
        #-- allocate initial plasma state file

        inps = Namelist()
        inps["inps"]["global_label"] = ['fastran_model_equilibrium']
        inps["inps"]["tokamak_id"  ] = [tokamak_id]
        inps["inps"]["runid"       ] = [run_id]
        inps["inps"]["shot_number" ] = [int(shot_number)]
        inps["inps"]["time"        ] = [0.0]
        inps["inps"]["nspec_ion"   ] = instate["n_ion"]
        inps["inps"]["nspec_imp"   ] = instate["n_imp"]
        inps["inps"]["nspec_beam"  ] = [1]
        inps["inps"]["nspec_fusion"] = [1]
        inps["inps"]["nspec_rfmin" ] = [0]
        inps["inps"]["nspec_gas"   ] = instate["n_ion"]
        inps["inps"]["z_ion"       ] = instate["z_ion"]
        inps["inps"]["a_ion"       ] = instate["a_ion"]
        inps["inps"]["z_imp"       ] = instate["z_imp"]
        inps["inps"]["a_imp"       ] = instate["a_imp"]
        inps["inps"]["z_beam"      ] = [1]
        inps["inps"]["a_beam"      ] = [2]
        inps["inps"]["z_fusion"    ] = [2]
        inps["inps"]["a_fusion"    ] = [4]
        inps["inps"]["z_gas"       ] = instate["z_ion"]
        inps["inps"]["a_gas"       ] = instate["a_ion"]
        inps["inps"]["nrho"        ] = [nrho]
        inps["inps"]["nrho_eq"     ] = [nrho]
        inps["inps"]["nth_eq"      ] = [101]
        inps["inps"]["nrho_eq_geo" ] = [nrho]
        inps["inps"]["nrho_gas"    ] = [nrho]
        inps["inps"]["nrho_nbi"    ] = [nrho]
        inps["inps"]["nrho_ecrf"   ] = [nrho]
        inps["inps"]["nrho_icrf"   ] = [nrho]
        inps["inps"]["nrho_fus"    ] = [nrho]
        inps["inps"]["nrho_anom"   ] = [nrho]
        inps.write("inps")

        print 'pstool init'
        logfile=open("pstool_init.log","w")
        retcode = subprocess.call([pstool_bin, "init"],stdout=logfile,stderr=logfile)
        logfile.close()
        if (retcode != 0):
           raise Exception('Error executing ', pstool_bin)

        #----------------------------------------------------------
        #-- run initial equilibrium

        inefit = Namelist()
        inefit["inefit"]["rho"  ] = rho
        inefit["inefit"]["press"] = zeros(nrho)
        inefit["inefit"]["jpar" ] = zeros(nrho)
        inefit["inefit"]['ip'   ] = [instate["ip"][0]*1.0e6]
        inefit["inefit"]['b0'   ] = instate["b0"]
        inefit["inefit"]['r0'   ] = instate["r0"]
        inefit["inefit"]["nlim" ] = instate["nlim" ]
        inefit["inefit"]["rlim" ] = instate["rlim" ]
        inefit["inefit"]["zlim" ] = instate["zlim" ]
        inefit["inefit"]["nbdry"] = instate["nbdry"]
        inefit["inefit"]["rbdry"] = instate["rbdry"]
        inefit["inefit"]["zbdry"] = instate["zbdry"]
        inefit.write("inefit")

        shot = int(shot_number)
        time = 0
        fixbdry_kfile_init(shot,time,"inefit")

        try:
            efit_bin = os.path.join(self.BIN_PATH, 'efitd90 129 129')
        except:
            efit_bin = 'efitd90 129 129'

        print "run efit"

        kfile = "k%06d.%05d"%(shot,time)
        args = "2\n 1\n "+kfile
        command = 'echo \"%s\"'%args + ' | ' + efit_bin

        f=open("xefit","w")
        f.write(command)
        f.close()

        cwd = services.get_working_dir()
        task_id = services.launch_task(1, cwd, "sh xefit", logfile='efit.log')
        retcode = services.wait_task(task_id)
        if (retcode != 0):
           print 'Error executing ', 'efit'
           raise

        shutil.copyfile("g%06d.%05d"%(shot,time), cur_eqdsk_file)

        #----------------------------------------------------------
        #-- load geqdsk to plasma state

        print 'pstool load geqdsk'

        shutil.copyfile("g%06d.%05d"%(shot,time), "geqdsk")
        logfile = open('pstool_geqdsk.log', 'w')
        retcode = subprocess.call([pstool_bin, "load","geqdsk", "1.0d-6"],
                      stdout=logfile,stderr=logfile)
        if (retcode != 0):
           print 'Error executing ', pstool_bin
           raise
        logfile.close()

        #----------------------------------------------------------
        #-- load innubeam to plasma state

        if "innubeam" in self.INPUT_FILES:
            print 'pstool load innubeam'
            logfile=open("pstool_innubeam.log","w")
            retcode = subprocess.call([pstool_bin, "load", "innubeam"],
                          stdout=logfile,stderr=logfile)
            logfile.close()
            if (retcode != 0):
               print 'Error executing ', pstool_bin
               raise

        #----------------------------------------------------------
        #-- compose inital plasma state profile

        print 'initial plasma state'

        init_plasmastate("instate","ps.nc" )
        shutil.copyfile("ps.nc",cur_state_file)
        #r0  = instate["r0"][0]
        #b0  = abs(instate["b0"][0])
        #ip  = instate["ip"][0]*1.0e6
        #ps = plasma_state_file("ps.nc",r0=r0,b0=b0,ip=ip)
        #instate2ps(instate,ps)
        #shutil.copyfile("ps.nc",cur_state_file)

        #----------------------------------------------------------
        #-- boundary condition plasma state file

        inbc = Namelist()
        inbc["inbc"]["r0"   ] = instate["r0"]
        inbc["inbc"]["b0"   ] = instate["b0"]
        inbc["inbc"]["ip"   ] = instate["ip"]
        inbc["inbc"]["nbdry"] = instate["nbdry"]
        inbc["inbc"]["rbdry"] = instate["rbdry"]
        inbc["inbc"]["zbdry"] = instate["zbdry"]
        inbc.write(cur_bc_file)

        # --------------------------------------------------------------
        # Update plasma state

        services.update_plasma_state()

        # --------------------------------------------------------------
        # archive output files

        services.stage_output_files(timeStamp, self.OUTPUT_FILES)

        #io_update_state(cur_state_file,cur_bc_file)

    def checkpoint(self, timeStamp=0):

        print 'modeleq_init.checkpoint() called'

        services = self.services
        services.stage_plasma_state()
        services.save_restart_files(timeStamp, self.RESTART_FILES)

    def finalize(self, timeStamp=0):

        print 'modeleq_init.finalize() called'

#------------
# Local

def init_plasmastate(f_instate, f_ps):

    #------------------------------------------------------------------
    # read instate and compose profiles

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
    omega = 150.0e3*(1.0-rho**1.5)**1.5
    zeff = zeff*ones(nrho)
    jtot = (jpar_axis-jpar_sep)*(1.0-rho**jpar_alpha)**jpar_beta+jpar_sep
    density_beam = (nbeam_axis-nbeam_sep)*(1.0-rho**nbeam_alpha)**nbeam_beta+nbeam_sep
    density_alpha = zeros(nrho)
    tbeami = tbeami*ones(nrho)

    #------------------------------------------------------------------
    # load ps

    r0  = instate["r0"][0]
    b0  = abs(instate["b0"][0])
    ip  = instate["ip"][0]*1.0e6

    ps = plasma_state_file(f_ps,r0=r0,b0=b0,ip=ip)

    #------------------------------------------------------------------
    # density

    a = sum(f_ion*z_ion)
    b = sum(f_imp*z_imp)
    c = sum(f_ion*z_ion)
    d = sum(f_imp*z_imp*z_imp)

    zne_adj = ne
    zzne_adj = ne*zeff

    # depletion due to beam ions

    zne_adj = zne_adj - 1.0*density_beam
    zzne_adj = zzne_adj - 1.0**2*density_beam

    # effective main ion and impurity densities

    nion = (zne_adj *d-zzne_adj*b)/(a*d-b*c)
    nimp = (zzne_adj*a-zne_adj *c)/(a*d-b*c)

    # density ion

    density_ion = array([f_ion[k]*nion for k in range(n_ion)])
    density_imp = array([f_imp[k]*nimp for k in range(n_imp)])

    density_th = array([sum(tmp) for tmp in density_ion.transpose()]) 
    density_th+= array([sum(tmp) for tmp in density_imp.transpose()]) 

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
    # beam

    ps["nbeami"][0] = 1.0e19*ps.node2cell(density_beam)
    ps["eperp_beami"][0] = 2.0*ps.node2cell(tbeami)
    ps["epll_beami"][0] = ps.node2cell(tbeami)

    #--------------------------------------------------------------
    # current

    ps.load_j_parallel(1.0e6*jtot)

    ps.load_j_parallel_CD(rho,zeros(nrho),"nb")
    ps.load_j_parallel_CD(rho,zeros(nrho),"ec")
    ps.load_j_parallel_CD(rho,zeros(nrho),"ic")
    ps.load_j_parallel_CD(rho,zeros(nrho),"bs")
    ps.load_j_parallel_CD(rho,zeros(nrho),"oh")

    #-------------------------------------------------------------
    # heating

    ps.load_profile (rho,zeros(nrho),"pbe","vol")
    ps.load_profile (rho,zeros(nrho),"pbi","vol")
    ps.load_profile (rho,zeros(nrho),"peech","vol")
    ps.load_profile (rho,zeros(nrho),"picrf_totals","vol",k=0)
    ps.load_profile (rho,zeros(nrho),"picth","vol")

    #-------------------------------------------------------------
    # update plasma state

    ps.close()
