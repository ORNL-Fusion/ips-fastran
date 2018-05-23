#! /usr/bin/env python

"""
 -----------------------------------------------------------------------
 fastran init component
 -----------------------------------------------------------------------
"""

import os
import shutil
from numpy import *

from component import Component

from Namelist import Namelist
import efit_io
from plasmastate import plasmastate
import instate_model
import instate_io
import cesol_io
import efit_io
import subprocess

class fastran_init (Component):

    def __init__(self, services, config):

        Component.__init__(self, services, config)
        print 'Created %s' % (self.__class__)

    def init(self, timestamp=0.0):

        print ('fastran_init.init() called')

    def step(self, timeStamp):

        #-- entry

        print ('fastran_init.step() called')
        services = self.services

        #-- run identifiers

        tokamak_id = services.get_config_param('TOKAMAK_ID')
        shot_number = services.get_config_param('SHOT_NUMBER')

        print 'tokamak_id =', tokamak_id
        print 'shot_number =',shot_number

        #-- stage input files

        input_dir_id = getattr(self,"INPUT_DIR_ID","")
        if input_dir_id:
            services.component_ref.config['INPUT_DIR'] = services.component_ref.config['INPUT_DIR']+"_%d"%int(float(input_dir_id))
        print 'INPUT_DIR =',services.component_ref.config['INPUT_DIR']

        services.stage_input_files(self.INPUT_FILES)

        #-- plasma state file names

        cur_state_file = services.get_config_param('CURRENT_STATE')
        cur_eqdsk_file = services.get_config_param('CURRENT_EQDSK')
        cur_bc_file = services.get_config_param('CURRENT_BC')
        cur_instate_file = services.get_config_param('CURRENT_INSTATE')
        try:
            cur_fastran_file = services.get_config_param('CURRENT_FASTRAN')
        except:
            cur_fastran_file = ''

        #-- set default

        init_method = getattr(self,'INIT_METHOD','instate')
        instate_method = getattr(self,'INSTATE_METHOD','')
        f_instate = getattr(self,'INSTATE','')
        f_inps = getattr(self,'INPS','')
        f_ingeqdsk = getattr(self,'INGEQDSK','')
        f_eped = getattr(self,'EPED','')

        #-- instate / dakota binding

        instate = Namelist(f_instate)

        var_list = [ var.upper() for var in instate["instate"].keys()]
        for key in var_list:

            if hasattr(self,key):
                instate["instate"][key][0] = float(getattr(self, key))
                print key, 'updated'

            if hasattr(self, "SCALE_"+key):
                scale = float(getattr(self, "SCALE_"+key))
                instate["instate"][key] = scale*array(instate["instate"][key])
                print key,'scaled',scale

        var_list = [ var.upper() for var in instate["flux"].keys()]
        for key in var_list:

            if hasattr(self,key):
                instate["flux"][key][0] = float(getattr(self, key))
                print key, 'updated'

        instate.write(f_instate)

        #-- expand instate

        if instate_method != '':

            instate_model.instate_model(f_instate,instate_method)

        #-- alloc plasma state file

        ps = plasmastate('ips',1)

        instate = Namelist(f_instate)

        print 'init from instate:', f_instate

        ps.init(
            cur_state_file,
            global_label = ['ips'],
            runid = ['ips'],
            shot_number = [int(shot_number)],
            nspec_th = [instate["instate"]["n_ion"][0] +  instate["instate"]["n_imp"][0]] ,
            nspec_beam = [ instate["instate"]["n_beam"][0] ],
            nspec_fusion = [ instate["instate"]["n_fusion"][0] ],
            nspec_rfmin = [ instate["instate"]["n_min"][0] ],
            nspec_gas = [ instate["instate"]["n_ion"][0] ],
            z_ion = instate["instate"]["z_ion"] + instate["instate"]["z_imp"],
            a_ion = instate["instate"]["a_ion"] + instate["instate"]["a_imp"],
            z_beam = instate["instate"]["z_beam"],
            a_beam = instate["instate"]["a_beam"],
            z_fusion = instate["instate"]["z_fusion"],
            a_fusion = instate["instate"]["a_fusion"],
            z_rfmin = instate["instate"]["z_min"],
            a_rfmin = instate["instate"]["a_min"],
            z_gas = instate["instate"]["z_ion"],
            a_gas = instate["instate"]["a_ion"],
            nrho = [101],
            time = [0.0],
            nlim = instate["instate"]["nlim"]
        )

        #-- input equlibrium

        if f_ingeqdsk:

            ps.load_geqdsk(f_ingeqdsk,keep_cur_profile=False)
            shutil.copyfile(f_ingeqdsk,cur_eqdsk_file)

        else:

            # r0 = instate["instate"]["r0"][0]
            # b0 = instate["instate"]["b0"][0]
            #
            # rb = array(instate["instate"]["rbdry"])
            # zb = array(instate["instate"]["zbdry"])
            #
            # R = 0.5*( max(rb) + min(rb) )
            # Z = 0.5*( max(zb) + min(zb) )
            # a = 0.5*( max(rb) - min(rb) )
            # kappa = 0.5*( ( max(zb) - min(zb) )/ a )
            # delta_u = ( R - rb[argmax(zb)] )/a
            # delta_l = ( R - rb[argmin(zb)] )/a
            # delta = 0.5*( delta_u + delta_l )
            #
            # ps.analytic_volume(b0, r0, a, kappa, delta)
            #
            # open(cur_eqdsk_file,"w").close()

            efit_bin = os.path.join(self.BIN_PATH, self.BIN)

            ishot = int(services.get_config_param('SHOT_NUMBER'))
            itime = int(services.get_config_param('TIME_ID'))

            efit_init(ishot, itime, f_instate, efit_bin, R0_scale=0, B0_scale=0)

            if cur_eqdsk_file != "g%06d.%05d"%(ishot,itime):
                shutil.copyfile("g%06d.%05d"%(self.ishot,self.itime), cur_eqdsk_file)

            ps.load_geqdsk(cur_eqdsk_file,keep_cur_profile=False)

        #-- load instate to plasma state

        print 'instate to ps:',f_instate
        instate_io.instate_to_ps(f_instate,ps)

        #-- write plasma state file

        ps.store(cur_state_file)

        #-- boundary condition state file

        instate = Namelist(f_instate)["inbc"]
        if not instate: instate = Namelist(f_instate)["instate"]

        inbc = Namelist()

        inbc["inbc"]["r0"] = instate["r0"]
        inbc["inbc"]["b0"] = instate["b0"]
        inbc["inbc"]["ip"] = instate["ip"]
        inbc["inbc"]["nbdry"] = instate["nbdry"]
        inbc["inbc"]["rbdry"] = instate["rbdry"]
        inbc["inbc"]["zbdry"] = instate["zbdry"]

        inbc["inbc"]["nlim"] = instate["nlim"]
        inbc["inbc"]["rlim"] = instate["rlim"]
        inbc["inbc"]["zlim"] = instate["zlim"]

        inbc.write(cur_bc_file)

        #-- instate, fastran nc file

        if f_instate: shutil.copyfile(f_instate,cur_instate_file)

        if cur_fastran_file: open(cur_fastran_file,"w").close()

        #-- eped

        if f_eped:
            cesol_io.eped_to_plasmastate(cur_state_file, f_eped, cur_instate_file)

        #-- update plasma state

        services.update_plasma_state()

        #-- archive output files

        services.stage_output_files(timeStamp, self.OUTPUT_FILES)

    def finalize(self, timeStamp=0.0):

        print 'fastran_init.finalize() called'

    #efit_bin = os.path.join(self.BIN_PATH, self.BIN)
    #if cur_eqdsk_file != "g%06d.%05d"%(self.ishot,self.itime):
    #    shutil.copyfile("g%06d.%05d"%(self.ishot,self.itime), cur_eqdsk_file)

def efit_init(ishot, itime, f_instate, efit_bin, R0_scale=0, B0_scale=0):

    #--- scaled GS

    scaled_gs = 0
    if R0_scale > 0 and B0_scale > 0: scaled_gs = 1

    #--- write efit run script

    kfile = "k%06d.%05d"%(ishot,itime)
    if scaled_gs: kfile = "k%06d.%05d_s"%(ishot,itime)

    args = "2\n 1\n "+kfile
    command = 'echo \"%s\"'%args + ' | ' + efit_bin

    f=open("xefit","w")
    f.write(command)
    f.close()

    #--- initial force free equilibrium

    print 'INIT EFIT, force free'

    instate = Namelist(f_instate)["instate"]

    nrho = 101
    inefit = Namelist()
    inefit["inefit"]["ip"   ] = [instate["ip"][0]*1.0e6]
    inefit["inefit"]["r0"   ] = instate["r0"]
    inefit["inefit"]["b0"   ] = instate["b0"]
    inefit["inefit"]["nrho" ] = [nrho]
    inefit["inefit"]["rho"  ] = linspace(0,1.0,nrho)
    inefit["inefit"]["press"] = nrho*[0.0]
    inefit["inefit"]["jpar" ] = nrho*[0.0]
    inefit["inefit"]["nlim" ] = instate["nlim" ]
    inefit["inefit"]["rlim" ] = instate["rlim" ]
    inefit["inefit"]["zlim" ] = instate["zlim" ]
    inefit["inefit"]["nbdry"] = instate["nbdry"]
    inefit["inefit"]["rbdry"] = instate["rbdry"]
    inefit["inefit"]["zbdry"] = instate["zbdry"]
    inefit.write("inefit")

    efit_io.fixbdry_kfile_init(ishot,itime,f_inefit="inefit")

    if scaled_gs:
        Rs = instate["r0"][0]/R0_scale
        Bs = instate["b0"][0]/B0_scale
        efit_io.scale_kfile(ishot,itime,Rs=Rs,Bs=Bs)

    f=open("log.efit","w")
    subprocess.call(["sh", "xefit"],stdout=f)

    if scaled_gs:
        g = readg("g%06d.%05d"%(ishot,itime))
        scaleg(g,R0=Rs,B0=Bs)
        writeg(g,ishot,itime,129)
