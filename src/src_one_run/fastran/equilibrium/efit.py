"""
 -----------------------------------------------------------------------
 efit component
 -----------------------------------------------------------------------
"""

import os
import re
import shutil
import glob
import subprocess
from numpy import *
import time as timer
from ipsframework import Component
from Namelist import Namelist
from fastran.equilibrium import efit_io
from fastran.plasmastate.plasmastate import plasmastate
from fastran.solver.inmetric_io import ps_to_inmetric
from fastran.equilibrium.efit_eqdsk import readg, writeg, scaleg
from fastran.util.input_default import input_default
from fastran.util.fastranutil import freeze

class efit(Component):
    def __init__(self, services, config):
        Component.__init__(self, services, config)
        print('Created %s' % (self.__class__))

    def _get_shot_time(self):
        ishot = int(self.services.get_config_param('SHOT_NUMBER'))
        itime = int(float(self.TIME_ID)) if hasattr(self, 'TIME_ID') else int(self.services.get_config_param('TIME_ID'))
        return ishot, itime

    def _write_xefit(self, ishot, itime):
        efit_bin = os.path.join(self.BIN_PATH, self.BIN)
        kfile = 'k%06d.%05d' % (ishot, itime)
        if self.user_inputs['scale'] == 'enabled':
            kfile = kfile + "_s"
        args = "2\n 1\n " + kfile
        command = 'echo \"%s\"' % args + ' | ' + efit_bin
        f = open('xefit', "w")
        f.write(command)
        f.close()

    def init(self, timeid=0):
        print('>>> efit.init() started')

        # -- user input parameters
        self.user_inputs = {}
        self.user_inputs['scale'] = input_default(
            self, key='SCALE', default='disabled', alias=['SCALE_GS'], maps={'0':'disabled', '1':'enabled'})
        self.user_inputs['r0_scale'] = float(input_default(
            self, key='R0_SCALE', default='1.7', alias=['r0_scale'], maps={})) # default to D3D EFIT
        self.user_inputs['b0_scale'] = float(input_default(
            self, key='B0_SCALE', default='2.0', alias=['b0_scale'], maps={})) # default to D3D EFIT

        # -- get shot and time
        ishot, itime = self._get_shot_time()
        if itime < 0:
            itime = timeid

        print('efit time = {}'.format(itime))

        # -- write efit run script
        self._write_xefit(ishot, itime)

        # -- initial force free equilibrium
        init_run = int(getattr(self, "INIT_RUN", 0))
        print('init_run = ', init_run)

        if init_run:
            print('force free')
            self.services.stage_state()
            self.initial_equilibrium(ishot, itime)
            self.services.update_state()

    def step(self, timeid=0):
        print('>>> efit.step() started')

        # -- freeze/resume
        if freeze(self, timeid, 'efit'): return None

        # -- set efit time
        ishot, itime = self._get_shot_time()
        if itime < 0:
            itime = timeid

        # -- stage plasma state files
        self.services.stage_state()

        # -- get plasma state file names
        ps_backend = getattr(self, 'PS_BACKEND', 'PS')

        cur_eqdsk_file = self.services.get_config_param('CURRENT_EQDSK')
        cur_state_file = self.services.get_config_param('CURRENT_STATE')
        cur_instate_file = self.services.get_config_param('CURRENT_INSTATE')

        use_instate = int(getattr(self, 'USE_INSTATE', '0'))
        instate = Namelist(cur_instate_file, 'r')

        # -- generate inefit
        f_inefit = 'inefit'
        f_inefit_prev = 'inefit0'

        model_pressure = input_default(
            self, key='PRESSURE', default='kinetic', alias=[], maps={})

        relax_profile = input_default(
            self, key='RELAX_PROFILE', default='0', alias=[], maps={}, type='float')

        relax_metric = input_default(
            self, key='RELAX_METRIC', default='0.5', alias=[], maps={}, type='float')

        betan_target = input_default(
            self, key='BETAN_TARGET', default='-1.0', alias=[], maps={}, type='float')

        if relax_profile > 0 and timeid >= 0:
            shutil.copyfile(f_inefit, f_inefit+'0')

        if ps_backend == 'PS' and use_instate == 0:
            efit_io.io_input_from_state(
                f_ps=cur_state_file, f_instate=cur_instate_file, f_inefit=f_inefit,
                model_pressure=model_pressure, betan_target=betan_target)
            nrho = 257
        elif ps_backend == 'INSTATE' or use_instate == 1:
            efit_io.io_input_from_instate(f_instate=cur_instate_file, f_inefit=f_inefit, model_pressure=model_pressure)
            nrho = instate['instate']['nrho'][0]

        if relax_profile > 0 and timeid > 0:
            print(f'relax_profile = {relax_profile}')
            inefit0 = Namelist(f_inefit_prev)
            inefit = Namelist(f_inefit)
            for key in ['press', 'jpar']:
                inefit['inefit'][key] = 0.5 * (array(inefit0['inefit'][key]) + array(inefit['inefit'][key]))
            inefit.write(f_inefit)

        # -- scaled GS
        if self.user_inputs['scale'] == 'enabled':
            Rs = instate['instate']['r0'][0] / self.user_inputs['r0_scale']
            Bs = instate['instate']['b0'][0] / self.user_inputs['b0_scale']
        else:
            Rs = 1.
            Bs = 1.

        # -- user inputs
        topology = input_default(self, key='TOPOLOGY', default='')
        error = float(input_default(self, key='ERROR', default='1.0e-4'))
        max_error_efit = float(input_default(self, key='MAX_EFIT_ERROR', default='0.1'))

        # -- write efit run script
        self._write_xefit(ishot, itime)

        # -- run efit, force-free
        init_run_step = int(getattr(self, 'INIT_RUN_STEP', '0'))
        print('init_run_step =', init_run_step)

        cwd = self.services.get_working_dir()

        niter = int(getattr(self, 'NITER', '5'))
        print('niter = ', niter)

        pat = re.compile('\s*it=')
        for k in range(niter):
            print('*** generate kfile %d'% k)
            t0 = timer.time()

            if k == 0 and init_run_step:
                print('force free')
                efit_io.fixbdry_kfile_init(ishot, itime, f_inefit=f_inefit)
                if self.user_inputs['scale'] == 'enabled':
                    efit_io.scale_kfile(ishot, itime, Rs=Rs, Bs=Bs)
                iconv = 0

            else:
                if k > 0:
                    relax = 1
                    print('apply underrelaxation')
                else:
                    relax = 0

                ps = plasmastate('ips', 1)

                geqdsk = readg(cur_eqdsk_file)

                r0 = geqdsk['rzero']
                b0 = abs(geqdsk['bcentr'])
                ip = geqdsk['cpasma']
                print('r0 =', r0)
                print('b0 =', b0)
                print('ip =', ip)

                try:
                    ps.init_from_geqdsk(cur_eqdsk_file, nrho=nrho, nth=101)
                except:
                    raise Exception('error in init_from_geqdsk')
                if k > 0:
                    inmetric_prev = inmetric
                inmetric = ps_to_inmetric(ps, r0, b0, ip)
                if k > 0 and relax_metric > 0:
                    print(f'inmtric under-relax: {relax_metric}')
                    inmetric_keys = ['psi', 'rhob', 'ipol', 'vol', 'gr2i']
                    for key in inmetric_keys:
                        inmetric['inmetric'][key] = 0.5*(array(inmetric_prev['inmetric'][key]) + array(inmetric['inmetric'][key]))

                iconv = efit_io.fixbdry_kfile(ishot, itime, Namelist(f_inefit, 'r'),
                            inmetric['inmetric'], relax=relax, topology=topology, error=error)
                if self.user_inputs['scale'] == 'enabled':
                    efit_io.scale_kfile(ishot, itime, Rs=Rs, Bs=Bs)

            if iconv:
                print('converged')
                break

            print('run efit')
            t1 = timer.time()

            if int(getattr(self, 'SERIAL', '0')) == 1:
                print('efit, subprocess')
                logfile = open('efit.log', 'w')
                retcode = subprocess.call(['sh xefit'], stdout=logfile, stderr=logfile, shell=True)
                logfile.close()
            else:
                task_id = self.services.launch_task(1, cwd, 'sh xefit', logfile='efit.log')
                retcode = self.services.wait_task(task_id)

            if (retcode != 0):
                raise Exception('Error executing efit')

            t2 = timer.time()
            print('elapsed time = %6.3f %6.3f' % (t1 - t0, t2 - t1))

            error_efit = []
            logfile = open('efit.log', 'r')
            lines = logfile.readlines()
            for line in lines:
                if pat.search(line):
                    try:
                        e = float(line.split('err=')[1].split()[0])
                    except:
                        e = 0.
                    error_efit.append(e)
            logfile.close()
            print('error_efit =', error_efit[-1])
            if abs(error_efit[-1]) > max_error_efit:
                raise Exception('EFIT not converged')

            if self.user_inputs['scale'] == 'enabled':
                shutil.copyfile('g%06d.%05d' % (ishot, itime), 'g%06d.%05d_s' % (ishot, itime))
                g = readg('g%06d.%05d' % (ishot, itime))
                scaleg(g, R0=Rs, B0=Bs)
                writeg(g, ishot, itime, 129)

            # -- update local geqdsk state
            if cur_eqdsk_file != 'g%06d.%05d' % (ishot, itime):
                shutil.copyfile('g%06d.%05d' % (ishot, itime), cur_eqdsk_file)

        # -- load geqdsk to plasma state file
        if ps_backend == 'PS':
            ps = plasmastate('ips', 1)
            ps.read(cur_state_file)
            ps.load_geqdsk(cur_eqdsk_file)
            ps.store(cur_state_file)

        elif ps_backend == 'INSTATE':
            instate['inmetric'] = inmetric['inmetric']
            instate.write(cur_instate_file)

        # -- update plasma state files
        self.services.update_state()

        # -- archive output files
        if self.OUTPUT_FILES:
            self.services.stage_output_files(timeid, self.OUTPUT_FILES, save_plasma_state=False)
        else:
            output_files = ' '.join(["{0}{1:06d}.{2:05d}".format(header, ishot, itime) for header in ['g', 'a', 'k', 'm']])
            output_files += " inefit"
            print('stage_output_files: ', output_files)
            self.services.stage_output_files(timeid,  output_files, save_plasma_state=False)

    def finalize(self, timeid):
        print('>>> efit.finalize() started')

        ishot, itime = self._get_shot_time()
        if itime < 0:
            itime = timeid

        dir_state = self.services.get_config_param('STATE_WORK_DIR')
        for header in ['k', 'a']:
            f = '%s%06d.%05d' % (header, ishot, itime)
            shutil.copyfile(f, os.path.join(dir_state, f))

    def initial_equilibrium(self, ishot, itime):
        # -- initial force free equilibrium
        ps_backend = getattr(self, 'PS_BACKEND', 'PS')
        cur_eqdsk_file = self.services.get_config_param('CURRENT_EQDSK')
        cur_state_file = self.services.get_config_param('CURRENT_STATE')
        cur_instate_file = self.services.get_config_param('CURRENT_INSTATE')

        instate = Namelist(cur_instate_file)['instate']
        nrho = 101
        inefit = Namelist()
        inefit['inefit']['ip'] = [instate['ip'][0]*1.e6]
        inefit['inefit']['r0'] = instate['r0']
        inefit['inefit']['b0'] = instate['b0']
        inefit['inefit']['nrho'] = [nrho]
        inefit['inefit']['rho'] = linspace(0, 1., nrho)
        inefit['inefit']['press'] = nrho*[0.]
        inefit['inefit']['jpar'] = nrho*[0.]
        inefit['inefit']['nlim'] = instate['nlim']
        inefit['inefit']['rlim'] = instate['rlim']
        inefit['inefit']['zlim'] = instate['zlim']
        inefit['inefit']['nbdry'] = instate['nbdry']
        inefit['inefit']['rbdry'] = instate['rbdry']
        inefit['inefit']['zbdry'] = instate['zbdry']
        inefit.write('inefit')

        efit_io.fixbdry_kfile_init(ishot, itime, f_inefit='inefit')

        if self.user_inputs['scale'] == 'enabled':
            print('scaled ', self.user_inputs['r0_scale'], self.user_inputs['b0_scale'])
            Rs = instate['r0'][0]/self.user_inputs['r0_scale']
            Bs = instate['b0'][0]/self.user_inputs['b0_scale']
            efit_io.scale_kfile(ishot, itime, Rs=Rs, Bs=Bs)

        f = open('log.efit', 'w')
        subprocess.call(['sh', 'xefit'], stdout=f)

        if self.user_inputs['scale'] == 'enabled':
            g = readg('g%06d.%05d' % (ishot, itime))
            scaleg(g, R0=Rs, B0=Bs)
            writeg(g, ishot, itime, 129)

        if cur_eqdsk_file != 'g%06d.%05d' % (ishot, itime):
            shutil.copyfile('g%06d.%05d' % (ishot, itime), cur_eqdsk_file)

        ps_update = getattr(self, 'PS_UPDATE', 'disabled')
        print('ps_update', ps_update)
        if ps_backend == 'PS' or ps_update == 'enabled':
            ps = plasmastate('ips', 1)
            ps.read(cur_state_file)
            ps.load_geqdsk(cur_eqdsk_file)
            ps.store(cur_state_file)

        elif ps_backend == 'INSTATE':
            ps = plasmastate('ips', 1)
            geqdsk = readg(cur_eqdsk_file)
            r0 = geqdsk['rzero']
            b0 = abs(geqdsk['bcentr'])
            ip = geqdsk['cpasma']
            ps.init_from_geqdsk(cur_eqdsk_file, nrho=nrho, nth=101)
            inmetric = ps_to_inmetric(ps, r0, b0, ip)

            instate = Namelist(cur_instate_file)
            instate['inmetric'] = inmetric['inmetric']
            instate.write(cur_instate_file)
