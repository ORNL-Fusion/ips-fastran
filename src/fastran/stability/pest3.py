"""
 -----------------------------------------------------------------------
 pest3 component
 -----------------------------------------------------------------------
"""

import os
import shutil
import re
from numpy import *
from component import Component
from plasmastate import plasmastate
from Namelist import Namelist

class pest3(Component):
    def __init__(self, services, config):
        Component.__init__(self, services, config)
        print ('Created %s' % (self.__class__))

    def init(self, timeid=0):
        print('pest3.init() called')

    def step(self, timeid=0):
        print('pest3.step() started')

        #-- entry
        services = self.services

        shot_number = services.get_config_param('SHOT_NUMBER')
        time_id = services.get_config_param('TIME_ID')

        #-- excutable
        corsica_bin = 'echo 0 | '+os.path.join(self.BIN_PATH_CALTRANS, self.BIN_CALTRANS)
        print('corsica_bin = ', corsica_bin)

        pest3_bin = os.path.join(self.BIN_PATH, self.BIN)
        print('pest3_bin = ', pest3_bin)

        #-- stage plasma state files
        services.stage_state()

        #-- get plasma state file names
        cur_instate_file = services.get_config_param('CURRENT_INSTATE')
        cur_eqdsk_file = services.get_config_param('CURRENT_EQDSK')

        #-- stage input files
        services.stage_input_files(self.INPUT_FILES)

        #-- prepare run
        shutil.copyfile(cur_eqdsk_file, 'eqdsk')
        f_inbas = getattr(self,'INBAS')

        #-- run teq
        print('run corsica')
        try:
            open("xcorsica", "w").write(corsica_bin+" "+f_inbas)
            cwd = services.get_working_dir()
            task_id = services.launch_task(1, cwd, 'sh xcorsica', logfile = 'xcorsica.log')
            retcode = services.wait_task(task_id)
        except Exception:
            raise Exception('...in launch_task, xcorsica')

        #-- get output
        teqdsk = "i%06d.%05d_test"%(int(shot_number), int(time_id))
        print('TEQDSK = ', teqdsk)
        shutil.copy(teqdsk, "ieqdsk")

        try:
            cmd = pest3_bin + ' -i8 -l16 -f\"ieqdsk\" -n1 -m\"xxxx\" -r\"xxxx\" -k\"1000 1100 1200\" '
            cmd+= '-s0 -D0.25 -W\"0 0.26 0.5 0.8 0.9 1.0\" -p\"2.0 2.0\" -M2 -b0.40 -K\"1\" -X '
            open("xpest3", "w").write(cmd)
            logfile = 'pest%06d.%05d'%(int(shot_number), int(time_id))
            cwd = services.get_working_dir()
            task_id = services.launch_task(1, cwd, 'sh xpest3', logfile = logfile)
            retcode = services.wait_task(task_id)
        except Exception:
            raise Exception('...in launch_task, pest3')

        ideal, qlist, deltap, rlist, llist = self.collect(logfile)
#       ideal, qlist, deltap = self.collect(logfile)

        instate = Namelist(cur_instate_file)
        instate["stab"]["ideal"] = [ideal]
        instate["stab"]["qlist"] = qlist
        instate["stab"]["deltap"] = deltap
        instate["stab"]["rlist"] = rlist
        instate["stab"]["llist"] = llist
        instate.write(cur_instate_file)

        #-- update plasma state files
        services.update_state()

        #-- archive output files
        services.stage_output_files(timeid, self.OUTPUT_FILES)

    def finalize(self, timeid=0):
        print('pest3.finalize() called')

    def collect(self, outfile):
        pat_w = re.compile("ideal instabilitie(s)")
        pat_q = re.compile("\s*safety factor")
        pat_r = re.compile("\s<rs d psi/d r>")
        #pat_l = re.compile("\s*Lambdas [psi_s norm]")
        pat_l = re.compile("Lambdas")

        lines = open(outfile,"r").readlines()
        ideal = 0
        for k, line in enumerate(lines):
            if pat_w.search(line):
                print(line)
                ideal = 1
        deltap = []
        qlist = []
        rlist = []
        llist = []
        for k, line in enumerate(lines):
            if pat_q.search(line):
                q = int(line.split("=")[-1])
                qlist.append(q)
                i_surface = int(lines[k-1].split()[0])
                print(q, i_surface)
                for line_r in lines[k+3:]:
                    if pat_r.search(line_r):
                        print(line_r)
                        rlist.append(float(line_r.split("=")[-1].split()[0]))
                        break
                for line_l in lines[k+4:]:
                    if pat_l.search(line_l):
                        print(line_l)
                        llist.append(float(line_l.split("=")[-1]))
                        break
                for line_D in lines[k+5:]:
                    pat = re.compile("%d %d"%(i_surface, i_surface))
                    if pat.search(line_D):
                        print(line_D)
                        deltap.append(float(line_D.split("=")[-1].split()[0]))
                        break

        print('Q =', qlist)
        print('DELTAP = ', deltap)

        return ideal, qlist, deltap, rlist, llist
#       return ideal, qlist, deltap
