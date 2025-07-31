from numpy import *
from glob import glob
import re
from fastran.plasmastate.plasmastate import plasmastate
from Namelist import Namelist


def read_nubeam_log(fname):
    pat_bt = re.compile("\s*beam-target neutrons")
    pat_bb = re.compile("\s*beam-beam neutrons")

    for line in open(fname,"r").readlines():
        if pat_bt.search(line):
            neutrons_bt = float( line.split()[-1].replace('D', 'E') )
        if pat_bb.search(line):
            neutrons_bb = float( line.split()[-1].replace('D', 'E') )

    print("%8.3e %8.3e"%(neutrons_bt, neutrons_bb))
    return neutrons_bt, neutrons_bb


def read_nubeam_log_avg(fnames):
    out_bt = []
    out_bb = []
    for fname in fnames:
        bt, bb = read_nubeam_log(fname)
        out_bt.append(bt)
        out_bb.append(bb)
    return average(array(out_bt)), average(array(out_bb))


def update_instate(f_state, f_instate):
    ps = plasmastate('ips', 1)
    ps.read(f_state)

    n0norm = ps['n0norm'][0, 0, :]
    sc0 = ps['sc0'][0]
    n0 = ps.cell2node_bdry(n0norm) * sc0 * 1.e-19

    instate = Namelist(f_instate)
    instate['instate']['n0'] = n0
    instate.write(f_instate)


if __name__ == "__main__":
    files = glob("log.nubeam_*")
    sorted(files)
    bt, bb = read_nubeam_log_avg(files)
    print(bt,bb)
