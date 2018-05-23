from numpy import *
from glob import glob
import re

def read_nubeam_log(fname):

    pat_bt = re.compile("\s*beam-target neutrons")
    pat_bb = re.compile("\s*beam-beam neutrons")

    for line in open(fname,"r").readlines():
        if pat_bt.search(line):
            # print line.split()[-1].replace('D','E')
            neutrons_bt = float( line.split()[-1].replace('D','E') )
        if pat_bb.search(line):
            # print line.split()[-1].replace('D','E')
            neutrons_bb = float( line.split()[-1].replace('D','E') )

    print "%8.3e %8.3e"%(neutrons_bt, neutrons_bb)
    return neutrons_bt, neutrons_bb

def read_nubeam_log_avg(fnames):

    out_bt = []
    out_bb = [] 
    for fname in fnames:
        bt, bb = read_nubeam_log(fname)
        out_bt.append(bt)
        out_bb.append(bb)
    return average(array(out_bt)), average(array(out_bb))

files = glob("/Users/jmpark/samples/nubeamlog/log.nubeam_*")
sorted(files)
bt, bb = read_nubeam_log_avg(files)
print bt,bb

#bt, bb = read_nubeam_log("log.nubeam")

