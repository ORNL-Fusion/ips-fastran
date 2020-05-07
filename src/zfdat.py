"""
 -----------------------------------------------------------------------
 io for fastran formatted data
 -----------------------------------------------------------------------
"""

import sys
import re
from numpy import *

def write_f(name, fvec, unit, file):
    file.write("@R %s %d [%s]\n"%(name.upper(), len(fvec), unit))
    for i in range(len(fvec)):
        file.write(" %16.9e"%fvec[i])
        if (i+1)%4 == 0 :
           file.write("\n")
    file.write("\n")

def write_i(name, ivec, unit, file):
    file.write("@I %s %d []\n"%(name.upper(), len(ivec)))
    for i in range(len(ivec)):
        file.write(" %16d"%ivec[i])
        if (i+1)%4 == 0 :
           file.write("\n")
    file.write("\n")

def write_formatted(fdata, filename, vars=None):
    f = open(filename, "w")
    if vars:
        for var in vars:
            if type(fdata[var]) == type(1):
                write_i(var,fdata[var],'',f)
            else:
                write_f(var,fdata[var],'',f)
    else:
        for var in fdata:
            if type(fdata[var][0]) == type(1):
                write_i(var,fdata[var],'',f)
            else:
                write_f(var,fdata[var],'',f)
    f.close()

def read_formatted(file, iverb=False):
    fdata = {}

#---read formatted data
    if iverb:
        print ("####################################")
        print ("READ ",file)

    f = open(file,"r")
    lines = f.readlines()
    f.close()
    nlines = len(lines)

    pat_R = re.compile("@")

#---parsing
    for k in range(nlines):
        if pat_R.search(lines[k]):
           p = lines[k].split()
           vname = p[1]
           ndata = int(p[2])
           #print (vname,ndata)
           fdata[vname] = [] # zeros(ndata)
           nread = 0
           for j in range(ndata):
               if pat_R.search(lines[k+j+1]):
                  print ('incompatible file format')
               if p[0][1] == 'R':
                   d = [ float(x) for x in lines[k+j+1].split() ]
               elif p[0][1] == 'I':
                   d = [ int(x) for x in lines[k+j+1].split() ]
               else:
                   raise Exception('type mismatch')
               nread = nread+len(d)
               fdata[vname] = fdata[vname] + d
               if nread == ndata: break
           if p[0][1]=='R':
              fdata[vname] = array(fdata[vname])

#    for vname in fdata:
#        fdata[vname] = array(fdata[vname])
    return fdata

def scale(fname, vname, factor):
    fdata = read_formatted(fname)
    fdata[vname] = factor*fdata[vname]
    write_formatted(fdata, fname+"_new")

def from_dict(indict, filename):
    f = open(filename,"w")
    for var in indict.keys():
        if type(indict[var][0]) == type(1):
            write_i(var, indict[var], '', f)
        else:
            write_f(var, indict[var], '', f)
    f.close()

def checkout(fname):
    pass

def compare(fname,fname1):
    pass

if __name__=="__main__":
    if len(sys.argv) < 2:
        raise Exception('usage: zdata.py mode')
    mode = sys.argv[1]
    if mode == 'checkout':
       checkout(sys.argv[2])
    elif mode == 'compare':
       compare(sys.argv[2], sys.argv[3])
