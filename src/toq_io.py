#!/usr/bin/env python

import string
from numpy import *
import re

def __line2vec(line, itype=0):

    if itype==1:
        ipos = [0,6,14,22,30,38,46,54,61,68,75]
        vec = [float(line[ipos[k]:ipos[k+1]]) for k in range(len(ipos)-1)]
    elif itype==2:
        ipos = [0,8,15,23,31,39,47,55,62,69,77]
        vec = [float(line[ipos[k]:ipos[k+1]]) for k in range(len(ipos)-1)]
    else:
        str = string.split(string.strip(line))
        vec = [ float(v) for v in str ]

    return vec

def read_peddata_profile(filename='peddata'):

    print "read",filename
    lines = open(filename,"r").readlines()

    data = []
    for k in range(12,len(lines),2):
        try:
            data.append(__line2vec(lines[k],1)+__line2vec(lines[k+1],2))
        except:
            print 'skip data from %d to %d'%(k/2,len(lines)/2)
            break
    data = array(data).transpose()
    peddata = {
        "psin"  : data[0],
        "ptot"  : data[2]*1.0e3,
        "Te"    : data[3],
        "ne"    : data[4],
        "pp"    : data[5],
        "jtot"  : data[6],
        "beta"  : data[7],
        "betan" : data[8],
        "nus"   : data[9],
        "rho"   : data[10],
        "q"     : data[12],
        "s"     : data[13],
        "alpha" : data[14],
        "betap" : data[17],
        "Ti"    : data[19]
        }

    return peddata

def read_peddata(filename='peddata'):

    print "read",filename
    file = open(filename,"r")
    lines = file.readlines()
    file.close()

    data = {"ped":{},"top":{} }

    line_1 = __line2vec(lines[2])
    line_2 = __line2vec(lines[3])

    data["ped"]["psin" ] = line_1[0]
    data["ped"]["p"    ] = line_1[2]
    data["ped"]["te"   ] = line_1[3]
    data["ped"]["ne"   ] = line_1[4]
    data["ped"]["beta" ] = line_1[7]
    data["ped"]["betan"] = line_1[8]
    data["ped"]["nustar"]= line_1[9]
    data["ped"]["rho"  ] = line_2[0]
    data["ped"]["q"    ] = line_2[2]
    data["ped"]["betap"] = line_2[7]

    line_1 = __line2vec(lines[6])
    line_2 = __line2vec(lines[7])

    data["top"]["psin" ] = line_1[0]
    data["top"]["p"    ] = line_1[2]
    data["top"]["te"   ] = line_1[3]
    data["top"]["ne"   ] = line_1[4]
    data["top"]["beta" ] = line_1[7]
    data["top"]["betan"] = line_1[8]
    data["top"]["nustar"]= line_1[9]
    data["top"]["rho"  ] = line_2[0]
    data["top"]["q"    ] = line_2[2]
    data["top"]["betap"] = line_2[7]

    return data

def read_dskbal(filename='dskbal'):

    print "read",filename
    file = open(filename,"r")
    lines = []
    for k in range(100):
        lines.append(file.readline())
    # lines = file.readlines()
    file.close()

    data = {}

    for line in lines:
       if re.search(r'.*q\(95\)',line):
          data['q95'] = float(line.split()[2])
       if re.search(r'.*q\(axis\)',line):
          data['q0'] = float(line.split()[2])
       if re.search(r'.*q\(min\)',line):
          data['qmin'] = float(line.split()[2])
       if re.search(r'.*beta poloidal',line):
          data['betap'] = float(line.split('=')[1].split()[0])
       if re.search(r'li.*=.*li\(3\).*=',line):
          data['li'] = float(line.split()[2])
          data['li3'] = float(line.split()[-1])
       if re.search(r'.*end',line): break

    return data

def read_dskfixb(filename='dskfixb'):

    print "read",filename
    lines = open(filename,"r").readlines()

    R = []
    Z = []
    k = 0
    while k<len(lines):
       if re.search(r'npsi,ngeom',lines[k]):
          npsi, ngeom = [ int(v) for v in lines[k+1].split() ]
          k+=2
          print npsi, ngeom
       elif re.search(r'R of boundary points',lines[k]):
          nt = 0
          while len(R) < ngeom:
              k += 1
              R += [ float(v) for v in lines[k].split() ]
       elif re.search(r'Z of boundary points',lines[k]):
          nt = 0
          while len(Z) < ngeom:
              k += 1
              try:
                  Z += [ float(v) for v in lines[k].split() ]
              except:
                  Z += [ float(lines[k][i*19:i*19+19]) for i in range(4) ]
       else:
           k+=1

    return R, Z
