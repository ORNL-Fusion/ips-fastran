#!/usr/bin/env python

"""
 -----------------------------------------------------------------------
  misc. utils
  JM
 -----------------------------------------------------------------------
"""

import os,sys,pickle

def fname_d3data(rdir,shot,time,prefix,ext=''):
    return rdir+'/'+'%s%06d.%05d'%(prefix,shot,time)+ext 

def get_shot_time(fname):
    import os
    shot = int(os.path.basename(fname).split('.')[0][-6:])
    time = int(os.path.basename(fname).split('.')[-1].split('_')[0])
    return shot,time

def check_file(fnames):
    for fname in fnames:
        if not os.path.exists(fname):
            print 'CANNOT FIND FILE:',fname
            sys.exit()

def load_savefile(fnames):
    data = [] 
    for fname in fnames:
        f = open(fname,'r')
        d = pickle.load(f)
        f.close()
        data.append(d)
    return data

def outscreen(s):
    print 72*'#'
    print s+":"
    print 72*'#'

def screen_out(s):
    print 72*'#'
    print s+":"
    print 72*'#'

def error(s):
    print 30*'-'
    print s+"...exit"
    sys.exit()

def linear(x,xgrid,y):

    ngrid = len(xgrid)
    if x <= xgrid[0]:
       return(y[0])
    elif x >= xgrid[-1]:
       return(y[-1])
    else:
       klo=0;
       khi=ngrid-1;
       while (khi-klo > 1):
          k=(khi+klo) >> 1
          if xgrid[k] > x :
              khi=k
          else:
              klo=k
    h = xgrid[khi]-xgrid[klo]
    a=(xgrid[khi]-x)/h;
    b=(x-xgrid[klo])/h;
 
    return a*y[klo]+b*y[khi]
