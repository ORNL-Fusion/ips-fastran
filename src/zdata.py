#!/usr/bin/env python

"""
 -----------------------------------------------------------------------
 io for fastran formatted data
 -----------------------------------------------------------------------
"""

import sys,os,re
from numpy import *
from collections import OrderedDict

class zdata():

    def __init__(self):

        self.fdata = OrderedDict()
        self.ftype = OrderedDict()

    def __call__(self,key):
        return self.fdata[key]

    def __getitem__(self,key):
        return self.fdata[key.upper()]

    def __setitem__(self,key,value):
        self.fdata[key.upper()] = array(value)
        self.ftype[key.upper()] = 'R'

    def read(self,fname,iverb=False):
    
        self.fdata = OrderedDict()
        self.ftype = OrderedDict()

        if iverb:
            print "####################################"
            print "READ ",file
    
        lines = open(fname,"r").readlines()
        nlines = len(lines)
    
        pat_R = re.compile("@")
    
        for k in range(nlines):
    
            if pat_R.search(lines[k]):
               p = lines[k].split()
               vname = p[1].upper()
               ndata = int(p[2])
               ftype = p[0][1]
               #print vname,ndata
               self.fdata[vname] = [] 
               self.ftype[vname] = ftype
               nread = 0
               for j in range(ndata):
                   if pat_R.search(lines[k+j+1]): 
                      print 'incompatible file format'
                   if ftype == 'R':
                       d = [ float(x) for x in lines[k+j+1].split() ]
                   elif ftype == 'I':
                       d = [ int(x) for x in lines[k+j+1].split() ]
                   else:
                       print 'type mismatch'
                       sys.exit()
                   nread = nread+len(d)
                   self.fdata[vname] = self.fdata[vname] + d
                   if nread == ndata: break
               if p[0][1]=='R':
                   self.fdata[vname] = array(self.fdata[vname])   

    def write(self,filename,vars=None):
    
        f = open(filename,"w")

        if vars: keys = vars
        else: keys = self.fdata.keys()

        for var in keys:
            if self.ftype[var] == "R":
                self.write_f(var,'',f)
            elif self.ftype[var] == "I":
                self.write_i(var,'',f)
            else:
                print 'type mismatch'
                sys.exit()

        f.close()

    def write_f(self,name,unit,file):
    
        fvec = self.fdata[name] 
        file.write("@R %s %d [%s]\n"%(name.upper(),len(fvec),unit))
        for i in range(len(fvec)):
            file.write(" %16.9e"%fvec[i])
            if (i+1)%4 == 0 :
               file.write("\n")
        file.write("\n")
    
    def write_i(self,name,unit,file):
    
        ivec = self.fdata[name]
        file.write("@I %s %d []\n"%(name.upper(),len(ivec)))
        for i in range(len(ivec)):
            file.write(" %16d"%ivec[i])
            if (i+1)%4 == 0 :
               file.write("\n")
        file.write("\n")
    
    
if __name__=="__main__":

   data = zdata()
   data.read("inmetric")

   print data["PMHD"]
   data["PMHD"] = [1.0]
   data.write("aigo")

   pass
