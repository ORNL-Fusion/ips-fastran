from numpy import *

class readblock():

    def __init__(self,fname,nhead):
        self.inputdata = loadtxt(fname,skiprows=nhead).transpose()
        self.keys = [ v for v in open(fname,"r").readlines()[nhead-1].split() ]
        try:
            self.ndata = len(self.inputdata[0])
        except:
            self.ndata = 1
            self.inputdata = array([self.inputdata]).transpose()
        print self.inputdata
    def __getitem__(self,key):
        return self.inputdata[self.keys.index(key)]        
