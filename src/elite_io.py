import os
from glob import glob
from numpy import *

class elite_data():

    def __init__(self, slef.nmodes = [5,10,15,20,25] ):

        keylist_eq = [
            "alpha_max",
            "pprime_max",
            "jedge_n1_max",
            "jedge_n1_sep",
            "loc_alpha_max",
            "loc_pprime_max",
            "loc_jedgen1_max",
        ]

        keylist = [
           "gamma0",
           "gamma1",
           "gamsq",
           "del",
           "gamr",
           "gamr2",
           "tmatasym",
        ]

        self.data = {}

        for key in keylist_eq: self.data[key] = 0.0
        for key in keylist: self.data[key] = zeros(len(self.nmodes))

        self.runid = 'serdel2an'

    def read(self,rdir):

        for k in range(len(self.nmodes)):
            fname = os.path.join(rdir,self.runid+"%d.gamma"%self.nmodes[k])
            lines = open(fname,"r").readlines()
            tmp = [ float(v) for v in lines[1].split() ]
            self.data["gamma0"][k] = tmp[1]
            self.data["gamma1"][k] = tmp[2]
            self.data["gamsq"][k] = tmp[3]
            self.data["del"][k] = tmp[4]
            self.data["gamr"][k] = tmp[5]
            self.data["gamr2"][k] = tmp[6]
            self.data["tmatasym"][k] = tmp[7]

        fname = os.path.join(rdir,self.runid+"%d.jedge"%self.nmodes[-1])
        lines = open(fname,"r").readlines()
        tmp = [ float(v) for v in lines[7].split() ]
        self.data["alpha_max"] = tmp[0]
        self.data["pprime_max"] = tmp[1]
        self.data["jedge_n1_max"] = tmp[2]
        self.data["jedge_n1_sep"] = tmp[3]
        tmp = [ float(v) for v in lines[9].split() ]
        self.data["loc_alpha_max"] = tmp[0]
        self.data["loc_pprime_max"] = tmp[1]
        self.data["loc_jedgen1_max"] = tmp[2]

def collect():

    pass

if __name__ == "__main__":

    rdir = 'n00020'
    
    elite = elite_data()
    elite.read(rdir)
    print elite.data
