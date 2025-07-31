from numpy import *
from Namelist import Namelist
from fastran.util.zinterp import zinterp

class pdata():
    def __init__(self):
        self.keylist = [
           "ne",
           "te",
           "ni",
           "ti",
           "omeg",
        ]
        self.data = {}
        self.units = {}

    def read(self,fname):
        lines = open(fname,"r").readlines()
        k = 0
        while k < len(lines):
            desc = lines[k].split()
            if desc[1] == "psinorm":
                npsi = int(desc[0])
                key = desc[2].split("(")[0]
                self.data[key] = {"x":zeros(npsi), "y":zeros(npsi), "y1":zeros(npsi)}
                print(key, npsi)
                for k2 in range(npsi):
                    k = k+1
                    tmp = lines[k].split()
                    self.data[key]['x' ][k2] = float(tmp[0])
                    self.data[key]['y' ][k2] = float(tmp[1])
                    self.data[key]['y1'][k2] = float(tmp[2])
            k = k+1

    def write(self,fname):
        f = open(fname,"w")
        for key in self.keylist:
            npsi = len(self.data[key]["x"])
            print(key, npsi)
            f.write("%d psinorm %s(%s) d%s/dpsi\n" %(npsi, key, self.units[key], key))
            for k in range(npsi):
                f.write("%f %f %f\n"%(self.data[key]["x"][k], self.data[key]["y"][k], self.data[key]["y1"][k]))

    def load_instate(self,fname, npsi=256):
         instate = Namelist(fname)
         psi_in = array(instate["inmetric"]["psi"])
         psi_in = psi_in/psi_in[-1]
         rho_in = array(instate["inmetric"]["rho"])

         psi = linspace(0,1,npsi)
         rho = zinterp(psi_in,rho_in)(psi)

         ne = zinterp(rho_in, array(instate["instate"]["ne"])*0.1 )
         ni = zinterp(rho_in, array(instate["instate"]["ni"])*0.1 )
         te = zinterp(rho_in, array(instate["instate"]["te"]) )
         ti = zinterp(rho_in, array(instate["instate"]["ti"]) )
         omeg = zinterp(rho_in, array(instate["instate"]["omega"])*1.0e-3 )

         self.data["ne"] = { "x":psi, "y":ne(rho), "y1":ne(rho,der=1) }
         self.data["ni"] = { "x":psi, "y":ni(rho), "y1":ni(rho,der=1) }
         self.data["te"] = { "x":psi, "y":te(rho), "y1":te(rho,der=1) }
         self.data["ti"] = { "x":psi, "y":ti(rho), "y1":ti(rho,der=1) }
         self.data["omeg"] = { "x":psi, "y":omeg(rho), "y1":omeg(rho,der=1) }

         self.units["ne"] = "10^20/m^3"
         self.units["ni"] = "10^20/m^3"
         self.units["te"] = "KeV"
         self.units["ti"] = "KeV"
         self.units["omeg"] = "kRad/s"

def collect():
    pass

if __name__ == "__main__":
    p = pdata()
    p.load_instate('i220010.00020')
    print(p.data)
    p.write('p220010.0020')
