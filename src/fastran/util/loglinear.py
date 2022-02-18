import numpy as np
from Namelist import Namelist

class loglinear():
    def __init__(self,fname,ykey):
        input = Namelist(fname)
        self.data = {}
        self.keys = []
        for key in input[ykey].keys():
             self.data[key.lower()] = input[ykey][key][0]
             if key.lower() != 'const':
                 self.keys.append(key.lower())

    def get(self,**kwargs):
        rval = self.data['const']
        for key in kwargs:
            rval *= kwargs[key.lower()]**self.data[key.lower()]
        return rval

    def get_(self,ps):
        rval = self.data['const']
        for key in self.keys:
            rval *= ps[key.lower()][0]**self.data[key.lower()]
        return rval

    def latex(self):
        rval = self.data['const']
        str = "%5.2f"%rval
        for key in self.keys:
            str += key.lower()+"^%.2f"%self.data[key.lower()]
        return str

class linear():
    def __init__(self,fname,ykey):
        input = Namelist(fname)
        self.data = {}
        self.keys = []
        for key in input[ykey].keys():
             self.data[key.lower()] = input[ykey][key][0]
             if key.lower() != 'const':
                 self.keys.append(key.lower())

    def get(self,**kwargs):
        rval = self.data['const']
        for key in kwargs:
            rval += kwargs[key.lower()]*self.data[key.lower()]
        return rval

    def get_(self,ps):
        rval = self.data['const']
        for key in self.keys:
            rval += ps[key.lower()]*self.data[key.lower()]
        return rval

    def latex(self):
        return ""
