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

