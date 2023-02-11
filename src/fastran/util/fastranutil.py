from Namelist import Namelist

def namelist_default(nml, key1, key2, val):
    if nml[key1][key2]:
        return nml[key1][key2]
    else:
        return val

def config_default(sim_conf, key1, key2, val):
     if key1 in sim_conf:
         if key2 in sim_conf[key1]:
             return_val = sim_conf[key1][key2]
         else:
             return_val = val
     else:
         return_val = val
     return return_val

if __name__=="__main__":
    test = Namelist()
    test["instate"]["ip"] = [1.0]
    print (namelist_default(test,"instate","ip",[2.0]))
    print (namelist_default(test,"instate","xx",[0.0]))
    print (namelist_default(test,"instate0","xx",[0.0]))
