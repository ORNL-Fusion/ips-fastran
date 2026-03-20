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

def int_timeid(timeid):
    return int(str(timeid).split('_')[-1])

def freeze(component, timeid, component_name):
    itime = int_timeid(timeid)
    ifreeze = int(getattr(component, 'FREEZE', -1))
    iresume = int(getattr(component, 'RESUME', -1))
    if ifreeze >= 0 and itime >= ifreeze:
        if iresume < 0 or itime < iresume:
            print(f'{component_name} skipped, FREEZE = {ifreeze}, RESUME = {iresume}, TIMEID = {itime}')
            return True
    else:
        return False

if __name__=='__main__':
    test = Namelist()
    test['instate']['ip'] = [1.0]
    print (namelist_default(test, 'instate', 'ip', [2.0]))
    print (namelist_default(test, 'instate', 'xx', [0.0]))
    print (namelist_default(test, 'instate0', 'xx', [0.0]))
