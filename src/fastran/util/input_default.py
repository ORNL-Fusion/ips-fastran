from Namelist import Namelist

def input_default(self, key='', default='', alias=[], maps={}, type='', verb=True):
    val = ''
    for _key in [key] + alias:
       if hasattr(self, _key):
           val = getattr(self, _key)
           if verb:
               print('INPUT:', _key, 'found', val, ', mapped value = ', maps.get(val, val))
           val =  maps.get(val, val)
           break
    if val == '':
       if verb:
           print('INPUT:', key, 'not found, use default value:', default)
       val = default
    if type == 'int':
        val = int(val)
    elif type == 'float':
        val = float(val)
    elif type == 'bool':
        val = bool(val)
    return val

if __name__=='__main__':
    class test(): pass

    check = test()
    check.b = '10'
    print( input_default(check, key='a', default='enabled', alias=['b'], maps={'10':'enabled'}, verb=True ) )
