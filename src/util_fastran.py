from Namelist import Namelist

def namelist_default(nml, key1, key2, val):

    if nml[key1][key2]:
        return nml[key1][key2]
    else:
        return val


if __name__=="__main__":

    test = Namelist()
    test["instate"]["ip"] = [1.0]
    print namelist_default(test,"instate","ip",[2.0])
    print namelist_default(test,"instate","xx",[0.0])
    print namelist_default(test,"instate0","xx",[0.0])
