#-----------------------------------------------------------------------
# read outone
# rewritten for speedup, Sep 2009
# last modified 2010 Jan
#-----------------------------------------------------------------------

import sys
import re
import pickle
import numpy
import Namelist

__outoneBlock = []
__outoneBlockPattern = {
      "current": [ ['j','r','curden','curohm','curboot'],
                   ["curden","curohm","curboot","curbeam","curbe","curbet", "currf","curpar"],
                   [(2,2),(2,3),(2,4),(2,5),(2,6),(2,7),(2,8),(2,10)] ],
      "energy":  [ ['j','r','r/a','te','ti','we'],
                   ["we","wi","wbeam","walp","te","ti" ],
                   [(2,5),(2,6),(2,7),(2,8),(2,3),(2,4)] ],
      "psi":     [ ['j','rho','psi','rmaj'],
                   ["rho","psi"],
                   [(3,1), (3,2)] ],
      "q":       [ ['j','r','te','ti','tn1'],
                   ["q"],
                   [(2,9)] ],
      "density" :[ ['j','r','r/a','zeff','ene'],
                   ["zeff","ene","enbeam","enalp","end","enh","enc"],
                   [(2,3),(2,4),(2,5),(2,6),(2,7),(2,8),(2,9) ] ],
      "rotation":[ ['j','r','omega'],
                   ["omega"],
                   [(3,2)] ],
      #"qe1":     [ ['j','r','r/a','qbeame','qrfe'],
      #             ["qbeame","qrfe","qtfuse"],
      #             [(2,3), (2,4), (2,6)] ],
      "qe1":     [ ['j','r','r/a','qbeame','qrfe'],
                   ["qbeame","qrfe"],
                   [(2,3), (2,4)] ],
      "qe2":     [ ['j','r','r/a','1.*\s+ qconde'],
                   ["qconde","qconve","qdelt","qohm","qione","qrad"],
                   [(2,4),(2,5),(2,6),(2,8),(2,9),(2,10)]],
      #"qi1":     [ ['j','r','r/a','qbeami'],
      #             ["qbeami","qrfi","qtfusi"],
      #             [(2,3),(2,4),(2,6)]],
      "qi1":     [ ['j','r','r/a','qbeami'],
                   ["qbeami","qrfi"],
                   [(2,3),(2,4)]],
      "qi2":     [ ['j','r','r/a','1.*\s+ qcondi'],
                   ["qcondi","qconvi","qioni","qcx"],
                   [(2,4),(2,5),(2,8), (2,9)]],
      "torque":  [ ['j','r.*','spbolt'],
                   ["storqueb","storque"],
                   [(1,3),(1,4)]],
#     "particle":[ ['particle .* name: d'],
#                   ["sion","srecom","scx","sbeam","sbcx","sfusion","ssum"],
#                   [(5,3),(5,4),(5,5),(5,6),(5,7),(5,8),(5,10) ]],
      "particle":[ ['particle sources .* for electrons'],
                   ["sion","srecom","scx","sbeam","sbcx","sfusion","ssum"],
                   [(4,3),(4,4),(4,5),(4,6),(4,7),(4,8),(4,9) ]],
      "sigma":   [ ['j','r','ftrap'],
                   ["ftrap","eta"],
                   [(2,2),(2,3)]],
      "neutral":  [['j','r','tn1','ennw1'],
                   ["tn1","ennw1","ennv1","volsn1"],
                   [(2,2),(2,3),(2,4),(2,5)]],
      "flux":     [['j','r','cond', 'conv'],
                   ["fconde","fconve","ffluxe","fcondi","fconvi","ffluxi"],
                   [(2,2),(2,3),(2,4),(2,5),(2,6),(2,9)]],
      "chi":      [['selected transport coefficients'],
                   ["chiineo","chie","chii"],
                   [(5,2),(5,3),(5,4)] ]
}

def mkpattern(str):
    p = "\s*"
    if len(str)==1:
        p=p+' '+str[0]
    else:
        for s in str: p=p+' '+s+' \s+'
    return p

def set_outoneBlock(nj, namep, namei, namen):
    global __outoneBlock

    __outoneBlock = []

    j = nj+3 
    __outoneBlockPattern["ip"] = \
      [ ['j','r','curden','curohm','curboot'],
        ["Ip","Iohm","Ibs","Inb","Irf"],
        [(j,1),(j,2),(j,3),(j,4),(j,7)] ] 

    __outoneBlockPattern["density"] = \
       [ ['j','r','r/a','zeff','ene'],
         ["zeff","ene","enbeam","enalp"],
         [(2,3),(2,4),(2,5),(2,6)],
       ]
    vnames =  __outoneBlockPattern["density"][1]
    offset =  __outoneBlockPattern["density"][2] 
    k = offset[-1][-1]+1
    for v in namep: 
        vnames.append("en"+v)
        offset.append((2,k))
        k=k+1
    for v in namei: 
        vnames.append("en"+v)
        offset.append((2,k))
        k=k+1
    for v in namen: 
        vnames.append("enn"+v)
        offset.append((2,k))
        k=k+1

    for b in list(__outoneBlockPattern.keys()):
        tmp = {}
        tmp["groupname"] = b
        tmp["pattern"] = re.compile(mkpattern(__outoneBlockPattern[b][0]))
        tmp["vnames"] = __outoneBlockPattern[b][1]
        tmp["offset"] = __outoneBlockPattern[b][2]
        tmp["locator"] = []
        tmp["time"] = []
        __outoneBlock.append(tmp)

def readoutone(foutone='outone', isummary=False, iverb=False):
#---read outone
    print("####################################")
    print("# READ OUTONE")

    filename = foutone

    f = open(filename,"r")
    lines = f.readlines()
    f.close()
    n = len(lines)

    onetwo_ver = lines[0].split(":")[-1].strip()
    onetwo_nw = int(lines[2].split()[5])
    onetwo_nh = int(lines[2].split()[8])
    onetwo_kj = int(lines[3].split()[5])
    print(("onetwo version :",onetwo_ver))
    print(("      nw,nh,kj : %d %d %d"%(onetwo_nw,onetwo_nh,onetwo_kj)))

#---read outone/namelis1
    print("####################################")
    print("# READ INONE NAMELIS1 in OUTONE")

    inone = Namelist.Namelist()
    inone.read(filename,only=["namelis1"])
    #if "NJ" in inone["namelis1"]:
    if "NJ" in inone["namelis1"].keys():
        nj = inone["namelis1"]["nj"][0]
        print("nj taken form inone/namelist1",nj)
    else:
        nj = onetwo_kj
        print("nj taken form kj",nj)
    namep = inone['namelis1']['namep']
    namei = inone['namelis1']['namei']
    namen = inone['namelis1']['namen']
    nprim = inone['namelis1']['nprim']
    nimp  = inone['namelis1']['nimp' ]
    nnew  = inone['namelis1']['nneu' ]
    print("namep :",namep)
    print("namei :",namei)
    print("namen :",namen)

#---init outoneBlock
    set_outoneBlock(nj,namep,namei,namen)

#---init return data 
    outone = {}
    for p in __outoneBlock:
      p["locator"] = []
      p["time"] = []

    pat0 = re.compile("1time =")

#---parsing outone
    print("####################################")
    print("# PARSING OUTONE: LOCATE")

    for k in range(n):
        if pat0.search(lines[k]):
           time = (lines[k].split())[2]

        for p in __outoneBlock:
            if p["pattern"].search(lines[k]):
                p["locator"].append(k)
                p["time"].append(time)
                if iverb: print(p["groupname"],time)

    print("####################################")
    print("# PARSING OUTONE: PUT")

    for p in __outoneBlock:
        ntime = len(p["time"])
        nvname = len(p["vnames"])

        if iverb: print(p["vnames"])

        for ivname in range(nvname):

            vname = p["vnames"][ivname]
            outone[vname] = {}
            outone[vname]["t"] = []
            outone[vname]["y"] = []

        for itime in range(ntime):

            loc = p["locator"][itime]
            time = p["time"][itime]

            for ivname in range(nvname):

                 vname   = p["vnames"][ivname]
                 offset  = p["offset"][ivname]
                 print(vname)
                 y = []
                 if p["groupname"] in ['sigma','flux','chi']:
                     y.append(0.0)
                     for j in range(loc+offset[0],loc+offset[0]+nj-1):
                         y.append( float(lines[j].split()[offset[1]]) )
                 elif p["groupname"] in ['ip']:
                     for j in range(loc+offset[0],loc+offset[0]+1):
                         y.append( float(lines[j].split()[offset[1]]) )
                 else:
                     for j in range(loc+offset[0],loc+offset[0]+nj):
                         y.append( float(lines[j].split()[offset[1]]) )

                 outone[vname]["t"].append(float(time))
                 outone[vname]["y"].append(numpy.array(y))

    #print len(outone["te"]["t"]) 
    #print len(outone["te"]["y"] )

#---global parameter
    #if isummary:
    #    summary = zsummary.getsummary(foutone)
    #    for var in summary.keys():
    #        outone[var] = summary[var]

    outone["shot"] = 0
    outone["namep"] = namep
    outone["namei"] = namei
    outone["namen"] = namen
    outone["inone"] = inone['namelis1']

    print("####################################")
    print("RETURN OUTONE")

    return outone

#-----------------------------------------------------------------------
# utils

def timeavg(data, atime, scale=1.0):
    t = numpy.array(data["t"])
    i = (t>=atime[0]) & (t<atime[1])

    y = scale*numpy.array(data["y"])[i] 
    y_avg = numpy.average(y,axis=0)
    y_err = numpy.std(y,axis=0)

    return {"y":y,"y_avg":y_avg,"y_err":y_err}  

def timeslice(data, time, iprn='OUTONE'):
    t = numpy.array(data["t"])
    if time < 0:
        y = numpy.array(data["y"])[-1]
        if iprn: print(iprn+":time slice taken at",t[-1],"(last)")
    elif time == 0:
        y = numpy.array(data["y"])[0]
        if iprn: print(iprn+":time slice taken at",t[0],"(last)")
    else:
        i = t >= time 
        y = numpy.array(data["y"])[i][0]
        if iprn: print(iprn+":time slice taken at",t[i][0],"(%d)"%time)
    return y 

#-----------------------------------------------------------------------
# standalone

if __name__=="__main__":
    from optparse import OptionParser

    parser = OptionParser()
    parser.add_option("--shot",
        action="store",type="int",dest="shot",default=-1)
    parser.add_option("--onetwodir",
        action="store",type="string",dest="onetwodir",default=".")
    parser.add_option("--prefix",
        action="store",type="string",dest="prefix",default="")
    parser.add_option("--isummary",
        action="store_true",dest="isummary",default=False)
    parser.add_option("--iverb",
        action="store_true",dest="iverb",default=False)

    (options,args) = parser.parse_args(sys.argv[1:])

    foutone = options.onetwodir+'/'+'outone'

    outone = readoutone(
              foutone=foutone
             ,isummary=options.isummary
             ,iverb=options.iverb)

    f = open('outone.save'+options.prefix,'w')
    pickle.dump(outone,f)
    f.close() 

