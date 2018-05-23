#!/usr/bin/env Python

#-----------------------------------------------------------------------
# fastran statefile 
# last modified : Jan 2013
# JM
#------------------------------------------------------------------------

import sys,re,os,os.path,pickle
from numpy import *
from optparse import OptionParser
import Namelist
import zoutone,zefitdata,zefitutil,zinterp,zfdat

#-----------------------------------------------------------------------
# main routine

def write(
        shot,
        time,
        dtavg,
        outone,
        geq,
        nrho_instate=101,
        tokamak_id='D3D',
        iDT=False,
        iscale=False,
        ):

    #------------------------------------------------------------------ 
    # time window

    if dtavg > 0:
       timew = [time-dtavg,time+dtavg]
    else:
       timew = None

    #------------------------------------------------------------------ 
    # configuration from geq

    ip  = geq['cpasma']
    r0  = geq["rzero" ]
    b0  = geq["bcentr"]

    rmajor = geq["ps"]["Rmajor_mean"][-1]
    aminor = geq["ps"]["rMinor_mean"][-1]
    kappa  = geq["ps"]["elong"][-1]
    delta  = geq["ps"]["triang"][-1]
    pmhd   = geq["ps"]["P_eq"]
    qmhd   = geq["ps"]["q_eq"]

    bdry = zefitutil.getbdry(geq)

    #------------------------------------------------------------------ 
    # instaet namelist

    instate = Namelist.Namelist()

    #------------------------------------------------------------------ 
    # configuration

    instate["instate"]["tokamak_id"] = [tokamak_id]

    instate["instate"]["ip"    ] = [ip*1.0e-6]
    instate["instate"]["b0"    ] = [b0]
    instate["instate"]["r0"    ] = [r0]
    instate["instate"]["rmajor"] = [rmajor]
    instate["instate"]["aminor"] = [aminor]
    instate["instate"]["kappa" ] = [kappa]
    instate["instate"]["delta" ] = [delta] 

    #------------------------------------------------------------------ 
    # composition

    atable = {"d":2,"h":1,"t":3,"c":12,"he":4}
    ztable = {"d":1,"h":1,"t":1,"c":6 ,"he":2}

    if iDT:
        nion = 2
        nimp = 2
        ion  = ['d','t']
        imp  = ['c','he']
        zfrac_ion = [0.5,0.5]
        zfrac_imp = [1.0,0.0]
    else: 
        nion = 1
        nimp = 1
        ion  = ['d']
        imp  = ['c']
        zfrac_ion = [1.0]
        zfrac_imp = [1.0]

    print "nion = %d"%nion,"==>",ion,"zfrac_ion=",zfrac_ion
    print "nimp = %d"%nimp,"==>",imp,"zfrac_imp=",zfrac_imp

    instate["instate"]["n_ion"] = [nion]
    instate["instate"]["z_ion"] = [ztable[s] for s in ion]
    instate["instate"]["a_ion"] = [atable[s] for s in ion]
    instate["instate"]["f_ion"] = zfrac_ion

    instate["instate"]["n_imp"] = [nimp]
    instate["instate"]["z_imp"] = [ztable[s] for s in imp]
    instate["instate"]["a_imp"] = [atable[s] for s in imp]
    instate["instate"]["f_imp"] = zfrac_imp

    #------------------------------------------------------------------ 
    # profile
    
    rho = arange(len(outone['te']['y'][0]))
    nrho = len(rho)
    rho = rho/(nrho-1.0)

    _fastranData = {
        "ne"       :{"map": "ene"     , "scale":1.0e-13} ,
        "ns_d"     :{"map": "end"     , "scale":1.0e-13} ,
        "ns_h"     :{"map": "enh"     , "scale":1.0e-13} ,
        "ns_t"     :{"map": "ent"     , "scale":1.0e-13} ,
        "ns_c"     :{"map": "enc"     , "scale":1.0e-13} ,
        "ns_he"    :{"map": "enhe"    , "scale":1.0e-13} ,
        "zeff"     :{"map": "zeff"    , "scale":1.0    } ,
        "te"       :{"map": "te"      , "scale":1.0    } ,
        "ti"       :{"map": "ti"      , "scale":1.0    } ,
        "omega"    :{"map": "omega"   , "scale":1.0    } ,
        "j_tot"    :{"map": "curpar"  , "scale":1.0e-2 } ,
        "j_oh"     :{"map": "curohm"  , "scale":1.0e-2 } ,
        "j_nb"     :{"map": "curbeam" , "scale":1.0e-2 } ,
        "j_bs"     :{"map": "curboot" , "scale":1.0e-2 } ,
        "j_ec"     :{"map": "currf"   , "scale":1.0e-2 } ,
        "j_ic"     :{"map": None      , "scale":1.0e-2 } ,
        "pe_nb"    :{"map": "qbeame"  , "scale":1.0    } ,
        "pe_ec"    :{"map": "qrfe"    , "scale":1.0    } ,
        "pe_ic"    :{"map": None      , "scale":1.0    } ,
        "p_ohm"    :{"map": "qohm"    , "scale":1.0    } ,
        "p_ei"     :{"map": None      , "scale":1.0    } ,
        "p_rad"    :{"map": "qrad"    , "scale":1.0    } ,
        "pe_ionization" 
                   :{"map": "qione"   , "scale":1.0    } ,
        "pe_fus"   :{"map": "qtfuse"  , "scale":1.0    } ,
        "pi_nb"    :{"map": "qbeami"  , "scale":1.0    } ,
        "pi_ec"    :{"map": None      , "scale":1.0    } ,
        "pi_ic"    :{"map": "qrfi"    , "scale":1.0    } ,
        "pi_ionization" 
                   :{"map": "qioni"   , "scale":1.0    } ,
        "pi_cx"    :{"map": "qcx"     , "scale":1.0    } ,
        "pi_fus"   :{"map": "qtfusi"  , "scale":1.0    } ,
        "torque_nb":{"map": "storqueb", "scale":0.1    } ,
        "torque_in":{"map": None      , "scale":0.1    } ,
        "se_nb"    :{"map": None     , "scale":1.0e-13 } ,
        "se_ionization" 
                   :{"map": None     , "scale":1.0e-13 } ,
        "si_nb"    :{"map": None     , "scale":1.0e-13 } ,
        "si_ionization" 
                   :{"map": None     , "scale":1.0e-13 } ,

        "q"        :{"map": "q"      , "scale":1.0     } , 
        "psipol"   :{"map": None     , "scale":1.0     } , 
        "density_beam"
                   :{"map": "enbeam"  , "scale":1.0e-13} ,
        "wbeam"    :{"map": "wbeam"   , "scale":1.0    } ,
        "density_alpha" 
                   :{"map": "enalp"   , "scale":1.0e-13} ,
        "walpha"   :{"map": "walp"    , "scale":1.0    } ,
        "chie"     :{"map": None      , "scale":1.0    } ,
        "chii"     :{"map": None      , "scale":1.0    } ,
    }

    profile_list = [
        "ne","te","ti","zeff","omega",
        "j_tot","j_oh","j_bs","j_nb","j_ec","j_ic",
        "pe_nb","pe_ec","pe_ic","pe_fus","pe_ionization","p_rad",
        "pi_nb","pi_ec","pi_ic", "pi_fus","pi_ionization", "pi_cx",
        "p_rad","p_ohm","p_ei",
        "torque_nb","torque_in",
        "se_nb","se_ionization","si_nb","si_ionization",
        "q","psipol",
        "density_beam","wbeam",
        "density_alpha","walpha",
        "chie","chii",
        ]

    output = {}
    for p in profile_list:
        p12 = _fastranData[p]["map"]
        print 'process profile '+ p
        if outone.has_key(p12):
            if timew:
                y = timeavg(outone[p12],timew) ['y_avg']
                yerr =timeavg(outone[p12],timew) ['y_err']
            elif time < 0:
                y = timeslice(outone[p12],-1,iprn=None)
            elif time == 0:
                y = timeslice(outone[p12],0,iprn=None)
            else:
                y = timeslice(outone[p12],time,iprn=None)
        else:
            print 'no variable ',p,p12
            output[p] = zeros(nrho)
        output[p] = {'y':y*_fastranData[p]["scale"],'units':''}
    output['q']['y'] = abs(output['q']['y']) 

    instate["instate"]["nrho"] = [nrho]
    instate["instate"]["rho"] = rho

    for key in profile_list:
        instate["instate"][key] = output[key]["y"] 


    zero_list = [
        "j_oh","j_bs","j_nb","j_ec","j_ic",
        "pe_nb","pe_ec","pe_ic","pe_fus","pe_ionization","p_rad",
        "pi_nb","pi_ec","pi_ic", "pi_fus","pi_ionization", "pi_cx",
        "p_rad","p_ohm","p_ei",
        "torque_nb","torque_in",
        "se_nb","se_ionization","si_nb","si_ionization",
        "q","psipol",
        "density_beam","wbeam",
        "density_alpha","walpha",
        "chie","chii",
        ]
    
    for key in zero_list:
        instate["instate"][key] = nrho*[0.0]

    instate["instate"]["nbdry"] = [bdry["nbdry"]]
    instate["instate"]["rbdry"] = bdry["rbdry"]
    instate["instate"]["zbdry"] = bdry["zbdry"]
    instate["instate"]["nlim" ] = [bdry["nlim" ]]
    instate["instate"]["rlim" ] = bdry["rlim" ]
    instate["instate"]["zlim" ] = bdry["zlim" ]

    #------------------------------------------------------------------ 
    # write

    instate.write("instate")

#-----------------------------------------------------------------------
# utils

def timeavg(data,atime,scale=1.0):

    t = array(data["t"])
    i = (t>=atime[0]) & (t<atime[1])

    y = scale*array(data["y"])[i] 
    y_avg = average(y,axis=0)
    y_err = std(y,axis=0)

    return {"y":y,"y_avg":y_avg,"y_err":y_err}  

def timeslice(data,time,iprn=None):
    iprn = 'tmp'
    t = array(data["t"])
    if time < 0:
        y = array(data["y"])[-1]
        if iprn: print "***"+iprn+":time slice taken at",t[-1],"(last)"
    elif time == 0:
        y = array(data["y"])[0]
        if iprn: print "***"+iprn+":time slice taken at",t[0],"(last)"
    else:
#       i = t >= time
        i = t > time-2
        y = array(data["y"])[i][0]
        if iprn: print "***"+iprn+":time slice taken at",t[i][0]
        if abs(time-t[i][0]) > 2.0:
           print "****** Warning:",time,t[i][0]
    return y 


def scale(f,const):

    for i in range(len(f)):
        f[i] *= const

#-----------------------------------------------------------------------
# standalone

if __name__ == "__main__":

    f = open("outone.save")
    outone = pickle.load(f)
    f.close()

    gfile = os.path.join('.','geqdsk.fastran')
    geq = zefitdata.efitdata(gfile,"ps")

    #--------------------------------------------------------------
    # write instate 

    write(  shot=123456
          , time=-1
          , dtavg=-1
          , outone=outone
          , geq=geq
          , iDT=True
          , iscale=False)

