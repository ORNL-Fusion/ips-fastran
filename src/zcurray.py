#! /usr/bin/env python

"""
 -----------------------------------------------------------------------
 utils for curray IO
 JM
 -----------------------------------------------------------------------
"""

import sys,os,shutil
import subprocess
from numpy import *

#-------------------
#--- zcode libraries
from zinterp import *
import zefitutil
from zplasmastate import plasma_state_file
from Namelist import Namelist

def write_f(fvec,num,ncol,file):

    for i in range(num):
        file.write(" %16.9e"%fvec[i])
        if (i+1)%ncol == 0 :
           file.write("\n")
    file.write("\n")

def line2vec(line):
    vec = []
    tmp = line.split()
    for k in range(len(tmp)):
        vec.append(float(tmp[k]))
    return array(vec) 

def readvec(lines,k0,nread):
    tmp = []
    for k in range(k0,k0+nread):
        vec = line2vec(lines[k])
        for i in range(len(vec)): tmp.append(vec[i])
    return array(tmp),k0+nread

def get_nread(n,ncol=5):
    nr= n/ncol
    if n%ncol>0: nr+=1
    return nr

def wrt_curray_input(f_state,f_incurray,f_geqdsk):

    # =================================================================
    # read plasma state file
    
    ps  = plasma_state_file(f_state)

    ps_xe  = 1.6022e-19
    ps_mp  = 1.6726e-27

    z_ion = [round(x) for x in ps["qatom_S"][1:]/ps_xe ]
    a_ion = [round(x) for x in ps["m_S"][1:]/ps_mp ]
    n_ion = len(z_ion)

    nrho  = ps.nrho
    rho   = ps["rho"][:]
    ne    = ps["ns"][0,:]
    te    = ps["Ts"][0,:]
    ti    = ps["Ti"][:]
    ne    = ps.cell2node(ne)
    te    = ps.cell2node(te)
    ti    = ps.cell2node(ti)

    ni = {}
    for k in range(n_ion):
        ns = ps["ns"][k+1,:]
        ni[k] = ps.cell2node(ns)

    n_imp =  len(ps.data.dimensions["dim_nspec_imp0"])

    print "n_ion = ",n_ion
    print "z_ion = ",z_ion
    print "a_ion = ",a_ion
    print "n_imp = ",n_imp

    # =================================================================
    # write trxplt.out
    # [m^-3, keV]

    f = open("trxpl.out","w")

    nspec = n_ion+1 #<=========

    f.write("%13d\n"%nrho)
    f.write("%13d\n"%nspec)
    f.write("\n")
    write_f(rho,nrho,5,f)

    # -------------
    # electron

    f.write("\n")
    f.write("%f\n"%0.0)
    f.write("%f\n"%0.0)
    f.write("\n")
    write_f(ne,nrho,5,f)
    f.write("\n")
    write_f(te,nrho,5,f)

    # -------------
    # main ions + impurities

    for k in range(n_ion):
        f.write("\n")
        f.write("%f\n"%z_ion[k])
        f.write("%f\n"%a_ion[k])
        f.write("\n")
        write_f(ni[k],nrho,5,f)
        f.write("\n")
        write_f(ti,nrho,5,f)

    # -------------
    # beam ions

    f.write("\n")
    write_f(zeros(nrho),nrho,5,f)

    f.close()

    # -------------
    # write curray.in

    incurray = Namelist(f_incurray,case="lower")
    incurray["incurray"]["eqdsk_name"] = [f_geqdsk]
    incurray["incurray"]["nprim"     ] = [n_ion-n_imp]
    incurray["incurray"]["nspec"     ] = [n_ion]
    incurray["incurray"]["nimp"      ] = [n_imp]
    incurray["incurray"]["nminor"    ] = [0]
    incurray["incurray"]["atm"       ] = a_ion
    incurray["incurray"]["azi"       ] = z_ion
    incurray["incurray"]["nprofs"    ] = [nrho]

    curray_in = Namelist()
    for key in incurray["incurray"].keys():
        curray_in["input"][key] = incurray["incurray"][key]
    curray_in.write("curray_in")


def read_curray_output():

    f=open("raytrout","r")
    lines = f.readlines()
    f.close()

    nrho  = int(lines[0].split()[0])
    nspec = int(lines[0].split()[1])

    print 'nrho = ',nrho
    print 'nspec = ',nspec

    k = 1

    nread = get_nread(nrho,ncol=5)

    rho,k = readvec(lines,k,nread)
    pe_ic,k = readvec(lines,k,nread) # [MW/m^3]
    pi_ic = []
    for ii in range(nspec):
        tmp,k = readvec(lines,k,nread)
        pi_ic.append(tmp)
    pi_ic = sum(array(pi_ic),axis=0)
    j_ic,k = readvec(lines,k,nread)
    j_ic = j_ic*0.01 # [MA/m^2]

    #print rho
    #print pe_ic
    #print pi_ic
    #print j_ic

    outcurray = Namelist()
    outcurray["outcurray"]["rho"  ] = rho
    outcurray["outcurray"]["nrho" ] = [nrho]
    outcurray["outcurray"]["pe_ic"] = pe_ic
    outcurray["outcurray"]["pi_ic"] = pi_ic
    outcurray["outcurray"]["j_ic" ] = j_ic
    outcurray.write("outcurray")

    return outcurray
    
###############################################################################
# io

def io_update_instate(f_instate,f_outcurray,f_incurray):

    # read incurray

    incurray = Namelist(f_incurray)

    try:
       j_multi = incurray["adjust"]["j_multi"][0]
    except:
       j_multi = 1.0
       print 'no j_multi inuput, set 1.0'

    # read outcurray

    outcurray = Namelist(f_outcurray)
    rho   = outcurray["outcurray"]["rho"  ]
    pe_ic = outcurray["outcurray"]["pe_ic"]
    pi_ic = outcurray["outcurray"]["pi_ic"]
    j_ic  = outcurray["outcurray"]["j_ic" ]

    j_ic = j_multi*array(j_ic)

    # update local infastran

    instate = Namelist(f_instate)
    rho_in = instate["instate"]["rho"]
    instate["instate"]["j_ic" ] = interp1d(rho,j_ic ,kind='cubic')(rho_in)
    instate["instate"]["pe_ic"] = interp1d(rho,pe_ic,kind='cubic')(rho_in)
    instate["instate"]["pi_ic"] = interp1d(rho,pi_ic,kind='cubic')(rho_in)
    instate["instate"]["j_ic" ] = interp1d(rho,j_ic ,kind='cubic')(rho_in)
    instate.write(f_instate)

def io_update_state(f_state,f_geq,f_outcurray,f_incurray):

    # read plasma state

    #ps  = plasma_state_file(f_state)

    geq = zefitutil.readg(f_geq) 
    r0  = geq["rzero" ]
    b0  = abs(geq["bcentr"])
    ip  = geq['cpasma']
    ps  = plasma_state_file(f_state,r0=r0,b0=b0,ip=ip)

    # read incurray

    incurray = Namelist(f_incurray)

    try:
       j_multi = incurray["adjust"]["j_multi"][0]
    except:
       j_multi = 1.0
       print 'no j_multi inuput, set 1.0'

    # read outcurray

    outcurray = Namelist(f_outcurray)

    rho = outcurray["outcurray"]["rho"  ]
    jp  = outcurray["outcurray"]["j_ic" ]
    pe  = outcurray["outcurray"]["pe_ic"]
    pi  = outcurray["outcurray"]["pi_ic"]

    jp = j_multi*array(jp)

    # update local infastran

    rho_ps = ps["rho"][:]
    jp_ps = 1.0e6*interp1d(rho,jp,kind='cubic')(rho_ps)
    pe_ps = 1.0e6*interp1d(rho,pe,kind='cubic')(rho_ps)
    pi_ps = 1.0e6*interp1d(rho,pi,kind='cubic')(rho_ps)

    ps.load_j_parallel_CD(rho_ps,jp_ps,"ic")
    ps.load_profile(rho_ps,pe_ps,"pmine","vol")
    ps.load_profile(rho_ps,pi_ps,"pmini","vol")

    ps.close()

###############################################################################
# check

if __name__ == "__main__":

    wrt_curray_input("instate")
    outcurray = read_curray_output()
    io_update_instate("instate","outcurray")


