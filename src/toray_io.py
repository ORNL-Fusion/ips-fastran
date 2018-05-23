#!/usr/bin/env python

"""
 ------------------------------------------------------------------------------
  utils for toray
 ------------------------------------------------------------------------------
"""

import sys,os,shutil
from numpy import *
import netCDF4
from scipy.interpolate import interp1d
from scipy import optimize

from Namelist import Namelist
from zinterp import zinterp

def write_toray_input(geq,prof,intoray):

    #geq["__psipol__" ] = ps["psipol"][:]
    #geq["__phit__"   ] = ps["phit"  ][:]
    #geq["__hfact__"  ] = ps["gncfh" ][:]
    #geq["__gb1__"    ] = ps["gb1"   ][:]
    #geq["__gb2__"    ] = ps["gb2"   ][:]
    #geq["__gri__"    ] = ps["gri"   ][:]

    #------------------------------------------------------------------
    # set misc

    lcentr = 2
    gasep  = 1.0e-4

    #------------------------------------------------------------------
    # psiin

    ledge  = prof["nrho"]

    f = open("psiin","w")

    f.write("%4d"%lcentr)
    f.write("%4d"%ledge)
    f.write("%4d"%geq["nw"])
    f.write("%4d"%geq["nh"])
    f.write("\n")

    f.write("%16.9e"%geq["xdim"])
    f.write("%16.9e"%geq["ydim"])
    f.write("%16.9e"%geq["rzero"])
    f.write("%16.9e"%geq["rgrid"])
    f.write("\n")

    f.write("%16.9e"%geq["rmaxis"])
    f.write("%16.9e"%geq["zmaxis"])
    f.write("%16.9e"%geq["ssimag"])
    f.write("%16.9e"%geq["ssibry"])
    f.write("%16.9e"%abs(geq["bcentr"]))
    f.write("\n")

    f.write("%16.9e"%gasep)
    f.write("\n")

    write_formatted(f,abs(geq["fpol"]),"%16.9e",5)
    write_formatted(f,geq["psirz"].reshape(geq["nw"]*geq["nh"]),"%16.9e",5)

    psin = geq["__psipol__"]/geq["__psipol__"][-1]
    psir = (geq["ssibry"]-geq["ssimag"])*psin+geq["ssimag"]
    psir = 1.0e5*psir  # cgs unit!

    tmp = psir
    nrho_m = len(tmp)
    rho_m = arange(nrho_m)/(nrho_m-1.0)
    tmp = interp1d(rho_m,tmp,kind='cubic')
    psir = tmp(prof["rho"])
    psir[0] = 0.0

    write_formatted(f,psir,"%16.9e",5)
    write_formatted(f,geq["qpsi"],"%16.9e",5)

    f.write("%5d"%geq["nbdry"])
    f.write("%5d"%geq["nlim"] )
    f.write("\n")

    tmp = []
    for k in range(geq["nbdry"]):
        tmp.append(geq["rbdry"][k])
        tmp.append(geq["zbdry"][k])
    tmp = array(tmp)
    write_formatted(f,tmp,"%16.9e",5)

    tmp = []
    for k in range(geq["nlim"]):
        tmp.append(geq["rlim"][k])
        tmp.append(geq["zlim"][k])
    tmp = array(tmp)
    write_formatted(f,tmp,"%16.9e",5)

    nrho = len(psin)
    rho = arange(nrho)/(nrho-1.0)

    rhopsi = interp1d(psin,rho,kind='cubic')

    psin_geq = arange(geq["nw"])/(geq["nw"]-1.0)
    rho_eval = rhopsi(psin_geq)
    rho_eval[ 0]=0.0
    rho_eval[-1]=1.0

    tmp = geq["__hfact__"]/abs(geq["bcentr"])
    nrho_m = len(tmp)
    rho_m = arange(nrho_m)/(nrho_m-1.0)
    tmp = interp1d(rho_m,tmp,kind='cubic')
    tmp = tmp(rho_eval)
    write_formatted(f,tmp,"%16.9e",5)

    tmp = geq["__gb2__"]/geq["bcentr"]**2
    nrho_m = len(tmp)
    rho_m = arange(nrho_m)/(nrho_m-1.0)
    tmp = interp1d(rho_m,tmp,kind='cubic')
    tmp = tmp(rho_eval)
    write_formatted(f,tmp,"%16.9e",5)

    tmp = geq["__gb1__"]/abs(geq["bcentr"])
    nrho_m = len(tmp)
    rho_m = arange(nrho_m)/(nrho_m-1.0)
    tmp = interp1d(rho_m,tmp,kind='cubic')
    tmp = tmp(rho_eval)
    write_formatted(f,tmp,"%16.9e",5)

    tmp = geq["__gri__"]*geq["rzero"]  #<=======
   #tmp = geq["__gri__"]*geq["__R_axis__"]
    nrho_m = len(tmp)
    rho_m = arange(nrho_m)/(nrho_m-1.0)
    tmp = interp1d(rho_m,tmp,kind='cubic')
    tmp = tmp(rho_eval)
    write_formatted(f,tmp,"%16.9e",5)

    f.close()

    #------------------------------------------------------------------
    # echin

    time    = 0.0
    nrho    = prof["nrho"]

    r0      = geq["rzero" ]*100.0  # [cm]
    bt0     = abs(geq["bcentr"])*1.0e4  # [Gauss]
    ra      = 1.0e2*(geq["__phit__"][-1]/pi/abs(geq["bcentr"]))**0.5 # [cm]

    f = open("echin","w")

    f.write("%16.9e" % time)
    f.write("\n")
    f.write("%4d"    % intoray["idamp"])
    f.write("%4d"    % prof   ["nrho" ])
    f.write("%4d"    % intoray["nray" ])
    f.write("%4d"    % intoray["nbfld"])
    f.write("\n")
    f.write("%16.9e" % intoray["freq" ])
    f.write("%16.9e" % intoray["wrfo" ])
    f.write("%16.9e" % intoray["x"    ])
    f.write("%16.9e" % intoray["z"    ])
    f.write("%16.9e" % intoray["thet" ])
    f.write("\n")
    f.write("%16.9e" % intoray["phai" ])
    f.write("%16.9e" % intoray["hlw"  ])
    f.write("%16.9e" % intoray["ratw" ])
    f.write("%16.9e" % r0              )
    f.write("%16.9e" % bt0             )
    f.write("\n")
    f.write("%16.9e" % ra              )
    f.write("\n")

    psin = interp1d(rho,psin,kind='cubic')(prof["rho"])
    psin[0]=0.0
    psin[-1]=0.0

    write_formatted(f,psin,"%16.9e",5)
   #write_formatted(f,psin_eval,"%16.9e",5)
    write_formatted(f,prof["zeff"],"%16.9e",5)
    write_formatted(f,prof["ne"]*1.0e13,"%16.9e",5) # cm^-3
    write_formatted(f,prof["te"],"%16.9e",5) # kev

    f.close()

def read_toray_output(f_ncfile="toray.nc"):

    #------------------------------------------------------------------
    # read

    nc = netCDF4.Dataset(f_ncfile,'r',format='NETCDF4')
    nrho = nc.variables["ledge"][0]
    rho  = nc.variables["xrho"  ][:]  #normalized sqrt(tor flux), edge
    jec  = nc.variables["currf" ][:]  #A/cm^2 per incident W, cell
    qec  = nc.variables["weecrh"][:]  #W/cm^3 per incident W, cell
    iec  = nc.variables["tidept"][:]  #A per incident W, cell
    pec  = nc.variables["tpowde"][:]  #W per incident W, cell
    nc.close()

    Iec = iec[-1]
    Pec = pec[-1]

    jec_out = zeros(nrho)
    qec_out = zeros(nrho)
    for i in range(1,nrho-1):
        jec_out[i] = 0.5*(jec[i]+jec[i-1])*1.0e4
        qec_out[i] = 0.5*(qec[i]+qec[i-1])*1.0e6
    jec_out[ 0] = 0.0
    jec_out[-1] = 0.0
    qec_out[ 0] = 0.0
    qec_out[-1] = 0.0

    print "Pec abs. fraction = %5.3f"%Pec
    print "Iec (kA/MW) = %5.1f"%(Iec*1.0e3)

    #------------------------------------------------------------------
    # return

    d = Namelist()
    d["outtoray"]["nrho"] = [nrho]
    d["outtoray"]["rho" ] = rho
    d["outtoray"]["jec" ] = jec_out # A/m^2/W
    d["outtoray"]["qec" ] = qec_out # W/m^3/W
    d["outtoray"]["Iec" ] = [Iec]   # [A/W]
    d["outtoray"]["Pec" ] = [Pec]   # []
    d.write("outtoray")

    return d

###############################################################################
# utils read/write

def write_formatted(f,vec,desc,ncol):

    n = len(vec)
    for k in range(n):
        f.write(desc%vec[k])
        if (k+1)%ncol == 0: f.write("\n")
    if n%ncol != 0: f.write("\n")

def line2vec(line,nlen=16):
    ncol = len(line)/nlen
    vec = zeros(ncol)
    for k in range(ncol):
        tmp = line[k*nlen:(k+1)*nlen]
        vec[k] = float(tmp)
    return vec

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

###############################################################################
# utils adjust

def cur_integral(rho,jpar,geq,r0,b0):

    vol  = zinterp(geq["rho"][:],geq["vol"][:])(rho)
    g_eq = zinterp(geq["rho"][:],geq["g_eq"][:])(rho)
    #r0 = geq.r0
    #b0 = abs(geq.b0)

    I = 0.0
    nrho = len(rho)
    for i in range(1,nrho):
        jpar_m = 0.5*(jpar[i]+jpar[i-1])
        dvol   = (vol[i]-vol[i-1])
        ipol_m = 0.5*(g_eq[i]+g_eq[i-1])
        ipol_m/= r0*b0
        I += jpar_m*dvol/ipol_m**2/(2.0*pi*r0)
    return I

def vol_integral(rho,f,geq):

    val = 0.0
    nrho = len(rho)
    vol  = zinterp(geq["rho"][:],geq["vol"][:])(rho)
    for i in range(nrho-1):
        f_m = 0.5*(f[i+1]+f[i])
        dvol = vol[i+1]-vol[i]
        val += f_m*dvol
    return val

def fit_gaussian(rho,f):

    fitfunc = lambda p, x: p[0]*exp(-(x-p[1])**2/(2*p[2]**2))
    errfunc = lambda p, x, y: fitfunc(p, x)-y

    p0 = [1.0,1.0,0.1]
    try:
        p, success = optimize.leastsq(errfunc, p0[:], args=(rho,f))
    except:
        p = [-1.0,-1.0,-1.0]

    return p

def adjust_jec(width,xmid,fmax,Iec,nrho,geq,r0,b0):

    gauss = lambda p, x: p[0]*exp(-(x-p[1])**2/(2*p[2]**2))

    def func(f):
        x = arange(nrho)/(nrho-1.0)
        y = gauss( [f,xmid,width], x )
        return cur_integral(x,y,geq,r0,abs(b0))-Iec

    tmp = optimize.bisect(func,0.0,10*fmax,xtol=1.0e-6)

    x = arange(nrho)/(nrho-1.0)
    y = gauss( [tmp,xmid,width], x )

    return y

def adjust_qec(width,xmid,fmax,Pec,nrho,geq):

    gauss = lambda p, x: p[0]*exp(-(x-p[1])**2/(2*p[2]**2))

    def func(f):
        x = arange(nrho)/(nrho-1.0)
        y = gauss( [f,xmid,width], x )
        return vol_integral(x,y,geq)-Pec

    tmp = optimize.bisect(func,0.0,10.0*fmax,xtol=1.0e-6)

    x = arange(nrho)/(nrho-1.0)
    y = gauss( [tmp,xmid,width], x )

    return y

def gauss_asym(p,x):

    rvec = zeros(len(x))
    for k in range(len(x)):
        if x[k]>=p[1]: rvec[k] = p[0]*exp(-(x[k]-p[1])**2/(2*p[2]**2))
        if x[k]< p[1]: rvec[k] = p[0]*exp(-(x[k]-p[1])**2/(2*p[3]**2))
    return rvec

def adjust_jec_asym(width,width2,xmid,fmax,Iec,nrho,geq,r0,b0):

    def func(f):
        x = arange(nrho)/(nrho-1.0)
        y = gauss_asym( [f,xmid,width,width2], x )
        return cur_integral(x,y,geq,r0,abs(b0))-Iec

    tmp = optimize.bisect(func,0.0,10*fmax,xtol=1.0e-6)

    x = arange(nrho)/(nrho-1.0)
    y = gauss_asym( [tmp,xmid,width,width2], x )

    return y

def adjust_qec_asym(width,width2,xmid,fmax,Pec,nrho,geq):

    def func(f):
        x = arange(nrho)/(nrho-1.0)
        y = gauss_asym( [f,xmid,width,width2], x )
        return vol_integral(x,y,geq)-Pec

    tmp = optimize.bisect(func,0.0,10.0*fmax,xtol=1.0e-6)

    x = arange(nrho)/(nrho-1.0)
    y = gauss_asym( [tmp,xmid,width,width2], x )

    return y

###############################################################################
# io plasma state

def input_from_state(geq,ps,intoray,k,nrho=101):

    rho_in   = ps["rho"][:]
    ne_in    = ps["ns"][0,:]*1.0e-19
    te_in    = ps["Ts"][0,:]
    zeff_in  = ps["Zeff"][:]
    ne_in    = ps.cell2node(ne_in)
    te_in    = ps.cell2node(te_in)
    zeff_in  = ps.cell2node(zeff_in)

    rho = arange(nrho)/(nrho-1.0)
    prof = {
        "nrho" : nrho,
        "rho"  : rho,
        "ne"   : interp1d(rho_in,ne_in,kind='cubic')(rho),
        "te"   : interp1d(rho_in,te_in,kind='cubic')(rho),
        "zeff" : interp1d(rho_in,zeff_in,kind='cubic')(rho)
    }

    intoray_k = {}
    for v in intoray["intoray"].keys():
        if v not in ["ntoray"]:
            intoray_k[v] = intoray["intoray"][v][k]

    geq["__psipol__" ] = ps["psipol"][:]
    geq["__phit__"   ] = ps["phit"  ][:]
    geq["__hfact__"  ] = ps["gncfh" ][:]
    geq["__gb1__"    ] = ps["gb1"   ][:]
    geq["__gb2__"    ] = ps["gb2"   ][:]
    geq["__gri__"    ] = ps["gri"   ][:]

    write_toray_input(geq,prof,intoray_k)

def update_state(geq,ps,intoray):

    ntoray = intoray["intoray"]["ntoray"][0]
    r0  = geq["rzero" ]
    b0  = abs(geq["bcentr"])

    for k in range(ntoray):

        rfpow = intoray["intoray"]["rfpow"][k]*1.0e-6 #[MW]

        f_ncfile = "toray_%d.nc"%k
        outtoray = read_toray_output(f_ncfile)

        nrho = outtoray["outtoray"]["nrho"][0]
        rho  = outtoray["outtoray"]["rho" ]     #[]
        jec  = outtoray["outtoray"]["jec" ]     #[A/m^2/W]
        qec  = outtoray["outtoray"]["qec" ]     #[W/m^3/W]
        Pec  = outtoray["outtoray"]["Pec" ][0]  #[]
        Iec  = outtoray["outtoray"]["Iec" ][0]  #[A/W]

        try:
           j_multi = intoray["adjust"]["j_multi"][k]
        except:
           j_multi = 1.0
           print 'no j_multi inuput, set 1.0'

        try:
           width = intoray["adjust"]["width"][k]
        except:
           width = 0.0
           print 'no width inuput, set 0'
        try:
           width2 = intoray["adjust"]["width2"][k]
        except:
           width2 = 0.0
           print 'no width2 inuput, set 0'
        print 'width: ',width,width2

        if width > 0.0:

            print 'adjust width'

            jec_max=0.0
            for j in range(len(rho)):
                if jec[j] > jec_max:
                   jec_max = jec[j]
                   j_find = j
            print 'rho_peak, jec_peak :',rho[j_find],jec_max

            if width2 > 0.0:
                print 'call asym'
                fit = adjust_jec_asym(
                    0.5*width,0.5*width2,rho[j_find],jec_max,Iec,nrho,ps,r0,b0)
                jec = fit

                fit = adjust_qec_asym(
                    0.5*width,0.5*width2,rho[j_find],qec[j_find],Pec,nrho,ps)
                qec = fit

            else:
                fit = adjust_jec(
                    0.5*width,rho[j_find],jec_max,Iec,nrho,ps,r0,b0)
                jec = fit

                fit = adjust_qec(
                    0.5*width,rho[j_find],qec[j_find],Pec,nrho,ps)
                qec = fit


        # scale_pow = vol_integral(rho,qec,ps)
        # print 'scale_pow =',scale_pow, rfpow/scale_pow
        #
        # scale_cur = cur_integral(rho,jec,ps,r0,abs(b0))
        # print 'scale_cur =',scale_cur, Iec/scale_cur
        #
        # if k==0:
        #    jec_sum = jec*rfpow*j_multi*Iec/scale_cur
        #    qec_sum = qec*rfpow/scale_pow
        # else:
        #    jec_sum+= jec*rfpow*j_multi*Iec/scale_cur
        #    qec_sum+= qec*rfpow/scale_pow

        if k==0:
           jec_sum = jec*rfpow*j_multi
           qec_sum = qec*rfpow
        else:
           jec_sum+= jec*rfpow*j_multi
           qec_sum+= qec*rfpow

    try:
       j0_seed = intoray["seed"]["j0"][0]
       x0_seed = intoray["seed"]["x0"][0]
       drho_seed = intoray["seed"]["drho"][0]
    except:
       j0_seed = 0.0
       drho_seed = 0.0
    print '***** seed current: j0,drho',j0_seed,drho_seed
    if j0_seed > 0.0:
        jec_seed = j0_seed*exp(-(rho-x0_seed)**2/(2*drho_seed**2))
        jec_sum += jec_seed

    # update plasma state file

    rho_ps = ps["rho"][:]
    jec_ps = 1.0e6*interp1d(rho,jec_sum,kind='cubic')(rho_ps)
    qec_ps = 1.0e6*interp1d(rho,qec_sum,kind='cubic')(rho_ps)

    ps.load_j_parallel(rho_ps,jec_ps,"rho_ecrf","curech",r0,b0)
    ps.load_vol_profile(rho_ps,qec_ps,"rho_ecrf","peech")

    # out = Namelist()
    # out["check"]["rho"] = rho
    # out["check"]["qec"] = qec_sum
    # out["check"]["rho_ps"] = rho_ps
    # out["check"]["qec_ps"] = qec_ps
    # out.write("check.dat")


###############################################################################
# check

if __name__ == "__main__":

    pass

#------------------------------------------------------------------
# NOTE

    """

    *** psiin ***

    read (2, '(4i4)') lcentr,ledge,nr,nz
    read (2, 200 )   rdim,zdim,rcenter,rinside
    read (2, 200 )   rmaxis,zmaxis,psimax,psilim,b0
    read (2, 200 )   gasep
    read (2, 200 )   (fpsi(i),i=1,nr)                # same as geqdsk
    read (2, 200 )   ((psi(i,j),i=1,nr),j=1,nz)      # same as geqdsk
    read (2, 200 )   (psir(i),i=1,ledge)             # normalized psi on rho grid
    read (2, 200 )   (qpsi(i),i=1,nr)                # same as geqdsk
    read (2,'(2i5)') npbdry,nlimtr
    read (2, 200 )   (rpbdry(i),zpbdry(i),i=1,npbdry)
    read (2, 200 )   (rlim(i),zlim(i),i=1,nlimtr)
    read (2, 200 ) (h_factr_rf(i), i=1,nr)  ! <SQRT(1-B/Bmax)> # on equi normalized psi grid, onetwo pass it as <SQRT(1-B/Bmax)>/Btor
    read (2, 200 ) (bsq_avg_rf(i),i=1,nr)   ! <B^2/B0^2>       # on equi normalized psi grid
    read (2, 200 ) (b_avg_rf(i),i=1,nr)     ! <B/B0>, B0==Btor # on equi normalized psi grid
    read (2, 200 ) (r0rinv_rf(i),i=1,nr)    ! <R0/R>           # on equi normalized psi grid, onetwo pass it using R0 = R(magnetic asis)
    200 format (5e16.9)

    *** echin ***

    read (2, 1000) idamp,j12,nray,nbfld
    read (2, 1001) fmu0,rfmod,x00,z00,thet, phai,bhalf,bsratio,rmajs, b0,rmins
    read (2, 1001) (psinrmr(j), j=1,j12)
    read (2, 1001) (zef(j),j=1,j12)
    read (2, 1001) (ene(j),j=1,j12)
    read (2, 1001) (ete(j),j=1,j12)

    *** echout ***

    write (3, 800) ledge
800 format (i4)
    write (3, 790) (xbouni(i),i=2,ledge+1)
    voltord2=voltor/2._p_
    write (3, 790) voltord2
    do 785 i=2,ledge
785 work12(i)=2._p_*dx2i(i)
    write (3, 790) (work12(i),i=2,ledge)
    write (3, 790) (weecrh(i),i=2,ledge)
    write (3, 790) (tpowde(i),i=2,ledge)
    write (3, 790) (wiecrt(i),i=2,ledge)
    write (3, 790) (tidept(i),i=2,ledge)
    !added 08/05/2005
    write (3, 790) (h_factr(i),i=2,ledge+1)
    write (3, 790) (bsq_avg(i),i=2,ledge+1)
    write (3, 790) (b_avg(i),i=2,ledge+1)
    write (3, 790) (r0rinv(i),i=2,ledge+1)
    write (3, 790) (currf(i),i=2,ledge+1)
    write (3, 790) (rjpdrho(i),i=2,ledge+1)
790 format (5e16.9)
    """

    """
    time    = 1.8      ![sec]

    idamp   = 8
    nrho    = 101
    nray    = 30
    nbfld   = 3

    freq    = 1.1e11   ! [Hz]
    wrfo    = 0.0      !
    xec     = 240.0    ! [cm]
    zec     = 67.94    ! [cm]
    thetec  = 100.0    ! [deg]
    phaiec  = 195.0    ! [deg]
    hlwec   = 1.7
    ratwec  = 1.0
    r0      = 169.55   ! [cm]
    bt0     = 1.75e4   ! [Gauss]
    ra      = 82.14    ! [cm]
    """
