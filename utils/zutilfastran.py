import sys,os,re,glob,shutil
import pickle
import numpy
import Namelist
import trplot_fastran
import zfdat
import zutil

def vol_integral(y,volp,rho,rhob):

    x = numpy.arange(len(y))/(len(y)-1.0)
    nrho = len(rho)
    tmp = 0.0
    for i in range(1,nrho):
        drhoc = rhob*(rho[i]-rho[i-1])
        volpc = 0.5*(volp[i]+volp[i-1])
        xc = 0.5*(rho[i]+rho[i-1])
        yc = zutil.linear(xc,x,y)
        tmp += yc*volpc*drhoc

    return tmp

def lin_integral(x,y):

    nx = len(x)
    tmp = 0.0
    for i in range(1,nx):
        dx = x[i]-x[i-1]
        xc = 0.5*(x[i]+x[i-1])
        yc = 0.5*(y[i]+y[i-1])
        tmp += xc*yc*dx

    return tmp

def get_pressure(ncfile,mfile):

    p = trplot_fastran.trplot(ncfile)
    p.info()
    p.readall()

    metric = zfdat.read_formatted(mfile)

    ne = p.d["nee"  ]["y"][-1]
    ni = p.d["nii"  ]["y"][-1]
    te = p.d["te"   ]["y"][-1]
    ti = p.d["ti"   ]["y"][-1]
    wb = p.d["wbeam"]["y"][-1]

    press = 1.602e3*(ne*te+ni*ti) + 2.0/3.0*1.0e6*wb

    volp = metric["VOLP"]
    rho  = metric["RHO"]
    rhob = metric["RHOB"][0]
    ones = len(rho)*[1.0]

    pintg = vol_integral(press,volp,rho,rhob)
    vol   = vol_integral(ones,volp,rho,rhob)

    paxis = press[0]
    pavg  = pintg/vol
    ppeak = paxis/pavg
    print 'paxis,pavg,ppeak:',paxis,pavg,ppeak
    print 'vol:',vol

    return {"x":p.d["nee"]["x"],
            "y":press,
            "ppeak":ppeak }
                         
def get_tauR(ncfile,mfile):

    p = trplot_fastran.trplot(ncfile)
    p.info()
    p.readall()

    metric = zfdat.read_formatted(mfile)
    rho  = metric["RHO"]
    rhob = metric["RHOB"][0]

    sigma = p.d["sigma"]["y"][-1]
    R0 = 1.75
    c = 2.0*numpy.pi*rhob**2
    cond = c*lin_integral(rho,sigma)/(2.0*numpy.pi*R0)
    tauR = 0.17*R0*cond
    print 'tauR:', tauR

    return {"tauR":tauR}

def read_rlog(fname='rlog'):

    f = open(fname,"r")
    lines = f.readlines()

    d = {}
    for k in range(len(lines)):
        p = re.compile("\s*BETAN")
        if p.search(lines[k]):
            betan = lines[k].split('=')[-1].strip()
        p = re.compile("\s*fni")
        if p.search(lines[k]):
            fni = lines[k].split()[-2].strip()
        p = re.compile("\s*STORED ENERGY")
        if p.search(lines[k]):
            tmp = lines[k].split()
            we = float(tmp[-3])
            wi = float(tmp[-2])
        p = re.compile("\s*Ip,Ibs")
        if p.search(lines[k]):
            tmp = lines[k].split()
            Ip = float(tmp[-4])
            Ibs = float(tmp[-3])
            Inb = float(tmp[-2])
            Irf = float(tmp[-1])

    d['betan'] = float(betan)
    d['fni'] = float(fni)
    d['we'] = float(we)
    d['wi'] = float(wi)
    d['Ip'] = float(Ip)
    d['Ibs'] = float(Ibs)
    d['Inb'] = float(Inb)
    d['Irf'] = float(Irf)

    print d
    
    return d
