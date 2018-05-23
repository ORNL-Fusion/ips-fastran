from numpy import *

mu0 = 4.0*pi*1.0e-7

def wth(ne,ni,nz,te,ti):
    return 1.5*1.602e3*(ne*te+(ni+nz)*ti)

def betan(w,vol,ip,b0,a0):
    wsum = sum([(vol[i+1]-vol[i])*w[i] for i in range(len(w)-1)])
    betan = wsum/vol[-1]/1.5
    betan *= 2.0*mu0/b0**2
    betan /= fabs(ip/(a0*b0))
    betan *= 1.0e2
    return betan

def betan_local(w,vol,ip,b0,a0):
    mu0 = 4.0*pi*1.0e-7
    betan = wsum/1.5
    betan *= 2.0*mu0/b0**2
    betan /= fabs(ip/(a0*b0))
    betan *= 1.0e2
    return betan

#betan = 1.602e5*(ne*te+nth*ti)*2.*mu0/b0**2*a0*b0/ip

def qmhd(a0,r0,b0,ip,kappa,delta):
    eps = a0/r0
    return 5.*a0**2/r0*b0/ip*(1.+kappa**2*(1.+2.*delta**2-1.2*delta**3))/2.*(1.17-0.65*eps)/(1.-eps**2)**2

def calculate_ion_density(z_ion,z_imp,f_ion,f_imp,ne,zeff,density_beam):

    a = sum(f_ion*z_ion)
    b = sum(f_imp*z_imp)
    c = sum(f_ion*z_ion)
    d = sum(f_imp*z_imp*z_imp)

    zne_adj = ne
    zzne_adj = ne*zeff

    zne_adj = zne_adj - 1.0*density_beam
    zzne_adj = zzne_adj - 1.0**2*density_beam

    nion = (zne_adj *d-zzne_adj*b)/(a*d-b*c)
    nimp = (zzne_adj*a-zne_adj *c)/(a*d-b*c)

    density_ion = array([f_ion[k]*nion for k in range(n_ion)])
    density_imp = array([f_imp[k]*nimp for k in range(n_imp)])

    density_th = array([sum(tmp) for tmp in density_ion.transpose()]) 
    density_th+= array([sum(tmp) for tmp in density_imp.transpose()]) 

    return density_ion, density_imp, density_th
