import os
import sys
import glob
import netCDF4
import numpy as np
from fastran.equilibrium.efit_eqdsk import readg, reada, getbdry

def get_shot_time(fname):
    ishot = int(os.path.basename(fname).split('.')[0][-6:])
    itime = int(os.path.basename(fname).split('.')[-1].split('_')[0])
    return (ishot, itime)

def wrt_netcdf(ishot, efitdir, efit_id='', tmin=0, tmax=sys.maxsize, nbdry_max=101, ncfile='', bin_fluxsurf=''):
    gfiles = sorted(glob.glob(os.path.join(efitdir, 'g%06d.?????%s'%(ishot, efit_id))))
    
    times = []
    for gfile in gfiles:
        time = get_shot_time(gfile)[-1]
        if time < tmin or time > tmax: continue
        times.append(time)

    nrho = 101
    rho = np.arange(nrho)/(nrho - 1.)
    
    ntimes = len(times)
    ip = np.zeros(ntimes)
    bt = np.zeros(ntimes)
    r0 = np.zeros(ntimes)

    nbdry = np.zeros(ntimes, 'i4')
    rbdry = np.zeros((ntimes, nbdry_max))
    zbdry = np.zeros((ntimes, nbdry_max))
    jpar = np.zeros((ntimes, nrho))
    jtor = np.zeros((ntimes, nrho))
    p_eq = np.zeros((ntimes, nrho))
    q_eq = np.zeros((ntimes, nrho))
    
    for k in range(ntimes):
        gfile = os.path.join(efitdir, 'g%06d.%05d%s'%(ishot, times[k], efit_id))
        print ('processing', gfile)
        
        g  = readg(gfile)
        
        ip[k] = g['cpasma']
        bt[k] = g['bcentr']
        r0[k] = g['rzero']
        bdry = getbdry(g, nskip=1)
        nbdry[k] = bdry['nbdry']
        rbdry[k, 0:nbdry[k]] = bdry['rbdry'] 
        zbdry[k, 0:nbdry[k]] = bdry['zbdry']
    
        os.system(bin_fluxsurf + " %s"%gfile)
        current = np.loadtxt("current.dat").transpose()
        jpar[k, :] = current[1]
        jtor[k, :] = current[2]
        p_eq[k, :] = current[3]
        q_eq[k, :] = current[4]

        if k == 0:
            nlim = bdry["nlim"]
            rlim = bdry["rlim"]
            zlim = bdry["zlim"]

    if ncfile:
        print('append to ' + ncfile)
        nc = netCDF4.Dataset(ncfile, 'r+', format='NETCDF4_CLASSIC')
    else:
        ncfile = 'g%06d.nc'%ishot
        print('dump to ' + ncfile)
        nc = netCDF4.Dataset(ncfile, 'w', format='NETCDF4_CLASSIC')
    
    dim_time_eq = nc.createDimension('time_eq', ntimes)
    dim_rho_eq = nc.createDimension('rho_eq', nrho)
    dim_bdry = nc.createDimension('bdry', nbdry_max)
    dim_lim = nc.createDimension('lim', nlim)
#   dim_one = nc.createDimension('one', 1)
    
    var_time_eq = nc.createVariable('time_eq', 'i4', ('time_eq',))
    var_rho_eq = nc.createVariable('rho_eq', 'f8', ('rho_eq',))
    var_nbdry = nc.createVariable('nbdry', 'i4', ('time_eq',))
#   var_nlim = nc.createVariable('nlim', 'i4', ('one',))
    
    var_ip = nc.createVariable('ip', 'f8', ('time_eq',))
    var_bt = nc.createVariable('bt', 'f8', ('time_eq',))
    var_r0 = nc.createVariable('r0', 'f8', ('time_eq',))
    
    var_rbdry = nc.createVariable('rbdry', 'f8', ('time_eq', 'bdry'))
    var_zbdry = nc.createVariable('zbdry', 'f8', ('time_eq', 'bdry'))
    var_rlim = nc.createVariable('rlim', 'f8', ('time_eq', 'lim'))
    var_zlim = nc.createVariable('zlim', 'f8', ('time_eq', 'lim'))

    var_j_tot = nc.createVariable('j_tot', 'f8',('time_eq', 'rho_eq'))
    var_jtor = nc.createVariable('jtor', 'f8',('time_eq', 'rho_eq'))
    var_p_eq = nc.createVariable('p_eq', 'f8',('time_eq', 'rho_eq'))
    var_q_eq = nc.createVariable('q_eq', 'f8',('time_eq', 'rho_eq'))
    
    var_time_eq[:] = times
    var_rho_eq[:] = rho
    var_nbdry[:] = nbdry
#   var_nlim[:] = nlim
    
    var_rbdry[:, :] = rbdry
    var_zbdry[:, :] = zbdry
    var_rlim[:] = rlim
    var_zlim[:] = zlim

    var_j_tot[:, :] = jpar * 1.e-6
    var_jtor[:, :] = jtor * 1.e-6
    var_p_eq[:, :] = p_eq
    var_q_eq[:, :] = q_eq
    
    var_ip[:] = ip
    var_bt[:] = bt
    var_r0[:] = r0
    
    nc.close()

if __name__=="__main__":
    ishot = 153648
    efit_id = ''
    efitdir = '/fusion/pillar-archive/u/parkjm/013/20132204/153648/fit00'
    tmin = 3000
    tmax = 4000
    nbdry_max = 101

    wrt_netcdf(ishot, efitdir, efit_id, tmin, tmax, nbdry_max)


