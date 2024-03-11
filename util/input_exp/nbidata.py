'''
 ----------------------------------------------------------------------
 get nbi parameters
 ----------------------------------------------------------------------
'''

import os
import sys
import numpy as np
from scipy import interpolate
import Namelist
import pmds

import netCDF4
from fastran.util.zinterp import zinterp
from input_exp.d3d import shotinfo
from input_exp import d3dbeam

def boxsmooth(y, n):
    ys = np.zeros(len(y))
    for i in range(n, len(y) - n):
        tmp = 0
        for k in range(2*n + 1):
            tmp+=y[i - n - k]
        tmp = tmp/(2*n + 1.)
        ys[i] = ys
    return ys

def get_nbi(shot, tmin=0, tmax=10000, dt=20, ncfile=None):
    # connect mds and open nb tree
    print('connecting to mds server')
    mds_server = 'atlas.gat.com'
    pmds.mdsconnect(mds_server)

    print('open nb tree')
    pmds.mdsopen('NB',shot)

    # get data
    print('get data', shot)

    time_mds = pmds.mdsvalue('DIM_OF(\\PINJ)')
    pinj_mds = {}
    einj_mds = {}
    for src in ['30L','30R','15L','15R','21L','21R','33L','33R']: 
        pinj_mds[src] = pmds.mdsvalue(f'\\PINJ_{src}')
        einj_mds[src] = pmds.mdsvalue(f'.NB{src}:NBVAC_SCALAR')

    try:
        btilt_15 = pmds.mdsvalue('.NB15L.OANB:BLPTCH_CAD')
    except:
        try:
           btilt_15 = pmds.mdsvalue('.NB15R.OANB:BLPTCH_CAD')
        except:
           btilt_15 = 0.

    if btilt_15 > 0.:
        try:
            stilt_15L = pmds.mdsvalue('.NB15L.OANB:SRCPTCH')/60.0
        except:
            stilt_15L = 0.
        try:
            stilt_15R = pmds.mdsvalue('.NB15R.OANB:SRCPTCH')/60.0
        except:
            stilt_15R = 0.
    else:
        stilt_15L = 0.
        stilt_15R = 0.

    print('shot :', shot)
    print('150  beam   tilt [deg] :', btilt_15)
    print('150L source tilt [min] :', stilt_15L)
    print('150R source tilt [min] :', stilt_15R)
    
    LTO = shotinfo(shot)
    print('LTO :', LTO)
    if LTO >= 3:
        nb_list = ['30L','30R',
                   '15L_UP', '15L_MU', '15L_ML', '15L_LO',
                   '15R_UP', '15R_MU', '15R_ML', '15R_LO',
                   '21L_UP', '21L_MU', '21L_ML', '21L_LO',
                   '21R_UP', '21R_MU', '21R_ML', '21R_LO',
                   '33L','33R']
    elif LTO >= 2: 
        nb_list = ['30L','30R',
                   '15L_UP', '15L_MU', '15L_ML', '15L_LO',
                   '15R_UP', '15R_MU', '15R_ML', '15R_LO',
                   '21L','21R',
                   '33L','33R']
    else:
        nb_list = ['30L','30R',
                   '15L','15R',
                   '21L','21R',
                   '33L','33R']

    # process data
    times = np.arange(tmin, tmax, dt)
    num_times = len(times)

    pnbi = np.zeros((len(nb_list), num_times))
    einj = np.zeros(len(nb_list))
    ffull = np.zeros(len(nb_list))
    fhalf = np.zeros(len(nb_list))
    nbsrc = len(nb_list)*['']

    for k, src in enumerate(nb_list):
        src0 = src.split('_')[0]
        oanb = len(src.split('_')) == 2
        print(src0, oanb)
        for i in range(num_times):
            ind = (time_mds >= times[i]-0.5*dt ) & (time_mds <= times[i]+0.5*dt )
            pnbi[k, i] = np.average(pinj_mds[src0][ind])
            if oanb:
               scale = 0.25
            else:
               scale = 1.0
            pnbi[k, i] *= scale
            if pnbi[k, i] <= 3.5e4*scale: pnbi[k, i]=0.0
        einj[k] = einj_mds[src0] 
        nbsrc[k] = src
        mix = d3dbeam.beam_species_mix(einj[k] * 1.e-3)
        ffull[k] = mix['after'][0]
        fhalf[k] = mix['after'][1]

    print(einj)
    print(ffull)
    print(fhalf)

    # disconnect mds
    pmds.mdsdisconnect()

    print (times)

    # dump to netcdf
    if ncfile:
        print('append to ' + ncfile)
        nc = netCDF4.Dataset(ncfile, 'r+', format='NETCDF4_CLASSIC')
    else:
        ncfile = 'nbi_%06d.nc'%shot
        print('dump to ' + ncfile)
        nc = netCDF4.Dataset(ncfile, 'w', format='NETCDF4_CLASSIC')

    dim_nbi = nc.createDimension('nbi', len(nb_list))
    dim_time_nbi = nc.createDimension('time_nbi', len(times))
    
    nc_time_nbi = nc.createVariable('time_nbi', 'f8', ('time_nbi', ))

    nc_einj = nc.createVariable('einj', 'f8', ('nbi', ))  
    nc_ffull = nc.createVariable('ffull', 'f8', ('nbi', ))  
    nc_fhalf = nc.createVariable('fhalf', 'f8', ('nbi', ))  
    nc_pnbi = nc.createVariable('pnbi', 'f8', ('nbi', 'time_nbi' ,))        

    nc_time_nbi[:] = times[:]
    nc_einj[:]  = einj[:]          
    nc_ffull[:]  = ffull[:]          
    nc_fhalf[:]  = fhalf[:]          
    nc_pnbi[:, :] = pnbi[:, :]          
    
    nc.close()

    # output
    return { 
        'nb_list'     : nb_list,
        'time_nbi'    : times,      # [msec]
        'pnbi'        : pnbi,       # [W]
        'einj'        : einj,       # [eV]
        'btilt_15'    : btilt_15,   # [deg]
        'stilt_15L'   : stilt_15L,
        'stilt_15R'   : stilt_15R 
    }

if __name__ == '__main__':

    shot= 153648 #147634 
    nbi = get_nbi(shot)
