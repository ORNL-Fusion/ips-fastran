"""
 ----------------------------------------------------------------------
 get ech parameters
 ----------------------------------------------------------------------
"""

import os
import sys
from numpy import *
import Namelist
import pmds

import netCDF4
from fastran.util.zinterp import zinterp

def boxsmooth(y,n):
    ys = zeros(len(y))
    for i in range(n,nlen(y)-n):
        tmp = 0
        for k in range(2*n+1):
            tmp+=y[i-n-k]
        tmp = tmp/(2*n+1.0)
        ys[i] = ys
    return ys

def get_ech(shot, tmin=2000, tmax=6000, dt=20, ncfile=''):

    # connect mds and open rf tree
    print('connecting to mds server')
    mds_server = 'atlas.gat.com'
    pmds.mdsconnect(mds_server)

    print('open rf tree')
    pmds.mdsopen('RF',shot)

    # get data
    print('get data')
    num_antennas = pmds.mdsvalue('\TOP.ECH:NUM_ANTENNAS')
    num_systems = pmds.mdsvalue('\TOP.ECH:NUM_SYSTEMS')
    ectime = pmds.mdsvalue('\TOP.ECH:ECTIME')
    ectimes = pmds.mdsvalue('\TOP.ECH:ECTIMES')

    times = arange(tmin,tmax,dt)
    num_times = len(times)

    gname = num_systems*['']

    pech = zeros((num_systems, num_times))
    aziang = zeros((num_systems, num_times))
    polang = zeros((num_systems, num_times))
    frequency = zeros(num_systems)
    dispersion = zeros(num_systems)
    launch_r = zeros(num_systems)
    launch_z = zeros(num_systems)

    for k in range(num_systems):
        id = k+1
        gname[k] = pmds.mdsvalue('\TOP.ECH.SYSTEM_%d.GYROTRON:NAME'%id)
        frequency[k] = pmds.mdsvalue('\TOP.ECH.SYSTEM_%d.GYROTRON:FREQUENCY'%id)
        dispersion[k] = pmds.mdsvalue('\TOP.ECH.SYSTEM_%d.ANTENNA:DISPERSION'%id)
        launch_r[k] = pmds.mdsvalue('\TOP.ECH.SYSTEM_%d.ANTENNA:LAUNCH_R'%id)
        launch_z[k] = pmds.mdsvalue('\TOP.ECH.SYSTEM_%d.ANTENNA:LAUNCH_Z'%id)
        aziang_a = pmds.mdsvalue('\EC%sAZIANG'%gname[k][0:3])
        polang_a = pmds.mdsvalue('\EC%sPOLANG'%gname[k][0:3])
        pech_a = pmds.mdsvalue('\EC%sFPWRC'%gname[k][0:3])

        print('processing: %d, %s'%(k, gname[k]))

        for i in range(num_times):
            ind = (ectime >= times[i] - 0.5 * dt) & (ectime <= times[i] + 0.5 * dt)
            pech[k, i] = average(pech_a[ind])
            if pech[k, i] <= 3.5e4: pech[k, i] = 0.

        aziang_spl = zinterp(ectimes,aziang_a)
        polang_spl = zinterp(ectimes,polang_a)

        aziang[k] = aziang_spl(times)
        polang[k] = polang_spl(times)

    # disconnect mds
    pmds.mdsdisconnect()

    # output
    output = { 
        'num_systems' : num_systems,  
        'gname'       : gname,
        'num_times'   : num_times,  
        'times'       : times,      # [msec]
        'frequency'   : frequency,  # [Hz]
        'dispersion'  : dispersion, # []
        'launch_r'    : launch_r,   # [m]
        'launch_z'    : launch_z,   # [m]
        'aziang'      : aziang,     # [degree]
        'polang'      : polang,     # [degree]
        'pech'        : pech        # [w]
    }

    # dump to netcdf
    if ncfile:
        print('append to ' + ncfile)
        nc = netCDF4.Dataset(ncfile, 'r+', format='NETCDF4_CLASSIC')
    else:
        ncfile = 'ech_%06d.nc'%shot
        print('dump to ' + ncfile)
        nc = netCDF4.Dataset(ncfile, 'w', format='NETCDF4_CLASSIC')

    dim_ech = nc.createDimension('ech', output['num_systems'])
    dim_time = nc.createDimension('time', output['num_times'])
    
    ech = nc.createVariable('ech', 'i4', ('ech',))
    time = nc.createVariable('time_ech', 'f8', ('time',))

    frequency  = nc.createVariable('frequency', 'f8', ('ech',))  
    dispersion = nc.createVariable('dispersion', 'f8', ('ech',)) 
    launch_r   = nc.createVariable('launch_r', 'f8', ('ech',))      
    launch_z   = nc.createVariable('launch_z', 'f8', ('ech',))      
    aziang     = nc.createVariable('aziang', 'f8', ('ech', 'time',))      
    polang     = nc.createVariable('polang', 'f8', ('ech', 'time',))       
    pech       = nc.createVariable('pech', 'f8', ('ech', 'time',))        
    
    ech[:]   = range(output['num_systems'])
    time[:] = output['times']

    frequency[:]  = output['frequency'] 
    dispersion[:] = output['dispersion']
    launch_r[:]   = output['launch_r']  
    launch_z[:]   = output['launch_z']  
    aziang[:,:]   = output['aziang']         
    polang[:,:]   = output['polang']          
    pech[:,:]     = output['pech']           
    
    nc.close()

    i = logical_and( output['times'] >= tmin, output['times'] <= tmax)
    print(i)
    rfpow = zeros(num_systems)
    thet = zeros(num_systems)
    phai = zeros(num_systems)
    for k in range(num_systems):
        rfpow[k] = average(output['pech'][k, :])
        phai[k] = average(output['aziang'][k, :])
        thet[k] = average(output['polang'][k, :])
    print(rfpow)
    print(phai)
    print(thet)

    intoray = Namelist.Namelist()
    intoray['intoray']['ntoray'] = [num_systems]
    intoray['intoray']['idamp'] = num_systems * [8] 
    intoray['intoray']['nray'] = num_systems * [30]
    intoray['intoray']['nbfld'] = num_systems * [3]
    intoray['intoray']['freq'] = output['frequency']
    intoray['intoray']['wrfo'] = num_systems * [0.]
    intoray['intoray']['x'] = output['launch_r'] * 100.
    intoray['intoray']['z'] = output['launch_z'] * 100.
    intoray['intoray']['hlw'] = num_systems * [1.7]
    intoray['intoray']['ratw'] = num_systems * [1.0]
    intoray['intoray']['rfpow'] = rfpow
    intoray['intoray']['thet'] = thet
    intoray['intoray']['phai'] = phai
    print(intoray)
    intoray.write('intoray_%06d'%shot)

    #------------------------------------------------------------------
    # return

    return output

if __name__ == '__main__':
    from optparse import OptionParser

    parser = OptionParser()
    parser.add_option('--shot', action='store', type='int', dest='shot', default=0)
    (options, args) = parser.parse_args(sys.argv[1:])

    ech = get_ech(options.shot)




