#!/usr/bin/env Python

"""
 ======================================================================
 Last modified : Aug 2013
 driver for fastran/ips
 JM
"""

import os,sys,glob,shutil,pickle,re

from numpy import *
from scipy import interpolate
import Namelist
import netCDF4
import zechdata,znbidata,zgaprofile,plasmastate,efit_eqdsk
import zd3dutil,znubeaminput
import xfastran_env

#######################################################################
# write input files
#

def interp2d(x,y,z,x0,y0):
    return interpolate.RectBivariateSpline(x,y,z,kx=1,ky=1)(x0,y0)[0][0]

#----------------------------------------------------------------------
# intoray

def wrt_intoray(ncfile,shot,tmin,tmax,outdir):

    ech = netCDF4.Dataset(ncfile,'r',format='NETCDF4')

    launch_r   = ech.variables['launch_r']

    id         = ech.variables['id'][:]
    time       = ech.variables['time'][:]
    frequency  = ech.variables['frequency'][:]
    dispersion = ech.variables['dispersion'][:]
    launch_r   = ech.variables['launch_r'][:]
    launch_z   = ech.variables['launch_z'][:]
    aziang     = ech.variables['aziang'][:,:]
    polang     = ech.variables['polang'][:,:]
    pinj       = ech.variables['pinj'][:,:]

    nech = len(id)

    polang_avg = zeros(nech)
    aziang_avg = zeros(nech)
    pinj_avg = zeros(nech)

    ind = ((time >= tmin) & (time <= tmax))

    for k in range(nech):
        polang_avg[k] = average(polang[k,:][ind])
        aziang_avg[k] = average(aziang[k,:][ind])
        pinj_avg[k] = average(pinj[k,:][ind])

    intoray = Namelist.Namelist()

    intoray["intoray"]["ntoray"] = [nech]
    intoray["intoray"]["idamp" ] = nech*[8]
    intoray["intoray"]["nray"  ] = nech*[30]
    intoray["intoray"]["nbfld" ] = nech*[3]
    intoray["intoray"]["freq"  ] = frequency[:]
    intoray["intoray"]["wrfo"  ] = nech*[0.0]
    intoray["intoray"]["x"     ] = launch_r[:]*1.0e2
    intoray["intoray"]["z"     ] = launch_z[:]*1.0e2
    intoray["intoray"]["hlw"   ] = dispersion[:]
    intoray["intoray"]["ratw"  ] = nech*[1.0]
    intoray["intoray"]["rfpow" ] = pinj_avg[:]
    intoray["intoray"]["thet"  ] = polang_avg[:]
    intoray["intoray"]["phai"  ] = aziang_avg[:]
    
    intoray["adjust"]["width"  ] = nech*[0.0]
    intoray["adjust"]["j_multi"] = nech*[1.0]

    intoray["edata"]["igafit"]    = [1]
    intoray["edata"]["netcdfdat"] = [True]

    intoray.write(outdir+'/intoray')

    ech.close()

#----------------------------------------------------------------------
# innubeam

def wrt_innubeam(ncfile,shot,tmin,tmax,outdir,Db=0.0):

    nbi = netCDF4.Dataset(ncfile,'r',format='NETCDF4')

    id        = nbi.variables["id"][:,:]
    time      = nbi.variables["time"][:]     
    einj      = nbi.variables["einj"][:]     
    pinj      = nbi.variables["pinj"][:]     
    btilt_15  = nbi.variables["btilt_15"][0] 
    stilt_15L = nbi.variables["stilt_15L"][0]
    stilt_15R = nbi.variables["stilt_15R"][0]

    nnbi = len(id)

    pinj_avg = zeros(nnbi)

    ind = ((time >= tmin) & (time <= tmax))

    for k in range(nnbi):
        pinj_avg[k] = average(pinj[k,:][ind])

    # template file names

    f_oanb_parm_data = xfastran_env.dir_template+"/oanb_parm_data.dat"
    f_nubeam_tmpl = xfastran_env.dir_template+"/innubeam_tmpl_d3d.dat"

    # read beam centerline information

    f=open(f_oanb_parm_data,"r")
    nbsrc = pickle.load(f)
    f.close()

    # read nubeam template file

    nubeam_tmpl = Namelist.Namelist(f_nubeam_tmpl)

    # get power, injection energy, etc from MDS

    nb_list = nnbi*[""]
    for k in range(nnbi):
        for kk in range(len(id[k])):
            nb_list[k] += id[k][kk]

    pinja  = zeros(nnbi)
    einja  = zeros(nnbi)
    ffulla = zeros(nnbi)
    fhalfa = zeros(nnbi)

    for k in range(nnbi):
        pinja[k] = pinj_avg[k]
        einja[k] = einj[k]
        if pinja[k] < 1.0e3: einja[k]=80.0e3 #<==to prevent nubeam crash
        mix = zd3dutil.beam_species_mix(einja[k]*1.0e-3)
        ffulla[k] = mix["after"][0]
        fhalfa[k] = mix["after"][1]
       #print k,pinja[k],einja[k],ffulla[k],fhalfa[k]

    # -----------------------------------------------------------------
    # LTO 1: after 210 counter, 2: after 150 off-axis

    if shot > 143703: LTO = 2 
    elif shot > 125042: LTO = 1
    else: LTO = 0

    # -----------------------------------------------------------------
    # distrbute power to 4 segments for 150LR
    # apply clipping loss power at absolute collimaor and port

    segments = ["up","mu","ml","lo"]
    x = nbsrc["stilt"]
    y = nbsrc["btilt"]

    stilt = {"15L":stilt_15L,"15R":stilt_15R}
    pinj_150 = {"15L":pinja[2],"15R":pinja[3]}

    for src in ["15L","15R"]:
        dist = []; loss = []
        for segment in segments:
            id = src[-1].lower()+segment
            dist.append (interp2d(x,y,nbsrc["dist"][id],stilt[src],btilt_15))
            loss.append (interp2d(x,y,nbsrc["loss"][id],stilt[src],btilt_15))
        dist = array(dist)
        loss = array(loss)
        pinj_150[src] = pinj_150[src]*dist/sum(dist)*(1.0-loss)
    
    nubeam_nbi = {}
    for k, src in enumerate(nb_list):
       #print k,src
        if src in ['15L','15R'] and LTO==2:
            nubeam_nbi[src] = {
                "pinja" :[ pinj_150[src][i] for i in range(4) ],
                "einja" :4*[einja[k]],
                "ffulla":4*[ffulla[k]], 
                "fhalfa":4*[fhalfa[k]],
                "tbona" :4*[(10.0, 0.0)[pinja[k]>1.0e3]],
                "tboffa":4*[(10.1,10.0)[pinja[k]>1.0e3]]
            }
        else:
            nubeam_nbi[src] = {
                "pinja" :[pinja[k]], 
                "einja" :[einja[k]],
                "ffulla":[ffulla[k]], 
                "fhalfa":[fhalfa[k]],
                "tbona" :[(10.0, 0.0)[pinja[k]>1.0e3]],
                "tboffa":[(10.1,10.0)[pinja[k]>1.0e3]]
            }

    # get nubeam namelist for beam geomentry

    nubeam_geo = znubeaminput.get_nubeam_geo(
                     btilt_15,stilt_15L,stilt_15R,nbsrc,LTO)

    geo_varlist = nubeam_geo["30L"].keys()
    nbi_varlist = ["pinja","einja","ffulla","fhalfa"]

    # combine all

    nubeam_namelist = {}
    for v in geo_varlist: 
        nubeam_namelist[v] = nubeam_geo["30L"][v]
    for v in nbi_varlist: 
        nubeam_namelist[v] = nubeam_nbi["30L"][v]
    
    for s in ["30R","15L","15R","21L","21R","33L","33R"]:
        for v in geo_varlist: 
            nubeam_namelist[v] += nubeam_geo[s][v]
        for v in nbi_varlist: 
            nubeam_namelist[v] += nubeam_nbi[s][v]

    # write namelist file

    out_namelist = Namelist.Namelist()

    nbeam = len(nubeam_namelist["pinja"])

    out_namelist["nubeam_run"]["dt_nubeam"] = [0.020]
    out_namelist["nubeam_run"]["nstep"]     = [20]
    out_namelist["nubeam_run"]["navg"]     = [10]

    out_namelist["nbi_config"]["nbeam"] = [nbeam]
    out_namelist["nbi_config"]["abeama"] = nbeam*[2.0] # hard-coded
    out_namelist["nbi_config"]["xzbeama"] = nbeam*[1.0] # hard-coded

    for v in nubeam_tmpl["nbi_init"].keys():
        out_namelist["nbi_init"][v] = nubeam_tmpl["nbi_init"][v]
    for v in nubeam_tmpl["nbi_update"].keys():
        out_namelist["nbi_update"][v] = nubeam_tmpl["nbi_update"][v]
    for v in nubeam_namelist:
        out_namelist["nbi_config"][v] = nubeam_namelist[v]

    out_namelist["nbi_model"]["difb_0"  ] = [Db]
    out_namelist["nbi_model"]["difb_a"  ] = [Db]
    out_namelist["nbi_model"]["difb_in" ] = [2]
    out_namelist["nbi_model"]["difb_out"] = [2]
    out_namelist["nbi_model"]["nkdifb"  ] = [3] 

    out_namelist["nbi_init"]["wghta"] = [20.0]
    out_namelist["nbi_init"]["nzones"] = [20]
    out_namelist["nbi_init"]["nzone_fb"] = [10]
    out_namelist["nbi_init"]["nznbma"] = [50]
    out_namelist["nbi_init"]["nznbme"] = [100]
    out_namelist["nbi_init"]["nptcls"] = [5000]
    out_namelist["nbi_init"]["ndep0"] = [5000]
    out_namelist["nbi_init"]["nsigexc"] = [1]
    out_namelist["nbi_init"]["nmcurb"] = [4] 
    out_namelist["nbi_init"]["nseed"] = [410338673]
    out_namelist["nbi_init"]["nltest_output"] = [0]
    out_namelist["nbi_init"]["nsdbgb"] = [2]
    
    out_namelist["nbi_update"]["nltest_output"] = [0]
    out_namelist["nbi_update"]["nbbcx_bb"] = [0]

    f_nubeam_namelist = outdir+"/innubeam"
    out_namelist.write(f_nubeam_namelist)

#-----------------------------------------------------------------------
# instate

def wrt_instate(ncfile_prf,ncfile_efit,shot,tmin,tmax,time_plasma_shape,outdir,fn_instate="instate"):

    #------------------
    # read profile data

    nc_prf = netCDF4.Dataset(ncfile_prf,'r',format='NETCDF4')
    
    prf = {} 

    time_prf = nc_prf.variables["time"][:]
    rho_prf  = nc_prf.variables["rho"][:]

    prfdat = {} 
    prfdat["ne"   ] = nc_prf.variables["ne"][:,:]
    prfdat["te"   ] = nc_prf.variables["te"][:,:]
    prfdat["ti"   ] = nc_prf.variables["ti"][:,:]
    prfdat["omega"] = nc_prf.variables["omega"][:,:]
    prfdat["nz"   ] = nc_prf.variables["nz"][:,:]
    prfdat["rad"  ] = nc_prf.variables["rad"][:,:]

    ind = ((time_prf >= tmin) & (time_prf <= tmax))

    for v in prfdat.keys():
        avg = average(prfdat[v][ind],axis=0)
        err = std(prfdat[v][ind],axis=0)
        prf[v] = {'y':avg,'e':err}

    zeff = 30.0*prf['nz']['y']/prf['ne']['y'] + 1.0
    prf['zeff'] = {'y':zeff,'e':0.0}

    nc_prf.close()

    #------------------
    # read efit data

    nc_efit = netCDF4.Dataset(ncfile_efit,'r',format='NETCDF4')
    
    efit = {} 

    time_efit = nc_efit.variables["time"][:]
    rho_efit  = nc_efit.variables["rho"][:]

    efitdat = {}
    efitdat["ip"    ] = nc_efit.variables["ip"    ][:] 
    efitdat["b0"    ] = nc_efit.variables["b0"    ][:]
    efitdat["r0"    ] = nc_efit.variables["r0"    ][:]
    efitdat["rmajor"] = nc_efit.variables["rmajor"][:]
    efitdat["aminor"] = nc_efit.variables["aminor"][:]
    efitdat["kappa" ] = nc_efit.variables["kappa" ][:]
    efitdat["delta" ] = nc_efit.variables["delta" ][:]
    efitdat["jtot"  ] = nc_efit.variables["jpar"  ][:,:] # <====
    efitdat["p_eq"  ] = nc_efit.variables["p_eq"  ][:,:]

    nbdry = nc_efit.variables["nbdry"][:]
    nlim = nc_efit.variables["nlim"][:]
    rbdry = nc_efit.variables["rbdry"][:,:]
    zbdry = nc_efit.variables["zbdry"][:,:]
    rlim = nc_efit.variables["rlim"][:,:]
    zlim = nc_efit.variables["zlim"][:,:]

    nc_efit.close()

    ind = ((time_efit >= tmin) & (time_efit <= tmax))

    for v in ["ip","b0","r0","rmajor","aminor","kappa","delta"]:
        avg = average(efitdat[v][ind])
        err = std(efitdat[v][ind])
        efit[v] = {'y':avg,'e':err}

    for v in ["jtot","p_eq"]:
        avg = average(efitdat[v][ind],axis=0)
        err = std(efitdat[v][ind],axis=0)
        efit[v] = {'y':avg,'e':err}

    ind = (time_efit == time_plasma_shape)
    nbdry = nbdry[ind][0]
    rbdry = rbdry[ind,0:nbdry][0]
    zbdry = zbdry[ind,0:nbdry][0]
    nlim  = nlim[ind][0]
    rlim  = rlim[ind,0:nlim][0]
    zlim  = zlim[ind,0:nlim][0]

    rb0 = []
    zb0 = []
    nskip = nbdry/50
    for i in range(0,nbdry,nskip):
        rb0.append(rbdry[i])
        zb0.append(zbdry[i]) 
    nbdry = len(rb0)
    rbdry = rb0
    zbdry = zb0

    #------------------
    # write instate

    instate = Namelist.Namelist()

    instate["instate"]["tokamak_id"] = ["D3D"]
    instate["instate"]["density_model"] = [0]

    for v in ["ip","b0","r0","rmajor","aminor","kappa","delta"]:
        instate["instate"][v] = [efit[v]['y']]
    instate["instate"]["ip"] = [efit["ip"]['y']*1.0e-6]

    instate["instate"]["n_ion"] = [1]
    instate["instate"]["z_ion"] = [1]
    instate["instate"]["a_ion"] = [2]
    instate["instate"]["f_ion"] = [1.0]

    instate["instate"]["n_imp"] = [1]
    instate["instate"]["z_imp"] = [6]
    instate["instate"]["a_imp"] = [12]
    instate["instate"]["f_imp"] = [1.0]

    instate["instate"]["n_min"] = [0]
    instate["instate"]["z_min"] = [1]
    instate["instate"]["a_min"] = [1]
    instate["instate"]["n_beam"] = [1]
    instate["instate"]["z_beam"] = [1]
    instate["instate"]["a_beam"] = [2]
    instate["instate"]["n_fusion"] = [1]
    instate["instate"]["z_fusion"] = [2]
    instate["instate"]["a_fusion"] = [4]

    instate["instate"]["nrho"] = [len(rho_prf)]
    instate["instate"]["rho"] = rho_prf
    instate["instate"]["ne" ] = prf["ne"]["y"] 
    instate["instate"]["te" ] = prf["te"]["y"]
    instate["instate"]["ti" ] = prf["ti"]["y"]
    instate["instate"]["omega" ] = prf["omega"]["y"]
    instate["instate"]["zeff" ] = prf["zeff"]["y"]
    instate["instate"]["p_rad" ] = prf["rad"]["y"]

    instate["instate"]["j_tot"] = efit["jtot"]["y"]*1.0e-6
    instate["instate"]["p_eq" ] = efit["p_eq"]["y"]

    zero_list = [
        "j_oh","j_bs","j_nb","j_ec","j_ic",
        "pe_nb","pe_ec","pe_ic","pe_fus","pe_ionization",#"p_rad",
        "pi_nb","pi_ec","pi_ic", "pi_fus","pi_ionization", "pi_cx",
        "p_ohm","p_ei",
        "torque_nb","torque_in",
        "se_nb","se_ionization","si_nb","si_ionization",
        "q","psipol",
        "density_beam","wbeam",
        "density_alpha","walpha",
        "chie","chii",
        ]
    
    for key in zero_list:
        instate["instate"][key] = len(rho_prf)*[0.0]

    instate["instate"]["nbdry" ] = [nbdry]
    instate["instate"]["rbdry" ] = rbdry
    instate["instate"]["zbdry" ] = zbdry
    instate["instate"]["nlim"  ] = [nlim-1]
    instate["instate"]["rlim"  ] = rlim[:-1]
    instate["instate"]["zlim"  ] = zlim[:-1]

    instate.write(outdir+"/"+fn_instate)

#-----------------------------------------------------------------------
# infastran

def wrt_infastran(input, outdir):

    def log2int(vec):
        rvec = [ (0,1)[vec[k]] for k in range(len(vec)) ]
        return rvec

    infastran = Namelist.Namelist(xfastran_env.dir_template+"/infastran_tmpl")

    infastran["infastran"]["solve_mhd"] = [-1]

    infastran["infastran"]["solve_te"] = log2int(input["in_model"]["solve_te"]) 
    infastran["infastran"]["solve_ti"] = log2int(input["in_model"]["solve_ti"]) 
    infastran["infastran"]["solve_v" ] = log2int(input["in_model"]["solve_rotation"])
    infastran["infastran"]["solve_ne"] = log2int(input["in_model"]["solve_density"])
    infastran["infastran"]["solve_j" ] = [0]
    if input["in_model"]["solve_current"][0]:
        infastran["infastran"]["relax_j"] = [1]
    if input["in_model"]["solve_te"][0] or input["in_model"]["solve_ti"][0] or \
       input["in_model"]["solve_density"][0] or \
       input["in_model"]["solve_rotation"][0] :
        infastran["infastran"]["isolver"] = [1]

    infastran["infastran"]["model_chi"] = input["in_model"]["model_transport"]
    infastran["infastran"]["model_bootstrap"] = input["in_model"]["model_bootstrap"]
    infastran["infastran"]["model_neoclass"] = input["in_model"]["model_neoclass"]

    infastran.write(outdir+"/infastran")

    shutil.copyfile(xfastran_env.dir_template+"/intglf_tmpl",outdir+"/intglf")    
   
#-----------------------------------------------------------------------
# ingeqdsk: dummy

def wrt_ingeqdsk(outdir):

    f=open(outdir+"/ingeqdsk","w")
    f.close()

#-----------------------------------------------------------------------
# configuration files

def wrt_config(input,outdir):

    f=open(xfastran_env.dir_template+"/fastran_scenario.config")
    lines = f.readlines()
    f.close()

    shot_number = input["in_d3d"]["shot"][0]
    sim_root = os.path.realpath(input["in_d3d"]["run_dir"][0])
    input_dir_sim = os.path.realpath(input["in_d3d"]["input_dir"][0])

    names = "INIT DRIVER EQ TR "
    if input["in_model"]["call_nubeam"][0]:
       names += "NB "
    if input["in_model"]["call_toray"][0]:
       names += "EC "

    f=open(outdir+"/fastran_scenario.config","w")

    for k,line in enumerate(lines):

        p0=re.compile("\s*SHOT_NUMBER = ")
        p1=re.compile("\s*SIM_ROOT = ")
        p2=re.compile("\s*INPUT_DIR_SIM = ")
        p3=re.compile("\s*NAMES = ")
        p4=re.compile("\s*USER = ")

        if p0.search(line):
           f.write("SHOT_NUMBER = %d\n"%shot_number)
        elif p1.search(line):
           f.write("SIM_ROOT = %s\n"%sim_root)
        elif p2.search(line):
           f.write("INPUT_DIR_SIM = %s\n"%input_dir_sim)
        elif p3.search(line):
           f.write("NAMES = %s\n"%names)
        elif p4.search(line):
           f.write("USER = %s\n"%os.environ["USER"])
        else:
           f.write(line)

    f.close()

    shutil.copyfile(xfastran_env.dir_template+"/venus.conf"
        ,"./venus.conf")    

#######################################################################
# driver
#

if __name__=="__main__":

    pass
