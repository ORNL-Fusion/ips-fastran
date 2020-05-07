#!/usr/bin/env Python

"""
 ----------------------------------------------------------------------
 calculate nubeam namelist variables
 last modified Jun, 2013
 ----------------------------------------------------------------------
"""

import os,sys,pickle
from numpy import *
from scipy import interpolate
import Namelist
import pmds
import zd3dutil

#######################################################################
# geometry routines for beamline
#

def find_tangent(p0,p1,p2):
    x0 = p0[0]
    y0 = p0[1]
    x1 = p1[0]
    y1 = p1[1]
    x2 = p2[0]
    y2 = p2[1]
    s  = (x1-x0)*(x2-x1)+(y1-y0)*(y2-y1)
    s /= (x2-x1)*(x2-x1)+(y2-y1)*(y2-y1)
    s *= -1
    x = x1+(x2-x1)*s
    y = y1+(y2-y1)*s
    p = [x,y]
    return array(p)

def rzt2xyz(p):
    deg2rad = (2.0*pi)/360.
    r = p[0] 
    z = p[1]
    t = deg2rad*p[2]
    pp = [r*cos(t),r*sin(t),z]
    return array(pp)

def length(p1,p2):
    return ((p1[0]-p2[0])**2+(p1[1]-p2[1])**2)**0.5

def nbigeo(ps_rzt,pe_rzt,src_angle,iLR): 

    # constant

    deg2rad = (2.0*pi)/360.
    rad2deg = 360.0/(2.0*pi)

    # 150NB source geometry

    s_w = 5.0
    s_h = 6.0
    alpha0 = 0.55
    alpha1 = 1.50

    # beamlet divergence 

    div_w = 0.5
    if iLR == "L": div_h = 1.0
    elif iLR == "R": div_h = 1.3
    else:
        print "iLR = L or R", iLR
        sys.exit()

    # focal length

    focl_w = 1.0e33
    focl_h = 1.0e33

    # 150NB absolute collimator

    l_col = 483.0
    col_rect = [ (14.376-2.125)*cos(deg2rad*4.33),
                 (14.376+2.125)*cos(deg2rad*4.33),
                 20.599,                          
                 17.381 ]

    # 150 NB port bottom edge 

    rp  = 244.6
    zp  =-27.45
    p_w = 50.0
    p_h = 50.0

    #------------------------------------------------------------------
    # output structure

    n_module = 4
    varlist = ["rtcena","xlbtna","xybsca"  
              ,"nbshapa","bmwidra","bmwidza" 
              ,"nbapsha" 
              ,"xlbapa","xybapa"  
              ,"rapedga","xzpedga" 
              ,"xrapoffa","xzapoffa"
              ,"divra","divza","foclra","foclza"   
              ,"nbapsh2" 
              ,"rapedg2","xzpedg2" 
              ,"xlbapa2","xrapoff2","xzapoff2"
              ,"nlco" ]
    nubgeo={}
    for var in varlist: nubgeo[var]=[]

    #------------------------------------------------------------------
    # source/colliator geometry

    # source

    alpha0 = deg2rad*alpha0
    alpha1 = deg2rad*alpha1

    alpha  = [alpha1,alpha0,-alpha0,-alpha1]

    xoffset = [ 2.0*s_h*sin(alpha0)+s_h*sin(alpha1),
                s_h*sin(alpha0),
                s_h*sin(alpha0),
                2.0*s_h*sin(alpha0)+s_h*sin(alpha1) ] 
    yoffset = [ 2.0*s_h*cos(alpha0)+s_h*cos(alpha1),
                s_h*cos(alpha0),
               -s_h*cos(alpha0),
               -2.0*s_h*cos(alpha0)-s_h*cos(alpha1) ] 

    # collimator

    da0 = col_rect[0] 
    db0 = col_rect[1] 
    dc0 = col_rect[2] 
    dd0 = col_rect[3] 

    dh = l_col*tan(deg2rad*src_angle)

    #------------------------------------------------------------------
    # loop over each module

    for im in range(n_module):

        #--------------------------------------------------------------
        # beamline
    
        ps = rzt2xyz(ps_rzt[im])
        pe = rzt2xyz(pe_rzt[im])
        angle = arctan((ps[2]-pe[2])/length(ps[0:2],pe[0:2])) 
        pt = find_tangent([0.0,0.0],ps[0:2],pe[0:2]) 

        r_t = length([0.0,0.0],pt) 
        l_st = length(ps[0:2],pt)/cos(angle) 
        z_s = ps[2] 

        #--------------------------------------------------------------
        # aperture 1: collimator

        x0 = xoffset[im]
        y0 = yoffset[im] 
        x1 = l_col/cos(deg2rad*src_angle) 
        y1 = -tan(alpha[im])*(x1-x0)+y0
        l_sa = ((x1-x0)**2+(y1-y0)**2)**0.5
        dc = (dc0-dh+y1)*cos(alpha[im])
        dd = (dd0+dh-y1)*cos(alpha[im])

        a_h = 0.5*(dc+dd)
        da_h = a_h-dc
        a_w = 0.5*(da0+db0)
        da_w = a_w-da0
        z_a = z_s-l_sa*sin(angle) 
        
        #--------------------------------------------------------------
        # aperture 2: port

        ps = rzt2xyz(ps_rzt[im])
        pe = rzt2xyz(pe_rzt[im])

        x1 = ps[0]; y1 = ps[1]; z1 = ps[2]
        x2 = pe[0]; y2 = pe[1]; z2 = pe[2]

        a = (x2-x1)**2+(y2-y1)**2
        b = 2*( (x2-x1)*x1+(y2-y1)*y1 )
        c = x1**2+y1**2-rp**2
    
        l1 = (-b+(b**2-4.0*a*c)**0.5)/(2.0*a)
        l2 = (-b-(b**2-4.0*a*c)**0.5)/(2.0*a)
        l = min([l1,l2])

        x = x1 + l*(x2-x1)
        y = y1 + l*(y2-y1)
        z = z1 + l*(z2-z1)
        r = (x**2+y**2)**0.5
        L = ((x-x1)**2+(y-y1)**2+(z-z1)**2)**0.5

        l_sa2 = L
        a_w2 = p_w
        da_w2 = 0.0
        a_h2 = p_h
        da_h2 = a_h2-(z-zp)

        #--------------------------------------------------------------
        # output

        nubgeo["rtcena"  ].append( r_t )             
        nubgeo["xlbtna"  ].append( l_st )
        nubgeo["xybsca"  ].append( z_s )
        nubgeo["nbshapa" ].append( 1 )
        nubgeo["bmwidra" ].append( s_w )
        nubgeo["bmwidza" ].append( s_h )
        nubgeo["divra"   ].append( deg2rad*div_w )
        nubgeo["divza"   ].append( deg2rad*div_h )
        nubgeo["foclra"  ].append( focl_w )
        nubgeo["foclza"  ].append( focl_h )
        nubgeo["nbapsha" ].append( 1 )
        nubgeo["xlbapa"  ].append( l_sa )
        nubgeo["xybapa"  ].append( z_a )
        nubgeo["rapedga" ].append( a_w )
        nubgeo["xzpedga" ].append( a_h )
        nubgeo["xrapoffa"].append( da_w )
        nubgeo["xzapoffa"].append( da_h )
        nubgeo["nbapsh2" ].append( 1 )
        nubgeo["rapedg2" ].append( a_w2 )
        nubgeo["xzpedg2" ].append( a_h2 )
        nubgeo["xlbapa2" ].append( l_sa2 )
        nubgeo["xrapoff2"].append( da_w2 )
        nubgeo["xzapoff2"].append( da_h2 )
        nubgeo["nlco"    ].append( True )

    return nubgeo

def interp2d(x,y,z,x0,y0):
    return interpolate.RectBivariateSpline(x,y,z,kx=1,ky=1)(x0,y0)[0][0]

def get_nubeam_geo(btilt,stilt_L,stilt_R,nbsrc,LTO):

    # -----------------------------------------------------------------
    # return data

    nubeam_namelist = {}

    if LTO>=2:

        # -----------------------------------------------------------------
        # 15L - offaxis beam
    
        segments = ["up","mu","ml","lo"]
        x = nbsrc["stilt"]
        y = nbsrc["btilt"]
    
        ps = []; pa = []
        for segment in segments:
            id = "l"+segment
            ps.append ([ interp2d(x,y,nbsrc["rs"][id],stilt_L,btilt),
                         interp2d(x,y,nbsrc["zs"][id],stilt_L,btilt),
                         interp2d(x,y,nbsrc["ts"][id],stilt_L,btilt) ] )
            pa.append ([ interp2d(x,y,nbsrc["ra"][id],stilt_L,btilt),
                         interp2d(x,y,nbsrc["za"][id],stilt_L,btilt),
                         interp2d(x,y,nbsrc["ta"][id],stilt_L,btilt) ] )
        nubeam_namelist["15L"] = nbigeo(ps,pa,stilt_L,"L")
    
        # -----------------------------------------------------------------
        # 15R - offaxis beam
    
        ps = []; pa = []
        for segment in segments:
            id = "r"+segment
            ps.append ([ interp2d(x,y,nbsrc["rs"][id],stilt_R,btilt),
                         interp2d(x,y,nbsrc["zs"][id],stilt_R,btilt),
                         interp2d(x,y,nbsrc["ts"][id],stilt_R,btilt) ] )
            pa.append ([ interp2d(x,y,nbsrc["ra"][id],stilt_R,btilt),
                         interp2d(x,y,nbsrc["za"][id],stilt_R,btilt),
                         interp2d(x,y,nbsrc["ta"][id],stilt_R,btilt) ] )
        nubeam_namelist["15R"] = nbigeo(ps,pa,stilt_R,"R")

    else:

        # -----------------------------------------------------------------
        # 15L
    
        nubeam_namelist["15L"] = { 
           "nlco"    : [True   ],
           "nbshapa" : [1      ], "bmwidra" : [6.0   ], "bmwidza" : [24.0],
           "foclra"  : [1.0e33 ], "foclza"  : [1.0e3 ],
           "divra"   : [0.00873], "divza"   : [0.0227],
           "rtcena"  : [114.6  ], "xlbtna"  : [802.8 ], "xybsca"  : [0.0 ],
          #"xbzeta"  : [314.289],
           "xlbapa"  : [186.1  ], "xybapa"  : [0.0   ],
           "nbapsha" : [1      ], "rapedga" : [8.85  ], "xzpedga" : [24.0],
           "xrapoffa": [0.0    ], "xzapoffa": [0.0   ],
           "nbapsh2" : [0      ], "rapedg2" : [0.0   ], "xzpedg2" : [0.0 ], 
           "xlbapa2" : [0.0    ], "xrapoff2": [0.0   ], "xzapoff2": [0.0 ]
        } 
    
        # -----------------------------------------------------------------
        # 15R
    
        nubeam_namelist["15R"] = { 
           "nlco"    : [True   ],
           "nbshapa" : [1      ], "bmwidra" : [6.0   ], "bmwidza" : [24.0],
           "foclra"  : [1.0e33 ], "foclza"  : [1.0e3 ],
           "divra"   : [0.00873], "divza"   : [0.0227],
           "rtcena"  : [76.2   ], "xlbtna"  : [817.3 ], "xybsca"  : [0.0 ],
          #"xbzeta"  : [320.16 ],
           "xlbapa"  : [186.1  ], "xybapa"  : [0.0   ],
           "nbapsha" : [1      ], "rapedga" : [8.85  ], "xzpedga" : [24.0],
           "xrapoffa": [0.0    ], "xzapoffa": [0.0   ],
           "nbapsh2" : [0      ], "rapedg2" : [0.0   ], "xzpedg2" : [0.0 ], 
           "xlbapa2" : [0.0    ], "xrapoff2": [0.0   ], "xzapoff2": [0.0 ]
        } 

    # -----------------------------------------------------------------
    # 30L

    nubeam_namelist["30L"] = { 
       "nlco"    : [True   ],
       "nbshapa" : [1      ], "bmwidra" : [6.0   ], "bmwidza" : [24.0],
       "foclra"  : [1.0e33 ], "foclza"  : [1.0e3 ],
       "divra"   : [0.00873], "divza"   : [0.0227],
       "rtcena"  : [114.6  ], "xlbtna"  : [802.8 ], "xybsca"  : [0.0 ],
      #"xbzeta"  : [314.289],
       "xlbapa"  : [186.1  ], "xybapa"  : [0.0   ],
       "nbapsha" : [1      ], "rapedga" : [8.85  ], "xzpedga" : [24.0],
       "xrapoffa": [0.0    ], "xzapoffa": [0.0   ],
       "nbapsh2" : [0      ], "rapedg2" : [0.0   ], "xzpedg2" : [0.0 ], 
       "xlbapa2" : [0.0    ], "xrapoff2": [0.0   ], "xzapoff2": [0.0 ]
    } 

    # -----------------------------------------------------------------
    # 30R

    nubeam_namelist["30R"] = { 
       "nlco"    : [True   ],
       "nbshapa" : [1      ], "bmwidra" : [6.0   ], "bmwidza" : [24.0],
       "foclra"  : [1.0e33 ], "foclza"  : [1.0e3 ],
       "divra"   : [0.00873], "divza"   : [0.0227],
       "rtcena"  : [76.2   ], "xlbtna"  : [817.3 ], "xybsca"  : [0.0 ],
      #"xbzeta"  : [320.16 ],
       "xlbapa"  : [186.1  ], "xybapa"  : [0.0   ],
       "nbapsha" : [1      ], "rapedga" : [8.85  ], "xzpedga" : [24.0],
       "xrapoffa": [0.0    ], "xzapoffa": [0.0   ],
       "nbapsh2" : [0      ], "rapedg2" : [0.0   ], "xzpedg2" : [0.0 ], 
       "xlbapa2" : [0.0    ], "xrapoff2": [0.0   ], "xzapoff2": [0.0 ]
    } 

    if LTO>=1:

        # -----------------------------------------------------------------
        # 21L
    
        nubeam_namelist["21L"] = { 
           "nlco"    : [False  ],
           "nbshapa" : [1      ], "bmwidra" : [6.0   ], "bmwidza" : [24.0],
           "foclra"  : [1.0e33 ], "foclza"  : [1.0e3 ],
           "divra"   : [0.00873], "divza"   : [0.0227],
           "rtcena"  : [76.2   ], "xlbtna"  : [817.3 ], "xybsca"  : [0.0 ],
          #"xbzeta"  : [159.84 ],
           "xlbapa"  : [186.1  ], "xybapa"  : [0.0   ],
           "nbapsha" : [1      ], "rapedga" : [8.85  ], "xzpedga" : [24.0],
           "xrapoffa": [0.0    ], "xzapoffa": [0.0   ],
           "nbapsh2" : [0      ], "rapedg2" : [0.0   ], "xzpedg2" : [0.0 ], 
           "xlbapa2" : [0.0    ], "xrapoff2": [0.0   ], "xzapoff2": [0.0 ]
        } 
    
        # -----------------------------------------------------------------
        # 21R
    
        nubeam_namelist["21R"] = { 
           "nlco"    : [False  ],
           "nbshapa" : [1      ], "bmwidra" : [6.0   ], "bmwidza" : [24.0],
           "foclra"  : [1.0e33 ], "foclza"  : [1.0e3 ],
           "divra"   : [0.00873], "divza"   : [0.0227],
           "rtcena"  : [114.6  ], "xlbtna"  : [802.8 ], "xybsca"  : [0.0 ],
          #"xbzeta"  : [165.711],
           "xlbapa"  : [186.1  ], "xybapa"  : [0.0   ],
           "nbapsha" : [1      ], "rapedga" : [8.85  ], "xzpedga" : [24.0],
           "xrapoffa": [0.0    ], "xzapoffa": [0.0   ],
           "nbapsh2" : [0      ], "rapedg2" : [0.0   ], "xzpedg2" : [0.0 ], 
           "xlbapa2" : [0.0    ], "xrapoff2": [0.0   ], "xzapoff2": [0.0 ]
        } 
    
    else: 

        # -----------------------------------------------------------------
        # 21L
        nubeam_namelist["21L"] = { 
           "nlco"    : [True  ],
           "nbshapa" : [1      ], "bmwidra" : [6.0   ], "bmwidza" : [24.0],
           "foclra"  : [1.0e33 ], "foclza"  : [1.0e3 ],
           "divra"   : [0.00873], "divza"   : [0.0227],
           "rtcena"  : [114.6  ], "xlbtna"  : [802.8 ], "xybsca"  : [0.0 ],
          #"xbzeta"  : [165.711],
           "xlbapa"  : [186.1  ], "xybapa"  : [0.0   ],
           "nbapsha" : [1      ], "rapedga" : [8.85  ], "xzpedga" : [24.0],
           "xrapoffa": [0.0    ], "xzapoffa": [0.0   ],
           "nbapsh2" : [0      ], "rapedg2" : [0.0   ], "xzpedg2" : [0.0 ], 
           "xlbapa2" : [0.0    ], "xrapoff2": [0.0   ], "xzapoff2": [0.0 ]
        } 
    
    
        # -----------------------------------------------------------------
        # 21R
    
        nubeam_namelist["21R"] = { 
           "nlco"    : [True  ],
           "nbshapa" : [1      ], "bmwidra" : [6.0   ], "bmwidza" : [24.0],
           "foclra"  : [1.0e33 ], "foclza"  : [1.0e3 ],
           "divra"   : [0.00873], "divza"   : [0.0227],
           "rtcena"  : [76.2   ], "xlbtna"  : [817.3 ], "xybsca"  : [0.0 ],
          #"xbzeta"  : [159.84 ],
           "xlbapa"  : [186.1  ], "xybapa"  : [0.0   ],
           "nbapsha" : [1      ], "rapedga" : [8.85  ], "xzpedga" : [24.0],
           "xrapoffa": [0.0    ], "xzapoffa": [0.0   ],
           "nbapsh2" : [0      ], "rapedg2" : [0.0   ], "xzpedg2" : [0.0 ], 
           "xlbapa2" : [0.0    ], "xrapoff2": [0.0   ], "xzapoff2": [0.0 ]
        } 


    # -----------------------------------------------------------------
    # 33L

    nubeam_namelist["33L"] = { 
       "nlco"    : [True   ],
       "nbshapa" : [1      ], "bmwidra" : [6.0   ], "bmwidza" : [24.0],
       "foclra"  : [1.0e33 ], "foclza"  : [1.0e3 ],
       "divra"   : [0.00873], "divza"   : [0.0227],
       "rtcena"  : [114.6  ], "xlbtna"  : [802.8 ], "xybsca"  : [0.0 ],
      #"xbzeta"  : [14.2889],
       "xlbapa"  : [186.1  ], "xybapa"  : [0.0   ],
       "nbapsha" : [1      ], "rapedga" : [8.85  ], "xzpedga" : [24.0],
       "xrapoffa": [0.0    ], "xzapoffa": [0.0   ],
       "nbapsh2" : [0      ], "rapedg2" : [0.0   ], "xzpedg2" : [0.0 ], 
       "xlbapa2" : [0.0    ], "xrapoff2": [0.0   ], "xzapoff2": [0.0 ]
    } 

    # -----------------------------------------------------------------
    # 33R

    nubeam_namelist["33R"] = { 
       "nlco"    : [True   ],
       "nbshapa" : [1      ], "bmwidra" : [6.0   ], "bmwidza" : [24.0],
       "foclra"  : [1.0e33 ], "foclza"  : [1.0e3 ],
       "divra"   : [0.00873], "divza"   : [0.0227],
       "rtcena"  : [76.2   ], "xlbtna"  : [817.3 ], "xybsca"  : [0.0 ],
      #"xbzeta"  : [20.1602],
       "xlbapa"  : [186.1  ], "xybapa"  : [0.0   ],
       "nbapsha" : [1      ], "rapedga" : [8.85  ], "xzpedga" : [24.0],
       "xrapoffa": [0.0    ], "xzapoffa": [0.0   ],
       "nbapsh2" : [0      ], "rapedg2" : [0.0   ], "xzpedg2" : [0.0 ], 
       "xlbapa2" : [0.0    ], "xrapoff2": [0.0   ], "xzapoff2": [0.0 ]
    }

    # -----------------------------------------------------------------
    # return

    return nubeam_namelist

#######################################################################
# beam power, evergy, ...
#

def get_nbi_param_mds(shot,tmin,tmax,nbsrc):

    print 'shot  = ',shot
    print 'time average between t = [%4d, %4d] msec'%(tmin,tmax)

    # -----------------------------------------------------------------
    # set

    segments = ["up","mu","ml","lo"]
    x = nbsrc["stilt"]
    y = nbsrc["btilt"]

    # -----------------------------------------------------------------
    # get nbi data from mds

    print 'connecting to mds server'
    mds_server = 'atlas.gat.com'
    pmds.mdsconnect(mds_server)

    print 'open nb tree'
    pmds.mdsopen('NB',shot)

    print 'get data'
    time_beam = pmds.mdsvalue('DIM_OF(\\PINJ)')
    pinj_mds = {}
    einj_mds = {}
    for src in ['30L','30R','15L','15R','21L','21R','33L','33R']:
        pinj_mds[src] = pmds.mdsvalue('\\PINJ_'+src)
        einj_mds[src] = pmds.mdsvalue('.NB'+src+':NBVAC_SCALAR')

    btilt   = pmds.mdsvalue('.NB15L.OANB:BLPTCH_CAD')
    stilt_L = pmds.mdsvalue('.NB15L.OANB:SRCPTCH')/60.0
    stilt_R = pmds.mdsvalue('.NB15R.OANB:SRCPTCH')/60.0

    stilt = {"15L":stilt_L,"15R":stilt_R}

    print 'btilt   [deg] = ',btilt
    print 'stilt_L [deg] = ',stilt["15L"]
    print 'stilt_R [deg] = ',stilt["15R"]

    ind = (time_beam>=tmin) & (time_beam<=tmax)
    pinja  = {}
    einja  = {}
    ffulla = {}
    fhalfa = {}
    for src in ['30L','30R','15L','15R','21L','21R','33L','33R']:
        pinja[src] = average(pinj_mds[src][ind])
        einja[src] = einj_mds[src] 
        if pinja[src] < 1.0e3: 
            einja[src]=80.0e3 #<==to prevent nubeam crash
        mix = zd3dutil.beam_species_mix(einja[src]*1.0e-3)
        ffulla[src] = mix["after"][0]
        fhalfa[src] = mix["after"][1]
        print src,pinja[src],einja[src],ffulla[src],fhalfa[src]

    pmds.mdsdisconnect()

    # -----------------------------------------------------------------
    # distrbute power to 4 segments for 150LR
    # apply clipping loss power at absolute collimaor and port

    pinj_150={}
    for src in ["15L","15R"]:
        dist = []; loss = []
        for segment in segments:
            id = src[-1].lower()+segment
            dist.append (interp2d(x,y,nbsrc["dist"][id],stilt[src],btilt))
            loss.append (interp2d(x,y,nbsrc["loss"][id],stilt[src],btilt))
        dist = array(dist)
        loss = array(loss)
        pinj_150[src] = pinja[src]*dist/sum(dist)*(1.0-loss)

    # -----------------------------------------------------------------
    # return
    
    retval = {}
    for src in ['30L','30R','21L','21R','33L','33R']:
        retval[src] = {
            "pinja" :[pinja[src]], 
            "einja" :[einja[src]],
            "ffulla":[ffulla[src]], 
            "fhalfa":[fhalfa[src]],
           #"tbona" :[(1.0e3,(tmin-300)*1.0e-3)[pinja[src]>1.0e3]],
           #"tboffa":[(2.0e3,(tmax+300)*1.0e-3)[pinja[src]>1.0e3]]
            "tbona" :[(10.0, 0.0)[pinja[src]>1.0e3]],
            "tboffa":[(10.1,10.0)[pinja[src]>1.0e3]]
        }
    for src in ['15L','15R']:
        retval[src] = {
            "pinja" :[ pinj_150[src][k] for k in range(4) ],
            "einja" :4*[einja[src]],
            "ffulla":4*[ffulla[src]], 
            "fhalfa":4*[fhalfa[src]],
           #"tbona" :4*[(1.0e3,(tmin-300)*1.0e-3)[pinja[src]>1.0e3]],
           #"tboffa":4*[(2.0e3,(tmax+300)*1.0e-3)[pinja[src]>1.0e3]]
            "tbona" :4*[(10.0, 0.0)[pinja[src]>1.0e3]],
            "tboffa":4*[(10.1,10.0)[pinja[src]>1.0e3]]
        }

    return retval,btilt,stilt

def get_nbi_param_user(inpow):

    pinja = {}
    einja = {}
    ffulla = {}
    fhalfa = {}

    srclist = ['30L','30R','15L','15R','21L','21R','33L','33R']

    for src in srclist:
        pinja[src] = inpow["inpow"]["pinj_"+src][0]*1.0e6
        einja[src] = inpow["inpow"]["einj_"+src][0]*1.0e3
        mix = zd3dutil.beam_species_mix(einja[src]*1.0e-3)
        ffulla[src] = mix["after"][0]
        fhalfa[src] = mix["after"][1]

    retval = {}
    for src in ['30L','30R','21L','21R','33L','33R']:
        retval[src] = {
            "pinja" :[pinja[src]], 
            "einja" :[einja[src]],
            "ffulla":[ffulla[src]], 
            "fhalfa":[fhalfa[src]],
            "tbona" :[(10.0, 0.0)[pinja[src]>1.0e3]],
            "tboffa":[(10.1,10.0)[pinja[src]>1.0e3]]
        }
    for src in ['15L','15R']:
        retval[src] = {
            "pinja" :[ pinja[src]/4.0 for k in range(4) ],
            "einja" :4*[einja[src]],
            "ffulla":4*[ffulla[src]], 
            "fhalfa":4*[fhalfa[src]],
            "tbona" :4*[(10.0, 0.0)[pinja[src]>1.0e3]],
            "tboffa":4*[(10.1,10.0)[pinja[src]>1.0e3]]
        }

    btilt = inpow["inpow"]["beam_tilt"][0]
    stilt = { "15L":inpow["inpow"]["src_tilt_15L"][0],
              "15R":inpow["inpow"]["src_tilt_15R"][0]}

    return retval,btilt,stilt

def get_nubeam_namelist(shot,tmin,tmax
        ,f_oanb_parm_data="oanb_parm_data.dat"
        ,f_nubeam_tmpl="nubeam_tmpl.dat"
        ,Db=0.0,inpow=None):

    f_oanb_parm_data = f_oanb_parm_data 

    # read beam centerline information

    f=open(f_oanb_parm_data,"r")
    nbsrc = pickle.load(f)
    f.close()

    # read nubeam template file

    nubeam_tmpl = Namelist.Namelist(f_nubeam_tmpl)

    # get power, injection energy, etc from MDS

    if inpow:
        nubeam_nbi,btilt,stilt = get_nbi_param_user(inpow)
    else:
        nubeam_nbi,btilt,stilt = get_nbi_param_mds(shot,tmin,tmax,nbsrc)

    # get nubeam namelist for beam geomentry

    nubeam_geo = get_nubeam_geo(
                     btilt,stilt["15L"],stilt["15R"],nbsrc=nbsrc)

    geo_varlist = nubeam_geo["30L"].keys()
    nbi_varlist = ["pinja","einja","tbona","tboffa","ffulla","fhalfa"]

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

    out_namelist["nbdrive_naml"]["nbeam"] = [nbeam]
    out_namelist["nbdrive_naml"]["abeama"] = nbeam*[2.0] # hard-coded
    out_namelist["nbdrive_naml"]["xzbeama"] = nbeam*[1.0] # hard-coded

    for v in nubeam_tmpl["nbdrive_naml"].keys():
        out_namelist["nbdrive_naml"][v] = nubeam_tmpl["nbdrive_naml"][v]
    for v in nubeam_namelist:
        out_namelist["nbdrive_naml"][v] = nubeam_namelist[v]

    out_namelist["nbdrive_naml"]["ebdmax"]=[max(array(nubeam_namelist["einja"]))]

    out_namelist["nbdrive_naml"]["adiff_a"    ] = 2*[Db]
    out_namelist["nbdrive_naml"]["adiff_0"    ] = 2*[Db]
    out_namelist["nbdrive_naml"]["adiff_ntime"] = [2]
    out_namelist["nbdrive_naml"]["adiff_time" ] = [0.0,10.0]
    out_namelist["nbdrive_naml"]["adiff_xpin" ] = [2.0,2.0] 
    out_namelist["nbdrive_naml"]["adiff_xpout"] = [1.0,1.0] 

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
    
    out_namelist["nbi_update"]["NLTEST_OUTPUT"] = [0]
    out_namelist["nbi_update"]["NBBCX_BB"] = [0]
    
    if shot > 0:
       f_nubeam_namelist = "nubeam_%06d.dat"%shot
    else:
       f_nubeam_namelist = "nubeam.dat"
    out_namelist.write(f_nubeam_namelist)

    return out_namelist

# ------------------------------------
# debug

if __name__ == "__main__":

    shot = 153648
    tmin = 3605
    tmax = 4805

    nubeam_namelist = get_nubeam_namelist(shot,tmin,tmax,Db=3000.0)

    #inpow = Namelist.Namelist("inpow")
    #nubeam_namelist = get_nubeam_namelist(-1,-1,-1,Db=3000.0,inpow=inpow)




