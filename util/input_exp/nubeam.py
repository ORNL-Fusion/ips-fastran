import os
import numpy as np
from scipy import interpolate
import netCDF4
import Namelist
from input_exp import d3d
from input_exp import nbidata
from input_exp import d3dbeam

def find_tangent(p0, p1, p2):
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
    return np.array(p)

def rzt2xyz(p):
    deg2rad = (2.0*np.pi)/360.
    r = p[0]
    z = p[1]
    t = deg2rad*p[2]
    pp = [r*np.cos(t),r*np.sin(t),z]
    return np.array(pp)

def length(p1, p2):
    return ((p1[0]-p2[0])**2+(p1[1]-p2[1])**2)**0.5

def nbigeo(ps_rzt, pe_rzt, src_angle, src, module):
    # constant
    deg2rad = (2.0*np.pi)/360.
    rad2deg = 360.0/(2.0*np.pi)

    # 150NB source geometry
    s_w = 5.0
    s_h = 6.0
    alpha0 = 0.55*deg2rad
    alpha1 = 1.50*deg2rad

    # beamlet divergence
    div_w = 0.5
    div_h = {"L":1.0, "R":1.3}[src]

    # focal length
    focl_w = 1.0e33
    focl_h = 1.0e33

    # 150NB absolute collimator
    l_col = 483.0
    col_rect = [ (14.376-2.125)*np.cos(deg2rad*4.33),
                 (14.376+2.125)*np.cos(deg2rad*4.33),
                 20.599,
                 17.381 ]

    # 150 NB port bottom edge
    rp = 244.6
    zp = -27.45
    p_w = 50.0
    p_h = 50.0

    # output structure
    n_module = 4
    varlist = [
        "rtcena", "xlbtna", "xybsca",
        "nbshapa", "bmwidra", "bmwidza",
        "nbapsha",
        "xlbapa", "xybapa",
        "rapedga", "xzpedga",
        "xrapoffa", "xzapoffa",
        "divra", "divza", "foclra", "foclza",
        "nbapsh2",
        "rapedg2", "xzpedg2",
        "xlbapa2", "xrapoff2", "xzapoff2",
        "nlco" ]
    nubgeo={}
    for var in varlist: nubgeo[var]=[]

    # source
    # alpha = [alpha1, alpha0, -alpha0, -alpha1]
    alpha = {
        "up":alpha1,
        "mu":alpha0,
        "ml":-alpha0,
        "lo":-alpha1}[module]

    xoffset = {
        "up":2.0*s_h*np.sin(alpha0)+s_h*np.sin(alpha1),
        "mu":s_h*np.sin(alpha0),
        "ml":s_h*np.sin(alpha0),
        "lo":2.0*s_h*np.sin(alpha0)+s_h*np.sin(alpha1)
    }[module]

    yoffset = {
        "up": 2.0*s_h*np.cos(alpha0)+s_h*np.cos(alpha1),
        "mu": s_h*np.cos(alpha0),
        "ml":-s_h*np.cos(alpha0),
        "lo":-2.0*s_h*np.cos(alpha0)-s_h*np.cos(alpha1)
    }[module]

    # collimator
    da0 = col_rect[0]
    db0 = col_rect[1]
    dc0 = col_rect[2]
    dd0 = col_rect[3]

    dh = l_col*np.tan(deg2rad*src_angle)

    # beamline
    ps = rzt2xyz(ps_rzt)
    pe = rzt2xyz(pe_rzt)
    angle = np.arctan((ps[2]-pe[2])/length(ps[0:2], pe[0:2]))
    pt = find_tangent([0.0, 0.0], ps[0:2], pe[0:2])

    r_t = length([0.0, 0.0], pt)
    l_st = length(ps[0:2], pt)/np.cos(angle)
    z_s = ps[2]

    # aperture 1: collimator
    x0 = xoffset
    y0 = yoffset
    x1 = l_col/np.cos(deg2rad*src_angle)
    y1 = -np.tan(alpha)*(x1-x0)+y0
    l_sa = ((x1-x0)**2+(y1-y0)**2)**0.5
    dc = (dc0-dh+y1)*np.cos(alpha)
    dd = (dd0+dh-y1)*np.cos(alpha)

    a_h = 0.5*(dc+dd)
    da_h = a_h-dc
    a_w = 0.5*(da0+db0)
    da_w = a_w-da0
    z_a = z_s-l_sa*np.sin(angle)

    # aperture 2: port
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

    # output
    nubgeo={}
    nubgeo["rtcena"  ] = [r_t]
    nubgeo["xlbtna"  ] = [l_st]
    nubgeo["xybsca"  ] = [z_s]
    nubgeo["nbshapa" ] = [1]
    nubgeo["bmwidra" ] = [s_w]
    nubgeo["bmwidza" ] = [s_h]
    nubgeo["divra"   ] = [deg2rad*div_w]
    nubgeo["divza"   ] = [deg2rad*div_h]
    nubgeo["foclra"  ] = [focl_w]
    nubgeo["foclza"  ] = [focl_h]
    nubgeo["nbapsha" ] = [1]
    nubgeo["xlbapa"  ] = [l_sa]
    nubgeo["xybapa"  ] = [z_a]
    nubgeo["rapedga" ] = [a_w]
    nubgeo["xzpedga" ] = [a_h]
    nubgeo["xrapoffa"] = [da_w]
    nubgeo["xzapoffa"] = [da_h]
    nubgeo["nbapsh2" ] = [1]
    nubgeo["rapedg2" ] = [a_w2]
    nubgeo["xzpedg2" ] = [a_h2]
    nubgeo["xlbapa2" ] = [l_sa2]
    nubgeo["xrapoff2"] = [da_w2]
    nubgeo["xzapoff2"] = [da_h2]
    nubgeo["nlco"    ] = [True]

    return nubgeo

def interp2d(x, y, z, x0, y0):
    return interpolate.RectBivariateSpline(x, y, z, kx=1, ky=1)(x0, y0)[0][0]

def nubeam_geo(btilt, stilt_L, stilt_R, LTO, f_oanb_parm_data='oanb_parm_data.nc'):
    nubeam_namelist = {}

    # -----------------------------------------------------------------
    # 15L, 15R
    if LTO >= 2:
        oanb = netCDF4.Dataset(f_oanb_parm_data, 'r', format='NETCDF4')
        x = oanb.variables["stilt"][:]
        y = oanb.variables["btilt"][:]
        for src in ["l", "r"]:
            for segment in ["up", "mu", "ml", "lo"]:
                id = "l"+segment
                rs = oanb.variables["rs_{}_{}".format(src, segment)][:, :]
                zs = oanb.variables["zs_{}_{}".format(src, segment)][:, :]
                ts = oanb.variables["ts_{}_{}".format(src, segment)][:, :]

                ra = oanb.variables["ra_{}_{}".format(src, segment)][:, :]
                za = oanb.variables["za_{}_{}".format(src, segment)][:, :]
                ta = oanb.variables["ta_{}_{}".format(src, segment)][:, :]

                ps = [ interp2d(x, y, rs, stilt_L, btilt),
                       interp2d(x, y, zs, stilt_L, btilt),
                       interp2d(x, y, ts, stilt_L, btilt) ]

                pa = [ interp2d(x, y, ra, stilt_L, btilt),
                       interp2d(x, y, za, stilt_L, btilt),
                       interp2d(x, y, ta, stilt_L, btilt) ]

                nubeam_namelist["15{}_{}".format(src.upper(), segment.upper())] = nbigeo(ps, pa, stilt_L, src.upper(), segment)
        oanb.close()
    else:
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

    # -----------------------------------------------------------------
    # 21L, 21R
    if LTO>=3:
        # 21L
        _nubeam_namelist = {
            "rtcena": [122.456917576, 122.39843872, 122.291059658, 122.243912327],
            "xlbtna": [850.214028782, 849.524257065, 848.148871449, 847.559317417],
            "xybsca": [208.02, 196.72, 185.34, 173.9],
            "nbshapa": [1, 1, 1, 1],
            "bmwidra": [6.0, 6.0, 6.0, 6.0],
            "bmwidza": [6.0, 6.0, 6.0, 6.0],
            "nbapsha": [1, 1, 1, 1],
            "rapedga": [12.1943563937, 12.1943563937, 12.1943563937, 12.1943563937],
            "xzpedga": [20.8428552257, 20.8490393773, 20.8490393773, 20.8428552257],
            "xlbapa": [453.683214255, 453.76331127, 453.76331127, 453.683214255],
            "xybapa": [52.8362001721, 48.5926639197, 45.4781981306, 41.2352837522],
            "xrapoffa": [4.04777542824, 4.04777542824, 4.04777542824, 4.04777542824],
            "xzapoffa": [-4.46981680614, 0.00602435873612, 3.29382360026, 7.76868597856],
            "divra": [0.00872664625997, 0.00872664625997, 0.00872664625997, 0.00872664625997],
            "divza": [0.0174532925199, 0.0174532925199, 0.0174532925199, 0.0174532925199],
            "foclra": [1e+33, 1e+33, 1e+33, 1e+33],
            "foclza": [1e+33, 1e+33, 1e+33, 1e+33],
            "nbapsh2": [1, 1, 1, 1],
            "rapedg2": [50.0, 50.0, 50.0, 50.0],
            "xzpedg2": [50.0, 50.0, 50.0, 50.0],
            "xlbapa2": [624.883056468, 625.477413129, 625.472461949, 626.013562484],
            "xrapoff2": [0.0, 0.0, 0.0, 0.0],
            "xzapoff2": [28.2732554343, 30.0120056491, 29.9970838718, 31.7070517593],
            "nlco": [True, True, True, True]
        }
        for k, segment in enumerate(["up", "mu", "ml", "lo"]):
            src = "21L_{}".format(segment.upper())
            nubeam_namelist[src] = {}
            for key in _nubeam_namelist.keys():
                nubeam_namelist[src][key] = [_nubeam_namelist[key][k]]

        # 21R
        _nubeam_namelist = {
            "rtcena": [77.6306180199, 77.6900136326, 77.8094778528, 78.8669484228],
            "xlbtna": [867.374512465, 866.495451488, 864.90199983, 864.032624085],
            "xybsca": [208.02, 196.72, 185.34, 173.9],
            "nbshapa": [1, 1, 1, 1],
            "bmwidra": [6.0, 6.0, 6.0, 6.0],
            "bmwidza": [6.0, 6.0, 6.0, 6.0],
            "nbapsha": [1, 1, 1, 1],
            "rapedga": [12.2071139209, 12.2071139209, 12.2071139209, 12.2071139209],
            "xzpedga": [20.8428552257, 20.8490393773, 20.8490393773, 20.8428552257],
            "xlbapa": [453.683214255, 453.76331127, 453.76331127, 453.683214255],
            "xybapa": [52.8453805851, 48.6034330998, 45.4887708784, 41.2492484256],
            "xrapoffa": [-4.08604800995, -4.08604800995, -4.08604800995, -4.08604800995],
            "xzapoffa": [-4.46981680614, 0.00602435873612, 3.29382360026, 7.76868597856],
            "divra": [0.00872664625997, 0.00872664625997, 0.00872664625997, 0.00872664625997],
            "divza": [0.0209439510239, 0.0209439510239, 0.0209439510239, 0.0209439510239],
            "foclra": [1e+33, 1e+33, 1e+33, 1e+33],
            "foclza": [1e+33, 1e+33, 1e+33, 1e+33],
            "nbapsh2": [1, 1, 1, 1],
            "rapedg2": [50.0, 50.0, 50.0, 50.0],
            "xzpedg2": [50.0, 50.0, 50.0, 50.0],
            "xlbapa2": [620.533019759, 621.121100267, 621.141864848, 621.915655011],
            "xrapoff2": [0.0, 0.0, 0.0, 0.0],
            "xzapoff2": [26.7727547458, 28.5751817186, 28.6478071129, 30.489610682],
            "nlco": [True, True, True, True]
        }
        for k, segment in enumerate(["up", "mu", "ml", "lo"]):
            src = "21R_{}".format(segment.upper())
            nubeam_namelist[src] = {}
            for key in _nubeam_namelist.keys():
                nubeam_namelist[src][key] = [_nubeam_namelist[key][k]]
    elif LTO>=1:
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


    # return
    return nubeam_namelist

def wrt_innubeam(
    shot,
    tmin,
    tmax,
    outdir='.',
    Db=0.,
    f_oanb_parm_data='oanb_parm_data.nc',
    f_nubeam_tmpl='innubeam_tmpl_d3d.dat'):
    # read nbi data from mdsplus
    nbi = nbidata.get_nbi(shot)

    nb_list = nbi['nb_list']
    print(nb_list)
    nnbi = len(nb_list)
    time = nbi["time_nbi"][:]
    einj = nbi["einj"][:]
    pinj = nbi["pnbi"][:]
    btilt_15 = nbi["btilt_15"]
    stilt_15L = nbi["stilt_15L"]
    stilt_15R = nbi["stilt_15R"]

    # time average of injection power, beam species mixes
    pinja  = np.zeros(nnbi)
    einja  = np.zeros(nnbi)
    ffulla = np.zeros(nnbi)
    fhalfa = np.zeros(nnbi)

    ind = ((time >= tmin) & (time <= tmax))
    for k in range(nnbi):
        pinja[k] = np.average(pinj[k,:][ind])
        einja[k] = einj[k]
        if pinja[k] < 1.e3: einja[k] = 80.e3 #<==to prevent nubeam crash
        mix = d3dbeam.beam_species_mix(einja[k] * 1.e-3)
        ffulla[k] = mix['after'][0]
        fhalfa[k] = mix['after'][1]
        print( f'{nb_list[k]:6} , pinja = {pinja[k]*1.e-6:6.2}, einja = {einja[k]*1.e-3:6.2f}, ffulla = {ffulla[k]:4.2}, fhalfa = {fhalfa[k]:4.2}')

    # tbona tboffa

    # innubeam namelist
    out_namelist = Namelist.Namelist()

    out_namelist["nubeam_run"]["dt_nubeam"] = [0.020]
    out_namelist["nubeam_run"]["nstep"]     = [20]
    out_namelist["nubeam_run"]["navg"]     = [10]

    out_namelist["nbi_config"]["nbeam"] = [nnbi]
    out_namelist["nbi_config"]["abeama"] = nnbi*[2.0] # hard-coded
    out_namelist["nbi_config"]["xzbeama"] = nnbi*[1.0] # hard-coded

    nubeam_tmpl = Namelist.Namelist(f_nubeam_tmpl)
    for v in nubeam_tmpl["nbi_init"].keys():
        out_namelist["nbi_init"][v] = nubeam_tmpl["nbi_init"][v]
    for v in nubeam_tmpl["nbi_update"].keys():
        out_namelist["nbi_update"][v] = nubeam_tmpl["nbi_update"][v]

    out_namelist['nbi_config']['pinja'] = pinja
    out_namelist['nbi_config']['einja'] = einja
    out_namelist['nbi_config']['ffulla'] = ffulla
    out_namelist['nbi_config']['fhalfa'] = fhalfa

    LTO = d3d.shotinfo(shot)
    print(f'shot = {shot}, LTS = {LTO}')
    nubeam_namelist =  nubeam_geo(btilt_15, stilt_15L, stilt_15R, LTO, f_oanb_parm_data=f_oanb_parm_data)

    var_list = nubeam_namelist[nb_list[0]].keys()
    for v in var_list:
        out_namelist['nbi_config'][v] = nubeam_namelist[nb_list[0]][v]
    for v in var_list:
        for key in nb_list[1:]:
            out_namelist['nbi_config'][v] += nubeam_namelist[key][v]

    out_namelist["nbi_model"]["difb_0"  ] = [Db]
    out_namelist["nbi_model"]["difb_a"  ] = [Db]
    out_namelist["nbi_model"]["difb_in" ] = [2]
    out_namelist["nbi_model"]["difb_out"] = [2]
    out_namelist["nbi_model"]["nkdifb"  ] = [3]

    f_nubeam_namelist = os.path.join(outdir, f'innubeam_{shot:06}')

    print(f_nubeam_namelist)
    print(out_namelist)
    out_namelist.write(f_nubeam_namelist)
