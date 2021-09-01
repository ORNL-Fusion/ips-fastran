"""
 -----------------------------------------------------------------------
 convert ps to instate
 -----------------------------------------------------------------------
"""

from numpy import *
from Namelist import Namelist
from plasmastate  import plasmastate
from efit_eqdsk import readg

def instate_to_ps(fn_instate, ps):
    #-----------------------------------------------------------
    #-- from instate
    instate = Namelist(fn_instate)["instate"]

    for key in instate.keys(): instate[key] = array(instate[key])

    nrho = instate["nrho"][0]
    rho = instate["rho"]

    n_ion = instate["n_ion"][0]
    z_ion = instate["z_ion"]
    a_ion = instate["a_ion"]
    f_ion = instate["f_ion"]

    n_imp = instate["n_imp"][0]

    z_imp = instate["z_imp"]
    a_imp = instate["a_imp"]
    f_imp = instate["f_imp"]

    ne = instate["ne"]
    te = instate["te"]
    ti = instate["ti"]
    omega = instate["omega"]
    zeff = instate["zeff"]

    ne[0] = ne[1]
    te[0] = te[1]
    ti[0] = ti[1]
    omega[0] = omega[1]
    zeff[0] = zeff[1]

    r0 = instate["r0"][0]
    b0 = abs(instate["b0"][0])

    density_model = instate["density_model"][0]

    #-- put zeros if not defined
    for key in [
        "j_oh", "j_bs", "j_nb", "j_ec", "j_ic", \
        "pe_nb", "pe_ec", "pe_ic", "pe_fus", "pe_ionization", "p_rad", \
        "pi_nb", "pi_ec", "pi_ic", "pi_fus", "pi_ionization", "pi_cx", "p_ohm", "p_ei", \
        "torque_nb", "torque_in", "se_nb", "se_ionization", "si_nb", "si_ionization", \
        "q", "psipol", \
        "density_beam", "wbeam", "density_alpha", "walpha", \
        "chie", "chii", "p_eq" ]:
        if key.upper() not in instate.keys(): instate[key] = zeros(nrho)

    density_beam = instate["density_beam"]
    density_alpha = instate["density_alpha"]
    wbeam = instate["wbeam"]
    torque = instate["torque_nb"]

    #-----------------------------------------------------------
    #-- update ps
    rho = ps["rho"][:]

    #-- density
    density_ion = {}
    for k in range(n_ion):
        density_ion[k] = zeros(nrho)

    density_imp = {}
    for k in range(n_imp):
        density_imp[k] = zeros(nrho)

    density_th = zeros(nrho)

    if density_model not in [0,1]:
        raise Exception("density_model error")

    if density_model == 0:
        print ('density_model = 0')

        a=0; b=0; c=0; d=0
        for k in range(n_imp):
            b = b+f_imp[k]*z_imp[k]
            d = d+f_imp[k]*z_imp[k]*z_imp[k]
        for k in range(n_ion):
            a = a+f_ion[k]*z_ion[k]
            c = c+f_ion[k]*z_ion[k]*z_ion[k]

        for i in range(nrho):
            zne_adj = ne[i]
            zzne_adj = ne[i]*zeff[i]

            # depletion due to beam ions
            zne_adj = zne_adj - 1.0*density_beam[i]
            zzne_adj = zzne_adj - 1.0**2*density_beam[i]

            # effective main ion and impurity densities
            nion = (zne_adj *d-zzne_adj*b)/(a*d-b*c)
            nimp = (zzne_adj*a-zne_adj *c)/(a*d-b*c)

            for k in range(n_ion):
                density_ion[k][i] = f_ion[k]*nion
            for k in range(n_imp):
                density_imp[k][i] = f_imp[k]*nimp

    elif density_model == 1:
        print ('density_model = 1: impurity = f*ne')

        for k in range(n_imp):
            density_imp[k] = ne*f_imp[k]
        nith = ne - density_beam - 2.0*density_alpha
        for k in range(n_imp):
            nith = nith - density_imp[k]
        for k in range(n_ion):
            density_ion[k] = nith*f_ion[k]
        zeff = nith + 2.0**2*density_alpha + density_beam
        for k in range(n_imp):
            zeff = zeff + z_imp[k]**2*density_imp[k]
        zeff = zeff/ne

    for k in range(n_ion):
        density_th = density_th + density_ion[k]
    for k in range(n_imp):
        density_th = density_th + density_imp[k]

    ps["ns"][0,:] = 1.0e19*ps.node2cell(ne)
    for k in range(n_ion):
        ps["ns"][k+1,:] = 1.0e19*ps.node2cell(density_ion[k])
    for k in range(n_imp):
        ps["ns"][k+n_ion+1,:] = 1.0e19*ps.node2cell(density_imp[k])
    ps["ni"][:] = 1.0e19*ps.node2cell(density_th)

    #-- beam
    ps["nbeami"][0][:] = 1.0e19*ps.node2cell(density_beam)
    ps["eperp_beami"][0][:] = 2.0*20.0*ones(nrho-1)
    ps["epll_beami"][0][:] = 20.0*ones(nrho-1)

    tbeam = wbeam/1.602e-3/(density_beam+1.0e-6)
    ps["eperp_beami"][0][:] = ps.node2cell(2.0*tbeam/3.0)
    ps["epll_beami"][0][:] = ps.node2cell(tbeam/3.0)

    #-- temperature
    ps["Ts"][0,:] = ps.node2cell(te)

    for k in range(n_ion):
        ps["Ts"][k+1,:] = ps.node2cell(ti)
    for k in range(n_imp):
        ps["Ts"][k+n_ion+1,:] = ps.node2cell(ti)

    ps["Ti"][:] = ps.node2cell(ti)

    #-- zeff
    ps["Zeff"][:] = ps.node2cell(zeff)
    ps["Zeff_th"][:] = ps.node2cell(zeff)

    #-- rotation
    ps["omegat"][:] = ps.node2cell(omega)

    #-- current
    j_tot = 1.e6*instate["j_tot"]
    j_tot[0] = j_tot[1]
    ps.load_j_parallel(rho, j_tot, "rho_eq", "curt", r0, b0, tot=True)

    for key in ["j_nb", "j_ec", "j_ic", "j_bs", "j_oh"]:
        if key.upper() not in instate.keys(): instate[key] = zeros(nrho)

    j_nb  = 1.e6*instate["j_nb"]
    j_ec  = 1.e6*instate["j_ec"]
    j_ic  = 1.e6*instate["j_ic"]
    j_bs  = 1.e6*instate["j_bs"]
    j_oh  = 1.e6*instate["j_oh"]

    ps.load_j_parallel(rho,j_nb, "rho_nbi", "curbeam", r0, b0)
    ps.load_j_parallel(rho,j_ec, "rho_ecrf", "curech", r0, b0)
    ps.load_j_parallel(rho,j_ic, "rho_icrf", "curich", r0, b0)
    ps.load_j_parallel(rho,j_bs, "rho", "curr_bootstrap", r0, b0)
    ps.load_j_parallel(rho,j_oh, "rho", "curr_ohmic", r0 ,b0)

    #-- MHD pressure
    ps["P_eq"][:] = instate["p_eq"]

    #-- heating
    for key in ["pe_nb", "pi_nb", "pe_ec", "pe_ic", "pi_ic", "pe_fus", "pi_fus"]:
        if key.upper() not in instate.keys(): instate[key] = zeros(nrho)

    pe_nb  = 1.e6*instate["pe_nb" ]
    pi_nb  = 1.e6*instate["pi_nb" ]
    pe_ec  = 1.e6*instate["pe_ec" ]
    pe_ic  = 1.e6*instate["pe_ic" ]
    pi_ic  = 1.e6*instate["pi_ic" ]
    pe_fus = 1.e6*instate["pe_fus"]
    pi_fus = 1.e6*instate["pi_fus"]

    ps.load_vol_profile (rho, pe_nb, "rho_nbi", "pbe")
    ps.load_vol_profile (rho, pi_nb, "rho_nbi", "pbi")
    ps.load_vol_profile (rho, pe_ec, "rho_ecrf", "peech")
    ps.load_vol_profile (rho, pe_ic, "rho_icrf", "picrf_totals", k=0)
    ps.load_vol_profile (rho, pi_ic, "rho_icrf", "picth")

    #-- particle source
    se_nb = 1.e19*instate["se_nb"]

    ps.load_vol_profile (rho, se_nb, "rho_nbi", "sbedep")

    #-- torque
    ps.load_vol_profile(rho, torque, "rho_nbi", "tqbe")

    #-- limiter
    # ps["rlim"][:] = instate["rlim"]
    # ps["zlim"][:] = instate["zlim"]

    #-- temp
    # ps["sc0"][:] = 2.0e21
    # ps["n0norm"][:] = 1.0e-10
    # ps["T0sc0"][:] = 0.01
    # ps["sc0_to_sgas"][:] = 1

def ps_to_instate(f_state, f_eqdsk, f_bc, f_instate, rdir='.'):
    geq = readg(f_eqdsk)
    r0  = geq["rzero" ]
    b0  = abs(geq["bcentr"])
    ip  = geq['cpasma']

    ps = plasmastate('ips', 1)
    ps.read(f_state)

    spec = ps.get_species()

    amain = spec["amain"]
    zmain = spec["zmain"]
    n_imp = spec["n_imp"]
    a_imp = spec["a_imp"]
    z_imp = spec["z_imp"]
    f_imp = spec["f_imp"]
    f_imp_0 = [f_imp[k][0] for k in range(n_imp)]

    n_ion = 1

    nrho  = len(ps["rho"])
    rho   = ps["rho"][:]
    ne    = ps["ns"][0,:]*1.0e-19
    te    = ps["Ts"][0,:]
    ti    = ps["Ti"][:]
    zeff  = ps["Zeff"][:]
    omega = ps["omegat"][:]
    ne    = ps.cell2node_bdry(ne)
    te    = ps.cell2node_bdry(te)
    ti    = ps.cell2node_bdry(ti)
    zeff  = ps.cell2node_bdry(zeff)
    omega = ps.cell2node_bdry(omega)
    q     = zeros(nrho)
    fp    = zeros(nrho)

    j_tot = 1.e-6*ps.dump_j_parallel(rho, "rho_eq", "curt", r0, b0, tot=True)
    j_nb  = 1.e-6*ps.dump_j_parallel(rho, "rho_nbi", "curbeam", r0, b0)
    j_ec  = 1.e-6*ps.dump_j_parallel(rho, "rho_ecrf", "curech", r0, b0)
    j_ic  = 1.e-6*ps.dump_j_parallel(rho, "rho_icrf", "curich", r0, b0)
    j_lh  = 1.e-6*ps.dump_j_parallel(rho, "rho_lhrf", "curlh", r0, b0)
    j_bs  = zeros(nrho)
    j_oh  = zeros(nrho)

    density_beam = ps.dump_profile(rho, "rho_nbi", "nbeami", k=0)*1.e-19
    wbeam = ps.dump_profile(rho, "rho_nbi", "eperp_beami", k=0) \
          + ps.dump_profile(rho, "rho_nbi", "epll_beami", k=0)
    wbeam = density_beam*wbeam*1.602e-3 #MJ/m**3

    pe_nb  = ps.dump_vol_profile(rho, "rho_nbi", "pbe" )*1.e-6
    pi_nb  = (ps.dump_vol_profile(rho, "rho_nbi", "pbi" )+ps.dump_vol_profile(rho, "rho_nbi", "pbth"))*1.e-6
    pth_nb = ps.dump_vol_profile(rho, "rho_nbi", "pbth")*1.e-6

    density_alpha = ps.dump_profile(rho, "rho_fus", "nfusi", k=0)*1.e-19
    walpha = ps.dump_profile(rho, "rho_fus", "eperp_fusi", k=0) \
         + ps.dump_profile(rho, "rho_fus", "epll_fusi", k=0)
    walpha = density_alpha*walpha*1.602e-3 #MJ/m**3

    pe_fus  = ps.dump_vol_profile(rho, "rho_fus", "pfuse" )*1.e-6
    pi_fus  = ps.dump_vol_profile(rho, "rho_fus", "pfusi" )*1.e-6
    pth_fus = ps.dump_vol_profile(rho, "rho_fus", "pfusth")*1.e-6

    pe_ec  = 1.e-6*ps.dump_vol_profile(rho, "rho_ecrf", "peech")
    pe_ic  = 1.e-6*ps.dump_vol_profile(rho, "rho_icrf", "picrf_totals", k=0)
    pe_lh  = 1.e-6*ps.dump_vol_profile(rho, "rho_lhrf", "pelh")

    pi_ec = zeros(nrho)
    pi_ic  = 1.0e-6*ps.dump_vol_profile(rho, "rho_icrf", "picrf_totals", k=1)
    pi_lh  = 1.e-6*ps.dump_vol_profile(rho, "rho_lhrf", "pilh")

    tqbe = ps.dump_vol_profile(rho, "rho_nbi", "tqbe")
    tqbi = ps.dump_vol_profile(rho, "rho_nbi", "tqbi")
    tqbjxb = ps.dump_vol_profile(rho, "rho_nbi", "tqbjxb")
    tqbth = ps.dump_vol_profile(rho, "rho_nbi", "tqbth")

    torque_nb= tqbe+tqbi+tqbjxb+tqbth
    torque_in= zeros(nrho)

    se_nb = 1.e-19*(ps.dump_vol_profile(rho, "rho_nbi", "sbedep")+ps.dump_vol_profile(rho, "rho_nbi", "sbehalo"))

    p_ei = zeros(nrho)
    p_rad = zeros(nrho)
    p_ohm = zeros(nrho)
    pe_ionization = zeros(nrho)
    pi_ionization = zeros(nrho)
    pi_cx = zeros(nrho)
    si_nb = zeros(nrho)
    chie = zeros(nrho)
    chii = zeros(nrho)

    se_ionization = zeros(nrho)
    si_ionization = zeros(nrho)

    #-------------------------------------------------------------------
    # write instate
    instate = Namelist(f_instate)
    TOKAMAK_ID = ['D3D']
    instate["instate"]["r0"    ] = [r0]
    instate["instate"]["b0"    ] = [b0]
    instate["instate"]["ip"    ] = [ip*1.0e-6]
    instate["instate"]["n_ion" ] = [n_ion]
    instate["instate"]["z_ion" ] = [zmain[0]]
    instate["instate"]["a_ion" ] = [amain[0]]
    instate["instate"]["f_ion" ] = [1.0]
    instate["instate"]["n_imp" ] = [n_imp]
    instate["instate"]["z_imp" ] = z_imp
    instate["instate"]["a_imp" ] = a_imp
    instate["instate"]["f_imp" ] = ones(n_imp)/n_imp
    instate["instate"]["nrho"  ] = [nrho]
    instate["instate"]["rho"   ] = rho
    instate["instate"]["ne"    ] = ne
    instate["instate"]["te"    ] = te
    instate["instate"]["ti"    ] = ti
    instate["instate"]["zeff"  ] = zeff
    instate["instate"]["omega" ] = omega
    instate["instate"]["j_tot" ] = j_tot
    instate["instate"]["j_oh"  ] = j_oh
    instate["instate"]["j_bs"  ] = j_bs
    instate["instate"]["j_nb"  ] = j_nb
    instate["instate"]["j_ec"  ] = j_ec
    instate["instate"]["j_ic"  ] = j_ic
    instate["instate"]["j_lh"  ] = j_lh
    instate["instate"]["pe_nb" ] = pe_nb
    instate["instate"]["pe_ec" ] = pe_ec
    instate["instate"]["pe_ic" ] = pe_ic
    instate["instate"]["pe_lh" ] = pe_lh
    instate["instate"]["pe_fus"] = pe_fus
    instate["instate"]["pe_ionization"] = pe_ionization
    instate["instate"]["p_rad" ] = p_rad
    instate["instate"]["pi_nb" ] = pi_nb
    instate["instate"]["pi_ec" ] = pi_ec
    instate["instate"]["pi_ic" ] = pi_ic
    instate["instate"]["pi_lh" ] = pi_lh
    instate["instate"]["pi_fus"] = pi_fus
    instate["instate"]["pi_ionization"]  = pi_ionization
    instate["instate"]["pi_cx" ] = pi_cx
    instate["instate"]["p_ohm" ] = p_ohm
    instate["instate"]["p_ei"  ] = p_ei
    instate["instate"]["torque_nb"] = torque_nb
    instate["instate"]["torque_in"] = torque_in
    instate["instate"]["se_nb" ] = se_nb
    instate["instate"]["se_ionization"] = se_ionization
    instate["instate"]["si_nb" ] = si_nb
    instate["instate"]["si_ionization"] = si_ionization
    instate["instate"]["q"     ] = q
   #instate["instate"]["psipol"] = psipol
    instate["instate"]["density_beam"] = density_beam
    instate["instate"]["wbeam" ] = wbeam
    instate["instate"]["density_alpha"] = density_alpha
    instate["instate"]["walpha"] = walpha
    instate["instate"]["chie"  ] = chie
    instate["instate"]["chii"  ] = chii
    instate["instate"]["nbdry"]  = [geq["nbdry"]]
    instate["instate"]["rbdry"]  = geq["rbdry"]
    instate["instate"]["zbdry"]  = geq["zbdry"]
    instate["instate"]["nlim" ]  = [geq["nlim"]]
    instate["instate"]["rlim" ]  = geq["rlim"]
    instate["instate"]["zlim" ]  = geq["zlim"]
    instate.write(f_instate)
