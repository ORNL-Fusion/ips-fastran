import re
from numpy import *
from Namelist import Namelist
from plasmastate import plasmastate
import outone_io
from efit_eqdsk import readg

def set_default():
    inone = Namelist()
    inone["namelis1"]["btor"     ] = [0.0]
    inone["namelis1"]["rmajor"   ] = [0.0]
    inone["namelis1"]["rminor"   ] = [200.0]
    inone["namelis1"]["zax"      ] = [0.0]
    inone["namelis1"]["rin"      ] = [350.0]
    inone["namelis1"]["rout"     ] = [890.0]
    inone["namelis1"]["totcur"   ] = [0.0]
    inone["namelis1"]["nj"       ] = [201]
    inone["namelis1"]["imesh"    ] = [0]
    inone["namelis1"]["time0"    ] = [0.0]
    inone["namelis1"]["timmax"   ] = [-1.0]
    inone["namelis1"]["dt"       ] = [0.001]
    inone["namelis1"]["dtmax"    ] = [0.01]
    inone["namelis1"]["dtmin"    ] = [0.0001]
    inone["namelis1"]["ilimdt"   ] = [1]
    inone["namelis1"]["itmax"    ] = [10]
    inone["namelis1"]["nmax"     ] = [20000000]
    inone["namelis1"]["relit"    ] = [0.01]
    inone["namelis1"]["relmax"   ] = [0.3]
    inone["namelis1"]["relmin"   ] = [0.1]
    inone["namelis1"]["steady_state"] = [1.0]
    inone["namelis1"]["theta"    ] = [0.8]
    inone["namelis1"]["prtlst"   ] = [0.0]
    inone["namelis1"]["iangrot"  ] = [1]
    inone["namelis1"]["itangrot" ] = [0]
    inone["namelis1"]["itenp"    ] = [0, 0]
    inone["namelis1"]["itte"     ] = [0]
    inone["namelis1"]["itti"     ] = [0]
    inone["namelis1"]["itxj"     ] = [0]
    inone["namelis1"]["inenez"   ] = [0]
    inone["namelis1"]["nprim"    ] = [1]
    inone["namelis1"]["namep"    ] = ['d']
    inone["namelis1"]["zfrac"    ] = [1.0]
    inone["namelis1"]["nimp"     ] = [1]
    inone["namelis1"]["namei"    ] = ['c']
    inone["namelis1"]["nneu"     ] = [1]
    inone["namelis1"]["namen"    ] = ['d']
    inone["namelis1"]["no_te_convection"] = [0]
    inone["namelis1"]["no_ti_convection"] = [0]
    inone["namelis1"]["taupin"   ] = [12.5]
    inone["namelis1"]["iwangrot" ] = [-3 ]
    inone["namelis1"]["jhirsh"   ] = [112] #[112]
    inone["namelis1"]["resistive"] = ['hirshman']
    inone["namelis1"]["wneo(1,1)"] = [1.0, 0.0, 0.0, 1.0, 0.0]
    inone["namelis1"]["wneo(1,2)"] = [0.0, 1.0, 0.0, 1.0, 0.0]
    inone["namelis1"]["wneo(1,3)"] = [0.0, 0.0, 1.0, 1.0, 0.0]
    inone["namelis1"]["wneo(1,4)"] = [0.0, 0.0, 0.0, 1.0, 0.0]
    inone["namelis1"]["wneo(1,5)"] = [0.0, 0.0, 0.0, 0.0, 2.0]
    inone["namelis1"]["wneot"    ] = [1.0]
    inone["namelis1"]["adjzeff"  ] = [1]
    inone["namelis1"]["nbctim"   ] = [1]
    inone["namelis1"]["bctime"   ] = [0.0]
    inone["namelis1"]["splninpt" ] = ['new']
    inone["namelis1"]["tportvb"  ] = [1]
    inone["namelis1"]["namelistvb"] = [1]
    inone["namelis1"]["newtonvb" ] = [0]
    inone["namelis1"]["fiziksvb" ] = [0]
    inone["namelis1"]["fusionvb" ] = [0]
    inone["namelis1"]["gridgenvb"] = [0]
    inone["namelis1"]["neucgvb"  ] = [0]
    inone["namelis1"]["zenvb"    ] = [0]
    inone["namelis1"]["do_eqplot"] = [0]
    inone["namelis1"]["run_preplt"] = [False]
    inone["namelis1"]["jbal"     ] = [1]
    inone["namelis1"]["jcoef"    ] = [1]
    inone["namelis1"]["jflux"    ] = [1]
    inone["namelis1"]["jprt"     ] = [1]
    inone["namelis1"]["jsourc"   ] = [1]
    inone["namelis1"]["jsxr"     ] = [0]
    inone["namelis1"]["jtfus"    ] = [0]
    inone["namelis1"]["momtm_file"] = [0]
    inone["namelis1"]["steps_per_plot"] = [1]

    inone["namelis2"]["ishot"    ] = [0]
    inone["namelis2"]["itime"    ] = [0.0]
    inone["namelis2"]["iyoka"    ] = [3]

    inone["namelis3"]["create_gcnmp_input"] = [False]
    inone["namelis3"]["iterdb"] = [1]
    inone["namelis3"]["wrt_kinetic_efit"] = [False]
    inone["namelis3"]["eqdskin"] = ['geqdsk']
    inone["namelis3"]["ieqdsk"] = [1]
    inone["namelis3"]["ifixshap"] = [1]
    inone["namelis3"]["irguess"] = [-1]
    inone["namelis3"]["nlimiter"] = [-1]

    inone.setHead("dee   diii-d\n\n---------\n")

    return inone

def write_inputfiles(f_state, f_eqdsk, f_infreya, dir_data=''):
    ps = plasmastate('ips', 1)
    ps.read(f_state)

    nrho  = len(ps["rho"])
    rho   = ps["rho"][:]
    ne    = ps.cell2node_bdry(ps["ns"][0,:])
    te    = abs(ps.cell2node_bdry(ps["Ts"][0,:]))
    ti    = abs(ps.cell2node_bdry(ps["Ti"][:]))
    zeff  = ps.cell2node_bdry(ps["Zeff"][:])
    omega = ps.cell2node_bdry(ps["omegat"][:])

    ps_xe  = 1.6022e-19
    ps_mp  = 1.6726e-27

    z_spec = [round(x) for x in ps["qatom_S"][1:]/ps_xe ]
    a_spec = [round(x) for x in ps["m_S"][1:]/ps_mp ]
    n_spec = len(z_spec)

    z_ion = []
    a_ion = []
    z_imp = []
    a_imp = []
    np = []
    nz = []

    spec_list = [ key.strip() for key in ps["S_name"][1:] ]
    print(spec_list)
    namep = []
    namei = ['c']

    for k in range(n_spec):
        if spec_list[k] in ['H', 'D', 'T']:
           print('ion', spec_list[k])
           z_ion.append(z_spec[k])
           a_ion.append(a_spec[k])
           np.append(ps.cell2node_bdry(ps["ns"][k+1,:]))
           namep.append(spec_list[k].lower())
        elif spec_list[k] == 'He4':
           print('He4 found')
           nhe4 = ps.cell2node_bdry(ps["ns"][k+1,:])
           namei.append('he')
           if min(nhe4) == 0:
               print('nhe4 zero')
               nhe4 = ones(len(nhe4))
           nhe4[-1] = abs(nhe4[-1]) #<=====
        else:
           print('imp', spec_list[k])
           z_imp.append(z_spec[k])
           a_imp.append(a_spec[k])
           nz.append(ps.cell2node_bdry(ps["ns"][k+1,:]))
#          nhe4 = zeros(nrho)
           nhe4 = ones(nrho)

    n_ion = len(z_ion)
    n_imp = len(z_imp)

    print('z_ion', z_ion)
    print('a_ion', a_ion)
    print('z_imp', z_imp)
    print('a_imp', a_imp)
    print('namep', namep)
    print('namei', namei)

    Z_C = 6.0
    Z_imp = zeros(nrho)
    dZ_eff = zeros(nrho)
    tmp =  zeros(nrho)

    for k in range(n_imp):
        tmp  += nz[k]
        Z_imp += z_imp[k]*nz[k]
        dZ_eff += z_imp[k]*(z_imp[k]-1.0)*nz[k]
    tmp += nhe4
    Z_imp += 2.*nhe4
    dZ_eff += 2.*(2.-1.)*nhe4

    Z_imp = Z_imp/tmp
    dZ_eff = dZ_eff/ne
    f_C = dZ_eff/(Z_C*(Z_C-1.0))
    nc = f_C*ne

    inone = set_default()

    if dir_data:
       inone["namelis1"]["onetwo_xsct_ext"] = [dir_data]

    innfreya = Namelist(f_infreya)
    for key in innfreya["innfreya"].keys():
        if key.lower() in ["rmajor","rminor","rin","rout","inenez","zfrac"]:
            inone["namelis1"][key] = innfreya["innfreya"][key]
        else:
            inone["namelis2"][key] = innfreya["innfreya"][key]

    inone["namelis1"]["nprim"] = [len(namep)]
    inone["namelis1"]["namep"] = namep
    inone["namelis1"]["nimp"] = [len(namei)]
    inone["namelis1"]["namei"] = namei
    inone["namelis1"]["nneu"] = [len(namep)]
    inone["namelis1"]["namen"] = namep

    inone["namelis1"]["inenez"] = [0]
    inone["namelis1"]["nj"] = [nrho]
    inone["namelis1"]["njene"] = [nrho]
    inone["namelis1"]["enein"] = ne*1.0e-6
    inone["namelis1"]["renein"] = rho
    inone["namelis1"]["njte"] =  [nrho]
    inone["namelis1"]["rtein"] = rho
    inone["namelis1"]["tein"] = te
    inone["namelis1"]["njti"] = [nrho]
    inone["namelis1"]["rtiin"] = rho
    inone["namelis1"]["tiin"] = ti
    inone["namelis1"]["njzef"] =  [nrho]
    inone["namelis1"]["rzeffin"] = rho
    inone["namelis1"]["zeffin"] = zeff
    inone["namelis1"]["rangrot"] = rho
    inone["namelis1"]["angrotin"] = omega

    inone["namelis1"]["njenp"] =len(namep)*[nrho]
    for k in range(len(namep)):
        inone["namelis1"]["renpin(1,%d)"%(k+1)] = rho
        inone["namelis1"]["enp(1,%d)"%(k+1)] = np[k]*1.0e-6

    inone["namelis1"]["njeni"] = len(namei)*[nrho]
    inone["namelis1"]["reniin(1,1)"] = rho
    inone["namelis1"]["reniin(1,2)"] = rho
    inone["namelis1"]["eni(1,1)"] = nc*1.0e-6
    if len(namei) > 1:
        inone["namelis1"]["eni(1,2)"] = nhe4*1.0e-6

    inone["namelis3"]["EQDSKIN"] = [f_eqdsk]

    inone.write("inone")

def write_inputfiles_instate(f_instate, f_eqdsk, f_infreya, dir_data=''):
    instate = Namelist(f_instate)

    rho = array(instate['instate']['rho'])
    nrho = len(rho)

    ne = array(instate['instate']['ne'])
    te = array(instate['instate']['te'])
    ti = array(instate['instate']['ti'])
    zeff = array(instate['instate']['zeff'])
    omega = array(instate['instate']['omega'])

    inone = set_default()

    if dir_data:
       inone["namelis1"]["onetwo_xsct_ext"] = [dir_data]

    innfreya = Namelist(f_infreya)
    for key in innfreya["innfreya"].keys():
        inone["namelis2"][key] = innfreya["innfreya"][key]

    # inone["namelis1"]["inenez"] = [1]
    # inone["namelis1"]["nprim"] = [2]
    # inone["namelis1"]["namep"] = ['d', 't']
    # inone["namelis1"]["zfrac"] = [0.5]
    # inone["namelis1"]["nimp"] = [1]
    # inone["namelis1"]["namei"] = ['c']
    # inone["namelis1"]["nneu"] = [2]
    # inone["namelis1"]["namen"] = ['d', 't']

    inone["namelis1"]["inenez"] = [1]
    inone["namelis1"]["nprim"] = [1]
    inone["namelis1"]["namep"] = ['d']
    inone["namelis1"]["zfrac"] = [1.0]
    inone["namelis1"]["nimp"] = [1]
    inone["namelis1"]["namei"] = ['c']
    inone["namelis1"]["nneu"] = [1] #[2]
    inone["namelis1"]["namen"] = ['d'] #['d', 'c']

    inone["namelis1"]["nj"] = [nrho]
    inone["namelis1"]["njene"] = [nrho]
    inone["namelis1"]["enein"] = ne*1.0e13
    inone["namelis1"]["renein"] = rho
    inone["namelis1"]["njte"] = [nrho]
    inone["namelis1"]["rtein"] = rho
    inone["namelis1"]["tein"] = te
    inone["namelis1"]["njti"] = [nrho]
    inone["namelis1"]["rtiin"] = rho
    inone["namelis1"]["tiin"] = ti
    inone["namelis1"]["njzef"] = [nrho]
    inone["namelis1"]["rzeffin"] = rho
    inone["namelis1"]["zeffin"] = zeff
    inone["namelis1"]["rangrot"] = rho
    inone["namelis1"]["angrotin"] = omega

    # inone["namelis1"]["njenp"] = 2*[nrho]
    # inone["namelis1"]["renpin(1,1)"] = rho
    # inone["namelis1"]["renpin(1,2)"] = rho
    # inone["namelis1"]["enp(1,1)"] = np[0]*1.0e-6
    # inone["namelis1"]["enp(1,2)"] = np[1]*1.0e-6
    # inone["namelis1"]["njeni"] = 2*[nrho]
    # inone["namelis1"]["reniin(1,1)"] = rho
    # inone["namelis1"]["reniin(1,2)"] = rho
    # inone["namelis1"]["eni(1,1)"] = nc*1.0e-6
    # inone["namelis1"]["eni(1,2)"] = nhe4*1.0e-6

    inone["namelis3"]["EQDSKIN"] = [f_eqdsk]

    inone["namelis2"]["extqerf_id"] = ['rf_e']
    inone["namelis2"]["extqerf_watts"] = [0]
    inone["namelis2"]["extqerf"] = [1.0]
    inone["namelis2"]["extqerf_nj"] = [nrho]
    inone["namelis2"]["extqerf_rho"] = rho
    inone["namelis2"]["extqerf_qe"] = instate["instate"]["pe_ic"]

    inone["namelis2"]["extqirf_id"] = ['rf_i']
    inone["namelis2"]["extqirf_watts"] = [0]
    inone["namelis2"]["extqirf"] = [1.0]
    inone["namelis2"]["extqirf_nj"] = [nrho]
    inone["namelis2"]["extqirf_rho"] = rho
    inone["namelis2"]["extqirf_qi"] = instate["instate"]["pi_ic"]

    # inone["namelis2"]["extcurrf_id"] = ['rf_j']
    # inone["namelis2"]["extcurrf_amps"] = [0]
    # inone["namelis2"]["extcurrf"] = [1.0]
    # inone["namelis2"]["extcurrf_nj"] = [nrho]
    # inone["namelis2"]["extcurrf_rho"] = rho
    # inone["namelis2"]["extcurrf_curr"] = instate["instate"]["j_ic"]

    inone.write("inone")

def read_output():
    outone = outone_io.readoutone(foutone='outone', isummary=False, iverb=False)
    pnbe = outone["qbeame"]["y"][-1]
    pnbe = outone["qbeami"]["y"][-1]
    pnbe = outone["enbeam"]["y"][-1]*1.0e-13
    pnbe = outone["curbeam"]["y"][-1]*0.01

def update_state(f_state, f_eqdsk, scales={}):
    #-- read geqdsk and plasma state
    ps = plasmastate('ips',1)
    ps.read(f_state)

    geq = readg(f_eqdsk)
    r0  = geq["rzero" ]
    b0  = abs(geq["bcentr"])
    ip  = geq['cpasma']

    #-- read outone
    outone = outone_io.readoutone(foutone='outone', isummary=False, iverb=False)

    pnbe = outone["qbeame"]["y"][-1]
    pnbi = outone["qbeami"]["y"][-1]
    density_beam = outone["enbeam"]["y"][-1]*1.0e-13
    density_beam = array(density_beam)
    density_beam[:10] = density_beam[10]
    wbeam = outone["wbeam"]["y"][-1]
    wbeam = array(wbeam)
    wbeam[:10] = wbeam [10]
    curbeam = outone["curbeam"]["y"][-1]*0.01
    curbe = outone["curbe"]["y"][-1]*0.01
    curbet = outone["curbet"]["y"][-1]*0.01
    j_nb = array(curbeam)+array(curbe)+array(curbet)
    torque =  outone["storque"]["y"][-1]*0.1
    ssum =  outone["ssum"]["y"][-1]*1.0e6
    tbeami = 2.0/3.0*1.0e3*wbeam/(1.602*density_beam)

    #--- scale
    print('len pnbe', len(pnbe), len(density_beam))
    print('tbeami =', tbeami[0])
    print('density_beam =', density_beam[0])

    if 'current' in scales:
        j_nb = array(j_nb)*scales['current']
    if 'particle' in scales:
        print('particle source scaled ', scales['particle'])
        ssum = array(ssum)*scales['particle']
    if 'torque' in scales:
        j_nb = array(torque)*scales['torque']

    #--- grid
    rho = ps["rho"][:]
    nrho = len(rho)

    #--- beam density, energy
    ps["nbeami"][0] = 1.0e19*ps.node2cell(density_beam)
    ps["eperp_beami"][0] = 2.0*ps.node2cell(tbeami)
    ps["epll_beami"][0] = ps.node2cell(tbeami)

    #--- current
    ps.load_j_parallel(rho, 1.0e6*j_nb, "rho_nbi", "curbeam", r0, b0)

    #--- heating
    ps.load_vol_profile (rho, 1.0e6*pnbe, "rho_nbi", "pbe")
    ps.load_vol_profile (rho, 1.0e6*pnbi, "rho_nbi", "pbi")
    ps.load_vol_profile (rho, torque, "rho_nbi", "tqbi")
    ps.load_vol_profile (rho, ssum, "rho_nbi", "sbedep")

    #--- store plasma state
    ps.store(f_state)

def read_outone_neutron(f):
    pat_1 = re.compile("\s*P DD")
    pat_2 = re.compile("\s*D\(D,n\)")
    pdd, neutron = 0.0, 0.0
    for line in open(f,"r").readlines():
        if pat_1.search(line):
            pdd = float(line.split()[2][1:])
        if pat_2.search(line):
            neutron = float(line.split()[2])
            break
    return pdd, neutron

def update_instate(f_instate,f_eqdsk,scales={}):
    #-- read geqdsk and plasma state
    instate = Namelist(f_instate)

    geq = readg(f_eqdsk)
    r0  = geq["rzero" ]
    b0  = abs(geq["bcentr"])
    ip  = geq['cpasma']

    #-- read outone
    outone = outone_io.readoutone(foutone='outone', isummary=False, iverb=False)
    pnbe = outone["qbeame"]["y"][-1]
    pnbi = outone["qbeami"]["y"][-1]
    density_beam = outone["enbeam"]["y"][-1]*1.0e-13
    density_beam = array(density_beam)
    density_beam[:10] = density_beam[10]
    wbeam = outone["wbeam"]["y"][-1]
    wbeam = array(wbeam)
    wbeam[:10] = wbeam [10]
    curbeam = outone["curbeam"]["y"][-1]*0.01
    curbe = outone["curbe"]["y"][-1]*0.01
    curbet = outone["curbet"]["y"][-1]*0.01
    j_nb = array(curbeam)+array(curbe)+array(curbet)
    torque =  outone["storque"]["y"][-1]*0.1
    ssum =  outone["ssum"]["y"][-1]*1.0e6
    tbeami = 2.0/3.0*1.0e3*wbeam/(1.602*density_beam)

    #--- neutron
    pdd, neutron = read_outone_neutron('outone')

    #--- scale
    print('len pnbe',len(pnbe),len(density_beam))
    print('tbeami =', tbeami[0])
    print('density_beam =', density_beam[0])

    if 'current' in scales:
        j_nb = array(j_nb)*scales['current']
    if 'particle' in scales:
        print('particle source scaled ', scales['particle'])
        ssum = array(ssum)*scales['particle']

    #--- write instate
    rho = instate['instate']['rho']
    nrho = len(rho)

    instate['instate']['wbeam'] = wbeam
    instate['instate']['density_beam'] = density_beam
    instate['instate']['j_nb'] = j_nb

    instate['instate']['pdd'] = [pdd]
    instate['instate']['neutron'] = [neutron]

    instate.write(f_instate)
