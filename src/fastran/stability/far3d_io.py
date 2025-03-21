"""
-----------------------------------------------------------------------
far3D IO
-----------------------------------------------------------------------
"""
import numpy as np
from Namelist import Namelist
from fastran.plasmastate.plasmastate import plasmastate
from fastran.state.instate import Instate


class far3d_io_profile():
    def __init__(self):
        self.key_map = {
            'rho'    : [0, 'Rho', '(norml_sqrt_toroid_flux)'],
            'q'      : [1, 'q', ''],
            'nbeam'  : [2, 'BeamIonDensity', '(10^14cm^-3)'],
            'nion'   : [3, 'IonDensity', '(10^14cm^-3)'],
            'ne'     : [4, 'ElecDensity', '(10^14cm^-3)'],
            'nalpha' : [5, 'AlphaDensity', '(10^14cm^-3)'],
            'nz'     : [6, 'ImpurityDensity', '(10^14cm^-3)'],
            'tbeam'  : [7, 'BeamIonEffectiveTemp', '(keV)'],
            'ti'     : [8, 'IonTemp', '(keV)'],
            'te'     : [9, 'ElectronTemp', '(keV)'],
            'talpha' : [10, 'Effective_Alpha_Temp', '(keV)'],
            'pbeam'  : [11, 'BeamPressure', '(kPa)'],
            'pth'    : [12, 'ThermalPressure', '(kPa)'],
            'pmhd'   : [13, 'EquilPressure', '(kPa)'],
            'omega'  : [14, 'TorRot', '(kHz)'],
            'vpol'   : [15, 'PolRot', '(10^5m/s)']
            }
        self.data = {}  
        pass

    def __getitem__(self, key):
        return self.data[key]

    def __setitem__(self, key, value):
        self.data[key] = value

    def read_profile(self, f_profile):
        with open(f_profile, 'r') as f:
           lines  =  f.readlines()
           self.b0 = float(lines[2])
           self.r0 = float(lines[4])
           self.a = float(lines[6])
           self.kappa = float(lines[8])
           self.delta = float(lines[10])
           self.imp = lines[12]
           self.a_main = float(lines[14])
           self.beta0 = float(lines[15][25:35])
           self.raxis = float(lines[15][42:])
           prof = np.loadtxt(lines[18:]).transpose()
           self.data = {}
           for key in self.key_map:
               self.data[key] = prof[self.key_map[key][0]]

    def from_state(self, f_instate, f_state):
        # todos: need to add He ash

        instate = Namelist(f_instate)
        self.b0 = instate['instate']['b0'][0]
        self.r0 = instate['instate']['r0'][0] 
        self.a0 = instate['instate']['a0'][0]
        self.delta = instate['instate']['delta'][0]
        self.kappa = instate['instate']['kappa'][0]
     
        ps = plasmastate('ips', 1)
        ps.read(f_state)
        spec = ps.get_species()

        nrho = len(ps['rho'])
        rho = ps['rho'][:]
        ne = ps['ns'][0, :]
        te = ps['Ts'][0, :]
        ti = ps['Ti'][:]
        zeff = ps['Zeff'][:]
        omega = ps['omegat'][:]
        ne = ps.cell2node_bdry(ne)
        te = ps.cell2node_bdry(te)
        ti = ps.cell2node_bdry(ti)
        zeff = ps.cell2node_bdry(zeff)
        omega = ps.cell2node_bdry(omega)
        ni = spec['np']
        nz = spec['nz']

        f_ion = ni / np.sum(ni, axis=0)
        a_ion = np.sum( spec['a_ion'] * np.average(f_ion, axis=1))
        self.a_ion = a_ion
        print('a_ion = ', a_ion)
       
        f_imp = nz / np.sum(nz, axis=0)
        a_imp = np.sum( spec['a_imp'] * np.average(f_imp, axis=1))
        self.a_imp = a_imp
        print('a_imp = ', a_imp)

        density_beam = ps.dump_profile(rho, 'rho_nbi', 'nbeami', k=0) 
        tbeam = ps.dump_profile(rho, 'rho_nbi', 'eperp_beami', k=0) +  ps.dump_profile(rho, 'rho_nbi', 'epll_beami', k=0) # keV
        wbeam = density_beam * tbeam * 1.602e-19 * 1.e3 # J/m**3
        print(tbeam)
        print(wbeam)

        density_alpha = ps.dump_profile(rho, 'rho_fus', 'nfusi', k=0) 
        talpha = ps.dump_profile(rho, 'rho_fus', 'eperp_fusi', k=0) +  ps.dump_profile(rho, 'rho_fus', 'epll_fusi', k=0) # keV
        walpha = density_alpha * talpha * 1.602e-19 * 1.e3 # J/m**3
        print(talpha)
        print(walpha)

        # pbeam = 2. / 3. * (wbeam + walpha)
        pbeam = 2. / 3. * wbeam # check this
        nth = np.sum(ni, axis=0) + np.sum(nz, axis=0)
        pth = 1.602e-19 * 1.e3 * (ne * te + nth * ti) 

        self['rho'] = rho
        self['q'] = ps['q_eq'][:]  
        self['nbeam'] = density_beam * 1.e-20
        self['nion'] = np.sum(ni, axis=0) * 1.e-20
        self['ne'] = ne  * 1.e-20
        self['nalpha'] = density_alpha * 1.e-20
        self['nz'] = np.sum(nz, axis=0) * 1.e-20    
        self['tbeam'] = tbeam  
        self['ti'] = ti 
        self['te'] = te     
        self['talpha'] = talpha 
        self['pbeam'] = pbeam * 1.e-3  
        self['pth'] = pth * 1.e-3  
        self['pmhd'] = ps['P_eq'] * 1.e-3  
        self['omega'] = omega / (2. * np.pi) * 1.e-3
        self['vpol'] = np.zeros(nrho)   

        print('done')

    def write_profile(self, f_profile):
        bt = self.b0
        r0 = self.r0
        a = self.a0
        delta = self.delta
        kappa = self.kappa
        a_ion = self.a_ion
        a_imp = 12 # hard coded, not used in Far3D
        name_imp = 'C' # hard coded, not used in Far3D

        header = 'PLASMA GEOMETRY\n'
        header += f'Vacuum Toroidal magnetic field at R={r0:7.5f}m [Tesla]\n'
        header += f' {bt:7.5f}\n'
        header += 'Geometric Center Major radius [m]\n'
        header += f' {r0:7.5f}\n'
        header += 'Minor radius [m]\n'
        header += f' {a:7.5f}\n'
        header += 'Avg. Elongation\n'
        header += f' {kappa:7.5f}\n'
        header += 'Avg. Top/Bottom Triangularity\n'
        header += f' {delta:7.5f}\n'
        header += 'Main Contaminant Species\n'
        header += f'      {a_imp}{name_imp}\n'
        header += 'Main Ion Species mass/proton mass\n'
        header += f'{a_ion}\n'
        header += f'TRYING TO GET TO BETA(0)=0.0296000 , Rmax=6.0440\n'
        header += '0\n'
        header += ''.join( [f'{self.key_map[key][1]}{ self.key_map[key][2]}, ' for key in self.key_map] )
        header += '\n'
        print(header) 
        print(self['rho'])
        with open(f_profile, 'w') as f:
            f.write(header)
            for k in range(1, len(self['rho'])):
                for key in self.key_map:
                    f.write(f'{self[key][k]:15.5f}')
                f.write('\n')

