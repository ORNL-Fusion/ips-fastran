'''
 -----------------------------------------------------------------------
 freegs io
 -----------------------------------------------------------------------
'''
import numpy as np
from Namelist import Namelist
import freegs
import freegs.geqdsk
import freegs.machine
from freegs.shaped_coil import ShapedCoil
import json
from fastran.util.zinterp import zinterp
from fastran.state.instate import Instate
from fastran.equilibrium.mhdeq_io import mhdeq_io


def read_coil_data(f_coil_data, scale=1.0):
    with open(f_coil_data, 'r') as f:
        coil_data = json.load(f)
    coils = []
    for key in coil_data:
        print(key, coil_data[key])
        shape = []
        for p in coil_data[key]:
            shape.append([p[0] * scale, p[1] * scale])
        coils.append([key, ShapedCoil(shape)])
    return coils


class freegs_io():
    def __init__(self, topology='DN', model_pressure='total'):
        self.topology = topology
        self.mhdeq = mhdeq_io(model_pressure)

    def from_instate(self, f_instate):
        self.instate = Instate(f_instate)

    def define_tokamak(self, f_coil_data, domain, nx=129, ny=129):
        # read coils
        coils_shape = read_coil_data(f_coil_data)
        rwall = self.instate['rlim']
        zwall = self.instate['zlim']

        # define tokamak from coil and wall
        self.tokamak = freegs.machine.Machine(
                           coils_shape, 
                           freegs.machine.Wall(rwall, zwall))

        # define freegs equilibriumm
        self.eq = freegs.Equilibrium(
                      tokamak=self.tokamak,
                      Rmin=domain['rmin'],
                      Rmax=domain['rmax'],
                      Zmin=domain['zmin'],
                      Zmax=domain['zmax'],
                      nx=nx,
                      ny=ny,
                      boundary=freegs.boundary.freeBoundaryHagenow)

    def set_profiles_simple(self):
        self.mhdeq.from_instate(self.instate)
        press = self.mhdeq.get_p() 
        r0 = self.instate['r0'][0]
        b0 = self.instate['b0'][0]
        ip = self.instate['ip'][0]
        self.profiles = freegs.jtor.ConstrainPaxisIp(
                            self.eq, 
                            press[0],  # Plasma pressure on axis [Pascals]
                            ip * 1.e6 ,  # Plasma current [Amps]
                            r0 * b0)  # Vacuum f=R*Bt

    def set_profiles(self, f_geqdsk, relax=0, relax_metric=0, relax_profile=0):
        self.mhdeq.from_instate(self.instate)
        self.mhdeq.load_geqdsk(f_geqdsk)
        self.mhdeq.get_ffprime_pprime_relax(relax)
        # press = mhdeq.get_p() 
        # print(press) 

        psin = self.mhdeq.psin 
        pprim = self.mhdeq.pprim
        ffprim = self.mhdeq.ffprim
        p = self.mhdeq.p
        f = self.mhdeq.f
        # print(psin, pprim, ffprim, p, f)

        pprim = zinterp(psin, pprim)
        ffprim = zinterp(psin, ffprim)
        pressure = zinterp(psin, p)
        fpol = zinterp(psin, f)

        r0 = self.instate['r0'][0]
        b0 = self.instate['b0'][0]

        def pprime_func(psin):
            return pprim(psin)

        def ffprime_func(psin):
           return ffprim(psin)

        def p_func(psin):
           return pressure(psin)

        def f_func(psin):
           return fpol(psin)

        self.profiles = freegs.jtor.ProfilesPprimeFfprime(
            pprime_func,
            ffprime_func,
            r0*b0,
            p_func=p_func,
            f_func=f_func)

    def solve(self, p0=[0.1, 0.5, 1.0], verbose=True, method=0):
        # find control points
        rbdry = self.instate['rbdry']
        zbdry = self.instate['zbdry']

        k_x1 = np.argmin(zbdry)
        k_x2 = np.argmax(zbdry)
        k_in = np.argmin(rbdry)
        k_out = np.argmax(rbdry)

        Rx_1, Zx_1 = rbdry[k_x1], zbdry[k_x1]
        Rx_2, Zx_2 = rbdry[k_x2], zbdry[k_x2]
        R_in, Z_in = rbdry[k_in], zbdry[k_in]
        R_out, Z_out = rbdry[k_out], zbdry[k_out]
        a0 = 0.5 * (R_out - R_in)

        if verbose:
            print('Rx_1 =', Rx_1)
            print('Rx_2 =', Rx_2)
            print('Zx_1 =', Zx_1)
            print('Zx_2 =', Zx_2)
            print('R_in, R_out =', R_in, R_out)
            print('a0 =', a0)

        # x-points 
        if self.topology == 'DN':
            xpoints = [(Rx_1, Zx_1), (Rx_2, Zx_2)]
        elif self.topology == 'LSN':
            xpoints = [(Rx_1, Zx_1)]
        elif self.topology == 'USN':
            xpoints = [(Rx_2, Zx_2)]

        # isoflux
        isoflux = [(R_out, 0., Rx_1, Zx_1), # Outboard midplane, lower X-point
                   (R_out, 0., Rx_2, Zx_2), # Outboard midplane, upper X-point
                   (R_in, 0., R_out, 0.)  # R_in and R_out
                  ]

        # profiles
        # inmhd = mhdeq_io(instate=self.instate)
        # press = inmhd.get_p() # self.instate['p_eq']

        # constraint
        constrain = freegs.control.constrain(xpoints=xpoints, isoflux=isoflux, gamma=1.e-12)
        constrain(self.eq)

        # solve 
        r0 = self.instate['r0'][0]
        b0 = self.instate['b0'][0]
        ip = self.instate['ip'][0]
        print('r0, b0, ip = ', r0, b0, ip)

        # solve
        if method == 0:
            paxis = self.profiles.paxis
            for k in range(len(p0)):
                print(f'freegs iteration {k}')
                profiles = freegs.jtor.ConstrainPaxisIp(
                    self.eq, 
                    p0[k] * paxis,  # Plasma pressure on axis [Pascals]
                    ip * 1.e6 ,  # Plasma current [Amps]
                    r0 * b0)  # Vacuum f=R*Bt

                freegs.solve(
                    self.eq,  # The equilibrium to adjust
                    profiles,  # The toroidal current profile function
                    constrain)  # Constraint function to set coil currents
        else:
            freegs.solve(
                self.eq,  # The equilibrium to adjust
                self.profiles,  # The toroidal current profile function
                constrain)  # Constraint function to set coil currents

        # eq now contains the solution
        print('Done!')
        print('Plasma current: %e Amps' % (self.eq.plasmaCurrent()))
        print('Plasma pressure on axis: %e Pascals' % (self.eq.pressure(0.)))
        print('Poloidal beta: %e' % (self.eq.poloidalBeta()))

    def write_geqdsk(self, f_geqdsk):
        r0 = self.instate['r0'][0]
        with open(f_geqdsk, 'w') as f:
            freegs.geqdsk.write(self.eq, f, R0=r0)

