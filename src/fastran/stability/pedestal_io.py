"""
 -----------------------------------------------------------------------
 pedetal io
 -----------------------------------------------------------------------
"""

import numpy as np
from Namelist import Namelist
from fastran.util.modelprofile import profile_pedestal
from fastran.util.zinterp import zinterp
from fastran.util.loglinear import loglinear
from fastran.util.formula import get_ni

def patch_pedestal(rho, yold, rhob, wped, yped, ysep, alpha=1.1, beta=1.1):
    nrho = len(rho)
    ynew = np.zeros(nrho)

    rhotop = 1.0-1.5*wped
    y1 = profile_pedestal(nrho, 1.0-0.5*wped, wped, yped, ysep, yold[0], alpha, beta)
    y0 = zinterp(rho/rhob, yold)

    for k in range(nrho):
        if (rho[k]<=rhotop):
           ynew[k] = y0( rho[k]/rhotop ) + y1(rhotop) - y0(1.0)
        else:
           ynew[k] = y1(rho[k])

    return ynew

def update_instate_pedestal(fn_instate, fn_inpedestal):
    instate = Namelist(fn_instate)
    ip = instate["instate"]["ip"][0]
    neped = instate["instate"]["neped"][0]
    nesep = instate["instate"]["nesep"][0]
    zeffped = instate["instate"]["zeff0"][0]

    fits = {}
    fits["pped"] = loglinear(fn_inpedestal, "pped")
    fits["wped"] = loglinear(fn_inpedestal, "wped")

    wped = fits["wped"].get(ip=ip, neped=neped, nesep=nesep)
    pped = fits["pped"].get(ip=ip, neped=neped, nesep=nesep)

    niped, nzped = get_ni(ne=neped, zeff=zeffped)
    teped= 0.5*pped/(1.602*neped)
    tiped= 0.5*pped/(1.602*niped)

    instate["instate"]["xmid"] = [1.0-0.5*wped]
    instate["instate"]["xwid"] = [wped]

    instate["instate"]["teped"] = [teped]
    instate["instate"]["tiped"] = [tiped]

    instate.write(fn_instate)
