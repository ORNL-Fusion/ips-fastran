import os
backend =  os.environ['PS_BACKEND']
if os.environ['PS_BACKEND'] == 'pyps':
    from fastran.plasmastate.plasmastate_pyps import plasmastate
elif os.environ['PS_BACKEND'] == 'pstool':
    from fastran.plasmastate.plasmastate_pstool import plasmastate
