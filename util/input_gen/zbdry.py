from numpy import *
from efit_eqdsk import readg
from Namelist import Namelist

f_eqdsk = sys.argv[1]
geq = readg(f_eqdsk)
rbdry = array(geq["rbdry"])[::2]
zbdry = array(geq["zbdry"])[::2]
nbdry = len(rbdry)

instate = Namelist()
instate["instate"]["nbdry"]  = [nbdry]
instate["instate"]["rbdry"]  = rbdry
instate["instate"]["zbdry"]  = zbdry
instate["instate"]["nlim" ]  = [geq["nlim"]]
instate["instate"]["rlim" ]  = geq["rlim"]
instate["instate"]["zlim" ]  = geq["zlim"]
instate.write("instate")
