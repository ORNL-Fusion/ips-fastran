#!/usr/bin/env Python

import os

#######################################################################
# common data
#

verid = "*** xfastran.py ver1.0/Aug2013\n"
dir_fastran  = os.environ["FASTRAN_ROOT"]
dir_template = os.path.join(dir_fastran,'template')

netcdf_form = 'NETCDF3_CLASSIC'
