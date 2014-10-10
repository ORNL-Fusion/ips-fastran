#!/task/imd/local64/bin/python

import sys,os

try:
    config = sys.argv[1]
except:
    err = 'no input configuration file'
    raise Exception(err)

host = os.environ["HOST"].split('.')[0]

if host in [ 'node%02d'%k for k in range(1,7) ]: 
   host = 'linux'
elif host in [ 'venusa' ]:
   host = 'venus' 

if host == 'venus':
    job = """#!/bin/csh
#$ -q all.q
#$ -o fastran_ips.log
#$ -e fastran_ips.err
#$ -cwd
#$ -V
#$ -l mem_free=4G,h_vmem=16G
source $FASTRAN_ROOT/setup/fastran.cshrc.venus
/task/imd/local64/bin/python $IPS_ROOT/bin/ips --config=%s --platform=$FASTRAN_ROOT/setup/venus.conf
"""%config
    f=open("qsub.job","w")
    f.write(job)
    f.close()
    os.system("qsub qsub.job")

else:
    cmd = "source $FASTRAN_ROOT/setup/fastran.bashrc.venus;"
    cmd+= "/task/imd/local64/bin/python $IPS_ROOT/bin/ips --config=%s --platform=$FASTRAN_ROOT/setup/%s.conf"%(config,host)
    cmd+= ">&fastran_ips.log&"
    print 'submit job on '+host
    print cmd
    os.system(cmd)

