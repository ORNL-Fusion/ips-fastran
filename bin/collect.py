import sys
import os
import shutil
from glob import glob
import tarfile
from numpy import *
from distutils.dir_util import copy_tree
from optparse import OptionParser

def filter_function(tarinfo):
   fname = os.path.abspath(tarinfo.name).split("/")
   if fname[-2] == "plasma_state" and fname[-3] == "simulation_results":
        print(fname)
        return None
   delete_files = ["bpltfil", "eqpltfil", "namelists", "qikone", "isllog", "runlog", "test_trnspt_mhd.txt"]
   if fname[-1] in delete_files: 
        print(fname)
        return None
   return tarinfo

rdir = sys.argv[1]
ishot = int(sys.argv[2])
itime = int(sys.argv[3])

source = os.path.join(rdir, "%06d.%05d"%(ishot, itime)) 
target = os.path.join(rdir, "%06d.%05d.summary"%(ishot, itime)) 

parser = OptionParser()
parser.add_option("-i", action="store_true", dest="iteration", default=False)
parser.add_option("-a", action="store_true", dest="archive", default=False)
parser.add_option("-s", action="store_false", dest="dakota", default=True)
(options,args) = parser.parse_args(sys.argv[1:])

print('options.iteration:', options.iteration)
print('options.archive:', options.archive)

if not os.path.exists(target):
    print('create direcotry: ', target)
    os.makedirs(target)
setup_dir = os.path.join(target, "__setup__") 
if not os.path.exists(setup_dir):
    os.makedirs(setup_dir)
if options.iteration:
    target_sub = os.path.join(target, '__conv__')
    if not os.path.exists(target_sub):
        print('create direcotry __conv__: ')
        os.makedirs(target_sub)

if options.archive:
    print(source)
    print(target)
    tar = tarfile.open(os.path.join(target, "run.tgz"), mode='w:gz')
    tar.add(source, filter=filter_function, arcname="run")
    tar.close()

if options.dakota:
    sim = sorted(glob(os.path.join(source, "SIMULATION_LIST.*")))[0]
    id = int(sim.split(".")[-1])
    
    print(sim, id)
    
    scans = []
    shots = []
    times = []
    with open(sim, "r") as f:
        lines = f.readlines()
        vars = lines[0].split()
        print(vars)
        k_time = vars.index('efit:TIME_ID')
        print('time column:', k_time)
        for line in lines[1:]:
            p = line.split()
            scans.append(p[0])
            shots.append(ishot)
            times.append(int(float(p[k_time])))
    scans=array(scans)
    shots=array(shots)
    times=array(times)
    
    print(scans)
    print('shots=', shots)
    print('times=', times)
    
    scan_list = scans[ argsort(times) ]
    shots = shots[ argsort(times) ]
    times = times[ argsort(times) ]
    
    print(scan_list)
    print(shots)
    print(times)
    print(72*'-')

else:
    scan_list = ["."]
    shots = [ishot]
    times = [itime]
    print(72*'-')

collect_output_list = [ 
    ["fastran_tr_fastran", "fastran.nc", "f", False ],
    ["fastran_tr_fastran", "i??????.?????", "i", True ],
    ["fastran_eq_efit", "g??????.?????", "g", False ],
    ["fastran_eq_efit", "a??????.?????", "a", False ],
    ["fastran_nb_nfreya", "outone", "o", True  ],
    ["fastran_eq_efit", "s??????.?????", "s", False  ],
]

collect_input_list = [
    "fastran_scenario.config",
    "submitjob",
    "fastran.modulefile"
]
if options.dakota:
    collect_input_list += "dakota.in"

for file in collect_input_list:
    shutil.copyfile(os.path.join(source, file), os.path.join(setup_dir, file))
copy_tree(os.path.join(source, "input"), os.path.join(setup_dir, "input"))

def lastiter(files):
    timelist = []
    for file in files:
        timelist.append(int(file.split("/")[-4]))
    imax = max(timelist) 
    return imax

for k, which_scan in enumerate(scan_list):
    print('\n*', k, which_scan)
    for collect in collect_output_list:

        rdir = collect[0]
        var = collect[1]
        prefix = collect[2]
        iconv = collect[3]

        print(rdir)
        if iconv:
            files = sorted(glob(os.path.join(source,which_scan)+"/work/%s_*/%s"%(rdir,var)))
            if len(files)==0: continue
            _file = files[-1]
            imax=4
        else:
            files = sorted(glob(os.path.join(source,which_scan)+"/simulation_results/*/components/%s_*/%s"%(rdir,var)))
            if len(files)==0: continue
            imax = lastiter(files)
            _file = sorted(glob(os.path.join(source,which_scan)+"/simulation_results/%d/components/%s_*/%s"%(imax,rdir,var)))[-1]
            print('iteration at', imax)

        try:
            shutil.copyfile(_file,target+"/%s%06d.%05d"%(prefix,shots[k],times[k]))
        except:
            pass
    
        if options.iteration:
            for file in files:
                try:
                    i =  int(file.split("/")[-4])
                    print(file, i)
                    if i<=100: 
                        shutil.copyfile(file,target_sub+"/%s%06d.%05d_%02d"%(prefix,shots[k],times[k],i))
                except:
                    print(file, 'pass')

