import importlib as port
utils=port.import_module('utils.utes')
import simulations as sim
import numpy as np
import plistlib as pll
import sys
import os
from itertools import product

#! ==============================================================================
# JUST CHANGE THIS FILE NAME NOT THE PATHS
#^ ==============================================================================
#? ==============================================================================
runs_name = "yuank1"

runs_num = 10     
end_time = 2000
start_time = 0.0
mons=200
"""name, repeats

the real end and start times (arbitrary units or based on rate constants)

"""

#? ==============================================================================
#^ ==============================================================================
# UNLESS YOU REALLY DGAF IM NOT YOUR DAD
#! ==============================================================================
#\\TODO: dont write so many times
#! ==============================================================================


#~ ==============================================================================
#& Making directories
#~ ==============================================================================

app_dir=os.getcwd()+'/'
results_dir=str(app_dir)+"results/"
runs_dir=str(results_dir)+runs_name+'.brid/'
utils_dir=str(app_dir)+"utils/"
params=pll.readPlist(utils_dir+'input.data.2')
os.makedirs(runs_dir,exist_ok=True)
os.makedirs(utils_dir,exist_ok=True)
# print(runs_dir)
# print(utils_dir)

rates=params['rates']
params['simulation']['tf']=end_time
params['simulation']['t0']=start_time
rates['subtraction']=10**-5
rates['fragmentation']=rates['subtraction']
#rates['fragmentation']=10**-5
a=(-3,0,1)
asw=[10.0**i for i in [int(j) for j in range(a[0],a[1],a[2])]]
b=(-3,-2,1)
bsw=[10**i for i in [int(j) for j in range(b[0],b[1],b[2])]]
c=(0,7,1)
csw=[10**i for i in [int(j) for j in range(c[0],c[1],c[2])]]
f=(0,7,1)
fsw=[10**i for i in [int(j) for j in range(f[0],f[1],f[2])]]
k=(-8,-4,1)
params['proteins']['monomers']=mons
ksw=[10.0**i for i in [int(j) for j in range(k[0],k[1],k[2])]]
for aa,kk in product(asw,ksw):
    print(aa,kk)
    rates['addition']=aa
    rates['coagulation']=aa
    rates['subtraction']=10.0**-3
    rates['fragmentation']=10.0**-3
    ext="_aa_{}_bb_{}_kk_{}.strup".format(aa,10**-3,kk)
    #ext2="_aa_{}_kk_{}.up".format(aa,kk)
    fn=runs_name+ext
    #fn2=runs_name+ext2
    rates['nucleation']=kk
    params['rates']=rates
    params['simulation']['runs']=runs_num
    # params['simulation'][]
    os.makedirs(runs_dir+fn+"/",exist_ok=True)
    os.makedirs(runs_dir+fn+"/averaged/",exist_ok=True)
    with open(utils_dir+'input.data.2','wb') as f:
        pll.dump(params,f)
    with open(runs_dir+fn+"/input.data.2",'wb') as f:
        pll.dump(params,f)
    for i in range(1,runs_num+1):
        sim.simulation(runs_dir+fn+"/"+fn+"_"+str(i)+".json")


