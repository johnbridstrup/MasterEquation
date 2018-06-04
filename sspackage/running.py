import importlib as port
utils=port.import_module('utils.utes')
import simulations as sim
import numpy as np
import plistlib as pll
import sys
import os

#! ==============================================================================
# JUST CHANGE THIS FILE NAME NOT THE PATHS
#^ ==============================================================================
#? ==============================================================================
runs_name = "middle_sweep_1"
runs_num = 200
end_time = 500.0
start_time = 0.0
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
params=pll.readPlist(utils_dir+'input.data')
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
a=(-2,3,6)
asw=np.logspace(a[0],a[1],a[2])
k=(-5,1,6)
ksw=np.logspace(k[0],k[1],k[2])
for aa in asw:
    rates['addition']=aa
    rates['coagulation']=aa
    for kk in ksw:
        ext="_aa_{}_kk_{}.strup".format(aa,kk)
        #ext2="_aa_{}_kk_{}.up".format(aa,kk)
        fn=runs_name+ext
        #fn2=runs_name+ext2
        rates['nucleation']=kk
        params['rates']=rates
        params['simulation']['runs']=runs_num
        # params['simulation'][]
        os.makedirs(runs_dir+fn+"/",exist_ok=True)
        os.makedirs(runs_dir+fn+"/averaged/",exist_ok=True)
        with open(utils_dir+'input.data','wb') as f:
            pll.dump(params,f)
        with open(runs_dir+fn+"/input.data",'wb') as f:
            pll.dump(params,f)
        for i in range(1,runs_num+1):
            sim.simulation(runs_dir+fn+"/"+fn+"_"+str(i)+".json")


