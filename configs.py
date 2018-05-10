import plistlib as pll
from datetime import datetime
import importlib as port
utils=port.import_module('utils.utes')
import os
import sys



"""configs

monomers
nucles
rates
steps
runs
filenames
model name
date

"""


if __name__=='__main__':
    fn='input.data'
    runs=10
    counts=300
    conc = 5
    units = {'value':5,
             'unit':'molar',
             'prefix':'micro'
    }
    nc=3
    mons=500
    rates=[1,0.0001,1,0.0001,0.000001]
    model_name='smoluchowski'
    outputs=['mass.data','number.data','polymers.data','skew.data','kurtosis.data','histogram.data','t_steps.data','t.data','state.data']
    date=str(datetime.now().isoformat(timespec='minutes'))
    path='/Users/John1/Development/StateVectors/results/'
    parser={
    'system':{
        'filename':fn,
        'path':path,
        'model':model_name,
        'date':date},
    'simulation':{
        'runs':runs,
        'steps':counts,
        'outputs':outputs},
    'proteins':{
        'nucleus':nc,
        'monomers':mons},
    'rates':rates}

    with utils.config_context(fn,path=path) as f:
        pll.dump(parser,f)



