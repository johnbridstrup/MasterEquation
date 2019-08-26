import importlib as port
import utils.utes as utils
import statevector as sv
import KMC
import sys
import os
import numpy as np
from shutil import copyfile
from astropy.table import Table, Column
import copy

k = port.import_module('kernels')


# noinspection SpellCheckingInspection
def simulation(fn):
    newpath = fn

    # path=utils.wd()+'/results/'+fn+'/'
    params = utils.load_config(utils.wd() + '/utils/input.data')
    # os.makedirs(path,exist_ok=True)
    # copyfile(utils.wd()+'/utils/input.data',path+'_input.data')
    rates = params['rates']
    proteins = params['proteins']
    sim_params = params['simulation']
    nc = proteins['nucleus']
    # FINISH output_files=params['simulation']['outputs']
    MM = proteins['monomers']
    M = k.Monomers(MM)
    global x
    x = sv.StateVector(M(), delete_zeros=True)
    global model
    model = KMC.Model(x, nc)
    # t0=sim_params['t0']
    tf = sim_params['tf']
    # bins=sim_params['bins']
    # bins+=1

    add = k.MonAdd(rates['addition'], M=MM, c=proteins['concentration']['value'])
    sub = k.MonSub(rates['subtraction'])
    coag = k.Coag(rates['coagulation'], M=MM, c=proteins['concentration']['value'])
    frag = k.Frag(rates['fragmentation'])
    nuc = k.Nuc(rates['nucleation'], nc=nc, M=MM, c=proteins['concentration']['value'])

    # print(x)
    model.add_propensity(add)
    model.add_propensity(sub)
    model.add_propensity(nuc)
    model.add_propensity(frag)
    model.add_propensity(coag)
    model.add_mechanisms(sv.SmoluchowskiModel(x, nc))
    # print(model.mechanisms)
    X = []
    Y = []
    # last=binstogram[0]
    # model.calculate_probability()
    # model.choose()
    # model.time_step()
    # model.advance()
    # binstogram=np.zeros((bins,MM))
    # tmp=np.zeros(MM)
    # curr=binstogram[0]
    looping = True
    countr = 0
    current_time = 0.0
    # bincounter=0
    # while model.data['t'][-1]>=t[bincounter]:
    #     binstogram[bincounter][0]+=MM
    #     bincounter+=1
    while looping:
        X.append(x)
        # countr+=1
        model.calculate_probability()
        model.choose()
        model.time_step()
        model.advance()
        current_time += model.t_step
        if countr > 300 or current_time > tf:
            looping = False
        X.append(model.data)
    Y.append(X[:])

    # model.data.save(utils.wd()+'/results/'+fn+'/'+fn+str(i)+'.json',model.data_list)

    model.save(newpath)
    x = sv.StateVector([500], delete_zeros=True)
    model = KMC.Model(x, 3)
    model.add_propensity(add)
    model.add_propensity(sub)
    model.add_propensity(nuc)
    model.add_propensity(frag)
    model.add_propensity(coag)
    model.add_mechanisms(sv.SmoluchowskiModel(x, 3))
    return Y


if __name__ == '__main__':
    simulation('.')
