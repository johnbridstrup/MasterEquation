import statevector as sv
import kernels2 as k
import KMC
import sys
import os
from shutil import copyfile
import importlib as port
utils=port.import_module('utils.utes')

def simulation(f,runs): #f=sys.argv[1], runs=int(sys.argv[2])
    fn=f
    path=utils.wd()+'/results/'+fn+'/'
    params=utils.load_config(utils.wd()+'/utils/input.data')
    os.makedirs(path,exist_ok=True)
    copyfile(utils.wd()+'/utils/input.data',path+'input.data')
    rates=params['rates']
    proteins=params['proteins']
    nc=proteins['nucleus']
    #FINISH output_files=params['simulation']['outputs']
    MM=proteins['monomers']
    M=k.Monomers(MM)
    x=sv.StateVector(M(),delete_zeros=True)
    model=KMC.Model(x,nc)
    
    add=k.MonomerAddition(rates[0],nc)
    sub=k.MonomerSubtraction(rates[1],nc)
    nuc=k.Nucleation(rates[4],nc)
    frag=k.Fragmentation(rates[4],nc)
    coag=k.Coagulation(rates[3],nc)
    # print(x)
    model.add_propensity(add)
    model.add_propensity(sub)
    model.add_propensity(nuc)
    model.add_propensity(frag)
    model.add_propensity(coag)
    model.add_mechanisms(sv.SmoluchowskiModel(x,nc))
    # print(model.mechanisms)


    for i in range(runs):
        looping=True
        countr=0
        while(looping):
            countr+=1
            model.calculate_probability()
            model.choose()
            model.time_step()
            model.advance()
            # print(x)
            #inp=input("0 to quit: ")
            # if inp=="0":
            #     looping=False
            # if inp=="1":
            #     print(model.data)
            if countr>300:
                looping = False


        
        
        # model.data.save(utils.wd()+'/results/'+fn+'/'+fn+str(i)+'.json',model.data_list)
        model.save(path+fn+str(i)+'.json')
        x=sv.StateVector([500],delete_zeros=True)
        model=KMC.Model(x,3)
        model.add_propensity(add)
        model.add_propensity(sub)
        model.add_propensity(nuc)
        model.add_propensity(frag)
        model.add_propensity(coag)
        model.add_mechanisms(sv.SmoluchowskiModel(x,3))

if __name__=='__main__':
    simulation(sys.argv[1],int(sys.argv[2]))

