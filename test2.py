import statevector as sv
import kernels2 as k
import KMC
import sys
import os
import importlib as port
utils=port.import_module('utils.utes')
from astropy.table import Table, Column

M=k.Monomers(400)
x=sv.StateVector(M(),delete_zeros=True)
model=KMC.Model(x,3)

add=k.MonomerAddition(0.1,3)
sub=k.MonomerSubtraction(0.0001,3)
nuc=k.Nucleation(.001,3)
frag=k.Fragmentation(.1,3)
coag=k.Coagulation(1,3)
# print(x)
model.add_propensity(add)
model.add_propensity(sub)
model.add_propensity(nuc)
model.add_propensity(frag)
model.add_propensity(coag)
model.add_mechanisms(sv.SmoluchowskiModel(x,3))
# print(model.mechanisms)


for i in range(int(sys.argv[2])):
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


    fn=sys.argv[1]
    os.makedirs(utils.wd()+'/results/'+fn,exist_ok=True)
    # model.data.save(utils.wd()+'/results/'+fn+'/'+fn+str(i)+'.json',model.data_list)
    model.save(utils.wd()+'/results/'+fn+'/'+fn+str(i)+'.json')
    x=sv.StateVector([500],delete_zeros=True)
    model=KMC.Model(x,3)
    model.add_propensity(add)
    model.add_propensity(sub)
    model.add_propensity(nuc)
    model.add_propensity(frag)
    model.add_propensity(coag)
    model.add_mechanisms(sv.SmoluchowskiModel(x,3))



