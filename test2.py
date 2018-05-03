import statevector as sv
import kernels2 as k
import KMC
from astropy.table import Table, Column

x=sv.StateVector([500],delete_zeros=True)
add=k.MonomerAddition(0.1,3)
sub=k.MonomerSubtraction(0.0001,3)
nuc=k.Nucleation(.001,3)
frag=k.Fragmentation(.1,3)
coag=k.Coagulation(1,3)
print(x)
model=KMC.Model(x,3)
model.add_propensity(add)
model.add_propensity(sub)
model.add_propensity(nuc)
model.add_propensity(frag)
model.add_propensity(coag)
model.add_mechanisms(sv.SmoluchowskiModel(x,3))
print(model.mechanisms)

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
    
model.data.save("run1.json")




