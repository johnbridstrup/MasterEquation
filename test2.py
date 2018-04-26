import statevector as sv
import kernels2 as k
import KMC
from astropy.table import Table, Column
M=500
c=5
nc=3
x=sv.StateVector([M],delete_zeros=True)
add=k.MonomerAddition(10,nc, M=M, c=c)
sub=k.MonomerSubtraction(0.0001,nc)
nuc=k.Nucleation(.01,nc, M=M, c=c)
frag=k.Fragmentation(.00001,nc)
coag=k.Coagulation(.001,nc,M=M,c=c)
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
    model.advance()
    # print(x)
    #inp=input("0 to quit: ")
    # if inp=="0":
    #     looping=False
    # if inp=="1":
    #     print(model.data)
    if countr>400:
        looping = False
    
model.data.save("results/run500x400st_highk.json")




