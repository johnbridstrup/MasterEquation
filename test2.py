import statevector as sv
import kernels2 as k
import KMC
x=sv.StateVector([500],delete_zeros=True)
add=k.MonomerAddition(0.1,3)
sub=k.MonomerSubtraction(0.0001,3)
nuc=k.Nucleation(100,3)
frag=k.Fragmentation(.1,3)
coag=k.Coagulation(1,3)
print(x)
model=KMC.Model(x,3)
model.add_propensity(add)
model.add_propensity(sub)
model.add_propensity(nuc)
model.add_propensity(frag)
model.add_propensity(coag)


looping=True
while(looping):
    model.calculate_probability()
    model.choose()
    model.time_step()
    model.advance()
    print(x)
    inp=input("0 to quit: ")
    if inp=="0":
        looping=False




