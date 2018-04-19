import kernels
import KMC
import statevector as sv

x=sv.StateVector([50],delete_zeros=True)
Smoluchowski=KMC.Model(kernels.smol_new,x)
looping=True
while(looping):
    Smoluchowski.ready()
    Smoluchowski.BulkProbability()
    Smoluchowski.choose()
    Smoluchowski.time_step()
    Smoluchowski.advance()
    inp=input("0 to stop: ")
    if inp=="0":
        looping = False