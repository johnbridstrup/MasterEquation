import importlib as port
import statevector as sv
utils=port.import_module('utils.utes')
k=port.import_module('kernels2')
import KMC
import sys
import os
import plistlib as pll



class my_model:
    def __init__(self,state,*args,**kwargs):
        self.params=None
        self.load_config()
        rates=self.params['rates']
        proteins=self.params['proteins']
        nc=proteins['nucleus']
        M=k.Monomers(proteins['monomers'])
        self.output_files=self.params['simulation']['outputs']
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
        

    def load_config(self):
        self.params=load_config()
    
def load_config():
    params=utils.load_config(utils.wd()+'/utils/input.data')
    return params

class StateSimulation:
    def __init__(self):
        self.model=my_model
    def Simulate(self):
        for i in range(int(sys.argv[2])):
            looping=True
            countr=0
            while(looping):
                countr+=1
                self.model.calculate_probability()
                self.model.choose()
                self.model.time_step()
                self.model.advance()
                # print(x)
                #inp=input("0 to quit: ")
                # if inp=="0":
                #     looping=False
                # if inp=="1":
                #     print(self.model.data)
                if countr>300:
                    looping = False


            fn=sys.argv[1]
            os.makedirs(utils.wd()+'/results/'+fn, mode=0o777, exist_ok=True)
            # self.model.data.save(utils.wd()+'/results/'+fn+'/'+fn+str(i)+'.json',self.model.data_list)
            self.model.save(utils.wd()+'/results/'+fn+'/'+fn+str(i)+'.json')
            x=sv.StateVector([500],delete_zeros=True)
            self.model=KMC.Model(x,3)
            self.model.add_propensity(add)
            self.model.add_propensity(sub)
            self.model.add_propensity(nuc)
            self.model.add_propensity(frag)
            self.model.add_propensity(coag)
            self.model.add_mechanisms(sv.SmoluchowskiModel(x,3))

if __name__ == '__main__':
    pass