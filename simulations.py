import importlib as port
import statevector as sv
utils=port.import_module('utils.utes')
k=port.import_module('kernels2')
import KMC
import sys
import copy
from shutil import copyfile
import os
import plistlib as pll
from abc import ABC,abstractmethod, ABCMeta, abstractproperty, abstractproperty

class ModelInterface(ABC):
    def __init__(self):
        self.mechanisms=None
    @abstractmethod
    def model_mechanisms(self):
        print(self.mechanisms)
class my_model:
    def __init__(self,*args,**kwargs):
        self.params=None
        self.load_config()
        print("my model prms",self.params)
        self.rates=self.params['rates']
        print("my model rts",self.rates)
        self.proteins=self.params['proteins']
        self.nc=self.proteins['nucleus']
        self.output_files=self.params['simulation']['outputs']
        self.M=self.proteins['monomers']
        self.state=sv.StateVector(self.M,delete_zeros=True)
        self.model=KMC.Model(self.state,self.nc)

        add=k.MonomerAddition(self.rates[0],self.nc)
        sub=k.MonomerSubtraction(self.rates[1],self.nc)
        nuc=k.Nucleation(self.rates[4],self.nc)
        frag=k.Fragmentation(self.rates[4],self.nc)
        coag=k.Coagulation(self.rates[3],self.nc)
        # print(self.state)
        self.model.add_propensity(add)
        self.model.add_propensity(sub)
        self.model.add_propensity(nuc)
        self.model.add_propensity(frag)
        self.model.add_propensity(coag)
        self.model.add_mechanisms(sv.SmoluchowskiModel(self.state,self.nc))
        print("my model model",self.model)
    def getModel(self):
        return self.model    

    def load_config(self):
        self.params=load_config()
        
    
def load_config():
    params=utils.load_config(utils.wd()+'/utils/input.data')
    return params

class StateSimulation:
    def __init__(self):
        the_model=my_model()
        self.model=the_model.getModel()
        print("SS model",self.model)
        self.fn=sys.argv[1]
        self.path=utils.wd()+'/results/'+self.fn+'/'
    def simulate(self):
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

            path=copy.copy(self.path)
            fn=copy.copy(self.fn)
            
            os.makedirs(path, exist_ok=True)
            self.model.save(path+fn+str(i)+'.json')
            self.model=my_model().getModel()
    def __call__(self):
        self.simulate()

if __name__ == '__main__':
    sim=StateSimulation()
    #sim()
# x=sv.StateVector([500],delete_zeros=True)
# self.model=KMC.Model(x,3)
# self.model.add_propensity(add)
# self.model.add_propensity(sub)
# self.model.add_propensity(nuc)
# self.model.add_propensity(frag)
# self.model.add_propensity(coag)
# self.model.add_mechanisms(sv.SmoluchowskiModel(x,3))