#! kernels.py
import KMC
import statevector as sv

class MergeKernel:
    def __init__(self,frequency,**kwargs):
        self.frequency=frequency
        try:
            self.fractal_dimension=kwargs["fractal_dimension"]
        except:
            self.fractal_dimension=1.0/3.0
    def __call__(self,r,s):
        return 1.0*self.frequency*(r**(self.fractal_dimension)+s**(self.fractal_dimension))

    
@KMC.Propensity
def Coagulation(s,*args,**kwargs):
    if kwargs["polymer_number"]>1:
        return 1.0*kwargs["frequency"](s[0],s[1])
    else:
        return 0.0
m_args={"reactant_range":True,"degree":2,"frequency":MergeKernel(0.001)}
@KMC.Propensity
def Fragmentation(s,*args,**kwargs):
    print("FAGGGG KERNEL     ",s)
    if (s-3)>0 and kwargs["polymer_number"]>0:
        try:
            print("THE FRAG PROB", 1.0*kwargs["frequency"]*(s-3))
            return 1.0*kwargs["frequency"]*(s-3)
        except Exception as e1:
            try:
                return 1.0*kwargs["frequency"](s)
            except Exception as e2:
                print("break kernel must be function or constant")
                raise e2 
        except:
            raise e1
    else:
        print("00000000000")
        return 0.0
    
frag_args={"range_start":1,"Range_stop":True,"frequency":1}
@KMC.Propensity
def MonomerAddition(s,*args,**kwargs):
    print(kwargs["polymer_number"])
    if args[0]>0 and kwargs["polymer_number"]>0:
        try:
            return 1.0*args[0]*kwargs["frequency"]*s
        except:
            try:
                return 1.0*kwargs["frequency"](s,args[0])
            except:
                pass
    else:
        return 0.0
mono_args={"Range":True,"frequency":10}
@KMC.Propensity
def Nucleation(s,*args,**kwargs):
    print("THE MONOMER COUNT",s)
    if s>kwargs["nucleus"]:
        try:
            results=kwargs["frequency"]
            for i in range(kwargs["nucleus"]):
                results*=s-i
            return 1.0*kwargs["frequency"]*results
        except:
            try:
                return 1.0*kwargs["frequency"](s,kwargs["nucleus"])
            except:
                raise
    else:
        return 0.0
nuc_args={"nucleus":3,"frequency":0.00001}
@KMC.Propensity
def MonomerSubtraction(s,*args,**kwargs):
    if s>1:
        try:
            return 2.0*kwargs["frequency"]
        except:
            try:
                return 1.0*kwargs["frequency"](s)
            except:
                return 0.0
    else:
        return 0.0
sub_args={"range_start":1,"Range_stop":True,"frequency":.0001}
smol_args=(mono_args,sub_args,nuc_args,frag_args,m_args)
smol_prps=[MonomerAddition,MonomerSubtraction,Nucleation,Fragmentation,Coagulation]
smol_new=[(MonomerAddition,mono_args),(MonomerSubtraction,sub_args),(Nucleation,nuc_args),(Fragmentation,frag_args),(Coagulation,m_args)]

        