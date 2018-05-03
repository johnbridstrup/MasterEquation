import importlib as port
import statevector as sv
import KMC
import numpy as np
import numpy.random as npr
import Helpers as halp

x=sv.StateVector([3,11,15,19,22])
Pa=KMC.MonomerAddition(.01)
Pn=KMC.Nucleation(3,.001)
Ps=KMC.MonomerSubtraction(13)
Pf=KMC.Fragmentation(3,28)
Pc=KMC.Coagulation(200)
print(x,"ZERO")
x+=5
x+=5
x+=2

print(x,"THREE")
print(x[2])
Pa(x)
Pn(x)
Ps(x)
Pf(x)
Pc(x)
pa=Pa.probability_vector
pn=Pn.probability_vector
ps=Ps.probability_vector
pf=Pf.probability_vector
pc=Pc.probability_vector
print(pa, "add")
print(pn, "nucleate")
print(ps, "sub")
print(pf, "Fragmenttion")
print(pc,"coag")
print(x.kurtosis(1),"kurtosis")
# class wrapper:
#     def __init__(self,func):
#         self.func=func
#     def __call__(self,*args,**kwargs):
#         num,*args=args
#         return [self.func(i,*args,**kwargs) for i in range(num)]
@KMC.Propensity 
def sq(s,*args,**kwargs):
    return np.sin(np.pi*s*(args[0]+1)/kwargs["length"])

y=sq(x,length=3,range=len(x.state))
sq_kargs={"length":3,"Range":True}
def size_dep(r,s):
    return .1*(r**.3333 + s**.3333)
@KMC.Propensity
def mergin(s,*args,**kwargs):
    return kwargs["frequency"](s[0],s[1])
Z=mergin(x,reactant_range=True,degree=2,frequency=size_dep)
m_args={"Reactant_range":True,"degree":2,"frequency":size_dep}



@KMC.Propensity
def frag(s,*args,**kwargs):
    if s-3>0:
        return kwargs["frequency"]*(s-3)
    else:
        return 0.0

frag_args={"range_start":1,"Range_stop":True,"frequency":100}
mod_args=(sq_kargs,m_args,frag_args)

mod=KMC.Model([[sq,mergin,frag],mod_args],x)

