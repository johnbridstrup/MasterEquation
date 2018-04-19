import numpy as np
import numpy.matlib as mat
import numpy.random as npr
import statevector as sv
import kernels2 as k
from abc import ABCMeta, ABC, abstractmethod
import itertools as itt
import functools as funcs
from functools import reduce
from operator import mul
import copy


    
class MarkovChain:
    def __init__(self):
        self.current_state=None
        self.possible_states=None
        self.Transitions=None

class KMC:
    def __init__(self,statevec):
        self.stop_criteria=False
        self.steps=0
        assert(type(statevec)==sv.StateVector)
                                               
    def reaction_propensities(self):
        pass
    def random_numbers(self):
        pass
    def generate_time_step(self):
        pass
    def determine_reaction(self):
        pass
    def update(self):
        pass
    def calculate_cumulative_function(self):
        pass
    def normalize_probability_function(self):
        pass
    def generate(self):
        pass
class Probability(metaclass=ABCMeta):
    def __init__(self, frequency):
        self.frequency=frequency
        self.probability_vector=None
    @abstractmethod
    def __call__(self,state,*args,**kwargs):
        pass
class Propensity(metaclass=ABCMeta):
    """Decorator to convert function "P(state variable) = probability of single state or transition" into "P(state vector) = vector of (probabilities, index). This is for linear propensities or propensities which are functions of only one state variable."
    
    use **kwargs for:
    indices = [list of indices] applies the function to given indices
    range = [start, end] or single integer applies function to indices in range(start,end) or range(single integer)
    index = single index: applies function to given index
    stoichiometry = for dat stoich
    use *args for:
    first argument will be the value of the state vector index

    example:
    import KMC
    @KMC.Propensity
    def harmonic(x,*args,**kwargs):
        if kwargs["length"] is None:
            kwargs["length"]=1
        return np.sin(x*args[0]*np.pi/kwargs["length"])
    x=sv.StateVector([1,1,1,1])
    print(harmonic(x,length=3,range=len(x.state)))
    # prints [(0.0, 0), (0.8660254037844386, 1), (0.8660254037844388, 2), (1.2246467991473532e-16, 3)]
    # note that index 3 is non-zero due to rounding error
    """
    def __init__(self,func,*args,**kwargs):
        self.propensity_function=func
        self.args=args
        self.kwargs=kwargs
        self.choice=None
    def __call__(self,state_vector,*args,**kwargs):
        self.state=state_vector
        print(kwargs)
        try:
            r_start=kwargs["range_start"]
            r_stop=len(self.state)
            if r_stop>r_start:
                return ((self.propensity_function(self.state[i],self.state[0],*args,**kwargs),i) for i in range(r_start,r_stop) if i>0)
            else:
                return 0.0
        except:
            try:
                rindices=kwargs["reactants"] # list of lists or tuples of reactant indices with repeats
                #pindices=kwargs["products"] # ditto for products

                return ((r,self.propensity_function(self.state[r],r,state=self.state,reactants=rindices)) for r in rindices)
            except:
                try:
                    reactant_range=kwargs["reactant_range"]
                    reaction_degree=kwargs["degree"]
                    react_sets=list(itt.combinations_with_replacement(range(reactant_range),reaction_degree))
                    return ((self.propensity_function([self.state[j] for j in i],i,frequency=kwargs["frequency"]),i) for i in react_sets)
                except:
                    try:
                        indices=kwargs["indices"]
                        return ((self.propensity_function(self.state[i],i,*args,**kwargs),i) for i in indices)
                    except:
                        try:
                            index=kwargs["index"]
                            return (self.propensity_function(self.state[index],index,*args,**kwargs),index)
                        except:
                            try:
                                index=kwargs["indices"]
                                return (self.propensity_function(self.state[index],index,*args,**kwargs),index)
                            except:
                                try:
                                    index_range=kwargs["range"]
                                    return ((self.propensity_function(self.state[i],i,*args,**kwargs),i) for i in range(index_range[0],index_range[1]))
                                except:
                                    try:
                                        return ((self.propensity_function(self.state[i],i,*args,**kwargs),i) for i in range(index_range))
                                    except:
                                        try:
                                            return ((self.propensity_function(self.state[0],*args,**kwargs),0))
                                        except:
                                            print("nothin")
                    
            
        
class InteractionProbability(Probability):
    def P_Pre(self,element,*args,**kwargs):
        return self.P(element,args,kwargs)
    @abstractmethod
    def P(self,element,*args,**kwargs):
        pass
    def __call__(self,state,*args,**kwargs):
        self.probability_vector=state.vectorize_non_zero(self.P_Pre,*args,**kwargs)
class IndexedProbability(Probability):
    def __init__(self, frequency):
        super().__init__(frequency)
        self.counter=0
    def __call__(self,state):
        self.counter=0
        self.probability_vector=state.vectorize_function(self.P_pre)
    def P_pre(self,element):
        self.counter+=1
        self.P(element)
    @abstractmethod
    def P(self,element):
        pass
class MonomerAddition(Probability):
    def __call__(self,state):
        self.probability_vector=state.vectorize_nonzero(self.Padd)
    def Padd(self,state,indices,*args,**kwargs):
        if indices[0]==0:
            indices=indices[1:]
            M=state[0]
        else:
            return [(0,0)]
        return [(M*state[i]*self.frequency,i) for i in indices]
class Nucleation(Probability):
    def __init__(self,nc,frequency):
        super().__init__(frequency)
        self.nc=copy.copy(nc)
        self.M_factor=1
    def __call__(self,state):
        M=state.sum(1)
        self.M_factor=1
        for i in range(copy.copy(self.nc)):
            self.M_factor*=M-i
        self.probability_vector=[(self.M_factor*self.frequency,0)]
class MonomerSubtraction(Probability):
    def __call__(self,state):
        self.probability_vector = state.vectorize_nonzero(self.Psub)
    def Psub(self,state,indices,*args,**kwargs):
        if indices[0]==0:
            indices=indices[1:]
        return [(2*self.frequency*state[i],i) for i in indices]
class Fragmentation(MonomerSubtraction):
    def __init__(self,nc,frequency):
        super().__init__(frequency)
        self.nc=nc
    def __call__(self,state):
        self.probability_vector = state.vectorize_nonzero(self.Pfrag)
        # if self.probability_vector is not None:
        #     self.probability_vector[1][0]=0
    def Pfrag(self,state,indices,*args,**kwargs):
        print(state,indices,"fraggin")
        if copy.copy(self.nc)>3:
            lim=copy.copy(self.nc)
        else: 
            lim=3
        [print("hur",i,state[i],(i-2)*state[i]*self.frequency) for i in indices if i>lim]
        out = [((i-2)*state[i]*self.frequency,i) for i in indices[1:] if i>lim]
        return out
class Coagulation(Probability):
    def __init__(self,frequency):
        super().__init__(frequency)
        self.counter=1
    def __call__(self,state):
        self.probability_vector = state.vectorize_nonzero(self.Pcoag)
    def Pcoag(self,state,indices,*args,**kwargs):
        #for i,v in enumerate(indices):
        print(indices,"indices")
        q=set(itt.chain(indices[1:],indices[1:]))
        w=list(itt.combinations(q,2))
        ww=[(i[0]*i[1]*self.frequency,i) for i in w]
        return ww

class Model:
    def __init__(self,state,nc):
        self.t_step=None
        self.nc=nc
        self.propensities=[]
        self.frequencies=[]
        self.state=state
        print(self.state)
        self.P=[]
        self.summed_P_vector=[]
        self.indices=[]
    # @property
    # def nc(self):
    #     return self.nc
    def add_propensity(self,propensity):
        self.propensities.append(propensity)
    def calculate_probability(self):
        self.P=[]
        self.summed_P_vector=[]
        for i in self.propensities:
            self.P.append(i(self.state.state))
        for i in self.P:
            try:
                self.summed_P_vector.append(sum(i))
            except:
                self.summed_P_vector.append(i)
        print(self.summed_P_vector)
        norm=sum(self.summed_P_vector)
        print(norm)
        self.norm_probability=[i/norm for i in self.summed_P_vector]
        print(self.norm_probability)
    def ready(self):
        self.P=[]
        self.indices=[]
        print(self.propensities)
        for ind,i in enumerate(self.propensities):
            print(i)
            self.P.append([])
            prop_func=i[0]
            ar=i[1]
            try:
                if ar["Range_stop"]==True:
                    ar["range_stop"]=len(self.state.state)
            except:
                try:
                    if ar["Range"]==True:
                        ar["range"]=len(self.state.state)
                except:
                    try:
                        if ar["Reactant_range"]==True:
                            ar["reactant_range"]=len(self.state.state)
                    except:
                        pass
            num=len(self.state.state)-1
            print("FUCK",len(self.state.state)-1)
            try:
                outp=list(prop_func(self.state.state,self.state[0],polymer_number=num,**ar))
            except:
                outp=prop_func(self.state,self.state[0],polymer_number=num,**ar)
            print(outp)
            print(ind)
            print(ind)
            print(ind)
            print(ind)
            print(ind)
            print(ind)
            print(ind)
            print(ind)
            self.P[ind].append(outp)
        print(self.P)
    def BulkProbability(self):
        self.bulk_P=[]
        for i in self.P:
            msum=0
            if type(i[0])==float or type(i[0])==int:
                self.bulk_P.append(i[0])
            else:
                for j in i:
                    try:
                        msum+=j[0]
                    except:
                        pass
                self.bulk_P.append(msum)
        sss=sum(self.bulk_P)
        self.norm_probability=[i/sss for i in self.bulk_P]

           
    def choose(self):
        self.choice=npr.choice(range(len(self.summed_P_vector)), p=self.norm_probability)
        print(self.choice)

    def time_step(self):
        self.t_step=(1/sum(self.summed_P_vector))*np.log10(1/npr.random())
        print(self.t_step,"Time Step")
    def advance(self):
        print("CHOICE",self.choice)
        if self.choice==0:
            index=npr.choice(range(len(self.state.state[1:])))+1
            self.state.state[index]=self.state.state[index]+1
            self.state.state[0]=self.state.state[0]-1
            print("adding", index, self.state.state[index])
        if self.choice==1:
            print("subbing")
            index=npr.choice(range(len(self.state.state[1:])))+1
            self.state.state[index]=self.state.state[index]-1
            if self.state.state[index]<self.nc:
                self.state.state[0]=self.state.state[0]+self.state.state[index]
                self.state.state[index]=self.state.state[index]-self.state.state[index]
            self.state.state[0]=self.state.state[0]+1
        if self.choice==2:
            print("nucleating")
            print(self.state.state)
            print(self.nc)
            self.state.state.append(self.nc)
            self.state.state[0]=self.state.state[0]-self.propensities[0].nc
        if self.choice==3:
            print("break")
            print(self.state.state)
            cs=np.cumsum([j for j in self.state.state[1:] if j>3])
            print(cs,"cs")
            ss=sum([j for j in self.state.state if j>3])
            print(ss,"ss")
            pindex=npr.choice(range(ss))
            print(pindex,"pindex")
            i=0
            count=0
            while(i<pindex):
                count+=1
                if self.state.state[count]>3:
                    i+=self.state.state[count]
                if i>pindex:
                    r=i-pindex
                    print(r,"r")
                    print(count,"count")
                    self.state.state[count]=self.state.state[count]-r
                    if self.state.state[count]<self.nc:
                        self.state.state[0]=self.state.state[0]+self.state.state[count]
                        del self.state.state[count]
                    if r<self.nc:
                        self.state.state[0]=self.state.state[0]+r
                    else:
                        self.state.state.append(r)
        print("SUM",sum(self.state.state))
        if self.choice==4:
            print("merge")
            index=npr.choice(range(len(self.state.state[1:])))+1
            r=self.state.state[index]
            del self.state.state[index]
            index=npr.choice(range(len(self.state.state[1:])))+1
            self.state.state[index]=self.state.state[index]+r
        