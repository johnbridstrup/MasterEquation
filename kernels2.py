#! kernels.py
import KMC
import statevector as sv
from abc import ABC,abstractmethod
from functools import wraps

def avo():
    return 6.022*10**23
def RateError(Exception):
    pass
def NoRatesError(RateError):
    pass
def RateTransformError(RateError):
   pass
def call_later(f):
    @wraps(f)
    def wrapper(*args, **kwargs):
        #print("calling later")
        return f
    #print("called later")
    return wrapper
class Rates:
    def __init__(self,R=None,T=None,*args,**kwargs):
        self.errs={}
        self.R=R
        self.T=T
        if R is not None:
            if T is not None:
                try:
                    self.rates=R*T
                except:
                    try:
                        self.rates=T(R)
                    except:
                        msg="can't transform even though R and T are given. Will hold for later"
                        print(msg)
                        self.errs['transform']=(RateTransformError,msg)
                        self.hold_ToR=call_later(T)
            else:
                if callable(R):
                    try:
                        self.rates=R(*args,**kwargs)
                    except:
                        try:
                            self.rates=R()
                        except:
                            msg="can't call R despite callable(R). will try holding it for later"
                            print(msg)
                            self.errs['rates']=(NoRatesError,msg)
                            self.hold_R=call_later(R)
        else:
            try:
                self.rates=T(kwargs)
            except:
                try:
                    self.rates=T()
                except:
                    msg="Holding T for later. Cant execute alone and have no R"
                    print(msg)
                    self.errs['transform']=(RateTransformError,msg)
                    self.hold_T=call_later(T)
class TransformFunction(ABC):
    def __init__(self,*args,**kwargs):
        self.args=args
        self.kwargs=kwargs
    @abstractmethod
    def __call__(self,R):
        pass      
def default_rates():
    R={}
    R['a']=1
    R['b']=0.000001
    R['am']=R['a']
    R['bf']=R['b']
    R['kn']=0.00001
    return R
class DefaultRates(Rates):
    def __init__(self):
        super().__init__(default_rates())

def I(inp,*ins,**kwins):
    return inp
# class IRates(Rates):
#     def     



class Monomers:
    def __init__(self,M=500):
        self.M=500
    def __call__(self):
        return [self.M]
class RateTransform(ABC):
    def __init__(self,f,**kwargs):
        self.freq=None
        # self.Volume=None
        self.params=kwargs
        self.bulkrate=f
    def __call__(self,**kwargs):
        return self.freq
    @abstractmethod
    def transform(self):
        pass
 
class Bimolecular(RateTransform):
    def __init__(self,f,M=500,c=1,same=False):
        super().__init__(f,M=M,c=c,same=same)
        self.transform()
    def transform(self):
        self.f=self.bulkrate*self.params['c']/self.params['M']
    def __call__(self,**kwargs):
        outp=self.freq
        outp*=self.params['N1']
    # @property
    # def Volume(self):
    #     pass
class N_Molecular(RateTransform):
    def __init__(self, f,M=500,c=1,nc=3,same=False):
        super().__init__(f,M=500,c=1,nc=3,same=same)
    def transform(self):
        nn=self.bulkrate
        if self.params['same']:
            for i in range(self.params['nc']):
                nn*=(self.params['M']-i)
            self.freq=nn
        else:
            raise ValueError('uhhhhh idk havent implemented it')      
class default_transform(RateTransform):
    def transform(self):
        self.freq=self.bulkrate
class Unimolecular(default_transform):
    pass
class Kernel(ABC):
    def __init__(self, f, bulk_transform=default_transform):
        self.bulk_transform=bulk_transform
        self.freq=f
        self.propensities=[]
        try:
            self.freq=bulk_transform(bulk_transform.freq,**bulk_transform.params)
        except:
            pass
    def add_propensity(self,prop,label=None):
        if label is None:
            label = prop.__name__
        try:
            for key,val in prop.items:
                self.propensities.append(val)
        except:
            try:
                [self.propensities.append(item) for label,item in prop.items()]
            except:
                try:
                    self.propensities.append(prop)
                except:
                    print('propensity adding is broked')
                    raise                  
class Coagulation(Kernel):
    def __init__(self, f, nc, bulk_transform=None, *args, **kwargs):
        super().__init__(f, bulk_transform)
        self.index_pairs = []
        
    def __call__(self, s):
        try:
            outp = []
            print("THE SSSSS ", s)
            if len(s) > 2:
                for ij, j in enumerate(s[1:]):
                    for ii, i in enumerate(s[ij:]):
                        self.index_pairs.append((ii+1, ij+ii+1))
                        outp.append(self.freq)
                return outp
            else:
                return 0.0
        except:
            return 0.0
class Fragmentation(Kernel):
    def __init__(self, f, nc, bulk_transform=default_transform, *args, **kwargs):
        super().__init__(f,bulk_transform)

    def __call__(self, s):
        try:
            if len(s) > 1:
                for i in s[1:]:
                    if i > 3:
                        return [(j-3)*self.freq for j in range(4, i)]
                    else:
                        return 0.0
            else:
                return 0.0
        except:
            return 0.0
class MonomerAddition(Kernel):
    def __init__(self, f, nc, bulk_transform=Bimolecular,*args, **kwargs):
        super().__init__(f,bulk_transform)
        self.nc = nc

    def __call__(self, s):
        try:
            if s[0] > 0 and len(s) > 1:
                return [s[i]*s[0]*self.freq for i in range(1, len(s))]
            else:
                return 0.0
        except:
            return 0.0
class Nucleation(Kernel):
    def __init__(self, f, nc,bulk_transform=N_Molecular,**kwargs):
        super().__init__(f,bulk_transform)
        self.nc = nc
        try:
            self.freq = self.freq *(kwargs["c"]/kwargs["M"])**self.nc
        except:
            pass

    def __call__(self, s, *args, **kwargs):
        prod = self.freq
        print("WE IN DAT CALL",type(s))
        print("____",type(s[0]),s[0],self.nc)
        try:
            if s[0] > self.nc:
                print(3)
                for i in range(self.nc):
                    prod = 1.0*prod*(1.0*s[0]-1.0*i)
                print("daprooood",prod)
                return prod
            else:
                return 0.0
        except Exception as e:
            print("catchin some't ", e)
            return 0.0
class MonomerSubtraction(Kernel):
    def __init__(self, f, bulk_transform=Unimolecular,*args, **kwargs):
        super().__init__(f,bulk_transform)

    def __call__(self, s, **kwargs):
        try:
            if len(s) > 1:
                return [2*self.freq for i in range(len(s)-1)]
            else:
                return 0.0
        except:
            return 0.0
class Propensities:
    def __init__(self,props,names=None):
        self.propensities={}
        try:
            for key,val in props.items():
                self.propensities[key]=val
        except:
            try:
                for name,prop in zip(names,props):
                    self.propensities[name]=prop
            except:
                try: 
                    for index,prop in zip(range(len(props)),props):
                        self.propensities[index]=prop
                except:
                    try:
                        self.propensities[names]=props
                    except:
                        try:
                            assert(type(props) != list)
                            self.propensities[0]=props
                        except:
                            print("WTF is this?: ",props)