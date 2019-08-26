#!/usr/bin/env python3

import KMC
import statevector as sv
from abc import ABC, abstractmethod
from functools import wraps


def avo():
    return 6.022*10**23


def RateError(Exception):
    pass


def NoRatesError(RateError):
    pass


def RateTransformError(RateError):
    pass
# def call_later(f):
#     @wraps(f)
#     def wrapper(*args, **kwargs):
#         # print("calling later")
#         return f
#     # print("called later")
#     return wrapper
# class Rates:
#     def __init__(self,R=None,T=None,*args,**kwargs):
#         self.errs={}
#         self.R=R
#         self.T=T
#         if R is not None:
#             if T is not None:
#                 try:
#                     self.rates=R*T
#                 except:
#                     try:
#                         self.rates=T(R)
#                     except:
#                         msg="can't transform even though R and T are given. Will hold for later"
#                         # print(msg)
#                         self.errs['transform']=(RateTransformError,msg)
#                         self.hold_ToR=call_later(T)
#             else:
#                 if callable(R):
#                     try:
#                         self.rates=R(*args,**kwargs)
#                     except:
#                         try:
#                             self.rates=R()
#                         except:
#                             msg="can't call R despite callable(R). will try holding it for later"
#                             # print(msg)
#                             self.errs['rates']=(NoRatesError,msg)
#                             self.hold_R=call_later(R)
#         else:
#             try:
#                 self.rates=T(kwargs)
#             except:
#                 try:
#                     self.rates=T()
#                 except:
#                     msg="Holding T for later. Cant execute alone and have no R"
#                     # print(msg)
#                     self.errs['transform']=(RateTransformError,msg)
#                     self.hold_T=call_later(T)
# class TransformFunction(ABC):
#     def __init__(self,*args,**kwargs):
#         self.args=args
#         self.kwargs=kwargs
#     @abstractmethod
#     def __call__(self,R):
#         pass
# def default_rates():
#     R={}
#     R['a']=1
#     R['b']=0.000001
#     R['am']=R['a']
#     R['bf']=R['b']
#     R['kn']=0.00001
#     return R
# class DefaultRates(Rates):
#     def __init__(self):
#         super().__init__(default_rates())

# def I(inp,*ins,**kwins):
#     return inp
# # class IRates(Rates):
# #     def


def factorial(n):
    if n < 1:
        return 1
    else:
        out = n*factorial(n-1)
        return out


def choose(n, k):
    if n > k:
        return factorial(n)/(factorial(k)*factorial(n-k))
    else:
        return 0.0


class Monomers:
    def __init__(self, M=500):
        self.M = M

    def __call__(self):
        return [self.M]


class RateTransform(ABC):
    def __init__(self, f, **kwargs):
        self.freq = None
        # self.Volume=None
        self.params = kwargs
        self.bulkrate = f
        # print("here3")

    def __call__(self, **kwargs):
        return self.freq

    @abstractmethod
    def transform(self):
        pass


class Bimolecular(RateTransform):
    def __init__(self, f, same=False, **kwargs):
        super().__init__(f, **kwargs)
        # print("here2")
        self.transform()
        self.same = same

    def transform(self):
        self.freq = self.bulkrate*self.params['c']/self.params['M']

    def __call__(self, **kwargs):
        return self.freq
    # @property
    # def Volume(self):
    #     pass


class N_Molecular(RateTransform):
    def __init__(self, f, M=500, c=1, nc=3, same=False):
        super().__init__(f, M=500, c=1, nc=3, same=same)
        print("nmol")
        self.transform()

    def transform(self):
        nn = self.bulkrate
        # print(nn)
        if self.params['same']:
            for i in range(self.params['nc']):
                nn *= (self.params['M']-i)
            self.freq = nn
        else:
            raise ValueError('uhhhhh idk havent implemented it')

    def __call__(self):
        return self.freq


class default_transform(RateTransform):
    def transform(self):
        self.freq = self.bulkrate


class Unimolecular(default_transform):
    pass


class Kernel(ABC):
    def __init__(self, f, bulk_transform=default_transform, **kwargs):
        self.bulk_transform = bulk_transform
        self.freq = None
        self.bulk_rate = f
        self.propensities = []
        try:
            # print("here1")
            self.freq = bulk_transform(self.bulk_rate, **kwargs)
        except Exception:
            self.freq = self.bulk_rate

    def add_propensity(self, prop, label=None):
        if label is None:
            label = prop.__name__
        try:
            for key, val in prop.items:
                self.propensities.append(val)
        except:
            try:
                [self.propensities.append(item)
                 for label, item in prop.items()]
            except:
                try:
                    self.propensities.append(prop)
                except:
                    # print('propensity adding is broked')
                    raise


class MonomerAddition(Kernel):
    def __init__(self, f, nc, bulk_transform=Bimolecular, *args, **kwargs):
        super().__init__(f, bulk_transform, **kwargs)
        self.nc = nc

    def __call__(self, s):
        # print("calling add")
        try:
            if s[0] > 0 and len(s) > 1:
                return [s[i]*s[0]*self.freq for i in range(1, len(s))]
            else:
                return 0.0
        except:
            return 0.0


class Coagulation(Kernel):
    def __init__(self, f, nc, bulk_transform=Bimolecular, *args, **kwargs):
        super().__init__(f, bulk_transform, **kwargs)
        self.index_pairs = []

    def __call__(self, s):
        try:
            outp = []
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
    def __init__(self, f, nc, bulk_transform=None, *args, **kwargs):
        super().__init__(f, bulk_transform, **kwargs)

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


class Nucleation(Kernel):
    def __init__(self, f, nc, bulk_transform=N_Molecular, **kwargs):
        super().__init__(f, bulk_transform, M=kwargs['M'], nc=nc, **kwargs)
        self.nc = nc
        try:
            self.freq = self.freq * (kwargs["c"]/kwargs["M"])**self.nc
        except:
            pass

    def __call__(self, s, *args, **kwargs):
        prod = self.freq
        # print("WE IN DAT CALL",type(s))
        # print("____",type(s[0]),s[0],self.nc)
        try:
            if s[0] > self.nc:
                # print(3)
                for i in range(self.nc):
                    prod = 1.0*prod*(1.0*s[0]-1.0*i)
                return prod
            else:
                return 0.0
        except:
            return 0.0


class MonomerSubtraction(Kernel):
    def __init__(self, f, bulk_transform = Unimolecular, *args, **kwargs):
        super().__init__(f, bulk_transform)

    def __call__(self, s, **kwargs):
        try:
            if len(s) > 1:
                return [2*self.freq for i in range(len(s)-1)]
            else:
                return 0.0
        except:
            return 0.0


class Propensities:
    def __init__(self, props, names=None):
        self.propensities = {}
        try:
            for key, val in props.items():
                self.propensities[key] = val
        except:
            try:
                for name, prop in zip(names, props):
                    self.propensities[name] = prop
            except:
                try:
                    for index, prop in zip(range(len(props)), props):
                        self.propensities[index] = prop
                except:
                    try:
                        self.propensities[names] = props
                    except:
                        try:
                            assert(type(props) != list)
                            self.propensities[0] = props
                        except:
                            pass
                            # print("WTF is this?: ",props)


class betterKernel(ABC):
    def __init__(self, bulk_rate_constant, **kwargs):
        self.bulk_rate_constant = bulk_rate_constant
        self.params = kwargs
        self.propensity = 0.0
        self.probability = 0.0

    @abstractmethod
    def transform(self):
        pass

    def __repr__(self):
        return repr(self.propensity)

    def __str__(self):
        return str(self.propensity)

    def __call__(self, s):
        return self.propensity


class MonAdd(betterKernel):
    def __init__(self, bulk_rate_constant, **kwargs):
        super().__init__(bulk_rate_constant, **kwargs)
        self.transform()

    def transform(self):
        self.propensity = self.bulk_rate_constant * \
            self.params['c']/self.params['M']

    def __call__(self, s):
        if len(s) > 1:
            self.probability = sum([self.propensity*i*s[0] for i in s[1:]])
            return self.probability
        else:
            return 0.0


class MonSub(betterKernel):
    def __init__(self, bulk_rate_constant, **kwargs):
        super().__init__(bulk_rate_constant, **kwargs)
        self.propensity = bulk_rate_constant

    def transform(self):
        pass

    def __call__(self, s):
        if len(s) > 1:
            self.probability = self.propensity*2*(len(s)-1)
            return self.probability
        else:
            return 0.0


class Coag(betterKernel):
    def __init__(self, bulk_rate_constant, **kwargs):
        super().__init__(bulk_rate_constant, **kwargs)
        self.transform()

    def transform(self):
        self.propensity = self.propensity = self.bulk_rate_constant * \
            self.params['c']/self.params['M']

    def __call__(self, s):
        L = len(s)
        if L > 2:
            self.probability = self.propensity*(factorial(L)/2)
            return self.probability
        else:
            return 0.0


class Nuc(betterKernel):
    def __init__(self, bulk_rate_constant, **kwargs):
        super().__init__(bulk_rate_constant, **kwargs)
        self.transform()

    def transform(self):
        self.propensity = self.bulk_rate_constant
        self.propensity *= (self.params['c'] /
                            (self.params['M']))**self.params['nc']

    def __call__(self, s):
        if s[0] > self.params['nc']:
            self.probability = self.propensity
            for i in range(self.params['nc']):
                self.probability *= (s[0]-i)
            return self.probability
        else:
            return 0.0


class Frag(betterKernel):
    def __init__(self, bulk_rate_constant, **kwargs):
        super().__init__(bulk_rate_constant)
        self.propensity = self.bulk_rate_constant

    def transform(self):
        self.propensity = self.bulk_rate_constant

    def __call__(self, s):
        if len(s) > 1:
            try:
                self.probability = [
                    (i-3)*self.propensity for i in s[1:] if i > 3]
                return self.probability
            except:
                return [0]
        else:
            return 0.0
