#! kernels.py
import KMC
import statevector as sv


class MergeKernel:
    def __init__(self, frequency, **kwargs):
        self.frequency = frequency
        try:
            self.fractal_dimension = kwargs["fractal_dimension"]
        except:
            self.fractal_dimension = 1.0/3.0

    def __call__(self, r, s):
        return 1.0*self.frequency*(r**(self.fractal_dimension)+s**(self.fractal_dimension))


class Coagulation:
    def __init__(self, f, nc, *args, **kwargs):
        self.freq = f
        self.index_pairs = []
        try:
            self.freq=self.freq*kwargs["c"]/kwargs["M"]
        except:
            pass
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


class Fragmentation:
    def __init__(self, f, nc, *args, **kwargs):
        self.freq = f

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


class MonomerAddition:
    def __init__(self, f, nc, *args, **kwargs):
        self.freq = f
        self.nc = nc
        try:
            self.freq = self.freq*kwargs["c"]/kwargs["M"]
        except:
            pass

    def __call__(self, s):
        try:
            if s[0] > 0 and len(s) > 1:
                return [s[i]*s[0]*self.freq for i in range(1, len(s))]
            else:
                return 0.0
        except:
            return 0.0


class Nucleation:
    def __init__(self, f, nc,**kwargs):
        self.freq = f
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


class MonomerSubtraction:
    def __init__(self, f, *args, **kwargs):
        self.freq = f

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