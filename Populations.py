import numpy as np
import numpy.random as npr
import scipy.stats as stats
from collections import deque
import copy


class PointError(ValueError):
    pass


class Point:
    def __init__(self, data, time=0.0):
        self.point = (data, time)
        self.value = data
        
        self.t = time
        
        self.nxt = None

    def __repr__(self):
        return repr(self.value)

    def __add__(self, other):
        try:
            return self.value+other.value
        except:
            return self.value+other
    def __iadd__(self, other):
        ## replace current with temp
        tmp1=Point(self.value,self.t)
        tmp2=Point(self.value,self.t+other.t)
        tmp1.nxt=self.nxt
        tmp2.nxt=tmp1
        self.value=self.value+other.value
        self.t=self.t+other.t
        self.nxt=tmp2
        return self
    def push_back(self,data,t):
        self.__iadd__(Point(data,t))
    def __radd__(self, other):
        return self.__add__(other)
    def __lt__(self, other):
        try:
            return self.value<other
        except:
            try:
                return self.value<other.value
            except:
                raise ValueError('Point compares with point, int or float')
    def __gt__(self, other):
        try:
            return self.value>other
        except:
            try:
                return self.value>other.value
            except:
                raise ValueError('Point compares with point, int or float')
    def __eq__(self, other):
        try:
            return self.value==other
        except:
            try:
                return self.value==other.value
            except:
                raise ValueError('Point compares with point, int or float')
    def __le__(self, other):
        try:
            return self.value<=other
        except:
            try:
                return self.value<=other.value
            except:
                raise ValueError('Point compares with point, int or float')
    def __ge__(self, other):
        try:
            return self.value>=other
        except:
            try:
                return self.value>=other.value
            except:
                raise ValueError('Point compares with point, int or float')
    def __ne__(self, other):
        try:
            return self.value!=other
        except:
            try:
                return self.value!=other.value
            except:
                raise ValueError('Point compares with point, int or float')
    
    def __getitem__(self, key):
        if key == 't' or key == 'time' or key == 'Time' or key == 'T' or key == 1:
            return self.t
        if key == 'value' or key == 'val' or key == 'Val' or key == 'Value' or key == 0:
            return self.value

    def __setitem__(self, key, val):
        if key == 't' or key == 'time' or key == 'Time' or key == 'T' or key == 1:
            self.t = val
        if key == 'value' or key == 'val' or key == 'Val' or key == 'Value' or key == 0:
            self.value = val
    def __iter__(self):
        curr=self
        while curr is not None:
            yield curr
            curr=curr.nxt
    def __mul__(self,other):
        try:
            return self.value*other.value
        except:
            try:
                return self.value * other
            except Exception as e:
                raise ValueError('multiplying by an unsupported type ::: ', e)
    def __rmul__(self,other):
        try:
            return self.__mul__(other)
        except Exception as e:
            raise ValueError('Unsupported Jawn',e)
    def get_all(self):
        return [(i, i.t) for i in self]
class PointUpdate(Point):
    def __init__(self,data,t):
        super().__init__(data,t)
    def __getitem__(self, key):
        if key == 't' or key == 'time' or key == 'Time' or key == 'T':
            return self.t
        elif key == 'value' or key == 'val' or key == 'Val' or key == 'Value':
            return self.value
        elif key == 'key':
            return 'vt'

class CompositionVector:
    _archive=[]
    @classmethod
    def archive(cls,point):
        cls._archive.append(point)
    def __init__(self, point=None, nc=0):
        self.__methods={'sum':self.__sum,'average':self.__avg, 'length':self.__len,}
        if point is not None:
            self.initialized=True
            try:
                self.state=np.array(point)
            except:
                print('gotta fit into a numpy array')
                raise
        else:
            self.initialized=False
    def add_state_variable(self,point):
        if not self.initialized:
            try:
                self.state=np.array(point)
                for i in self.state:
                    print(i.get_ID())
            except:
                print('gotta finnarray')
                raise
        else:
            self.state=np.append(self.state,point)
    def __getitem__(self, key):
        if isinstance(key,int):
            try:
                return self.state[key]
            except:
                try:
                    return self.state
                except:
                    raise IndexError("Index wrong dogg")
        elif key=='average':
            return sum(self.state)/len(self.state)
        elif key=='sum':
            return sum(self.state)
        elif key=='length':
            return len(self.state)
        elif isinstance(key,tuple):
            try:
                outp=self.state[key[0]]
                for i in range(key[1]):
                    outp=outp.next
                return outp
            except:
                print("that time step probably doesnt exist")
    def __repr__(self):
        return repr(self.state)
    def __str__(self):
        return str(self.state)
    def __setitem__(self,index,key):
        pass
    def __iter__(self):
        self.it=0
        return self
    def length(self):
        return len(self.state)
    def __next__(self):
        index=self.it
        self.it+=1
        while self.it<=len(self.state):
            return self.state[index]
        raise StopIteration
    def __sum(self):
        return sum(self.state)
    def __avg(self):
        return sum(self.state)/len(self.state)
    def __len(self):
        return len(self.state)
    def sub(self,index,num,**kwargs):
        try:
            self.state[index]+=Point(num,kwargs['t'])
        except:
            self.state[index]+=Point(num)
    def sum(self,shift=1):
        return sum(self.state[shift:])
    def skew(self,shift=1):
        return stats.skew(self.state[shift:])    
    def moment(self,n,shift=1):
        return stats.moment(self.state[shift:],moment=n)
    def histogram(self,shift=0):
        return np.histogram(self.state[shift:], bins=99, range=(1,100))
    def kurtosis(self,shift=1):
        return stats.kurtosis(self.state[shift:])
    def coefficient_of_variation(self,shift=1):
        return stats.variation(self.state[shift:])
    
        
def main():
    lst=[0,1,2,3,4,5]
    tlst=[0.0,0.1,0.2,0.3,0.4,0.5]
    points=[Point(i,j) for i,j in zip(lst,tlst)]
    state=CompositionVector(points)
    state[2]+=Point(11,11)
    state[4]+=Point(100,1000)
    for i in state:
        print(i.get_all())
    print(sum(state))

if __name__ == '__main__':
    main()
