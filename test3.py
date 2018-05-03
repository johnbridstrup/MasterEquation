import numpy as np
import numpy.random as npr
import scipy.stats as scps
import types
import copy
import itertools as itt
from abc import ABCMeta,ABC,abstractmethod


class testbase(ABC):
    @classmethod
    def __call__(cls, *args, **kwargs):
        print("call")
        cls.do_the_thing(*args,**kwargs)
    
    @classmethod
    @abstractmethod
    def do_the_thing(cls,*args,**kwargs):
        pass

class tester(testbase):
    def __init__(self,d):
        self.d=d
    def doit(self):
        print("d: {}".format(self.d))
    @classmethod
    def do_the_thing(cls,*args,**kwargs):
        print("doing the thing")

if __name__ == '__main__':
    tester.__call__()
    x=tester(3)
    y=x.doit
    x.d=5
    y()
