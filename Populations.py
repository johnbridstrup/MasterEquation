import numpy as np
import numpy.random as npr
import scipy.stats as stats
from collections import deque
import copy

class PointError(ValueError):
    pass
class Point:
    def __init__(self,data,time=0.0):
        self.point=(data,time)
        self.value=data
        self.t=time
        self.next=None
    def __repr__(self):
        return repr(self.value)
    def __add__(self,other):
        try:
            return self.value+other.value
        except:
            return self.value+other
    def __radd__(self,other):
        return self.__add__(other)
                  
class Population:
    _ID=0
    _counter=0
    _ID_index_map={}
    _archived=[]
    @classmethod
    def archive(cls,population):
        cls._archived.append(population)
    @classmethod
    def increment(cls):
        cls._counter+=1
        cls._ID+=1
    def __init__(self,point=None):
        self.head=point
        self._ID=Population._ID
        self._index=Population._counter
        Population._ID_index_map[self._index]=self._ID
    def __iter__(self):
        curr = self.head
        while curr is not None:
            yield curr
            curr=curr.next
    def push_back(self,point):
        point.next=self.head
        self.head=point
    def get_all(self):
        return [(i,i.t) for i in self]
    def __repr__(self):
        return repr(self.head.value)
    def __add__(self, other):
        return self.head.__add__(other)
    def __radd__(self,other):
        return self.head.__add__(other)

