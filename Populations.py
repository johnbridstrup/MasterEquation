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
        self.next = None

    def __repr__(self):
        return repr(self.value)

    def __add__(self, other):
        try:
            return self.value+other.value
        except:
            return self.value+other

    def __radd__(self, other):
        return self.__add__(other)

    def __getitem__(self, key):
        if key == 't' or key == 'time' or key == 'Time' or key == 'T':
            return self.t
        if key == 'value' or key == 'val' or key == 'Val' or key == 'Value':
            return self.value

    def __setitem__(self, key, val):
        if key == 't' or key == 'time' or key == 'Time' or key == 'T' or key == 1:
            self.t = val
        if key == 'value' or key == 'val' or key == 'Val' or key == 'Value' or key == 0:
            self.value = val


class Population:
    _ID = 0
    _counter = 0
    _ID_index_map = {}
    _archived = []

    @classmethod
    def archive(cls, population):
        cls._archived.append(population)

    @classmethod
    def increment(cls):
        cls._counter += 1
        cls._ID += 1

    def __init__(self, point=None):
        self.head = point
        self._ID = Population._ID
        self._index = Population._counter
        Population._ID_index_map[self._index] = self._ID
        Population.increment()

    def __iter__(self):
        curr = self.head
        while curr is not None:
            yield curr
            curr = curr.next

    def push_back(self, point):
        try:
            assert(isinstance(point, Point))
            point.next = self.head
            self.head = point
        except:
            try:
                tmp_next = Point(point[0], point[1])
                tmp_next.next = self.head
                self.head = point
            except:
                print("input to push_back must be Point or container of two elements")
                raise PointError

    def get_all(self):
        return [(i, i.t) for i in self]

    def __repr__(self):
        return repr(self.head.value)

    def __add__(self, other):
        return self.head.__add__(other)

    def __radd__(self, other):
        return self.head.__add__(other)


class CompositionVector:
    def __init__(self, population=None):
        if population is not None:
            try:
                assert(isinstance(population, Population))
                self.comp_vector = np.array(population)
            except:
                try:
                    assert(isinstance(population, Point))
                    self.comp_vector = np.array(Population(population))
                except:
                    print(
                        'Composition vector must contain instances of Population or Point')
                    raise PointError
        else:
            self.comp_vector = None
