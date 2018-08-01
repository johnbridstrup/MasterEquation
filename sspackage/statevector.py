#!/Users/John1/anaconda2/envs/py3/bin/python
# -*- coding: utf-8 -*-
#
# statevector.py
# @Author : John Bridstrup (john.bridstrup@gmail.com)
# @Link   :
# @Date   : 2018-4-12 15:12:08
import numpy as np
import numpy.random as npr
import scipy.stats as scps
import types
import copy
import itertools as itt
from abc import ABCMeta, ABC, abstractmethod

# ==============================================================================
# MODEL BUILDERS
# ==============================================================================
class ModelProcedure:
    def __init__(self):
        self.mechanisms={}
    
    def add_mechanism(self,mech,name):
        self.mechanisms[name]= mech
    def run(self,name):
        self.mechanisms[name]()
def SmoluchowskiModel(state,nc):
    model=ModelProcedure()
    model.add_mechanism(MonomerAddition(state,nc),"add")
    model.add_mechanism(MonomerSubtraction(state,nc),"subtract")
    model.add_mechanism(Nucleation(state,nc),"nucleate")
    model.add_mechanism(Fragmentation(state,nc),"break")
    model.add_mechanism(Coagulation(state,nc),"merge")
    return model
# ==============================================================================
# MODEL PRIMITIVES
# ==============================================================================
class StateOperation(ABC):
    """base class for state propagation

    algorithm must be overridden, but the others can be ignored. call them from within the previous function if it makes sense
    """

    def __init__(self, state, *args, **kwargs):
        self.state = state
        self.routines=[]
    def __call__(self, *args, **kwargs):
        # print("calling")
        self.algorithm(*args, **kwargs)
    
    @abstractmethod
    def algorithm(self, *args, **kwargs):
        pass


class ProteinAggregationBase(StateOperation, ABC):
    def __init__(self,state,nc):
        super().__init__(state)
        self.nc=nc
    def length(self):
        return self.state.length()

class MonomerAddition(ProteinAggregationBase):
    def algorithm(self, *args, **kwargs):
        index = npr.choice(range(self.length()-1))+1
        self.state.add(index, 1)
        self.state.sub(0, 1)


class MonomerSubtraction(ProteinAggregationBase):
    def algorithm(self, *args, **kwargs):
        index = npr.choice(range(self.length()-1))+1
        final = self.state.sub(index, 1)
        if final < self.nc:
            self.state.add(0, final)
            self.state.remove(index)
        self.state.add(0, 1)


class Nucleation(ProteinAggregationBase):
    def algorithm(self, *args, **kwargs):
        self.state.add_state_variable(self.nc)
        self.state.sub(0, self.nc)


class Fragmentation(ProteinAggregationBase):
    def algorithm(self, *args, **kwargs):
        
        ss = self.state.sum(condition=self.nc)
        print("ss",ss)
        print("state",self.state)
        print("sum",self.state.sum(condition=self.nc))
        pindex = npr.choice(range(ss))
        i = 0
        count = 0
        while(i < pindex):
            count += 1
            if self.state.get(count) > 3:
                i += self.state.get(count)
            if i > pindex:
                r = i-pindex
                # print(r, "r")
                # print(count, "count")
                brk = self.state.sub(count, r)
                if brk < self.nc:
                    self.state.add(0, brk)
                    self.state.remove(count)
                if r < self.nc:
                    self.state.add(0, r)
                else:
                    self.state.add_state_variable(r)


class Coagulation(ProteinAggregationBase):
    def algorithm(self, *args, **kwargs):
        index = npr.choice(range(self.length()-1))+1
        r = self.state.get(index)
        self.state.remove(index)
        index = npr.choice(range(self.length()-1))+1
        self.state.add(index, r)
# ==============================================================================
# exceptions
# ==============================================================================

class InputException(Exception):
    def __init__(self, node=None, msg=None):
        if msg is None:
            # Set some default useful error message
            msg = "Input is not an array or list"
        super(InputException, self).__init__(msg)

# ==============================================================================
# END RANDOM FUNCTIONS
# ==============================================================================


class StateVector:
    #` LITERALLY JUST MAKE IT A SLIGHTLY MORE SPECIALIZED ARRAY

    def __init__(self, the_input=None, initial_conditions=None, delete_zeros=False):
        if the_input is not None:
            try:
                self.state = np.array(the_input.state)
                self.initial_conditions = the_input.initial_conditions
                self.delete_zeros = the_input.delete_zeros
            except:
                try:
                    if type(the_input) == int:
                        raise InputException()
                    self.state = np.array(the_input)
                except (TypeError, IndexError, InputException) as e:
                    pass
                if type(the_input) == int:
                    self.state = np.zeros(the_input)
                else:
                    self.state = the_input
                    self.delete_zeros = delete_zeros
                if initial_conditions is not None:
                    try:
                        assert(any([i != 0 for i in self.state]))
                    except Exception as e:
                        pass
                    try:
                        for i, v in enumerate(initial_conditions):
                            self.state[i] = v
                    except ValueError as e:
                        self.state[0] = initial_conditions
                        pass

            else:
                if initial_conditions is not None:
                    try:
                        self.state = np.array(initial_conditions)
                    except:
                        # print("its some weird, not array-ifiable type")
                        raise
            self.delete_zeros = delete_zeros

    def __repr__(self):
        return repr(self.state)

    def __str__(self):
        return str(self.state)
    def __len__(self):
        return len(self.state)

    def __getitem__(self, index):
        try:
            return self.state[index]
        except (TypeError, IndexError) as e:
            # print("No string get implementiation, you must use an index")
            raise e

    def length(self):
        return len(self.state)

    def __setitem__(self, key, val):
        try:
            self.state[key] = val
        except (IndexError, KeyError) as e:
            raise e
# ==============================================================================
# Getting and setting elements
# ==============================================================================

    def get_index_by_value(self, val, random=False):
        if random == True:
            return npr.choice(np.argwhere(self.state == val).T[0])
        else:
            return np.argwhere(self.state == val.T[0])

    def set_state_variable(self, index, value):
        try:
            self.state[index] = value
        except:
            # print("Tried to access and set by index value")
            # print("index: ", index)
            # print("value: ", value)
            # print("state: ", self.state)
            raise

    def random_set_by_value(self, old_value, new_value, add=False):
        index = self.get_index_by_value(old_value, random=True)
        if add != True:
            try:
                self.set_state_variable(index, new_value)
            except:
                self.dev_info(state=self.state, old_val=old_value,
                              new_value=new_value, index=index)
        else:
            self.state[index] += new_value
            if self.delete_zeros and self.state[index] <= 0:
                try:
                    self.state.remove_state_variable(index)
                except:
                    self.dev_info(index=index, state=self.state)
                    raise

    def addOne(self, shift=1):
        index = npr.choice(range(len(self.state)))
        # print("INDEX", index)
        self.state[index] = self.state[index]+1
        if shift == 1:
            self.state[0] = self.state[0]-1

    def subOne(self, nc, shift=1):
        index = npr.choice(range(len(self.state[1:])))
        self.state[index] = self.state[index]-1
        if self.state[index] < nc:
            self.state[0] = self.state[0]+self.state[index]
            self.remove_state_variable(index)
        if shift == 1:
            self.state[0] = self.state[0]+1

    def breakOne(self, nc):
        index = npr.choice(
            [i for i, v in enumerate(self.state) if v > 3 and i > 0])
        # print("index", index)
        try:
            r = npr.choice(np.array(list(range(4, self.state[index]))))
            # print(r, "r")
            self.state[index] = self.state[index]-r
            if self.state[index] < nc:
                self.state[0] = self.state[0]+self.state[index]
                self.remove_state_variable(index)
            if r > nc:
                self.state = self.state+r
            else:
                self.state[0] = self.state[0]+r
        except:
            pass
            # print("WHAT THE FUCKING UCK YOU FUCK")

    def merge(self):
        r = npr.choice(self.state[1:])
        # print(r, "merrrge")
        self.random_remove_val(r)
        s = npr.choice(range(len(self.state[1:])))+1
        self.state[s] = self.state[s]+r

    def non_zero(self):
        return np.nonzero(self.state)[0]
# ==============================================================================
# Adding and removing elements
# ==============================================================================

    def remove_state_variable(self, inpt):
        try:
            self.state = np.delete(self.state, inpt)
        except:
            # print("Non-index delete not implemented yet")
            raise

    def random_remove_val(self, val):
        try:
            # print("here")
            index = npr.choice(np.argwhere(self.state == val).T[0])
            # print("random remove index: ", index)
            while(index != 0):
                index = npr.choice(np.argwhere(self.state == val).T[0])
            # print("random remove index: ", index)
            self.state = np.delete(self.state, index)
        except:
            print("heres what was going on")
            print("val: ", val)
            # print("index: ", index)
            # print("state: ", self.state)
            raise


# ==============================================================================
# Basic Operations
# ==============================================================================
    def get(self, index):
        return copy.deepcopy(self.state[index])

    def add(self, index, value):
        self.state[index] = self.state[index]+value

    def sub(self, index, value):
        self.state[index] = self.state[index]-value
        return self.state[index]

    def remove(self, index):
        self.state = np.delete(self.state, index)

    def mul(self, index, value):
        self.state[index] = self.state[index]*value

    def exp(self, index, value):
        self.state[index] = self.state[index]**value

    def func(self, index, fun, *args, **kwargs):
        if callable(fun):
            self.state[index] = fun(self.state[index], *args, **kwargs)

    def add_state_variable(self, val):
        try:
            self.state = np.append(self.state, val)
        except:
            print("add_state: failed to append")
# ==============================================================================
# Arithmetic and basic functions
# ==============================================================================

    def sum(self, shift=1, condition=None):
        if condition is not None:
            try:
                return sum([j for j in self.state[shift:] if j >= condition])
            except:
                try:
                    return sum([j for j in self.state[shift:] if condition(j)])
                except:
                    print("cant apply condition, trying regular sum")
        else:
            return np.sum(self.state[shift:])

# ==============================================================================
# functional jawns
# ==============================================================================
    def vectorize_function(self, func, *args, **kwargs):
        return np.apply_along_axis(func, 0, self.state, args, kwargs)

    def scalar_function(self, func, *args, **kwargs):
        return sum([func(i, *args, **kwargs) for i in self.state])

    def vectorize_non_zero(self, func, *args, **kwargs):
        indices = self.non_zero()
        return (indices, np.apply_along_axis(func, 0, self.state[indices], indices, args, kwargs))

    def vectorize_nonzero(self, func, *args, **kwargs):
        indices = self.non_zero()
        return func(self.state, indices, *args, **kwargs)

    def conditional_indices(self, func, condition, *args, **kwargs):
        pass
# ==============================================================================
# statistics
# ==============================================================================

    def moment(self, n, shift=0):
        return scps.moment(self.state[shift:], moment=n)

    def skew(self, shift=1):
        """Calculate skew

        returns skew from scipystats

        :param shift: If first shift-1 shouldn't count toward skew, start at shift; defaults to 0
        :param shift: int, optional
        :return: skew
        :rtype: float
        """

        return scps.skew(self.state[shift:])
    def histogram(self, shift=1):
        return np.histogram(self.state[shift:],bins=99,range=(1,100),)
    def kurtosis(self, shift=1):
        """Kurtosis

        scipy stats kurtosis

        :param shift: calculate starting at index 0+shift, defaults to 0
        :param shift: int, optional
        :return: Kurtosis
        :rtype: float
        """

        return scps.kurtosis(self.state[shift:])

    def coefficient_of_variation(self, shift=0):
        return scps.variation(self.state[shift:])
# ==============================================================================
# dev
# ==============================================================================

    def dev_info(self, **kwargs):
        for key, val in kwargs.items():
            print(str(key)+": "+str(val))


if __name__ == '__main__':
    sv_initial = StateVector(1, [100])
    print(sv_initial)
    sv_initial+1
    print(sv_initial)
    sv_initial[1] -= 1
    print(sv_initial)
