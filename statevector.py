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
import itertools as itt
from abc import ABCMeta,ABC,abstractmethod


# ==============================================================================
# RANDOM FUNCTIONS
# ==============================================================================

class InputException(Exception):
    def __init__(self,node=None,msg=None):
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
                self.state=np.array(the_input.state)
                self.initial_conditions=the_input.initial_conditions
                self.delete_zeros=the_input.delete_zeros
            except:
                try:
                    if type(the_input) == int:
                        raise InputException()
                    self.state = np.array(the_input)
                except (TypeError, IndexError, InputException) as e:
                    pass
                if type(the_input) == int:
                    self.state=np.zeros(the_input)
                else:
                    self.state=the_input
                    self.delete_zeros=delete_zeros
                if initial_conditions is not None:
                    try:
                        assert(any([i!=0 for i in self.state]))
                    except Exception as e:
                        pass
                    try:
                        for i,v in enumerate(initial_conditions):
                            self.state[i]=v
                    except ValueError as e:
                        self.state[0]=initial_conditions
                        pass
            
            else:
                if initial_conditions is not None:
                    try:
                        self.state=np.array(initial_conditions)
                    except:
                        print("its some weird, not array-ifiable type")
                        raise 
            self.delete_zeros=delete_zeros
    def __repr__(self):
        return repr(self.state)
    def __str__(self):
        return str(self.state)
    
    def __getitem__(self, index):
        try:
            return self.state[index]
        except (TypeError,IndexError) as e:
            print("No string get implementiation, you must use an index")
            raise e
    def __setitem__(self,key, val):
        try:
            self.state[key]=val
        except (IndexError,KeyError) as e:
            raise e
# ==============================================================================
# Getting and setting elements
# ==============================================================================
    def get_index_by_value(self, val, random=False):
        if random==True:
            return npr.choice(np.argwhere(self.state==val).T[0])
        else:
            return np.argwhere(self.state==val.T[0])

    def set_state_variable(self, index, value):
        try: 
            self.state[index]=value
        except:
            print("Tried to access and set by index value")
            print("index: ", index)
            print("value: ", value)
            print("state: ", self.state)
            raise
    def random_set_by_value(self,old_value, new_value,add=False):
        index=self.get_index_by_value(old_value,random=True)
        if add!=True:
            try:
                self.set_state_variable(index,new_value) 
            except:
                self.dev_info(state=self.state,old_val=old_value,new_value=new_value,index=index)
        else:
            self.state[index]+=new_value
            if self.delete_zeros and self.state[index]<=0:
                try:
                    self.state.remove_state_variable(index)
                except:
                    self.dev_info(index=index,state=self.state)
                    raise
    def addOne(self,shift=1):
        index=npr.choice(range(len(self.state)))
        print("INDEX",index)
        self.state[index]=self.state[index]+1
        if shift==1:
            self.state[0]=self.state[0]-1
    def subOne(self,nc,shift=1):
        index=npr.choice(range(len(self.state[1:])))
        self.state[index]=self.state[index]-1
        if self.state[index]<nc:
            self.state[0]=self.state[0]+self.state[index]
            self.remove_state_variable(index)
        if shift==1:
            self.state[0]=self.state[0]+1
    def breakOne(self,nc):
        index=npr.choice([i for i,v in enumerate(self.state) if v>3 and i>0])
        print("index",index)
        try:
            r=npr.choice(np.array(list(range(4,self.state[index]))))
            print(r,"r")
            self.state[index]=self.state[index]-r
            if self.state[index]<nc:
                self.state[0]=self.state[0]+self.state[index]
                self.remove_state_variable(index)    
            if r>nc:
                self.state=self.state+r
            else:
                self.state[0]=self.state[0]+r
        except:
            print("WHAT THE FUCKING UCK YOU FUCK")
        
        
    def merge(self):
        r=npr.choice(self.state[1:])
        print(r,"merrrge")
        self.random_remove_val(r)
        s=npr.choice(range(len(self.state[1:])))+1
        self.state[s]=self.state[s]+r
    def non_zero(self):
        return np.nonzero(self.state)[0]
# ==============================================================================
# Adding and removing elements
# ==============================================================================   
    def remove_state_variable(self, inpt):
        try:
            self.state = np.delete(self.state,inpt)
        except:
            print("Non-index delete not implemented yet")
            raise
    def random_remove_val(self,val):
        try:
            print("here")
            index=npr.choice(np.argwhere(self.state==val).T[0])
            print("random remove index: ",index)
            while(index!=0):
                index=npr.choice(np.argwhere(self.state==val).T[0])
            print("random remove index: ",index)
            self.state = np.delete(self.state,index)
        except:
            print("heres what was going on")
            print("val: ", val)
            # print("index: ", index)
            # print("state: ", self.state)
            raise
    
    def add_state_variable(self,val):
        try:
            self.state = np.append(self.state,val)
        except:
            print("add_state: failed to append")
# ==============================================================================
# Arithmetic and basic functions
# ==============================================================================
    def sum(self,shift=0):
        try:
            self.dev_info(state=self.state,tryit="now we in here")
            return np.sum(self.state[shift:])
        except Exception as e:
            print("sum: Cant take sum of state (test)")
            print([i for i in self.state[shift:]])
            try:
                pass
            except:
                raise e
# ==============================================================================
# functional jawns
# ==============================================================================
    def vectorize_function(self,func,*args,**kwargs):
        return np.apply_along_axis(func,0,self.state,args,kwargs)
    def scalar_function(self, func, *args, **kwargs):
        return sum([func(i,*args,**kwargs) for i in self.state])
    def vectorize_non_zero(self,func,*args,**kwargs):
        indices=self.non_zero()
        return (indices,np.apply_along_axis(func,0,self.state[indices],indices, args, kwargs))
    def vectorize_nonzero(self,func,*args,**kwargs):
        indices=self.non_zero()
        return func(self.state,indices,*args,**kwargs)
    def conditional_indices(self,func,condition,*args,**kwargs):
        pass
# ==============================================================================
# statistics
# ==============================================================================
    def moment(self,n, shift=0):
        return scps.moment(self.state[shift:],moment=n)
    def skew(self, shift=0):
        """Calculate skew
        
        returns skew from scipystats
        
        :param shift: If first shift-1 shouldn't count toward skew, start at shift; defaults to 0
        :param shift: int, optional
        :return: skew
        :rtype: float
        """

        return scps.skew(self.state[shift:])
    def kurtosis(self,shift=0):
        """Kurtosis
        
        scipy stats kurtosis
        
        :param shift: calculate starting at index 0+shift, defaults to 0
        :param shift: int, optional
        :return: Kurtosis
        :rtype: float
        """

        return scps.kurtosis(self.state[shift:])
    def coefficient_of_variation(self,shift=0):
        return scps.variation(self.state[shift:])
# ==============================================================================
# dev
# ==============================================================================
    def dev_info(self,**kwargs):
        for key,val in kwargs.items():
            print(str(key)+": "+str(val))

if __name__ == '__main__':
    sv_initial = StateVector(1,[100])
    print(sv_initial)
    sv_initial+1
    print(sv_initial)
    sv_initial[1]-=1
    print(sv_initial)


                
            
        