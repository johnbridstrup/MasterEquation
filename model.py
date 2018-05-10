#> model.py contains the machinery for building the model.
    #? What is a model?
    ##def:.. A Model :..
    #def: Describes what we observe (Aspect: State) 
    #def: Predicts what we might observe later (Aspect: Propensity)
    #def: Explains the observations  (Aspect: Transition)
        #`THEY WORK AT ANY RESOLUTION`#
        #   #^ System has: Model, Engine, Observables (Aspects)
        #   #   #/ State has: Variables, Units
        #   #   #   #~ Etc....
#> A model has mechanisms
#>  > A mechanism has aspects
#>  >   > Aspects have descriptions or procedures
from abc import ABC, ABCMeta, abstractmethod, abstractproperty

class Aspect(ABC):
    pass
    #~ the building blocks of mechanisms
class A
