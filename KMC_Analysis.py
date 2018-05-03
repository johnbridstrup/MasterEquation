import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import statevector as sv
import glob



if __name__=="__main__":
    file='run1.json'
    path='results/'
    





    fn=path+file
    data=pd.read_json(fn)

    

    print([i for i in data['polymers']])