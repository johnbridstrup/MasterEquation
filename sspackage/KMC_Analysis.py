import statevector as sv
import KMC
import kernels as k
import simulations
import numpy as np
import pandas as pd
import sys
import os
import importlib as port
import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot as plt
from glob import glob
utils = port.import_module('utils.utes')


def none_max(inp,maxx=0):
    return max(inp,maxx)
def ext_with_zeroes(l,mx):
    L=max([len(i) for i in l])
    return [i.extend([0] * (none_max(L,mx) - len(i))) for i in l]
def df_col2array(df,col,index_sort=True, max_list=10):
    tmp=[df[col][index] for index in range(len(df[col]))]
    ext_with_zeroes(tmp,max_list)
    return np.array(tmp)    
def array_ify(unordered_df, maxx=None):
    try:
        larr=[unordered_df[i] for i in range(len(unordered_df))]
       
        L=max([len(i) for i in larr])
       
        [i.extend([0] * (none_max(L,maxx) - len(i))) for i in larr]
        return larr
    except:
        larr=[unordered_df[i] for i in range(len(unordered_df))]
        
        return larr
def main(fn):
    pass

def plotter(z):
    pdat=array_ify(z['polymers'])
    t=array_ify(z['t'])
    plt.figure()
    for j in range(1,5): 
        plt.plot(t,[i[j] for i in pdat])
    plt.show()

def tmax(a,b,c):
    return max(max(max(list(a)),max(list(b))),max(list(c)))
def maxspace(a,b,c,N):
    return np.linspace(0,tmax(a,b,c),N)

def list_time(df):
    return array_ify(df['t'])
def prepare_data(path,tbins=50):
    output={}
    df=pd.read_json(path)
    output['data']=df
    output['polymers']=array_ify(df['polymers'])
    output['t']=list_time(df)
    return output
def mass(polies,nc=2):
    return np.sum(polies[nc:])
def funcs(x='tbd'):
    pass
def all_subdirs_of(b='.'):
    result = []
    for d in os.listdir(b):
        bd = os.path.join(b, d)
        if os.path.isdir(bd): result.append(bd)
    return result

class Analysis:
    MAX_LOADED_DATA = 30
    def __init__(self,path=None):
        self.data={}
        if path is None:
            try:
                self.path=max(all_subdirs_of(path),key=os.path.getmtime)
            except:
                print("no results in this piece, you can keep working tho")
        else:
            self.path=path
        self.files=glob(path+'/*.json')
    def filenames(self):
        
        return(self.files)
    def load(self,**kwargs):
        if kwargs:
            print("If it doesn't do you what want, just know the only keywords that do anything here are \"paths\" and load_funcs (returning tuple with (func,args")
            self.kwdata={}
            for path in kwargs["paths"]:
                try:
                    self.kwdata[os.path.basename(path)]=pd.read_json(path)
                except:
                    pass
            for func in kwargs["load_funcs"]:
                try:
                    func[0](func[1])
                except:
                    pass
        for num,f in enumerate(self.files):
            self.data[f]=pd.read_json(f)
            if num>=Analysis.MAX_LOADED_DATA:
                break
        
if __name__=="__main__":
    outs=prepare_data('results/burp/burp0.json')
    #print(outs['t'])
    outs['polymers']
    # result_files=glob('results/burp/*.json')
    # result_full_paths=[utils.wd()+'/'+i for i in result_files]
    # print(result_full_paths)
    # path='results/burp/burp0.json'
    # main(path)
    # print('argsssssss',sys.argv)
    # z=pd.DataFrame()
    # for i in result_full_paths:    
    #     z=pd.read_json(i)
    #     pdat=array_ify(z['polymers'])
    #     t=array_ify(z['t'])
