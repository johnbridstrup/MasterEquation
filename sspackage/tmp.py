import numpy as np
import pandas as pd
from bisect import bisect
import KMC_Analysis as kanal
import simulations as sim
import sys
import importlib as port
utils=port.import_module('utils.utes')
from functools import wraps
import matplotlib
matplotlib.use("agg")
import matplotlib.pyplot as plt
import glob



# def holder(*args):
#     print("placeholder until {} is implemented".format(str(args)))
# class Loop:
#     COMMANDS={"print":holder,"execute":holder,"info":holder,"add":holder}
#     def __init__(self,f):
#         self.funcs={"main":f}
#         self.names=[]
#     def add_func(self,f,name,*args,**kwargs):
#         if not any([(i==name) for i in self.names]):
#             self.funcs[name]=f
#             self.args[name]=args
#             self.kwargs[name]=kwargs
#     def call(self,names):
#         for name in names:
#             if name in self.names:
#                 return self.funcs[name](*self.args[name],**self.kwargs[name])
#     def __call__(self,*args,**kwargs):
#         print("")
        
# def plot(x,y):
#     fig=plt.figure()    


runs_name="freq_test2"
path=utils.wd()+'/results/'+runs_name+'/'
files=glob.glob(path+'**/*.json')
# print(files)
files_gen=(i for i in files)
for i in files_gen:
    j=i
    # print(j)
    df=pd.read_json(j)
    dff=kanal.df_col2array(df,'polymers')
    histo_df=pd.DataFrame()
    # print(dff)
    t=kanal.array_ify(df['t'])

    

    plt.figure()
    num=4
    plt.plot(t,[k[1:num] if len(k)>num else len(k) for k in dff])
    # print(j[-4:]+'.png')
    plt.savefig(j+'.png')
    # plt.close()
# A=kanal.Analysis(path)
# A.load()

# a_keys=list(A.data.keys())
# AA=kanal.df_col2array(A.data[a_keys[0]],'polymers',max_list=400)
# t=kanal.array_ify(A.data[a_keys[0]]['t'])




