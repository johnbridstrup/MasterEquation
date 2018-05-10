import numpy as np
import pandas as pd
import astropy as ap
import importlib as port
from copy import copy
import matplotlib
import matplotlib.pyplot as plt
utils = port.import_module('utils.utes')
kanal = port.import_module('KMC_Analysis')
def trymin(data):
  try:
    return min(data)
  except:
    try:
      return data.minx()
    except:
      try:
        return trymin([min(i) for i in data])
      except:
        raise TypeError("can't take the max of them jawns")
def trymax(data):
  try:
    return max(data)
  except:
    try:
      return data.max()
    except:
      try:
        return trymax([max(i) for i in data])
      except:
        raise TypeError("can't take the max of them jawns")
def recur_max(datmax):
  try:
    len(datmax)
    return recur_max(max(datmax))
  except:
    return datmax # recursion is fun

class DasBinner:
  @staticmethod
  def _t_shift(d):
    return d['t']-min(d['t'])
  def __init__(self, df1, *dfs, **labeled_args):
    self.tdata=[]
    self.dtdata=[]
    self.tsdata=[]
    self.ydata=[]
    self.monomer_data=[]
    self.sets=0
    try:
      self.tdata.append(kanal.array_ify(df1['t']))
      self.ydata.append([i[1:] for i in kanal.array_ify(df1['polymers'])])
      self.monomer_data.append([i[0] for i in kanal.array_ify(df1['polymers'])])
      self.dtdata.append(kanal.array_ify(df1['t_steps']))
      self.tsdata.append(kanal.array_ify(DasBinner._t_shift(df1)))
      self.sets+=1
    except:
      try:
        for df in df1:
          self.tdata.append(kanal.array_ify(df['t']))
          self.ydata.append([i[1:] for i in kanal.array_ify(df['polymers'])])
          self.monomer_data.append([i[0] for i in kanal.array_ify(df['polymers'])])
          self.dtdata.append(kanal.array_ify(df['t_steps']))
          self.tsdata.append(kanal.array_ify(DasBinner._t_shift(df)))
          self.sets+=1
      except:
        print("derp")
    if dfs:
      for df in dfs:
        self.tdata.append(kanal.array_ify(df['t']))
        self.ydata.append([i[1:] for i in kanal.array_ify(df['polymers'])])
        self.monomer_data.append([i[0] for i in kanal.array_ify(df['polymers'])])
        self.dtdata.append(kanal.array_ify(df['t_steps']))
        self.tsdata.append(kanal.array_ify(DasBinner._t_shift(df)))
        self.sets+=1
  def add_data(self,*dfs):
    try:
      for df in dfs:
        self.tdata.append(kanal.array_ify(df['t']))
        self.ydata.append([i[1:] for i in kanal.array_ify(df['polymers'])])
        self.monomer_data.append([i[0] for i in kanal.array_ify(df['polymers'])])
        self.dtdata.append(kanal.array_ify(df['t_steps']))
        self.tsdata.append(kanal.array_ify(DasBinner._t_shift(df)))
        self.sets+=1
    except:
      print("you had literaly one job")
  def bin(self):
    self.tbins=np.linspace(0,recur_max(self.tdata),100)
    self.ybinned=[0] * 100
    ybinned=[0] * 100
    bincount=1
    print(self.tbins)
    for i in self.ydata:
      bincount=1
      tr=self.tbins[bincount]
      for j in i:
        while j>tr:
          bincount+=1
          tr=self.tbins[bincount]
        ybinned[bincount-1]+=j
    self.ybinned=[jj/self.sets for jj in ybinned]
    print(self.ybinned)
  def time_series(self,name='sum'):
     
          
if __name__=='__main__':
  data=[]
  data.append(pd.read_json('results/freq_test2/freq_test2_aa_100.0_kk_1e-06/freq_test2_aa_100.0_kk_1e-06_1.json'))
  data.append(pd.read_json('results/freq_test2/freq_test2_aa_100.0_kk_1e-06/freq_test2_aa_100.0_kk_1e-06_2.json'))
  data.append(pd.read_json('results/freq_test2/freq_test2_aa_100.0_kk_1e-06/freq_test2_aa_100.0_kk_1e-06_3.json'))
  data.append(pd.read_json('results/freq_test2/freq_test2_aa_100.0_kk_1e-06/freq_test2_aa_100.0_kk_1e-06_4.json'))
  dater=DasBinner(data)
  plt.figure()  
  for x,y in zip(dater.tsdata,dater.ydata):
    plt.plot(x,[i[1] for i in y])
  plt.show()
  dater.bin()
        
