import pandas as pd
import numpy as np
from pathlib import Path
import os
from os.path import splitext
import plistlib as pll
import matplotlib
matplotlib.use("agg")
import matplotlib.pyplot as plt
# import collections.abc as cabc
# import astropy as ap
import importlib as port
# from contextlib import ExitStack as eStck
utils = port.import_module('utils.utes')
# kanal = port.import_module('KMC_Analysis')

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


class Data_Directory:
    def __init__(self,dire=None, **kargs):
        self._Paths={}
        self._results_directories={}
        self._directories={}
        self._filepaths={}
        self._paths={}
        self.CWD=Path.cwd()
        self.dirs=[]
        self.active_name=None
        self.filepaths=[]
        self.top=None
        self.run_names=[]
        self._params={}
        if dire is not None:
            
            try:
                self.path = os.getcwd() +'/'+ dire
                self.Path=Path(self.path)
                self.path_set=True
                self.Top=dire
            except:
                try:
                    self.Path=Path(dire)
                    self.path=dire
                    self.path_set=True
                except:
                    self.path_set=False
        else:
            self.path_set=False
        
        self.full_init(names=kargs['names'])
    def full_init(self,**args):
        try:
            self.set_resultspath_subdir(args['names'])
        except:
            for i in args['names']:
                self.new_run(i)
        
        self.get_resultspath_directories()
        self.add_dir_to_active()
        self.enum_param_sets()
        self.store_filepaths()
        
    def set_resultspath_subdir(self,res_sub='no_leading_forwardslash/'):
        try:
            self.active_name=splitext(res_sub)[0]
        except:
            self.active_name=res_sub
        self.run_names.append(self.active_name)
        if not self.path_set:
            self.Path=Path('./resultspath/').joinpath(res_sub)
            self.path+'/'+res_sub+'/'
        else:
            self.Path=self.Path.joinpath(res_sub)
            self.path+=res_sub
        self._Paths[self.active_name]= self.Path
        self.path=str(self.Path)
        self._paths[self.active_name]=str(self.Path)
        print("get_resultspath_directories to get all subdirs of resultspath/yourdir")
    def new_run(self,runname):
        self.active_name=runname
        self.Path=self.CWD.joinpath('results',runname)
        self._Paths[self.active_name]=self.Path
        self.path=str(self.Path)
        self._paths[self.active_name]=str(self.Path)
        self.run_names.append(runname)
    def get_resultspath_directories(self):
        self.resultspath_directories = [x for x in self.Path.iterdir() if x.is_dir()]
        self._results_directories[self.active_name]=self.resultspath_directories
        print("obj.enum_resultspath_directories() to see the subdirectory list")
    def print_resultspath_directories(self):
        for i in self.resultspath_directories:
            outp=i.parts[-1]
            print(outp)
    def add_dir_to_active(self,the_dir='all'):
        if the_dir=='all':
            self.dirs=[]
            for i in self.resultspath_directories:
                self.dirs.append(i)
            self._directories[self.active_name]=self.dirs
            self._filepaths[self.active_name]={}
    def enum_param_sets(self):
        #self.param_sets={}
        for i in self.dirs:
            """getting the params and info
            
            """
            # check = False
            tmpParams=pll.readPlist(str(i)+'/input.data')

            self._params[self.active_name]=tmpParams
            self._params[self.active_name]['path']=str(i)
    def store_filepaths(self,name=None):
        self._filepaths[self.active_name]={}
        for ind,i in enumerate(self.dirs):
            if name is None:
                self.filepaths.append([x for x in list(i.iterdir()) if '.json' in str(x)])
                self._filepaths[self.active_name][ind]=self.filepaths[-1]
    def _make_master_dict(self):
        self.master={'names':self.run_names}
        for i in self.run_names:
            self.master[i]=None
            self.master[i]={'top_directory':self._Paths[i],'set_directory':self._results_directories[i],'run_directory':self._directories[i],'filepath':self._filepaths[i],'parameters':self._params[i]}
                    
    #~ class now can return DICT of all filepaths.
    #~ Next is class that can open and handle files that they point to
    #~ Nexter is one that formats the data
    #~ Nextest is one that plots/calulates/whatev
class DataHandler:
    def __init__(self,data):
        self.the_data=data.master
        self.sets=data.master['names']
    def list_files(self):
        print(self.sets) 
    def calculate_volume(self):
        for i in self.sets:
            prms=self.the_data[i]['parameters']
            prfx=1
            if prms['proteins']['concentration']['prefix']=='micro':
                prfx==10**-6
            conc=prms['proteins']['concentration']['value']
            M=prms['proteins']['monomers']
            NA=6.022*10**23
            vol=M/(conc*NA)
            print(vol,'L')
            self.the_data[i]['parameters']['proteins']['volume']=vol
    def split_dict(self):
        self.parameters_dict={}
        self.proteins_dict={}
        self.simulation_dict={}
        self.system_dict={}
        self.rates_dict={}
        for i in self.sets:
            self.parameters_dict=self.the_data[i]['parameters']
            self.proteins_dict=parameters[i]['proteins']
            self.simulation_dict[i]=parameters[i]['simulation']
            self.system_dict[i]=parameters[i]['system']
            self.rates_dict[i]=parameters[i]['rates']
    def calculate_rate_constants(self):
        for i in self.sets:
            prfx=1
            if self.proteins_dict[i]['concentration']['prefix']=='micro':
                prfx==10**-6
            conc=prfx*self.proteins_dict[i]['concentration']['value']
            nc=self.proteins_dict[i]['nucleus']
            M=self.proteins_dict[i]['monomers']
            a=rates_dict[i][0]
            k=rates_dict[i][2]
            self.ka=a*conc/M
            print(self.ka)
            self.kn=k*(conc/M)**(nc-1)
            print(self.kn)
    def open_files(self):
        self.data={}
        for i in self.sets:
            self.data[i]={}
            for key,dic in self.the_data[i]['filepath'].items():
                self.data[i][key]=[]
                for file in dic:
                    #print(file)
                    self.data[i][key].append(pd.read_json(file))
    def get_lists(self,name=None,number=None,dat='polymers'):
        self.poly={}
        self.t={}
        for name in self.sets:
            self.poly[name]={}
            self.t[name]={}
            for key,df_list in self.data[name].items():
                self.poly[name][key]=[]
                self.t[name][key]=[]
                for df in df_list:
                    self.poly[name][key].append(df['polymers'])
                    self.t[name][key].append(df['t'])
    def fix_polies(self):
        for i in self.sets:
            for key,arr_list in self.poly[i].items():
                outv=[]
                outt=[]
                for ind,arr in enumerate(arr_list):
                    ext_with_zeroes(arr,5)
                    outv.append([arr[index] for index in range(len(arr))])
                    outt.append([self.t[i][key][ind][index] for index in range(len(arr))])
                self.poly[i][key]=outv
                self.t[i][key]=outt
    def plot_all(self):
        for i in self.sets:
            for key,arr in self.poly[i].items():
                t=self.t[i][key]
                #print(arr)
                plt.figure()
                for index,j in enumerate(arr):
                    plt.plot(t[index],[y[1] for y in j])
                    plt.savefig(str(self.the_data[i]['run_directory'][key])+'/'+i+'{}.png'.format(index))
                    plt.close()
        

if __name__=='__main__':
    data_files=Data_Directory('/results/',names='schoot1.brid')
    data_files._make_master_dict()
    data=DataHandler(data_files)
