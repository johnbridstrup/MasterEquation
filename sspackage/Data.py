import pandas as pd
import numpy as np
from pathlib import Path
import os
from os.path import splitext
import plistlib as pll
import matplotlib
#matplotlib.use('Agg', force=True)
import matplotlib.pyplot as plt
# import collections.abc as cabc
# import astropy as ap
import importlib as port
from copy import copy
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
        self._params['path']=[]
        self._params[self.active_name]=[]
        for i in self.dirs:
            
            """getting the params and info
            
            """
            # check = False
            tmpParams=pll.readPlist(str(i)+'/input.data')

            self._params[self.active_name].append(tmpParams)
            self._params['path'].append(str(i))
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
    def __init__(self,data,poly_length=5, bins=100):
        self.the_data=data.master
        self.bins=bins
        self._is_timeseries=False
        self.poly_length=poly_length
        self.sets=data.master['names']
    def list_files(self):
        print(self.sets) 
    def calculate_volume(self):
        self.vol={}
        for i in self.sets:
            self.vol[i]=[]
            for prms in self.the_data[i]['parameters']:
                prfx=1
                if prms['proteins']['concentration']['prefix']=='micro':
                    prfx==10**-6
                conc=prms['proteins']['concentration']['value']
                M=prms['proteins']['monomers']
                NA=6.022*10**23
                self.vol[i].append(M/(conc*NA))
                print(self.vol[i],'L')
                # self.the_data[i]['parameters']['proteins']['volume']=vol
    def split_dict(self):
        self.parameters_dict={}
        self.proteins_dict={}
        self.simulation_dict={}
        self.system_dict={}
        self.rates_dict={}
        for i in self.sets:
            self.parameters_dict[i]=self.the_data[i]['parameters']
            self.proteins_dict[i]=[j['proteins'] for j in self.parameters_dict[i]]
            self.simulation_dict[i]=[j['simulation'] for j in self.parameters_dict[i]]
            self.system_dict[i]=[j['system'] for j in self.parameters_dict[i]]
            self.rates_dict[i]=[j['rates'] for j in self.parameters_dict[i]]
    def calculate_rate_constants(self):
        self.kn={}
        self.ka={}
        for i in self.sets:
            self.kn[i]=[]
            self.ka[i]=[]
            for index,params in enumerate(self.proteins_dict[i]):
                prfx=1
                if params['concentration']['prefix']=='micro':
                    prfx==10**-6
                conc=prfx*params['concentration']['value']
                nc=params['nucleus']
                M=params['monomers']
                a=self.rates_dict[i][index][0]

                k=self.rates_dict[i][index][4]
                self.ka[i].append(a*conc/M)
                print(self.ka)
                self.kn[i].append(k*(conc/M)**(nc-1))
                print(self.kn)
    def open_files(self):
        self.data={}
        for i in self.sets:
            self.data[i]={}
            for key,dic in self.the_data[i]['filepath'].items():
                #print(key)
                self.data[i][key]=[]
                for file in dic:
                    #print(file)
                    self.data[i][key].append(pd.read_json(file))
                    #print(self.data[i][key][-1])
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
                    #ext_with_zeroes(arr,self.poly_length)
                    outv.append([arr[index] for index in range(len(arr))])
                    outt.append([self.t[i][key][ind][index] for index in range(len(arr))])
                self.poly[i][key]=outv
                self.t[i][key]=outt
    def make_master(self,nonan=True):
        self.master_data={}
        self.monomers={}
        if nonan:
            self.noNan_master_data={}
        for i in self.sets:
            self.master_data[i]={}
            self.monomers[i]={}
            if nonan:
                self.noNan_master_data[i]={}
            for j in self.data[i].keys():
                self.master_data[i][j]=[]
                self.monomers[i][j]=[]
                if nonan:
                    self.noNan_master_data[i][j]=[]
                for dff in self.data[i][j].copy():
                    dff['t']=pd.to_timedelta(dff['t'],unit='s')
                    dff.sort_values('t',inplace=True)
                    try:
                        dff.drop('t_steps',axis=1,inplace=True)
                    except:
                        print(dff['t_steps'])
                    ddff=pd.DataFrame(dff.polymers.values.tolist(),index=dff['t'].values).add_prefix('polymer_')
                    if nonan:
                        ddff_noNaN=ddff.fillna(0)
                    self.monomers[i][j].append(ddff['polymer_0'].copy())
                    ddff.drop('polymer_0',axis=1,inplace=True)
                    if nonan:
                        ddff_noNaN.drop('polymer_0',axis=1,inplace=True)
                    self.master_data[i][j].append(ddff)
                    if nonan:
                        self.noNan_master_data[i][j].append(ddff_noNaN)
    def make_series(self):
        self.series_list={}
        self.monomer_series={}
        for i in self.sets:
            self.series_list[i]={}
            self.monomer_series[i]={}
            for j in self.data[i].keys():
                self.series_list[i][j]=[]
                self.monomer_series[i][j]=[]
                for index,df in enumerate(self.data[i][j]):
                    sers=[]
                    sers_t=[]
                    df['t']=pd.to_timedelta(df['t'],unit='s')
                    df.sort_values('t',inplace=True)
                    for polies,tt in zip(df['polymers'].values,df['t'].values):
                        for indd,pol in enumerate(polies[:]):
                            try:
                                sers[indd].append(pol)
                                sers_t[indd].append(tt)
                            except:
                                sers.append([pol])
                                sers_t.append([tt])
                    self.series_list[i][j].append([pd.Series(data=pols,index=ts) for pols,ts in zip(copy(sers),copy(sers_t))])
                    self.monomer_series[i][j].append(self.series_list[i][j][index][0].copy())
                    del self.series_list[i][j][-1][0]
    def full_time(self):
        self.times={}
        for i in self.sets:
            self.times[i]={}
            for j in self.data[i].keys():
                self.times[i][j]=pd.concat([pd.DataFrame(ii.index.values) for ii in self.monomer_series[i][j]])
                self.times[i][j].sort_values(inplace=True,by=0)
                self.times[i][j].reset_index(inplace=True)
                self.times[i][j].drop('index',axis=1,inplace=True)
    def chunktime(self,n=100):
        self.chunk_time={}
        for i in self.sets:
            self.chunk_time[i]={}
            for j in self.data[i].keys():
                self.chunk_time[i][j]=[]
                for k in range(0,len(self.times[i][j].values),n):
                    self.chunk_time[i][j].append(self.times[i][j].values[k:k+n])
    def hist_generator(self,tbin=100,set=(0,0)):
        self.chunktime(tbin)
        run_name=self.sets[set[0]]
        paramset=set[1]
        L=self.proteins_dict[run_name][paramset]['monomers']
        
        for times in self.chunk_time[run_name][paramset]:
            hist_slice=np.zeros(L)
            pols=[]
            for t in times:
                
                for sers in self.series_list[run_name][paramset]:
                    for ser in sers:
                        try:
                            pols.append(ser.loc[t])
                            break
                        except:
                            pass
            hh=np.histogram(np.array(pols),bins=L-1,range=(1,L))
            for i in pols:
                print(i)
                hist_slice[i-1]+=1
            yield (hist_slice,hh)
                    
    def time_index(self):
        
        self.noNan_master_data={}
        self.series_data={}
        self.monomers={}
        self.series_list={}
        for i in self.sets:
            
            self.noNan_master_data[i]={}
            self.series_data[i]={}
            self.series_list[i]={}
            self.monomers[i]={}
            for j in self.data[i].keys():
                # self.master_data[i][j]=[]
                self.noNan_master_data[i][j]=[]
                self.monomers[i][j]=[]
                self.series_data[i][j]=[]
                self.series_list[i][j]=[]
                
                
                for index,df in enumerate(self.data[i][j]):
                    sers=[]
                    sers_t=[]
                    df['t']=pd.to_timedelta(df['t'],unit='s')
                    df.sort_values('t',inplace=True)
                    for polies,tt in zip(df['polymers'].values,df['t'].values):
                        for indd,pol in enumerate(polies[1:]):
                            try:
                                sers[indd].append(pol)
                                sers_t[indd].append(tt)
                            except:
                                sers.append([pol])
                                sers_t.append([tt])
                    self.series_list[i][j].append([pd.Series(data=pols,index=ts) for pols,ts in zip(sers,sers_t)])        
                    ser=pd.Series(data=df['polymers'].values,index=df['t'].values,name='polymers_{}'.format(index))
                    if ser.name=='polymers_0':
                        ser.rename('monomers')
                    #df.set_index('t',inplace=True)
                    try:
                        df.drop('t_steps',axis=1,inplace=True)
                    except:
                        print(df['t_steps'])
                    ddf=pd.DataFrame(df.polymers.values.tolist(),index=df['t'].values).add_prefix('polymer_')
                    # try:
                    #     df.rename(columns={'polymers':'polymers_{}'.format(index)},inplace=True)
                    # except Exception as e:
                    #     raise e
                    ddf_noNaN=ddf.fillna(0)
                    self.monomers[i][j].append(ddf['polymer_0'].copy())
                    ddf.drop('polymer_0',axis=1,inplace=True)
                    ddf_noNaN.drop('polymer_0',axis=1,inplace=True)
                    self.master_data[i][j].append(ddf)
                    self.noNan_master_data[i][j].append(ddf_noNaN)
                    self.series_data[i][j].append(ser)
    def check_t_overlap(self):
        for i in self.sets:
            for j in self.data[i].keys():
                pass
    #def concat_data(self):
  
        # for i in self.sets:
        #     self.master_data[i]={}
        #     for j in self.data_list[i]:
        #         print(j)
        #         self.master_data[i][j]=pd.concat(self.data_list[i][j],axis=1)
        #         self.master_data[i][j]['t']=self.master_data[i][j].index
                #self.master_data[i][j].set_axis(list(range(len(self.master_data[i][j]))),axis=0,inplace=True)
    # def re_index(self):
    #     for i in self.sets:
    #         for j in self.data_list[i]:
    #             self.master_data[i][j].set_axis(list(range(len(self.master_data[i][j]))),axis=0,inplace=True)
    def plot_all(self):
        for i in self.sets:
            for key,arr in self.poly[i].items():
                t=self.t[i][key]
                #print(arr)
                plt.figure()
                for index,j in enumerate(arr):
                    plt.plot(t[index],[y[2] for y in j])
                    plt.savefig(str(self.the_data[i]['run_directory'][key])+'/'+i+'{}_{}.png'.format(index,2))
                    plt.close()
    def min_maxt(self):
        self.max_t={}
        self.min_t={}
        for j in self.sets:
            self.max_t[j]={}
            self.min_t[j]={}
            for key, paramset in self.t[j].items():
                self.max_t[j][key]=max([i[-1] for i in paramset])
                self.min_t[j][key]=min([i[0] for i in paramset])

    # def bin_data(self, bins=None):
    #     if bins is None:
    #         bins=self.bins
    #     self.t_binned={}
    #     self.binned_data={}
    #     for i in self.sets:
    #         self.t_binned[i]={}
    #         self.binned_data[i]={}
    #         for key,paramset in self.poly[i].items():
                
    # def bin_average_data(self,bins=50):
    #     # self.histogram{}
    #     self.t_binned={}
    #     self.histogram={}
    #     for i in self.sets:
    #         self.histogram[i]={}
    #         self.t_binned[i]={}
            
    #         for key,paramset in self.poly[i].items():
    #             max_length=self.proteins_dict[i][int(key)]['monomers']
    #             self.t_binned[i][key]=np.linspace(self.min_t[i][key]*((bins-1)/bins),self.max_t[i][key],bins)
    #             self.histogram[i][key]=np.zeros((bins+1,max_length))
    #             for run_t,run in zip(self.t[i][key],paramset):
    #                 adder=np.zeros(max_length)
    #                 adder[0]=self.proteins_dict[i][int(key)]['monomers'] # number of monomers, number of polymers
    #                 bin_num=0
    #                 for t,pol_list in zip(run_t,run):
    #                     while self.t_binned[i][key][bin_num]<t:
    #                         if bin_num is 0:
    #                             self.histogram[i][key][0][0]+=self.proteins_dict[i][int(key)]['monomers']
    #                         else:
    #                             for iin,kval in enumerate(adder):
    #                                 self.histogram[i][key][bin_num][iin]+=kval
    #                                 #the first element is the number of monomers, the rest are the size of a single polymer
    #                         bin_num+=1
    #                     adder=np.zeros(max_length)
    #                     adder[0]=pol_list[0]
    #                     self.histogram[i][key][bin_num][0]+=pol_list[0]
    #                     for poly_index,poly_size in enumerate(pol_list[1:]):
    #                         self.histogram[i][key][bin_num][poly_size-1]+=1
    #                         adder[poly_index]+=1
    #             for bindex,hist_slice in enumerate(self.histogram[i][key]):
    #                 for poly_size,summed_value in enumerate(hist_slice):
    #                     self.histogram[i][key][bindex][poly_size]=summed_value/self.parameters_dict[i][int(key)]['simulation']['runs']
    def prep_data(self):
        self.split_dict()
        self.open_files()
        self.get_lists()
        self.fix_polies()
        self.make_master()
        self.make_series()
        self.full_time()
        self.chunktime()


class PolymerTimeSeries:
    def __init__(self,data):
        pass

class PolymerSet:
    def __init__(self,data=None):
        self.polymers=[]
        if data is not None:
            try:
                for i in data:
                    self.polymers.append(PolymerTimeSeries(i))
            except:
                self.polymers.append(PolymerTimeSeries(data))
            
            
                        
                    

if __name__=='__main__':
    data_files=Data_Directory('/results/',names='test1.brid')
    data_files._make_master_dict()
    data=DataHandler(data_files)
    data.prep_data()
    histgen=data.hist_generator()
