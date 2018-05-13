import numpy as np
import pandas as pd
import astropy as ap
import importlib as port
from copy import copy
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns
utils = port.import_module('utils.utes')
kanal = port.import_module('KMC_Analysis')

class DataAnalysis:
    def __init__(self, data_binner, *args, **kwargs):
        pass
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
        return datmax  # recursion is fun
def none_max(inp,maxx=0):
    return max(inp,maxx)
def ext_with_zeroes(l,mx=500):
    L=max([len(i) for i in l])
    return [i.extend([0] * (none_max(L,mx) - len(i))) for i in l]
class DasBinner:
    @staticmethod
    def _t_shift(d):
        return d['t']-min(d['t'])

    def __init__(self, df1, *dfs, **labeled_args):
        self.tdata = []
        self.dtdata = []
        self.tsdata = []
        self.ydata = []
        self.monomer_data = []
        self.sets = 0
        self.binned = False
        self._sorted=False
        try:
            self.tdata.append(kanal.array_ify(df1['t']))
            self.ydata.append([i[1:]
                               for i in kanal.array_ify(df1['polymers'])])
            self.monomer_data.append(
                [i[0] for i in kanal.array_ify(df1['polymers'])])
            self.dtdata.append(kanal.array_ify(df1['t_steps']))
            self.tsdata.append(kanal.array_ify(DasBinner._t_shift(df1)))
            self.sets += 1
        except:
            try:
                for df in df1:
                    self.tdata.append(kanal.array_ify(df['t']))
                    self.ydata.append([i[1:]
                                       for i in kanal.array_ify(df['polymers'])])
                    self.monomer_data.append(
                        [i[0] for i in kanal.array_ify(df['polymers'])])
                    self.dtdata.append(kanal.array_ify(df['t_steps']))
                    self.tsdata.append(kanal.array_ify(DasBinner._t_shift(df)))
                    self.sets += 1
            except:
                print("derp")
        if dfs:
            for df in dfs:
                self.tdata.append(kanal.array_ify(df['t']))
                self.ydata.append([i[1:]
                                   for i in kanal.array_ify(df['polymers'])])
                self.monomer_data.append(
                    [i[0] for i in kanal.array_ify(df['polymers'])])
                self.dtdata.append(kanal.array_ify(df['t_steps']))
                self.tsdata.append(kanal.array_ify(DasBinner._t_shift(df)))
                self.sets += 1

    def add_data(self, *dfs):
        try:
            for df in dfs:
                self.tdata.append(kanal.array_ify(df['t']))
                self.ydata.append([i[1:]
                                   for i in kanal.array_ify(df['polymers'])])
                self.monomer_data.append(
                    [i[0] for i in kanal.array_ify(df['polymers'])])
                self.dtdata.append(kanal.array_ify(df['t_steps']))
                self.tsdata.append(kanal.array_ify(DasBinner._t_shift(df)))
                self.sets += 1
                self.binned = False
        except:
            print("you had literaly one job")
    # def all_time(self):
    #     t_unsorted=[[] for _ in range(len(self.ydata))]
    #     [t_unsorted[set_index].append(t,data_index) for set_index,i in enumerate(self.tdata) for data_index,t in enumerate(i)]
    #     self.t_sorted=sorted(self.t_unsorted)
    #     self._sorted=True
    def all_time(self):
        self.t_unsorted=[(t,(set_index,data_index)) for set_index,i in enumerate(self.tdata) for data_index,t in enumerate(i)]
        self.t_sorted=sorted(self.t_unsorted)
        self._sorted=True
    def bin_time(self, binss='get',scale='linear'):
        if binss=='get':
            if len(self.tdata)<10:
                binss=int(len(self.tdata)*max([len(i) for i in self.tdata])/10)
            else:
                binss=200
        tbin=np.linspace(0,max(max(self.tdata)),binss)
        self.t_binned=np.linspace(0,max(max(self.tdata)),binss)
        self.data_binned=[[[] for _ in range(binss)] for __ in range(len(self.ydata))]
        if not self._sorted:
            self.all_time()
        bincount=0
        for i in self.t_sorted:
            while i[0] > tbin[bincount+1]:
                bincount+=1
            self.data_binned[i[1][0]][bincount].append((i[0],self.ydata[i[1][0]][i[1][1]]))
    @staticmethod
    def set_gen(data):
        yield from data
    def data_gen(self, data):
        self.y_gens=[DasBinner.set_gen(dset) for dset in DasBinner.set_gen(self.data_binned)]
    def series_gen(self):
        for a_gen in self.y_gens:
            y,t=zip(*[(i[1],i[0]) for j in a_gen for i in j if j!=[]])
            yield (y,t)
    
    def mass_gens(self,restart=True,beg=0):
        if restart:
            self.series_gen()
        y,t=next(self.series_gen())
        m=(sum(i) for i in y)
        yield (m,t)
    def length_gens(self):
        y,t=next(self.series_gen())
        P=(len(i) for i in y)
        yield (P,t)
    def mass_dev(self):
        
    def _dataSelector(self, name='ydata',seti=0):
        try:
            dat=self.__dict__[name][seti]
            yield from dat
            #? Yields a specific data set one value at a time
        except:
            if set is None:
                dat=self.__dict__[name]
                yield from dat
                #? Yields one full data set at a time
    def bin(self, binss=100):
        #self.tbins=np.linspace(0,recur_max(self.tdata), binss)
        self.bins=binss
        self.ymax = recur_max(self.ydata)
        ##
        #~ Stored variables
        ##
        self.tbins = np.linspace(0, recur_max(self.tdata), binss)
        self._mass = np.zeros(binss)
        self.stdev = np.zeros(binss)
        self._length = np.zeros(binss)
        self._number= np.zeros(binss)
        
        ##
        #~ Placeholder variables
        ##
        _mass = np.zeros(binss)
        _m2 = np.zeros(binss)
        _number = np.zeros(binss)
        _n2 = np.zeros(binss)
        _length = np.zeros(binss)
        _l2 = np.zeros(binss)
        sum_y2binned = np.zeros(binss)
        self.list_y2binned = np.zeros(binss)
        for i1,l1 in enumerate(self.tdata):
            for i2,l2 in enumerate(l1):
                bincount = 0
                tr = self.tbins[bincount]
                while l2 > tr:
                    bincount += 1
                    tr = self.tbins[bincount]
                _mass[bincount-1] += sum(self.ydata[i1][i2])
                _m2[bincount-1] += _mass[bincount-1]**2
                _number[bincount-1] += len(self.ydata[i1][i2])
                _n2[bincount-1] += _number[bincount-1]**2
                _length[bincount-1] += _mass[bincount-1]/_number[bincount-1]
                _l2[bincount-1] += _length[bincount-1]**2
                sum_y2binned[bincount-1] += sum([y**2 for y in self.ydata[i1][i2]])
                self.list_y2binned = [y**2 for y in self.ydata[i1][i2]]
                
        for jj, val in enumerate(_mass):
            if _number[jj]>0:
                self._mass[jj] = val/_number[jj]
                self.stdev[jj] = np.sqrt(sum_y2binned[jj]/_number[jj]-self._mass[jj]**2)
            else:
                self._mass[jj]=0
                self.stdev[jj]=0
        
        print(self.stdev)

        self.binned = True
    
    
    def time_series(self, name='sum'):
        types = {'mass': 0, 'number': 0, 'deviation': 0}
        if self.binned:
            print("do things here later")
    def time(self,binned=True,sorted=True):
        if binned:
            return self.data_gen(self.t_binned)
        elif sorted:
            self.data_gen(self.t_sorted)
    
            
            

if __name__ == '__main__':
    data = []
    data.append(pd.read_json(
        'results/runs_test2/runs_test2_aa_1e-05_kk_1e-09/runs_test2_aa_1e-05_kk_1e-09_1.json'))
    data.append(pd.read_json(
        'results/runs_test2/runs_test2_aa_1e-05_kk_1e-09/runs_test2_aa_1e-05_kk_1e-09_2.json'))
    data.append(pd.read_json(
        'results/runs_test2/runs_test2_aa_1e-05_kk_1e-09/runs_test2_aa_1e-05_kk_1e-09_3.json'))
    data.append(pd.read_json(
        'results/runs_test2/runs_test2_aa_1e-05_kk_1e-09/runs_test2_aa_1e-05_kk_1e-09_4.json'))
    dater = DasBinner(data)
    # plt.figure()
    # for x, y in zip(dater.tsdata, dater.ydata):
    #     plt.plot(x, [i[1] for i in y])
    # plt.show()
    dater.bin()
