from glob import glob
from collections import defaultdict, namedtuple
import biplist as bp
import pandas as pd
from math import log10
import matplotlib
from collections.abc import Sequence

from shutil import copyfile

rate_set = namedtuple('RateSet', 'a b k')

NA = 6.022 * 10 ** 23
PREFIXES = {'femto': 10.0 ** -15, 'pico': 10.0 ** -12, 'Angstrom': 10.0 ** -10,
            'nano': 10.0 ** -9, 'micro': 10.0 ** -6, 'milli': 10.0 ** -3, 'centi': 10.0 ** -2,
            'kilo': 10.0 ** 3, 'mega': 10.0 ** 6, 'giga': 10.0 ** 9, 'tera': 10.0 ** 12}


def yield_from_helper(itr):
    try:
        yield from itr
    except StopIteration:
        pass


class DataFiles:
    def __init__(self, sim_name):
        self.sim_name = sim_name
        self.filepaths = glob("results/" + self.sim_name + ".brid/**")
        self.run_dict = defaultdict(list)
        self.set_strings = [i.split(self.sim_name + "_")[1].split(".strup")[0] for i in self.filepaths]
        self.rate_sets = defaultdict(rate_set)
        for ind, v in enumerate(self.set_strings):
            _, a, _, b, _, k = v.split('_')
            print(v)
            self.rate_sets[v] = rate_set(a, b, k)
        for i, v in enumerate(self.filepaths):
            for j in glob(v + '/*.json'):
                self.run_dict[self.set_strings[i]].append(j)
        self.input_files = {}
        self.load_inputs()

    def __getitem__(self, item):
        if item in self.__dir__():
            return getattr(self, item)

    def load_inputs(self):
        for i, v in enumerate(self.filepaths):
            self.input_files[self.set_strings[i]] = bp.readPlist(v + '/input.data')


PROPERTY_LIST = ['rate_constants', 'stochastic_rate_constants', 'volume', 'system_properties', 'simulation_properties', \
                 'protein_properties']


class DataProperties(DataFiles):
    def __init__(self, sim_name):
        super().__init__(sim_name)
        self.simulation_properties = []
        self.protein_properties = []
        self.system_properties = []
        self.rate_constants = [rate_set(float(i.a), float(i.b), float(i.k)) for i in
                               self.rate_sets.values()]
        self.stochastic_rate_constants = []

        for dct in self.input_files.values():
            self.simulation_properties.append(dct['system'])
            self.protein_properties.append(dct['proteins'])
            self.system_properties.append(dct['simulation'])
        for index, (sysprops, rset) in enumerate(list(zip(self.protein_properties, self.rate_constants))):
            cprops = sysprops['concentration']
            c0 = cprops['value'] * PREFIXES[cprops['prefix']]
            num = sysprops['monomers']
            nc = sysprops['nucleus']
            self.system_properties[index]['volume'] = num / (NA * c0)
            self.stochastic_rate_constants.append(rate_set(rset.a * c0 / num, rset.b, rset.k * (c0 / num) ** (nc - 1)))

    # def __getattr__(self, item):
    #     if item in PROPERTY_LIST:
    #         return getattr(self,item)
    #     else:
    #         raise KeyError("Key not in {}".format(PROPERTY_LIST))

    def __getitem__(self, item):
        if item in PROPERTY_LIST:
            if item is 'volume':
                return [i['volume'] for i in self.system_properties]
            else:
                return getattr(self, item)
        else:
            raise KeyError("Key is not in {}".format(PROPERTY_LIST))


class FileKeyHandler:
    def __init__(self, filepath=None, data_file_handler=None, param_set_index=None, label=None, file_set=None):
        self.file = False
        if filepath is not None:
            self.file = filepath
            if label is None:
                self.label = filepath[-15:-5]
        if param_set_index is not None:
            if data_file_handler is not None:
                self.file = data_file_handler.run_dict[data_file_handler.set_strings[param_set_index]]
            else:
                raise AttributeError("If using indices, a DataFiles class must be passed as data_file_handler")
        if label is not None:
            self.label = label
        if file_set is not None:
            self.file = file_set
            try:
                self.labels = [i[-15:-5] for i in file_set]
            except IndexError:
                self.labels = [i for i in file_set]
            self.label = None

        if not self.file:
            raise AttributeError("file did not get set correctly")

    def __hash__(self):
        return hash(self.label)

    def __eq__(self, other):
        return (self.label == other) or (self.label == other.label)

    def __iter__(self):
        self.i = 0
        self.n = len(self.file)
        return self

    def __next__(self):
        if self.i < self.n:
            outp = self.file[self.i]
            self.label = self.labels[self.i]
            self.i += 1
            return outp
        else:
            raise StopIteration


class DataHandler(DataProperties):
    def __init__(self, sim_name):
        super().__init__(sim_name)
        self.raw_data = {}
        self.file_paths = {}
        self._file_paths_built = False
        self._zeros_appended = False
        self._agg = False
        self.loaded_indices = []
        self.loaded_runs = {}
        self.labels = {}

        self.M = defaultdict(list)
        self.P = defaultdict(list)
        self.L = defaultdict(list)
        self.monomers = defaultdict(list)
        self.logt = defaultdict(list)
        self.logM = defaultdict(list)
        self.logP = defaultdict(list)
        self.logL = defaultdict(list)

    def load_data_from_index(self, index, label=None, run_number=0):
        if label is None:
            label = str(index)+'_'+str(run_number)
        if index in self.loaded_indices and run_number in self.loaded_runs:
            print('Already loaded. Label is "{}"'.format(self.labels[index][run_number]))
        else:
            if index not in self.loaded_indices:
                self.loaded_indices.append(index)
                self.loaded_runs[index] = [run_number]
                self.labels[index] = {run_number: label}
            elif run_number not in self.loaded_runs:
                self.loaded_runs[index][run_number] = label

            self.raw_data[label] = pd.read_json(self.run_dict[self.set_strings[index]][run_number]).sort_index()
            self._append_zeros(label)
            self.raw_data[label+"_props_index"] = self.protein_properties[index]

    def _append_zeros(self, label, max_length=1):
        if not self._zeros_appended:
            for ind, i in enumerate(self.raw_data[label].values):
                if len(i[0]) > max_length:
                    return self._append_zeros(label, len(i[0]))
                elif len(i[0]) < max_length:
                    for ii in range(max_length - len(i[0])):
                        self.raw_data[label].values[ind][0].append(0)

    def _load_agg(self, agg = 'Agg'):
        if not self._agg:
            matplotlib.use(agg)

    def _t(self, label):
        return self.raw_data[label]['t']

    def _logt(self, label):
        if self.logt[label] == []:
            self.logt[label] = [log10(i) for i in self.raw_data[label]['t']]
        return self.logt[label]

    def _mass(self, label):
        if self.M[label] == []:
            self.monomers[label] = self.raw_data[label+"_props_index"]['monomers']
            self.M[label] = [self.monomers[label]-i[0] for i in list(self.raw_data[label]['polymers'])]
        return self.M[label]

    def _logmass(self, label):
        return [log10(i) for i in self._mass(label)]

    def plot_mass(self, label, logx=False, logy=False, xrange=None, **kwargs):
        if logx:
            t = self._logt(label)
            if xrange is not None:
                xrange = [log10(i) for i in xrange]
        else:
            t = self._t(label)

        if logy:
            y = self._logmass(label)
        else:
            y = self._mass(label)
        self.plot_series(label, t, y, x_range=xrange, **kwargs)

    def plot_series(self, label, t, y, show=False, save=True, agg='Agg', name='default', x_range=None, x_index_range=None):
        self._load_agg(agg)
        import matplotlib.pyplot as plt

        if name == 'default':
            name = '{}_mass.png'.format(label)
        else:
            name += '.png'

        plt.figure()
        if x_range is not None:
            plt.xlim(x_range[0], x_range[1])
        elif x_index_range is not None:
            plt.xlim(t[x_index_range[0]], t[x_index_range[1]])

        plt.plot(t, y)
        if save:
            plt.savefig(name)
        if show:
            plt.show()
        plt.close()

    def show_parameter_sets(self):
        print("PARAMETER HIERARCHY\n")
        print("||=======PARAMETER==HIERARCHY=======||")
        if not self._file_paths_built:
            print("||{:<7}||{:<7}||{:<7}||{:>7}||".format("== a ==", "== b ==", "== k ==", "=index="))
            for ind, i in enumerate(self.rate_sets.values()):
                if i.a not in self.file_paths:
                    self.file_paths[i.a] = {i.b: {i.k: ind}}
                    print("||{:<7}||{:<7}||{:<7}||{:>7}||".format(i.a, i.b, i.k, ind))
                else:
                    if i.b not in self.file_paths[i.a]:
                        self.file_paths[i.a][i.b] = {i.k: ind}
                        print("||{:<7}||{:<7}||{:<7}||{:>7}||".format(i.a, i.b, i.k, ind))
                    else:
                        if i.k not in self.file_paths[i.a][i.b]:
                            self.file_paths[i.a][i.b][i.k] = ind
                            print("||{:<7}||{:<7}||{:<7}||{:>7}||".format(i.a, i.b, i.k, ind))
                        else:
                            counter = 1
                            in_dict = True
                            while in_dict:
                                counter += 1
                                ky = str(i.k) + "_{}".format(counter)
                                if ky not in self.file_paths:
                                    self.file_paths[i.a][i.b][ky] = ind
                                    print("||{:<7}||{:<7}||{:<7}||{:>7}||".format(i.a, i.b, i.k, ind))
                                    in_dict = True
            self._file_paths_built = True
        else:
            print("||{:<7}||{:<7}||{:<7}||{:>7}||".format("== a ==", "== b ==", "== k ==", "=index="))
            for a_key, aa in self.file_paths.items():
                new_a = True
                for b_key, bb in aa.items():
                    new_b = True
                    for k_key, kk in bb.items():
                        if new_a:

                            print("||{:<7}||{:<7}||{:<7}||{:>7}||".format(a_key, b_key, k_key, self.file_paths[a_key][b_key][k_key]))
                            new_a = False
                            new_b = False
                        elif new_b:
                            print("||{:<7}||{:<7}||{:<7}||{:>7}||".format("", b_key, k_key, self.file_paths[a_key][b_key][k_key]))
                            new_b = False
                        else:
                            print("||{:<7}||{:<7}||{:<7}||{:>7}||".format("", "", k_key, self.file_paths[a_key][b_key][k_key]))
        print("||==================================||")

