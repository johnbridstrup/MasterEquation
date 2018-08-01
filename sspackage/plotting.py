import importlib as port
dat=port.import_module('Data')
import matplotlib
import matplotlib.pyplot as plt


name=dat.get_input("What is the folder name? ")
bins=100
t_range=(0,2000)

datafiles,data = dat.prep(name)

## Histograms
plot_histograms=True
paramsets=[(i,the_dict['rates']) for i,the_dict in list(enumerate(datafiles.master[name]['parameters']))]
if plot_histograms:
    for i,j in paramsets:
        print("parameter set {} has rate constants:  \n".format(i))
        [print("\t{} : {}".format(ii,jj)) for ii,jj in j.items()]
    which_sets=dat.get_input("Which histogram slices (int or 'all')? ")
    histgens=[]
    if which_sets=='all':
        (hist_gens,polymer_names) = ([data.histogram_factory(name,k) for k in range(len(paramsets))],{str(i):[j for j in data.noNan_master_data[name][0][i].keys()] for i in range(len(paramsets))})

    else:
        hist_gens = [data.histogram_factory(name,k) for k in range(len(paramsets))]
    poly_names=[data.noNan_master_data[name][0][i[0]].keys()[j] for i in paramsets for j in range(datafiles.master[name]['parameters'][0]['simulation']['runs']) ]
hists=get_hists(which_sets,data)
        





if plot_histograms:
    for i in which_slices:
        
        pass
