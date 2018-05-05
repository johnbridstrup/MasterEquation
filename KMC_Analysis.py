import numpy as np
import pandas as pd
import pickle
import matplotlib
import matplotlib.pyplot as plt
import statevector as sv
import glob
import importlib as port
import sys
utils=port.import_module('utils.utes')


def none_max(inp,maxx):
    if maxx is not None:
        return max(inp,maxx)
    else:
        return inp
    
def array_ify(unordered_df, maxx=None):
    try:
        larr=[unordered_df[i] for i in range(len(unordered_df))]
        print(type(larr))
        L=max([len(i) for i in larr])
        print(L)
        [i.extend([0] * (none_max(L,maxx) - len(i))) for i in larr]
        return larr
    except:
        larr=[unordered_df[i] for i in range(len(unordered_df))]
        print(type(larr))
        return larr
def main(fn):
    pass
    # dat=pd.read_json(sys.argv[1])
    # jawn_indices={column:index for index,column in enumerate(list(dat.columns))}
    # jawn=array_ify(dat)
    # print(jawn)
    # print(jawn_indices)
    # print (sys.argv)
    # print(jawn[jawn_indices['polymers']])
    
if __name__=="__main__":
    main('results/run1.json')
    print('argsssssss',sys.argv)
    z=pd.read_json('results/terps/terps0.json')
    pdat=array_ify(z['polymers'])
    t=array_ify(z['t'])
    plt.figure()
    for j in range(1,5): 
        plt.plot(t,[i[j] for i in pdat])
    plt.show()

    # polydat=[data['polymers'][i] for i in range(len(data['polymers']))]
    # pdat=np.array(polydat)
    # massdat=[data['mass'][i] for i in range(len(data['polymers']))]
    # numdat=[data['number'][i] for i in range(len(data['polymers']))]
    # stepdat=[data['t_steps'][i] for i in range(len(data['polymers']))]
    # tdat=[data['t'][i] for i in range(len(data['polymers']))]
    # kurtdat=[data['kurtosis'][i] for i in range(len(data['polymers']))]
    # skewdat=[data['skew'][i] for i in range(len(data['polymers']))]
    # histdat=[data['histogram'][i] for i in range(len(data['polymers']))]
    # statedat=[data['state'][i] for i in range(len(data['polymers']))]

    # print(polydat)

    # print(polydat[:][0:8])
    # print(tdat)
    # plt.plot(tdat,polydat[:][1])
    