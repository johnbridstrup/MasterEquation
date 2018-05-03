import numpy as np
import pandas as pd
import statevector as sv
import matplotlib.pyplot as plt

def getdata(fn):
    pass


if __name__ == '__main__':
    fn='run1000.json'
    data=pd.read_json(fn)
    hist = data['histogram']
    print(hist)
    print(len(hist))
    print(hist[0][0])
    plt.figure(1)
    plt.subplot(221)
    plt.plot(list(range(1,100)),hist[1][0])
    plt.subplot(222)
    plt.plot(list(range(1,100)),hist[10][0])
    plt.subplot(223)
    plt.plot(list(range(1,100)),hist[50][0])
    plt.subplot(224)
    plt.plot(list(range(1,100)),hist[88][0])
    plt.show()
    