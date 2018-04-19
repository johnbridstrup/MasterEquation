import os
import csv
import matplotlib
import matplotlib.animation as manimation
import matplotlib.pyplot as plt
import numpy as np
import configparser as cp

# ==============================================================================
# Filepath and directory utilities
# ==============================================================================
def formatFilename(name,folder,var_list,val_list,*args):
    for i,j in zip(var_list,val_list):
        name=name+"_"+i+"_"+str(j) 
    dirpath="/"
    for i in args:
        if type(i) is str:
            dirpath+=i+"/" 
    cur_dir=os.getcwd()+dirpath
    return cur_dir+folder+"/"+name
def formatDirname(name,folder,var_list,val_list,*args):
    for i,j in zip(var_list,val_list):
        name=name+"_"+i+"_"+str(j) 
    dirpath="/"
    for i in args:
        if type(i) is str:
            dirpath+=i+"/"
    return dirpath+folder+"/"+name
def formatDataFilename(name,var_list,val_list,*args):
    for i,j in zip(var_list,val_list):
        name=name+"_"+i+"_"+str(j) 
    dirpath="/"
    for i in args:
        if type(i) is str:
            dirpath+=i+"/" 
    cur_dir=os.getcwd()+dirpath
    return cur_dir+"Data/"+name
def formatMovieFilename(name,var_list,val_list,*args):
    for i,j in zip(var_list,val_list):
        name=name+"_"+i+"_"+str(j) 
    dirpath="/"
    for i in args:
        if type(i) is str:
            dirpath+=i+"/" 
    cur_dir=os.getcwd()+dirpath
    return cur_dir+"Movies/"+name
def formatPlotFilename(name,var_list,val_list,*args):
    for i,j in zip(var_list,val_list):
        name=name+"_"+i+"_"+str(j) 
    dirpath="/"
    for i in args:
        if type(i) is str:
            dirpath+=i+"/" 
    cur_dir=os.getcwd()+dirpath
    return cur_dir+"Plots/"+name

def makeDir(dirname,ext=None):
    if not os.path.exists("./"+dirname+"/"):
        os.makedirs("./"+dirname+"/")
def makeAbsDir(dirname,ext=None):
    if not ext:
        if not os.path.exists(dirname):
            os.makedirs(dirname)
    elif type(ext)==str:
        if not os.path.exists(dirname+"."+ext):
            os.makedirs(dirname+"."+ext)
    else:
        print("WHACK DOG EXT SHOULD BE A STRING")
def makeDirp(dirname,path):
    if not os.path.exists("./"+path+dirname):
        os.makedirs("./"+path+dirname)

# ==============================================================================
# Program specific data utilities
# ==============================================================================
def formatPolymerData(d,t, ddo=False):
    if ddo==True:
        x=[]
        tt=[]
        for i,v in enumerate(d):
            if i<len(d)-1:
                x.append(v)
            x.append(v)
        for i,v in enumerate(t):
            if i>0:
                tt.append(v)
            tt.append(v)
        return (x,tt)
    else:
        return(d,t)
def input_loop(callbacks,*arguments):
    looping=True
    while(looping):
        cmd_input=input("0 to end loop")
        if cmd_input=="0":
            looping=False
        for index,callback in enumerate(callbacks):
            callback(arguments[index])  
# ==============================================================================
# Movie utilities
# ==============================================================================
def movieWriter(xdata,ydata,filepath, frames): 
    FFMpegWriter = manimation.writers['ffmpeg']
    moviewriter = FFMpegWriter(fps=30)
    fig = plt.figure()
    l, = plt.plot([], [])
    plt.xlim(min(xdata),max(xdata))
    yy=0
    for i in ydata:
        for jj in i:
            if yy<jj:
                yy=jj
    plt.ylim(0,yy)
    # print(filepath)
    with moviewriter.saving(fig,filepath+'.mp4',dpi=100):
        for j in range(frames):
            l.set_data(xdata,ydata[j])
            moviewriter.grab_frame()
    plt.close(fig)
    plt.gca()
    plt.gcf()
    plt.cla()
# ==============================================================================
# CSV utilities
# ==============================================================================
def csvWriter(data,filepath, tag = None):
    if type(tag)==str:
        filepath+="_"
        filepath+=tag
    with open(filepath+'.csv','w',newline='') as f:
        
        writer = csv.writer(f,delimiter=',')
        # print(data,"data")
        data1=[[i] for i in data]
        writer.writerows(data1)
def csvLoad(fn):
    rdr=open(fn,'r')
    csvrdr=csv.reader(rdr,quoting=csv.QUOTE_MINIMAL)
    return csvrdr
def getColumn(filename, column):
    results = csv.reader(open(filename), delimiter=",",quoting=csv.QUOTE_NONNUMERIC)
    return [result[column] for result in results]
# ==============================================================================
# configparser utilities
# ==============================================================================
def inputReader(filename):
    config = cp.ConfigParser()
    config.read(filename)
    return config
def configWriter(filename, sections, op_val_pair_dict):
    parser=cp.ConfigParser()
    for i in sections:
        print(i)
        parser.add_section(i)
    for key,val in op_val_pair_dict.items():
        print(key,val[0],val[1])
        parser.set(key,val[0],val[1])
def model_config_writer(filename,strings,reactions):
    sections=["model"]
    op_val={}
    op_val["model"]=[("name",strings[0]),("type",strings[1]),("algorithm",strings[2]),("reactions",reactions)]
    print(op_val)
    configWriter(filename,sections,op_val)

    
# ==============================================================================
# List Utilities
# ==============================================================================
def ElementDivision(list1, list2,):
    return[i/j for i,j in zip(list1,list2)]
def list_over_scalar(list, scalar):
    return [i/scalar for i in list]
def normalize(list):
    magnitude=sum(list)
    return [j/magnitude for j in list]
class store_list:
    def __init__(self, data,filename):
         self.data=data
         self.storage=[]
         self.filename=filename
    def append(self):
        self.storage.append(self.data[:]) 
    def write(self):
        csvWriter(self.storage,self.filename)

