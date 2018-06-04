import sys
import os
import glob
from pathlib import Path
import importlib as port 
from functools import wraps
from numpy import ndarray
from datetime import datetime
import plistlib as pll

# class includes:
#     def __init__(self,path,includes=None):
#         self.path=path
#         self.incl={}
#         if includes is not None:
#             try:
#                 self.incl={key:item for key,item in includes.items()}
#             except:
#                 try:
#                     self.incl={str(index):item for index,item in enumerate(includes)}
#                 except:
#                     try:
#                         if includes==str
#                             self.incl=includes
#                     except:
#                         print("stop breaking things")
#                         raise
def path_leaf(path):
    head, tail = os.path.split(path)
    return tail or os.path.basename(head)  
"""
class includer:
    """"""include utility
    
    has a dictionary of paths to include files and semantics for importing
    a particular *.py file or function from within it.
""""""
    #& One dictionary for paths: {labels[i]:path.to.files}
    #& one dictionary for librarys {label from labels[i]:file or object}
    
    """
"""
    def __init__(self, libs=None,paths=None, labels=None):
        self._paths=paths
        self._libs=libs
        try:
            if labels is None:
                labels = [path_leaf(path) for path in paths]
                self.paths = {lab: path for lab,
                    path in zip(labels, paths)}
        except:
            pass
        self._labels=labels
        if paths is not None:
            if isinstance(paths,dict): #{libname : path to lib + libname}
                self.paths=paths
            elif isinstance(paths,list) and labels is not None:
                try:
                    self.paths={lab:path for lab,path in zip(labels,paths)}
                    try: 
                        self.libs={lab:lib for lab,lib in zip(labels,libs)}
                    except:
                        pass
                except:
                    try:
                        self.paths={}
                        if all(isinstance(i,str) for i in paths):
                            self.paths[labels]=paths
                            try:
                                self.libs={}
                                self.libs[labels]=libs
                            except:
                                pass
                        else:
                            raise TypeError('they arent all strings')
                    except:
                        try:
                            if isinstance(labels,str) and isinstance (paths,str):
                                self.paths={labels:paths}
                                try:
                                    self.libs={labels:paths}
                                except:
                                    pass
                            elif isinstance(labels,str) and isinstance (paths,list):
                                if all(isinstance(path,str) for path in paths):
                                    self.paths={labels:paths}
                                    try:
                                        if all(isinstance(lib,str) for lib in libs):
                                            self.libs={labels:libs}
                                        else:
                                            raise TypeError('not all strings')
                                    except:
                                        pass
                        except:
                            pass 
        else:
            
    def addpath(self,path={"label":"path/"}):
        self._labels.append(list(path.keys())[0])
    def addlib(self,lib={}):
        self.libs[(list(lib.keys())[0])]
    def __call__(self,libs):
        try:
            for key,path in self.paths.items():
                return {self._labels[key]:[path+'/'+lib for lib in libs[key]]}
        except:
            print("broked")
                               
    """
def wd():
    return os.getcwd()
def mdirs(fn,dirs="/results/"):
    os.makedirs(wd()+dirs+fn,exist_ok=True)
def file_list(dr=None,ext=None):
    if dr is not None:
        try:
            if ext is not None:
                try:
                    return glob.glob(dr+'*.'+ext)
                except:
                    try:
                        return [glob.glob(dr+'*.'+i) for i in ext]
                    except:
                        print("{} is not a valid path probably or maybe {} is some weird object".format(dr,ext))
                        raise
            else:
                try:
                    return glob.glob(dr+'*')
                except:
                    print("{} is not a valid path probably or maybe {} is some weird object".format(dr,ext))
                    raise
        except:
            try:
                return [file_list(i,ext) for i in dr]
            except:
                print("{} (or \"the directory\" if i'm to believe your nonsense) smells funky".format(dr))
                raise
def newdirs(dr, path=None, rel_path=None, Dirs=None):
    if path is None:
        path=wd()
    if rel_path is not None:
        path = path + '/' + rel_path
    print(path)
    if Dirs is None:
        Dirs=os.listdir(path)
    if not os.path.isdir(path):
        try:
            os.mkdir(path+'/'+dr)
        except:
            try:
                for i in dr:
                    newdirs(i,path,Dirs)
            except:
                print('{} is already a directory'.format(i))


class config_context:
    mode='wb'
    @classmethod
    def set_mode(cls,mode='wb'):
        cls.mode=mode
    @classmethod
    def __enter__(cls):
        cls.open_file = open(wd()+'/inputs/', cls.mode)
        return cls.open_file
    @classmethod
    def __exit__(cls, *args):
        cls.open_file.close()
class conf_context():
    def __init__(self, filename, mode='wb',write=True):
        self.filename = filename
        self.mode = mode
        if write:
            mkdirs(filename)

    def __enter__(self):
        self.open_file = open(self.filename, self.mode)
        return self.open_file

    def __exit__(self, *args):
        self.open_file.close()
def mkdirs(fn):
    os.makedirs(fn,exist_ok=True)
def load_config(fn):
    config_context.set_mode(mode='r')
    oup=pll.readPlist(fn)
    return oup
if __name__ == '__main__':
    newdirs("test")
    print(file_list(wd()))            
      
