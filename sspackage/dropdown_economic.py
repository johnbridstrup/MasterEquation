
# coding: utf-8

# Most examples work across multiple plotting backends, this example is also available for:
# 
# * [Matplotlib - dropdown_economic](../matplotlib/dropdown_economic.ipynb)

# In[1]:

import numpy as np
import Data
import pandas as pd
import holoviews as hv
hv.extension('bokeh')


# ## Declaring data

# In[2]:
def none_max(inp,maxx=0):
    return max(inp,maxx)
def ext_with_zeroes(l,mx):
    L=max([len(i) for i in l])
    return [i.extend([0] * (none_max(L,mx) - len(i))) for i in l]
def df_col2array(df,col,index_sort=True, max_list=10):
    tmp=[df[col][index] for index in range(len(df[col]))]
    ext_with_zeroes(tmp,max_list)
    return np.array(tmp)

macro_df = pd.read_csv('http://assets.holoviews.org/macro.csv', '\t')
my_df = pd.read_json('results/yuantest2.brid/yuantest2_aa_0.001_bb_0.001_kk_0.0001.strup/yuantest2_aa_0.001_bb_0.001_kk_0.0001.strup_1.json')
My_df=df_col2array(my_df,'polymers',max_list=8)
mm=np.array([sum(i[1:]) for i in My_df])
logt=np.array([np.log10(i) for i in my_df['t']])
my_df['t']=logt
my_df['polymers']=mm
key_dimensions   = [('t', 'time'), ('t', 'time')]
value_dimensions = [('polymers', 'Polymers'), ('t_steps', 'Time Steps')]
#macro = hv.Table(macro_df, key_dimensions, value_dimensions)
macro = hv.Table(my_df, key_dimensions, value_dimensions)


# ## Plot

# In[3]:


#get_ipython().run_cell_magic('opts', 'Overlay [width=700 height=400 show_frame=False]', "%%opts Curve (color='k') Scatter [color_index=2 size_index=2 scaling_factor=1.4] (cmap='Blues' line_color='k')\n%%opts VLine (color='k' line_width=1)\n%%opts Text (text_font_size='13px')\
#ngdp_curves = macro.to.curve('Year', 'GDP Growth')
ngdp_curves = macro.to.curve('t', 'polymers')

ngdp_unem_scatter = macro.to.scatter('t', ['polymers', 't_steps'])
#nannotations = hv.Arrow(1973, 8, 'Oil Crisis', 'v') * hv.Arrow(1975, 6, 'Stagflation', 'v') * hv.Arrow(1979, 8, 'Energy Crisis', 'v') * hv.Arrow(1981.9, 5, 'Early Eighties\\n Recession', 'v')

rend=hv.renderer('bokeh')
rend.save(ngdp_curves,'test2.html')