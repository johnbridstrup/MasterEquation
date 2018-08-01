
# coding: utf-8

# <div class="contentcontainer med left" style="margin-left: -50px;">
# <dl class="dl-horizontal">
#   <dt>Title</dt> <dd> Histogram Element</dd>
#   <dt>Dependencies</dt> <dd>Matplotlib</dd>
#   <dt>Backends</dt> <dd><a href='./Histogram.ipynb'>Matplotlib</a></dd> <dd><a href='../bokeh/Histogram.ipynb'>Bokeh</a></dd>
# </dl>
# </div>

# In[1]:


import numpy as np
import holoviews as hv
import matplotlib as mpl
mpl.use('TkAgg')
hv.extension('matplotlib')


# ``Histogram``s partition the `x` axis into discrete (but not necessarily regular) bins, showing counts in each as a bar. A ``Histogram`` accepts the output of ``np.histogram`` as input, which consists of a tuple of the histogram values with a shape of ``N`` and bin edges with a shape of ``N+1``. As a simple example we will generate a histogram of a normal distribution with 20 bins.

# In[ ]:


np.random.seed(1)
data = np.random.randn(10000)
frequencies, edges = np.histogram(data, 20)
print('Values: %s, Edges: %s' % (frequencies.shape[0], edges.shape[0]))
hv.Histogram((edges, frequencies))


# The ``Histogram`` Element will also expand evenly sampled bin centers, therefore we can easily cast between a linearly sampled Curve or Scatter and a Histogram.

# In[ ]:


xs = np.linspace(0, np.pi*2)
ys = np.sin(xs)
curve = hv.Curve((xs, ys))
curve + hv.Histogram(curve)


# The ``.hist`` method is an easy way to compute a histogram from an existing Element:

# In[ ]:


points = hv.Points(np.random.randn(100,2))

points.hist(dimension=['x','y'])

renderer=hv.renderer('matplotlib')
plt=renderer.save(curve,'html_test.html')

# The ``.hist`` method is just a convenient wrapper around the ``histogram`` operation that computes a histogram from an Element, and then adjoins the resulting histogram to the main plot. You can also do this process manually; here we create an additional set of ``Points``, compute a ``Histogram`` for the 'x' and 'y' dimension on each, and then overlay them and adjoin to the plot.

# In[ ]:


#get_ipython().run_cell_magic('opts', 'Histogram (alpha=0.3)', "from holoviews.operation import histogram\npoints2 = hv.Points(np.random.randn(100,2)*2+1).redim.range(x=(-5, 5), y=(-5, 5))\n\nxhist, yhist = (histogram(points2, bin_range=(-5, 5), dimension=dim) *\n                histogram(points,  bin_range=(-5, 5), dimension=dim) \n                for dim in 'xy')\n(points2 * points) << yhist(plot=dict(width=125)) << xhist(plot=dict(height=125))")


# For full documentation and the available style and plot options, use ``hv.help(hv.Histogram).``
