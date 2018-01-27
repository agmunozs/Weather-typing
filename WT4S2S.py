
# coding: utf-8

# In[20]:


import xarray as xr
from netCDF4 import Dataset
import pydap as dap
import numpy as np
import sklearn.cluster as cl
from eofs.xarray import Eof
import pandas as pd
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature
from PlotContourQuiver import *


# In[ ]:


get_ipython().system(u'rm -Rf hgt_ecmwf.nc')
get_ipython().system(u"curl -k -b '__dlauth_id=017a28e8531cac13efd89be8a7612c4c0754a83606f8f90270d14d84f62c28d7ff7fe8fbfb04c0495ddf938392d0bf3d9617e8b7' 'https://iridl.ldeo.columbia.edu/SOURCES/.ECMWF/.S2S/.ECMF/.reforecast/.control/.pressure_level_gh/.gh/P/500/VALUE/P/removeGRID/S/(last)/VALUE/S/removeGRID/hdate/(1997)/(2016)/RANGE/X/(-150)/(-50)/RANGE/Y/(20)/(80)/RANGE/data.nc' >hgt_ecmwf.nc")


# In[4]:


hgt = xr.open_dataset('./hgt_ecmwf.nc',decode_times=False,engine='scipy')['gh']
hgt=hgt.stack(time=('hdate', 'L')).T


# In[5]:


class KMean_XR_CI_One:
    '''A single simulation'''
    def __init__(self, pc, n_cluster,  n_init=1):
        self.centroid, self.cluster, _ = cl.k_means(pc, n_clusters=n_cluster, n_init=n_init)
        self.ReSort() # re-sort cluster labels
    
    def ReSort(self):
        # a helper function to re-sort the cluster labels
        x = np.int_(self.cluster)
        x_from = np.unique(x).argsort()
        counts = np.array([np.sum(x == xi) for xi in x_from])
        orders = (-counts).argsort().argsort()
        self.cluster = orders[x_from[x]]


# In[6]:


class KMean_XR_CI:
    
    def __init__(self, ds, n_cluster, prop_variance, var_name,  n_init=50):
        '''Initialize the object and perform calculations'''
        
        # input data
        assert len(ds.shape) is 3, ""
        self.rawdata = ds
        assert type(n_cluster) is int, "n_cluster must be int"
        self.n_cluster = n_cluster
        assert 0 < prop_variance < 1, "prop_variance must be between 0 and 1"
        self.prop_variance = prop_variance
        assert type(var_name) is str, "var_name must be a string"
        self.var_name = var_name
        assert type(n_init) is int, "n_init must be an integer"
        self.n_init = n_init
        
        # perform calculations
        self.Decompose() # create PC time series -- EOF space
        self.Cluster() # run clustering algorithm
        self.CalcCI()
    
    def Decompose(self):
        '''perform EOF decomposition'''
        solver = Eof(self.rawdata)
        var_frac = solver.varianceFraction()
        cumvar = np.cumsum(var_frac.values)
        self.npcs = np.where(cumvar >= self.prop_variance)[0].min()
        self.pc = solver.pcs(npcs = self.npcs) # time series of PC
        
    def Cluster(self):
        '''Run the k-means algorithm to get centroids & assignments'''
        k_means = [KMean_XR_CI_One(pc=self.pc.values, n_cluster=self.n_cluster, n_init=1) for i in range(self.n_init)]
        self.centroids = np.array([k.centroid for k in k_means])
        self.cluster = np.array([k.cluster for k in k_means])

    def CalcCI(self):
        def CalcACC(P,Q):
            '''Anomaly Correlation Coefficient'''
            acc = np.dot(P,Q) / np.sqrt(np.dot(P,P) * np.dot(Q,Q))
            return(acc)
        
        # initialize
        Aij = np.ones([self.n_init, self.n_cluster, self.n_cluster])
        Aprime = np.ones([self.n_init, self.n_cluster])
        Abest = np.ones(self.n_init)
        # compute classifiability index
        for n in range(0, self.n_init):
            for i in range(0, self.n_cluster):
                for j in range(0, self.n_cluster):
                    if i == j:
                        Aij[n, i, j] = np.nan
                    else:
                        P = self.centroids[n, i, :]
                        Q = self.centroids[n, j, :]
                        Aij[n, i, j] = CalcACC(P,Q)
                Aprime[n, :] = np.nanmax(Aij[n,:,:], axis=0)
            Abest = Aprime.min(axis=1)
        # extract useful
        self.classifiability = Abest.mean()
        self.best_sim = np.where(Abest == Abest.max())[0][0]
    
    def GetBestSim(self):
        return self.centroids[self.best_sim], self.cluster[self.best_sim], self.classifiability


# In[7]:


sol=KMean_XR_CI(ds = hgt, n_cluster=5, prop_variance=0.9, var_name='hgt', n_init=100)


# In[8]:


centroids, wt, classifiability = sol.GetBestSim()


# In[9]:


print('Classifiability Index:', classifiability)


# In[10]:


print(wt+1)


# In[11]:


wt_unique = np.sort(pd.unique(wt))


# In[12]:


clusters = np.sort(np.unique(wt))
prop_days = [np.mean(wt==i) for i in clusters]
my_proj = ccrs.Mercator()


# In[13]:


fig, axes = SetupAxes(ncol = 2, nax = len(wt_unique), proj = my_proj, figsize = [15, 9])
for i in range(len(wt_unique)):
    ax = axes[GetRowCol(i, axes)]
    AxisContourQuiver(
        ds = SubSetMean(ds=reanalysis, i=i, cat_ts=wt), 
        contour_name = 'hgt', quiver_name = None, maxval = 25,
        ax = ax, n_level = 31)
    ax.set_title('Weather Type {}: {:.1%} of days'.format(i+1, prop_days[i]))
FormatAxes(axes, coast=True, grid=True, border=True, river = False, extent=None)
plt.show()


# In[17]:


get_ipython().system(u'pip3.6 install cartopy')


# In[18]:


wt_unique

