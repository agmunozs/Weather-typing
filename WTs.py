# classifiability.py
'''This Module contains functions required to produce classifiability index applicable to K-means
   Clustering Process. This algorithm is a python translation of MATLAB functions written by N. Fauchereau, modified by 
   Á.G. Muñoz (AOS-Princeton,IRI-Columbia)
   This classifiability index is clearly set out in Michelangeli et al. 1995
   as defined by Cheng and Wallace 1993.
   Functions:
   index_classify
   ref_partition
   kmeanscluster
   corrmat'''
import sys
from numpy import *
import pyclimate.mvarstatools as pyclim
from scipy.cluster.vq import *

def corrmat(x,y):
  ''' cor = corrmat(x,y)
      
      Function called by index_classify to calculate correlations between two matrices
      
      Input: x, y observations need be first dimension and x.shape[0] == y.shape[0] must be true
      Returns: correlation matrix'''
  nlig, ncol = x.shape
  nlig1,ncol2= y.shape
  if nlig != nlig1:
    print 'x and y do NOT have same number of observations!'; sys.exit(1)
  x= pyclim.standardize(x)
  y= pyclim.standardize(y)
  y=mat(y);x=mat(x)
  cor=(y.transpose()*x)/(nlig-1)
  return asarray(cor)


def index_classify(centros,data):
  '''clindex, pcmax, dimax = index_classify(centros,data)
     
     Function to calculate Cheng and Wallace 1993 / Michelangeli et al. 1995
     classifiability index
     
     Input: centros (array) dimensions MUST be npts x nclasses x nexpe
            data (array) the obs x points data which was clustered
     Returns: clindex (tuple) contains classifiability index of given no. of clusters 
              pcmax
              dimax'''
  ndata, npts = data.shape
  average = mean(data,0)
  vartotale =sum(sum(data*data,0) - ndata*average**2)
              
  npoints, nclasses, nexpe =centros.shape
  fac= 1/float((nexpe*(nexpe-1)))
  centros = centros - tile(average[:,newaxis,newaxis],(1,nclasses,nexpe))
  pcmax = zeros((nexpe,nexpe)); dimax=zeros((nexpe,nexpe))
  for ne1 in xrange(nexpe):
    for ne2 in xrange(nexpe):
         x=centros[:,:,ne1]; y=centros[:,:,ne2]
         #calculation of correlations and selection of the best correlation 
         r = corrmat(x,y)
         pc = 1 - r
         pc = pc.min(0).max();
         pcmax[ne1,ne2] = pc
         # distance calculation and selection of the optimal distance
         x = tile(x,(nclasses,1)); x = x.reshape(npoints,nclasses*nclasses,order='F')
         y = tile(y,(1,nclasses))
         di2 = sum( (x-y) * (x-y) ,0)
         di2 = di2.reshape(nclasses,nclasses,order='F')
         dimax[ne1,ne2] = di2.min(0).max()
  
  # Index of Classifiability
  clindex1 = 1 - sum(pcmax*fac)              # pcmax=clindex1;
  clindex2 = sum(dimax*fac)/vartotale;       # dimax=clindex2;
  clindex=[clindex1,clindex2];

  return clindex, pcmax, dimax
  
  
def ref_partition(pcmax,nexpe):
  ''' iclref, parsim = partref(pcmax,nexpe)
  
      Function that returns best partition of data from all nexpe experiments
      
      Input: pcmax (array) output of index_classify
             nexpe (int) number of experiments
     Returns: iclref (1-dim array) choice of best partition
               parsim (int) the experiment that returned this best assignment'''
  pcmax2 = pcmax / (nexpe-1);
  parsim = sum(pcmax2.transpose());
  [parsim,iclref]=min(parsim);
  parsim=1-(parsim);
  return iclref, parsim


def kmeanscluster(data,nexp=50,nk=[2,3,4,5,6,7,8,9,10,11,12,13]):
  '''kclasses, clindex = kmeanscluster(data,nexp=50,nk=[2,3,4,5,6,7,8,9,10,11,12,13,14,15])
  
     Function (to be implemented as a class at later stage) that calculates k-means clustering
     Requires modules: pyclimate, scipy.cluster.vq
     
     Input:   data (array) obs/time x points
     Returns: clindex (array) which gives indication of how many clusters to use
              more variables are to be returned in the future'''
  obs = whiten(data)
  clindex = [];
  for j in xrange(len(nk)):
    nclasses = nk[j]
    centros = zeros((data.shape[1],nclasses,nexp))
    for i in xrange(50):
      ne=i
      cbook, totaldist = kmeans(obs,nclasses,iter=50)
      try:
        centros[:,:,ne] = cbook.transpose()
      except:
	print "Experiment no. ",ne+1," Failed: nclasses=",nclasses," & cbook shape is ",cbook.shape
    clind, pcmax, dimax = index_classify(centros,data)
    clindex.append(clind)
  
  #kclass, obs_dists = vq(obs,cbook)
  
  return array(clindex)