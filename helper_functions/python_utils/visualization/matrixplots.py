import matplotlib.pyplot as plt
import numpy as np 
from scipy.spatial.distance import pdist, squareform
from scipy.cluster.hierarchy import *
import matplotlib.pyplot as plt 
from mpl_toolkits.axes_grid1 import make_axes_locatable

import sys
sys.path.append('../statistics/')
from subsample import unbiasedDistanceMatrix 



def plotDistanceMatrix(data, trl_labels, grp_labels = None, sorted = False, sorting_order = None, 
                               metric = 'euclidean', clustering_method = 'single', biasCorrection = False, 
                               save_dendrogram = False, clip_negatives = True):
  '''
  Given input data X, plot a distance matrix of its classes. Supports multiple metrics
  and hierarchical clustering for distance matrix sorting. Inputs are:
  
  data ( array)      - features x reps matrix; distances are between rows
  trl_labels (list)  - reps iterable specifying trial condition 
  grp_labels         - labels for distance matrix; default = integer labeling
  sorted (Bool)      - if True, sort matrix by hierarchical clustering before plotting
  metric (str)       - metric to use for distance calculations: default = 'euclidean'
  clustering_method (str) - algorithm for hierarchical clustering; 'single' (nearest neighbor) or 'average' (UPGMA)
  save_dendrogram (Bool)  - if true, save a copy of the dendrogram as an .eps file
  
  Guy Wilson, NPTL, 2018.
  
  Returns:
  
  ax (matplotlib axes object) - axes for image 
  ordering (iterable)         - ordering applied to matrix rows and columns 
  
  TODO:
    - make more efficient for sorted condition (currently calculates dist mat 2x)
    - add multi-metric support for bias correction
  '''
  
  if grp_labels is None:
    grp_labels, inverse = np.unique(trl_labels, return_inverse = True)
  
  # calculate distances:
  if biasCorrection: 
    row_dists, _        = unbiasedDistanceMatrix(data, trl_labels, metric = metric, subtractMean = False)
      
    if clip_negatives:
      row_dists[row_dists < 0] = 0
    
  else: 
    data_in   = np.zeros((data.shape[0], len(grp_labels)))
    for i in range(len(grp_labels)):
      data_in[:, i] = np.mean(data[:, np.where(inverse == i)[0]], axis= 1)
    row_dists = pdist(data_in.transpose(), metric = metric)
    row_dists = squareform(row_dists)
  
  # optional hierarchical clustering for making matrix pretty:
  if sorted:
    if sorting_order is not None: 
      ordering       = sorting_order
    else:
      row_dists_dendrogram = row_dists.copy()
      row_dists_dendrogram[row_dists_dendrogram < 0] = 0   # small negative values indicate close to 0
    
      # linkage expects either compressed distance matrix or positions matrix 
      dists_compressed = np.ravel(row_dists_dendrogram[np.triu_indices(row_dists.shape[0], k = 1)])
      Z                = linkage(dists_compressed, optimal_ordering = True, method = clustering_method)
      tree             = dendrogram(Z, labels = grp_labels)
      ordering         = np.asarray(tree["leaves"]).astype(int)
      
  if save_dendrogram:
    plt.savefig('dendrogram.eps', format = 'eps')
    
  row_dists   = row_dists[ordering, :]
  row_dists   = row_dists[:, ordering]
  grp_labels  = grp_labels[ordering]

  # plot the results: 
  fig     = plt.figure()
  ax      = plt.gca()
  im      = ax.imshow(row_dists, vmin = np.min(row_dists[np.nonzero(row_dists)]))
  divider = make_axes_locatable(ax)
  cax     = divider.append_axes("right", size="5%", pad=0.2)
  ax.figure.colorbar(im, cax=cax)
  
  ax.set(xticks= np.arange(len(grp_labels)), yticks= np.arange(len(grp_labels)),
         xticklabels= grp_labels, yticklabels= grp_labels)
  

  return plt.gca(), row_dists, ordering



