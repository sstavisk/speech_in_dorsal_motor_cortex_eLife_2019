import numpy as np
from numpy.random import choice 

# module for sumbsampling-based estimations 


def subsample(data1, data2, numsamples = None):
  '''
  Given two arrays of data, subsample one randomly to have
  same length as other. Optional argument for subsampling both 
  with a specific sample #. Inputs are: 
  
  data1, data2 (1D iterable) - data to be subsampled
  numsamples (int)           - if not None, select this many samples randomly 
                               from both 
                               
  TODO:
    - can probz make things v easy by "randomly subsampling" an entire matrix for the larger one; 
      should help collapse cases
  '''
  
  if numsamples is None:
    numsamples = min(len(data1), len(data2))

  if len(data1) > len(data2):
    data1_subsample = choice(data1, numsamples)
    data2_subsample = data2
    
  elif len(data1) < len(data2):
    data1_subsample = data1
    data2_subsample = choice(data2, numsamples)
    
  else:  # data1 and data2 are same length
    data1_subsample = choice(data1, numsamples)
    data2_subsample = choice(data2, numsamples)
    
    
  return data1_subsample, data2_subsample
    
def unbiasedDistance(data1, data2, subtractMean = False, metric = 'euclidean'):
  '''
  Use a median-unbiased estimate of distances between different groups. 
  Inputs are: 
  
    data1 (features x reps array) - data matrix holding values of group 1
    data2 (features x reps array) - data matrix holding values of group 2
    subtractMean (Bool)           - if True, center each vector before computing difference size
    metric (str)                  - metric to use; can be 'euclidean', 'correlation', 'mahalanobis'
  
  Guy's version of Sergey's version of Frank's code.
  Guy Wilson, NPTL, 2018.
  
    % If class 1 and class 2 are of differnet numbers of trials, it loops through all
    % combinations of trial_i (from class 1) and all trial
    '''
  
  assert data1.shape[0] == data2.shape[0], "Groups have different feature numbers. Cannot calculate distances."

  data1_reps = data1.shape[1]
  data2_reps = data2.shape[1]

  if data1_reps == data2_reps:
    squaredDistEstimates = np.zeros((data1_reps,))

    for trl in range(data1_reps):
      bigSetIdx         = np.setdiff1d(np.arange(data1_reps), [trl])
      leaveOutIdx       = trl 

      meanDiff_bigSet   = np.mean(data1[:, bigSetIdx] - data2[:, bigSetIdx], axis = 1)
      meanDiff_smallSet = data1[:, leaveOutIdx] - data2[:, leaveOutIdx] 

      if subtractMean:  
        squaredDistEstimates[trl] = (meanDiff_bigSet - np.mean(meanDiff_bigSet)).dot(meanDiff_smallSet - np.mean(meanDiff_smallSet))
      else:
        squaredDistEstimates[trl] = meanDiff_bigSet.dot(meanDiff_smallSet)

  # unequal trial counts: 
  else:
    squaredDistEstimates = np.zeros((data1_reps * data2_reps ,))
    iteration            = 0  # will increment through squaredDistEstimates

    for data1_trl in range(data1_reps):
      for data2_trl in range(data2_reps):
        bigSetIdx_c1   = np.setdiff1d(np.arange(data1_reps), [data1_trl]) # class 1
        bigSetIdx_c2   = np.setdiff1d(np.arange(data2_reps), [data2_trl]) # class 2

        leaveOutIdx_c1 = data1_trl
        leaveOutIdx_c2 = data2_trl

        meanDiff_bigSet   = np.mean(data1[:, bigSetIdx_c1], axis = 1) - np.mean(data2[:, bigSetIdx_c2], axis = 1)
        meanDiff_smallSet = data1[:, leaveOutIdx_c1] - data2[:, leaveOutIdx_c2] 

        if subtractMean:  
          print('Using subtract mean')
          squaredDistEstimates[iteration] = (meanDiff_bigSet - np.mean(meanDiff_bigSet)).dot(meanDiff_smallSet - np.mean(meanDiff_smallSet))
        else:
          squaredDistEstimates[iteration] = meanDiff_bigSet.dot(meanDiff_smallSet)

        iteration += 1

  meanOfSquares      = np.mean(squaredDistEstimates)
  lessBiasedEstimate = np.sign(meanOfSquares) * np.sqrt(np.abs(meanOfSquares))

  return lessBiasedEstimate, meanOfSquares, squaredDistEstimates 
    
  
def unbiasedDistanceMatrix(data, labels, subtractMean = False, metric = 'euclidean'):
  '''
  Use a median-unbiased estimate of distances between different groups. 
  Inputs are: 
  
    data (features x reps array) - data matrix holding positions of classes/groupings
    labels (reps iterable)       - holds class labels for reps in data
    subtractMean (Bool)          - if True, center each vector before computing difference size
    metric (str)                 - metric to use; can be 'euclidean', 'correlation', 'mahalanobis'
    
    
    TODO:
      - implement correlation, mahalanobis distances
      - optimize speed??
  
  Guy's version of Sergey's version of Frank's code. 
  Guy Wilson, NPTL, 2018.
  '''

  unique, inverse = np.unique(labels, return_inverse = True)
  num_classes     = len(unique)
  dists           = np.zeros((num_classes, num_classes))
  
  for i in range(num_classes):
    for j in range(num_classes): 
      data1 = data[:, np.where(inverse == i)[0]]
      data2 = data[:, np.where(inverse == j)[0]]
      if i <= j:
        dists[i, j], _, _ = unbiasedDistance(data1, data2, subtractMean = subtractMean, metric = metric)
      else:
        dists[i, j] = dists[j, i]

  return dists, unique


  
  
  
  
  
  
  
  
  
  