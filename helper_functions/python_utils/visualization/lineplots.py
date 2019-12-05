import matplotlib.pyplot as plt
import numpy as np
import sys
sys.path.append('../firingrate/')
from firingrate import getTrialActivity 


def plotsd(data, color, time_bins = None, toggleSE = False, alpha = 0.2):
  '''
  Plotting with standard error shading. Inputs are:
  
  data (2D array)      - reps x time; SDs are taken along 1st axis
  color                - color value 
  time_bins (1D array) - timestamp of bins from data 
  toggleSE (Boolean)   - if true, use SE shading 
  '''
  
  numreps     = data.shape[0]
  if time_bins is None:
    time_bins = np.arange(0, data.shape[1])
  
  mean_signal = np.mean(data, axis = 0)
  sd_signal   = np.std(data, axis= 0) 
  
  if toggleSE:
    sd_signal /= np.sqrt(numreps)
  
  plt.plot(time_bins, mean_signal, color= color)
  plt.fill_between(time_bins, mean_signal - sd_signal, mean_signal + sd_signal, color=color, alpha=alpha)
  
  
  
  
def grandAveragePlot(data_struct, neurons = None, alignment = None, numbins = None):
  '''
  Plot FR average across all trials (aligned to trial start). Inputs are:
  
  data_struct (DataStruct) - uses the .FR attribute for plotting 
  neuron (array)           - neurons to plot (default: all)
  alignment (array)        - FR timestamps to align data to (default: 0)
  numbins (int)            - if not None, plot alignment:(alignment + numbins) time segment            
  '''
  
  sample_period = data_struct.FRbinsize / 1000
  stacked_trls  = getTrialActivity(data_struct, neurons, alignment, numbins)
  
  figs          = list()
  
  for i in range(len(neurons)):
      #fig = plt.figure()
      plt.plot(np.mean(stacked_trls[i], axis = 0) / sample_period)
      plt.ylabel('FR (Hz)')
      plt.xlabel('Time (msec)')
      #plt.show()
      
     # figs.append(fig)
      
 # return figs


def grandAverageSummaryPlot(data_struct, neuron, d= 0.15, cue_shift = 300, num_bins = 500, toggleSave = False, path = None):
  '''
  Summary plot for a given neuron. Top row is GA total time course. Bottom three plots are (left to right)
  aligned to start, goCue, and timeSpeech. Inputs are:
  
  data_struct (DataStruct) - dataset source
  d (float)                - size of broken axis markers
  neuron (int)             - neuron to plot
  cue_shift (int)          - how much backward in time to plot from goCue, timeSpeech
  num_bins (int)           - total timelength (in bins) of bottom plots
  toggleSave (Bool)        - toggle figure download; if False, generates image instead 
  path (Str)               - path for image outputs 
  '''
  
  
  #plt_neuron = 7
  d          = .015 # how big to make the diagonal lines in axes coordinates
  
  plt.subplot(211)
  grandAveragePlot(data_struct, neurons = [neuron])
  plt.title('Grand-Averaged Timecourse')
  ymin, ymax = plt.ylim()
  
  mean = np.mean(data_struct.FRgoCue)
  std  = np.std(data_struct.FRgoCue)
  plt.fill_between(np.arange(mean - std, mean + std), ymin, ymax, color= 'k', alpha=0.2)

  ymax += 10

  plt.subplot(234)
  grandAveragePlot(data_struct, neurons = [neuron], numbins = 500)
  plt.axvline(0, color = 'k', linestyle = '--')
  plt.ylim(bottom = ymin, top = ymax)
  ax1 = plt.gca()
  ax1.spines['top'].set_visible(False)
  ax1.spines['right'].set_visible(False)
  plt.title('Speech onset alignment')
  plt.title('Trial start alignment')

  plt.subplot(235)
  grandAveragePlot(data_struct, neurons = [neuron], alignment = data_struct.goCue - cue_shift, numbins = num_bins)
  plt.axvline(cue_shift, color = 'k', linestyle = '--')
  plt.ylim(bottom = ymin, top = ymax)
  ax2       = plt.gca()
  
  ax2.get_yaxis().set_visible(False)
  ax2.spines['top'].set_visible(False)
  ax2.spines['right'].set_visible(False)
  ax2.spines['left'].set_visible(False)
  plt.title('Go cue alignment')

  plt.subplot(236)
  grandAveragePlot(data_struct, neurons = [neuron], alignment = data_struct.timeSpeech - cue_shift, numbins = num_bins)
  plt.axvline(cue_shift, color = 'k', linestyle = '--')
  plt.ylim(bottom = ymin, top = ymax)
  ax3 = plt.gca()
  ax3.get_yaxis().set_visible(False)
  ax3.spines['top'].set_visible(False)
  ax3.spines['right'].set_visible(False)
  ax3.spines['left'].set_visible(False)
  plt.title('Speech onset alignment')

  kwargs = dict(transform=ax1.transAxes, color='k', clip_on=False)
  ax1.plot((1-d,1+d), (-d,+d), **kwargs)

  kwargs.update(transform=ax2.transAxes) 
  ax2.plot((-d,+d), (-d,+d), **kwargs)
  ax2.plot((1-d,1+d), (-d,+d), **kwargs)

  kwargs.update(transform=ax3.transAxes) 
  ax3.plot((-d,+d), (-d,+d), **kwargs)
  fig = plt.gcf()
  
  numlabels = (num_bins/100) + 2
  labels    = np.linspace(-1 * cue_shift - 100, num_bins - cue_shift, numlabels).astype('int')
  labels    = labels.astype('str').tolist()
  ax2.set_xticklabels(labels)
  ax3.set_xticklabels(labels)
  plt.tight_layout()
  
  if toggleSave:
    plt.savefig(path + str(neuron) + '.png')
  else:
    plt.show()
  
  
  
  
  