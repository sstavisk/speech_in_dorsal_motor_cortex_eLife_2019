import numpy as np


def getAllPhonemes(struct, start_offset, stop_offset, addNoise = 0, source = 'FR'):
  '''
  Given R struct with segmentations, pull out neural/EMA data for all phonemes.
  Inputs are:
  
    struct (datastruct) - holds data to process
    start_offset (int)  - pull activity beginning this far from phoneme onset
    stop_offset (int)   - grab data from this much past phoneme onset 
    addNoise (float)    - if > 0, add gaussian noise to our phoneme onset measurements
    source (str)        - whether to use FR, raster, HFLP, or EMA signal (default: FR)
  
  
  Outputs is a list of entries:
  
    [0] - phoneme ID (string denoting which one)
    [1] - trial # frrom which phoneme is coming
    [2] - phoneme order within a word (1 = first, 2 = second, etc)
    [3] - neural/EMA data array of form (channels x time)
    
  where all entries are of length = (total # phonemes) 
  '''
  
  phonemes_data_table = list()

  phonemes_id    = list()  # list of strings tracking phoneme type
  phonemes_trl   = list()  # list of ints tracking trial ID
  phonemes_order = list()  # order within word 
  phonemes_FR    = list()  # list of neurons x time arrays containing onset-related activity across units 

  phonemes_id = list()
  for trl in range(len(struct.speechLabel)):
      for i, ph in enumerate(struct.phonemes[trl]):
          
          start = struct.starts[trl][i] + start_offset
          stop  = struct.starts[trl][i] + stop_offset
          
          if addNoise > 0:
            offset = int(np.random.normal(0, addNoise))
            start += offset
            stop  += offset

          phonemes_id.append(ph)
          phonemes_trl.append(trl)
          phonemes_order.append(i + 1)
          if source == 'FR':
            phonemes_FR.append(struct.FR[trl][:, start:stop].copy()) 
          if source == 'spike':
            phonemes_FR.append(struct.spikeRaster[trl][start:stop, :].T.copy()) 
          if source == 'HLFP':
            phonemes_FR.append(struct.HLFP[trl][start:stop, :].T.copy()) 
          if source == 'EMA':
            phonemes_FR.append(struct.AKT[trl][:, start:stop].copy())


  phonemes_data_table.extend([np.asarray(phonemes_id), np.asarray(phonemes_trl), np.asarray(phonemes_order), np.dstack(phonemes_FR)])
  
  return phonemes_data_table
