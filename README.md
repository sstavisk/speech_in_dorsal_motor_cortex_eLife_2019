# speech_in_dorsal_motor_cortex_eLife_2019
Generates the analyses and figures of "Neural ensemble dynamics in dorsal motor cortex during speech in people with paralysis", eLife 2019, by Sergey D. Stavisky, Francis R. Willett, Guy H. Wilson, Brian A. Murphy, Paymon Rezaii, Donald T. Avansino, William D. Memberg, Jonathan P. Miller, Robert F. Kirsch, Leigh R. Hochberg, A. Bolu Ajiboye, Shaul Druckmann, Krishna V. Shenoy, Jaimie M. Henderson.

The scripts in the top directory are named according to what figure in the paper they relate to. For example, figure1_and_2_firing_rates.m will generate firing rate plots for specific neurons or electrodes, as in Figure 1 of the paper.

Dependencies are listed at the top of each script. In general, you will need to put /helper_functions and its subdirectories to your MATLAB path. A few of the figure scripts require additional code packs (for example, for dPCA and jPCA); I've provided the URLs and paper references for that code. 

These scripts will not run without the underlying study datasets, which will be hosted by eLife. A list of all the requisite datasets is provided in list_of_datasets.txt.

Sergey D. Stavisky, December 2019, Stanford University
