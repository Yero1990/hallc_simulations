#This code reads the theory data files (grouped in th_nq bins) and divided by
#the corresponding Ksig_cc1 to get the reduced Xsec


import sys
import os
import numpy as np
import pandas as pd
import numpy.ma as ma
from sys import argv  
import matplotlib.colors as mcolors
import matplotlib.pyplot as plt

'''
example of selecting a pm bin and plotting the FSI+PWIA/PWIA vs. thrq for a file

df = pd.read_csv('csec_calc_thrq28_2_1_1_0_12.data', comment='#')


thrq_pm_bin = (df['th_nq_mc'])[(df['pr']>0.820) & (df['pr']<0.860) ]

ratio_pm_bin = (df['ratio'])[(df['pr']>0.820) & (df['pr']<0.860) ]

# convert to numpy arrays
thrq_pm_bin.to_numpy()
ratio_pm_bin.to_numpy()

plt.plot(thrq, ratio, marker='o', linestyle='None')
plt.show()

'''
#Conversion Factors
dtr = np.pi / 180.
nb2ub = 1./1000.    # 1ub = 1000 nb
GeV2MeV = 1000.  # 1 GeV = 1000 MeV


# the actual PWIA model is within the data file, in a column named:
# crs0 : PWIA
# crs12 : PWIA + FSI
# ratio: crs12 /  crs0 
model_set = ["2_1_1_0_12", "3_1_1_0_12"]   #V18, CD-Bonn [the '12' stands for PWIA+FSI]

# central recoil angle
thrq_set = [28, 49, 55, 60, 66, 72, 79]

# set the central pm bin +/- bin width condition for plotting angular distribution
# (this is esentially to select the angular distirbution for  pmiss slice)

# pm_c +/- pm_bw [GeV/c]
pm_c = 0.840 
pm_bw = 0.02 


# loop over each model (V18, CD-Bonn)
for model in model_set:

    # loop over central recoil angle kin. setting
    for ithrq in thrq_set: 
        print('ithrq: ', ithrq)
        basename = 'q4_sig_avkin_thnq_pm/csec_calc_thnq%d_%s.data' % (ithrq, model)

        df = pd.read_csv(basename, comment='#')

        pm_min = pm_c - pm_bw
        pm_max = pm_c + pm_bw
        
        
        thrq  = ((df['th_nq_mc'])[(df['pr']>pm_min) & (df['pr']<pm_max) ]).to_numpy()
        ratio = ((df['ratio'])[(df['pr']>pm_min) & (df['pr']<pm_max) ]).to_numpy()

        plt.plot(thrq, ratio, marker='o', linestyle='None')

    plt.show()
