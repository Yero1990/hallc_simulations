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
from scipy.interpolate import interp1d

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
model_name     = ["V18", "CD-Bonn"]
model_set = ["2_1_1_0_12", "3_1_1_0_12"]   #[V18, CD-Bonn] [the '12' stands for PWIA+FSI]
#model_set = ["2_1_1_0_123"]   #V18, CD-Bonn [the '12' stands for PWIA+FSI]
#model_set = ["3_1_1_0_12"]   #V18, CD-Bonn [the '12' stands for PWIA+FSI]

# central recoil angle
#thrq_set = [28, 49, 55, 60, 66, 72, 79]
thrq_set = [49, 60, 72]
pm_set = [520, 560, 600, 640, 680, 720, 760, 800, 840, 880, 920, 960]


# set the central pm bin +/- bin width condition for plotting angular distribution
# (this is esentially to select the angular distirbution for  pmiss slice)

# pm_c +/- pm_bw [GeV/c]
#pm_c = 0.760
pm_bw = 0.02 

fig, ax = plt.subplots(4, 3)

idx_plot = 0 
# loop over each pm bin
for ipm_c in pm_set:
    # --> 1 TAB
    
    # define a common x range (that includes full th_rq range for multiple interpolations)
    x = np.linspace(0, 90, num=200, endpoint=True)

    imod = 0
    
    # loop over each model (V18, CD-Bonn)
    for model in model_set:

        # --> 1 TAB 
        theory_calc_name = model_name[imod]
        imod = imod+1
        
        total_f = []

        col = ['b', 'g', 'r']
        icol = 0
        # loop over central recoil angle kin. setting
        for ithrq in thrq_set:

            # --> 1 TAB 
            print('ithrq: ', ithrq)
            basename = 'q4_sig_avkin_thnq_pm/csv/csec_calc_thrq%d_%s.data' % (ithrq, model)

            df = pd.read_csv(basename, comment='#')
            pm_c = ipm_c/1000.
            pm_min = pm_c - pm_bw
            pm_max = pm_c + pm_bw

            print('pm_c, pm_min, pm_max', pm_c, pm_min, pm_max)

            # for central kin thrq=72, need to cut out thrq<70, to make more smooth transition
            if ithrq==72:   
                thrq  = ((df['th_nq_mc'])[(df['pr']>pm_min) & (df['pr']<pm_max) & (df['th_nq_mc']>60.)]).to_numpy()
                ratio = ((df['ratio'])[(df['pr']>pm_min) & (df['pr']<pm_max) & (df['th_nq_mc']>60.) ]).to_numpy()
                
            else:
                thrq  = ((df['th_nq_mc'])[(df['pr']>pm_min) & (df['pr']<pm_max) ]).to_numpy()
                ratio = ((df['ratio'])[(df['pr']>pm_min) & (df['pr']<pm_max) ]).to_numpy()

            print('ithrq = ', ithrq)
            print('thrq = ', thrq)
            # interpolate data
            f_ratio = interp1d(thrq, ratio, kind='linear', fill_value=np.nan, bounds_error=False)
            
            # append the interpolated function array to a total array
            total_f.append(f_ratio(x))

            print('idx_plot = ', idx_plot)
            ax = plt.subplot(4, 3, idx_plot+1)

            #if(model=="2_1_1_0_12"):   # AV18
                # plot the different thrq calculations separately (before avergaing)
            #    ax.plot(x, f_ratio(x), marker='None', color=col[icol], linestyle='--', label=r'$\theta_{nq} = %d$'%ithrq)
            #if(model=="3_1_1_0_12"):   # CD-Bonn
                # plot the different thrq calculations separately (before avergaing)
            #    ax.plot(x, f_ratio(x), marker='None', color=col[icol], linestyle='-', label=r'$\theta_{nq} = %d$'%ithrq)
            icol = icol+1
            
    
        # <-- TAB
        # Define the common x-values for averaging
        
        
        total_f_m0 = np.ma.masked_array(total_f[0], np.isnan(total_f[0]))
        total_f_m1 = np.ma.masked_array(total_f[1], np.isnan(total_f[1]))
        total_f_m2 = np.ma.masked_array(total_f[2], np.isnan(total_f[2]))

        data = []
        data.append(total_f_m0)
        data.append(total_f_m1)
        data.append(total_f_m2)
        
        total_f_avg = np.ma.average(data, axis=0)

        print(total_f_avg)
        f_ratio_total = interp1d(x, total_f_avg, kind='cubic')

       
        
        print('f_ratio_total = ', f_ratio_total(x))
        # plot the avergaed calculations for overlapping thrq
        #if(model=="2_1_1_0_12"):
        #    plt.plot(x, f_ratio_total(x), marker='None', color='k', linestyle='--', label=r' %s'%theory_calc_name)
        if(model=="3_1_1_0_12"):
            plt.plot(x, f_ratio_total(x), marker='None', color='k', linestyle='-', label=r' %s'%theory_calc_name)

    idx_plot = idx_plot+1
        
    
plt.legend()
plt.show()




'''
#--------------------------------------------------------------------------------------
# Here's how to average two interpolated functions in Python using scipy.interpolate.
#--------------------------------------------------------------------------------------

import numpy as np
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt

# Define the data points for the two functions
x1 = np.array([1, 2, 3, 4, 5])
y1 = np.array([2, 4, 1, 3, 5])
x2 = np.array([1.5, 2.5, 3.5, 4.5])
y2 = np.array([3, 1, 4, 2])

# Create the interpolation functions
f1 = interp1d(x1, y1, kind='linear')
f2 = interp1d(x2, y2, kind='linear', fill_value="extrapolate")

# Define the common x-values for averaging
x_new = np.linspace(max(min(x1), min(x2)), max(max(x1), max(x2)), num=100)

# Calculate the average of the two functions
y_avg = (f1(x_new) + f2(x_new)) / 2

# Plot the results
plt.plot(x1, y1, 'o', label='Data 1')
plt.plot(x2, y2, 'x', label='Data 2')
plt.plot(x_new, f1(x_new), '-', label='Interpolation 1')
plt.plot(x_new, f2(x_new), '--', label='Interpolation 2')
plt.plot(x_new, y_avg, '-', label='Average Interpolation')
plt.legend()
plt.xlabel('x')
plt.ylabel('y')
plt.title('Average of Two Interpolated Functions')
plt.show()
'''
