import numpy as np
import numpy.ma as ma
import pandas as pd
import matplotlib.pyplot as plt
from uncertainties import unumpy

# This script reads both (pwia, fsi) for a given th_rq setting, and overlays the
# FSI / PWIA for all the central thrq_c, binned in missing momentum (pm_bin)

thrq_c = [28, 49, ] #deg
pm_c = [800] # MeV/c


fig, ax = plt.subplots(5, 8, sharex='col', sharey='row')
fig.text(0.5, 0.01, r'Recoil Angle $\theta_{rq}$ [deg]', ha='center')
fig.text(0.01, 0.5, r'R = FSI / PWIA', va='center', rotation='vertical')
fig.set_size_inches(14,8, forward=True)
    
# loop over central recoil angles files for a given central momentum
for i in thrq_c:
    
    # read dataframe
    df_fsi = pd.read_csv('d2_pm800_thrq%i_fsi_rad_yield_80.0uA_1.0hr.csv'%(i), comment='#')
    df_pwia = pd.read_csv('d2_pm800_thrq%i_pwia_rad_yield_80.0uA_1.0hr.csv'%(i), comment='#')

    thrq_bin = df_fsi.x0
    pm_bins = df_fsi.y0[thrq_bin==thrq_bin[0]] # only use a set of pm_bins (avoid duplicates)

    pm_binw = (df_fsi.yup[0] - df_fsi.ylow[0])
    
    fsi_N = df_fsi.zcont
    fsi_Nerr = df_fsi.zcont_err

    pwia_N = df_pwia.zcont
    pwia_Nerr = df_pwia.zcont_err

    ratio = fsi_N / pwia_N
    ratio_err = ratio * np.sqrt((fsi_Nerr/fsi_N)**2 + (pwia_Nerr/pwia_N)**2)

    rel_err = ratio_err / ratio

    rel_err_thrs = 0.3 # mask >30 % relative error
    ratio = ma.masked_where(rel_err>rel_err_thrs, ratio)
    ratio_err = ma.masked_where(rel_err>rel_err_thrs, ratio_err)


    for idx, pm_bin in enumerate(pm_bins):
        print(idx)

        ax = plt.subplot(5, 8, idx+1)
        ax.errorbar(thrq_bin[df_fsi.y0==pm_bin], ratio[df_fsi.y0==pm_bin], ratio_err[df_fsi.y0==pm_bin], marker='o', linestyle='None')

        if i==thrq_c[0]:
            ax.set_title(r"%.2f $\pm$ %.2f MeV"%(pm_bin, pm_binw/2.), fontsize=10)


plt.tight_layout()
plt.show()

  
        
