import numpy as np
import numpy.ma as ma
import pandas as pd
import matplotlib.pyplot as plt


thrq_c = [28] #deg
pm_c = [800] # MeV/c

# loop over central recoil angles files for a given central momentum
for i in thrq_c:
    
    # read dataframe
    df_fsi = pd.read_csv('d2_pm800_thrq%i_fsi_rad_yield_80.0uA_1.0hr.csv'%(i), comment='#')
    df_pwia = pd.read_csv('d2_pm800_thrq%i_pwia_rad_yield_80.0uA_1.0hr.csv'%(i), comment='#')

    # loop over central missing momentum bin for a given thrq_c file
    for pm_bin in  df_fsi.y0:
    
        # get data columns and put conditions
        # x0 -> thrq, y0 -> pmiss_central +/- (bin_width/2), depends on which histogram was used

        # angular distributions (x0) for a given pmiss_central (y0)
        thrq_fsi = df_fsi.x0[df_fsi.y0==pm_bin]
        thrq_pwia = df_pwia.x0[df_pwia.y0==pm_bin]
        
        # counts  corresponding to pmiss_central (y0)
        fsi_N = df_fsi.zcont[df_fsi.y0==pm_bin]
        pwia_N = df_pwia.zcont[df_pwia.y0==pm_bin]
        
        # error in counts  corresponding to pmiss_central (y0)
        #thrq_fsi_Nerr = df_fsi.zcont_err[df_fsi.y0==0.7]
        #thrq_pwia_Nerr = df_pwia.zcont_err[df_pwia.y0==0.7]
        
        # mask zero counts (masked values will not be plotted)
        #thrq_0p7_fsi_N_mask = ma.masked_where(thrq_0p7_fsi_N==0, thrq_0p7_fsi_N)
        #thrq_0p7_pwia_N_mask = ma.masked_where(thrq_0p7_pwia_N==0, thrq_0p7_pwia_N)
        
        thrq_ratio = fsi_N / pwia_N
        
        #thrq_0p7_fsi_N_mask_rel_err = np.sqrt(thrq_0p7_fsi_N_mask) / thrq_0p7_fsi_N_mask
        #thrq_0p7_pwia_N_mask_rel_err = np.sqrt(thrq_0p7_pwia_N_mask) / thrq_0p7_pwia_N_mask
        
        #thrq_0p7_ratio_err = thrq_0p7_ratio * np.sqrt(thrq_0p7_fsi_N_mask_rel_err**2 + thrq_0p7_pwia_N_mask_rel_err**2)


        
plt.plot(thrq_0p7_, thrq_0p7_ratio, marker='o', linestyle='None')
plt.show()
