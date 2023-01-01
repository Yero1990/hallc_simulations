import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# read dataframe
df_12 = pd.read_csv('kin_summary_Eb12.00.txt', comment='#')
df_14 = pd.read_csv('kin_summary_Eb14.00.txt', comment='#')
df_16 = pd.read_csv('kin_summary_Eb16.00.txt', comment='#')
df_18 = pd.read_csv('kin_summary_Eb18.00.txt', comment='#')
df_20 = pd.read_csv('kin_summary_Eb20.00.txt', comment='#')
df_22 = pd.read_csv('kin_summary_Eb22.00.txt', comment='#')


fig, ax  = plt.subplots(2, 2)

# th_e vs kf
ax[0,0].plot( df_12['th_e'] , df_12['kf'], linestyle='', marker='o', mec='r', mfc='w', alpha=0.8)
ax[0,0].plot( df_14['th_e'] , df_14['kf'], linestyle='', marker='o', mec='b', mfc='w', alpha=0.8)

#ax[0,0].plot( df_16['th_e'] , df_16['kf'], linestyle='', marker='o', mec='orange', mfc='w', alpha=0.8)
#ax[0,0].plot( df_18['th_e'] , df_18['kf'], linestyle='', marker='o', mec='r', mfc='w', alpha=0.8)
#ax[0,0].plot( df_20['th_e'] , df_20['kf'], linestyle='', marker='o', mec='g', mfc='w', alpha=0.8)
#ax[0,0].plot( df_22['th_e'] , df_22['kf'], linestyle='', marker='o', mec='b', mfc='w', alpha=0.8)
ax[0,0].set(xlabel='SHMS Angle [deg]', ylabel='SHMS Momentum [GeV/c]')


# Pf vs th_p
ax[0,1].plot( df_12['th_p'], df_12['Pf'], linestyle='', marker='o', mec='r', mfc='w', alpha=0.8)
ax[0,1].plot( df_14['th_p'], df_14['Pf'], linestyle='', marker='o', mec='b', mfc='w', alpha=0.8)

#ax[0,1].plot( df_16['th_p'], df_16['Pf'], linestyle='', marker='o', mec='orange', mfc='w', alpha=0.8)
#ax[0,1].plot( df_18['th_p'], df_18['Pf'], linestyle='', marker='o', mec='r', mfc='w', alpha=0.8)
#ax[0,1].plot( df_20['th_p'], df_20['Pf'], linestyle='', marker='o', mec='g', mfc='w', alpha=0.8)
#ax[0,1].plot( df_22['th_p'], df_22['Pf'], linestyle='', marker='o', mec='b', mfc='w', alpha=0.8)
ax[0,1].set(xlabel='HMS Angle [deg]', ylabel='HMS Momentum [GeV/c]')


# Q2 vs xbj
ax[1,0].plot( df_12['xbj'], df_12['Q2'], linestyle='', marker='o', mec='r', mfc='w', alpha=0.8)
ax[1,0].plot( df_14['xbj'], df_14['Q2'], linestyle='', marker='o', mec='b', mfc='w', alpha=0.8)

#ax[1,0].plot( df_16['xbj'], df_16['Q2'], linestyle='', marker='o', mec='orange', mfc='w', alpha=0.8)
#ax[1,0].plot( df_18['xbj'], df_18['Q2'], linestyle='', marker='o', mec='r', mfc='w', alpha=0.8)
#ax[1,0].plot( df_20['xbj'], df_20['Q2'], linestyle='', marker='o', mec='g', mfc='w', alpha=0.8)
#ax[1,0].plot( df_22['xbj'], df_22['Q2'], linestyle='', marker='o', mec='b', mfc='w', alpha=0.8)
ax[1,0].set(xlabel=r'$x_{bj}$', ylabel=r'$Q^{2}$ [GeV$^{2}$]')

# Pr vs th_rq
ax[1,1].plot( df_12['th_rq'], df_12['Pr'], linestyle='', marker='o', mec='r', mfc='w', alpha=0.8, label=r'E$_{b}$=12 GeV')
ax[1,1].plot( df_14['th_rq'], df_14['Pr'], linestyle='', marker='o', mec='b', mfc='w', alpha=0.8, label=r'E$_{b}$=14 GeV')

#ax[1,1].plot( df_16['th_rq'], df_16['Pr'], linestyle='', marker='o', mec='orange', mfc='w', alpha=0.8, label=r'E$_{b}$=16 GeV')
#ax[1,1].plot( df_18['th_rq'], df_18['Pr'], linestyle='', marker='o', mec='r', mfc='w', alpha=0.8, label=r'E$_{b}$=18 GeV')
#ax[1,1].plot( df_20['th_rq'], df_20['Pr'], linestyle='', marker='o', mec='g', mfc='w', alpha=0.8, label=r'E$_{b}$=20 GeV')
#ax[1,1].plot( df_22['th_rq'], df_22['Pr'], linestyle='', marker='o', mec='b', mfc='w', alpha=0.8, label=r'E$_{b}$=22 GeV')
ax[1,1].set(xlabel=r'$\theta_{rq}$ [deg]', ylabel=r'Recoil Momentum, $P_{r}$ [GeV/c]')



plt.tight_layout()
plt.legend()
plt.show()
