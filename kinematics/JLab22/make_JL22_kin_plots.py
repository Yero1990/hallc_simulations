import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# read dataframe
df_11 = pd.read_csv('kin_summary_Eb11.00.txt', comment='#')
#df_12 = pd.read_csv('kin_summary_Eb12.00.txt', comment='#')
#df_14 = pd.read_csv('kin_summary_Eb14.00.txt', comment='#')
#df_16 = pd.read_csv('kin_summary_Eb16.00.txt', comment='#')
#df_18 = pd.read_csv('kin_summary_Eb18.00.txt', comment='#')
#df_20 = pd.read_csv('kin_summary_Eb20.00.txt', comment='#')
df_22 = pd.read_csv('kin_summary_Eb22.00.txt', comment='#')


fig, ax  = plt.subplots(3, 3)

# kf vs. th_e
ax[0,0].plot( df_11['th_e'] , df_11['kf'], linestyle='', marker='o', mec='g', mfc='w', alpha=0.8)
ax[0,0].plot( df_22['th_e'] , df_22['kf'], linestyle='', marker='s', mec='r', mfc='w', alpha=0.8)

ax[0,0].set(xlabel='SHMS Angle [deg]', ylabel='SHMS Momentum [GeV/c]')

# Q2 vs. th_e
ax[0,1].plot( df_11['th_e'] , df_11['Q2'], linestyle='', marker='o', mec='g', mfc='w', alpha=0.8)
ax[0,1].plot( df_22['th_e'] , df_22['Q2'], linestyle='', marker='s', mec='r', mfc='w', alpha=0.8)
ax[0,1].set(xlabel='SHMS Angle [deg]', ylabel=r'$Q^{2}$ [GeV$^{2}$]')

# xbj vs. th_e
ax[0,2].plot( df_11['th_e'] , df_11['xbj'], linestyle='', marker='o', mec='g', mfc='w', alpha=0.8)
ax[0,2].plot( df_22['th_e'] , df_22['xbj'], linestyle='', marker='s', mec='r', mfc='w', alpha=0.8)
ax[0,2].set(xlabel='SHMS Angle [deg]', ylabel=r'$x_{Bj}$')

# Pr vs. th_e
ax[1,0].plot( df_11['th_e'] , df_11['Pr'], linestyle='', marker='o', mec='g', mfc='w', alpha=0.8)
ax[1,0].plot( df_22['th_e'] , df_22['Pr'], linestyle='', marker='s', mec='r', mfc='w', alpha=0.8)
ax[1,0].set(xlabel='SHMS Angle [deg]', ylabel=r'Recoil Momentum, $P_{r}$ [GeV/c]')

# th_rq vs. th_e
ax[1,1].plot( df_11['th_e'] , df_11['th_rq'], linestyle='', marker='o', mec='g', mfc='w', alpha=0.8)
ax[1,1].plot( df_22['th_e'] , df_22['th_rq'], linestyle='', marker='s', mec='r', mfc='w', alpha=0.8)
ax[1,1].set(xlabel='SHMS Angle [deg]', ylabel=r'$\theta_{rq}$ [deg]')

# Pf vs th_p
ax[1,2].plot( df_11['th_e'], df_11['Pf'], linestyle='', marker='o', mec='g', mfc='w', alpha=0.8)
ax[1,2].plot( df_22['th_e'], df_22['Pf'], linestyle='', marker='s', mec='r', mfc='w', alpha=0.8)
ax[1,2].set(xlabel='SHMS Angle [deg]', ylabel='HMS Momentum [GeV/c]')

# Pf vs th_p
ax[2,0].plot( df_11['th_p'], df_11['Pf'], linestyle='', marker='o', mec='g', mfc='w', alpha=0.8)
ax[2,0].plot( df_22['th_p'], df_22['Pf'], linestyle='', marker='s', mec='r', mfc='w', alpha=0.8)
ax[2,0].set(xlabel='HMS Angle [deg]', ylabel='HMS Momentum [GeV/c]')

# Q2 vs xbj
ax[2,1].plot( df_11['xbj'], df_11['Q2'], linestyle='', marker='o', mec='g', mfc='w', alpha=0.8)
ax[2,1].plot( df_22['xbj'], df_22['Q2'], linestyle='', marker='s', mec='r', mfc='w', alpha=0.8)
ax[2,1].set(xlabel=r'$x_{bj}$', ylabel=r'$Q^{2}$ [GeV$^{2}$]')

# Pr vs th_rq
ax[2,2].plot( df_11['th_rq'], df_11['Pr'], linestyle='', marker='o', mec='g', mfc='w', alpha=0.8, label=r'E$_{b}$=11 GeV')
ax[2,2].plot( df_22['th_rq'], df_22['Pr'], linestyle='', marker='s', mec='r', mfc='w', alpha=0.8, label=r'E$_{b}$=22 GeV')
ax[2,2].set(xlabel=r'$\theta_{rq}$ [deg]', ylabel=r'Recoil Momentum, $P_{r}$ [GeV/c]')



plt.tight_layout()
plt.legend()
plt.show()
