import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# read dataframe
df_11 = pd.read_csv('kin_summary_Eb11.00.txt', comment='#')
df_12 = pd.read_csv('kin_summary_Eb12.00.txt', comment='#')


fig, ax  = plt.subplots(1, 3)

# th_e vs kf
ax[0].plot( (df_11[df_11['Pr']==1])['th_e'] , (df_11[df_11['Pr']==1])['kf'], linestyle='', marker='o', mec='r', mfc='w')
ax[0].plot( (df_11[df_11['Pr']==1.2])['th_e'] , (df_11[df_11['Pr']==1.2])['kf'], linestyle='', marker='^', mec='r', mfc='w')

ax[0].plot( (df_12[df_12['Pr']==1])['th_e'] , (df_12[df_12['Pr']==1])['kf'], linestyle='', marker='o', mec='b', mfc='w')
ax[0].plot( (df_12[df_12['Pr']==1.2])['th_e'] , (df_12[df_12['Pr']==1.2])['kf'], linestyle='', marker='^', mec='b', mfc='w')

ax[0].set(xlabel='SHMS Angle [deg]', ylabel='SHMS Momentum [GeV/c]')


# Pf vs th_p
ax[1].plot( (df_11[df_11['Pr']==1])['th_p'] , (df_11[df_11['Pr']==1])['Pf'], linestyle='', marker='o', mec='r', mfc='w')
ax[1].plot( (df_11[df_11['Pr']==1.2])['th_p'] , (df_11[df_11['Pr']==1.2])['Pf'], linestyle='', marker='^', mec='r', mfc='w')

ax[1].plot( (df_12[df_12['Pr']==1])['th_p'] , (df_12[df_12['Pr']==1])['Pf'], linestyle='', marker='o', mec='b', mfc='w')
ax[1].plot( (df_12[df_12['Pr']==1.2])['th_p'] , (df_12[df_12['Pr']==1.2])['Pf'], linestyle='', marker='^', mec='b', mfc='w')


ax[1].set(xlabel='HMS Angle [deg]', ylabel='HMS Momentum [GeV/c]')


# xbj vs th_rq
ax[2].plot( (df_11[df_11['Pr']==1])['th_rq'] , (df_11[df_11['Pr']==1])['xbj'], linestyle='', marker='o', mec='r', mfc='w', label='Eb = 11 GeV, Pr=1 GeV/c')
ax[2].plot( (df_11[df_11['Pr']==1.2])['th_rq'] , (df_11[df_11['Pr']==1.2])['xbj'], linestyle='', marker='^', mec='r', mfc='w', label='Eb = 11 GeV, Pr=1.2 GeV/c')

ax[2].plot( (df_12[df_12['Pr']==1])['th_rq'] , (df_12[df_12['Pr']==1])['xbj'], linestyle='', marker='o', mec='b', mfc='w', label='Eb = 12 GeV, Pr=1 GeV/c')
ax[2].plot( (df_12[df_12['Pr']==1.2])['th_rq'] , (df_12[df_12['Pr']==1.2])['xbj'], linestyle='', marker='^', mec='b', mfc='w', label='Eb = 12 GeV, Pr=1.2 GeV/c')

ax[2].set(xlabel=r'$\theta_{rq}$ [deg]', ylabel=r'$x_{bj}$')


plt.tight_layout()
plt.legend()
plt.show()
