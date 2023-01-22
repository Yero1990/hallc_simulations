import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# read dataframe
#df_11 = pd.read_csv('kin_summary_Eb11.00.txt', comment='#')

#df_12 = pd.read_csv('kin_summary_Eb12.00.txt', comment='#')

#df_13 = pd.read_csv('kin_summary_Eb13.00.txt', comment='#')

#df_14 = pd.read_csv('kin_summary_Eb14.00.txt', comment='#')
#df_16 = pd.read_csv('kin_summary_Eb16.00.txt', comment='#')
#df_18 = pd.read_csv('kin_summary_Eb18.00.txt', comment='#')
#df_20 = pd.read_csv('kin_summary_Eb20.00.txt', comment='#')
#df_22 = pd.read_csv('kin_summary_Eb22.00.txt', comment='#')


'''
# single file read
df = df_14

fig, ax  = plt.subplots(3, 3)

# kf vs. th_e
ax[0,0].plot( df['th_e'] , df['kf'], linestyle='', marker='^', mec='b', mfc='w', alpha=0.8) 
ax[0,0].set(xlabel='SHMS Angle [deg]', ylabel='SHMS Momentum [GeV/c]')

# Q2 vs. th_e
ax[0,1].plot( df['th_e'] , df['Q2'], linestyle='', marker='^', mec='b', mfc='w', alpha=0.8)
ax[0,1].set(xlabel='SHMS Angle [deg]', ylabel=r'$Q^{2}$ [GeV$^{2}$]')


# xbj vs. th_e
ax[0,2].plot( df['th_e'] , df['xbj'], linestyle='', marker='^', mec='b', mfc='w', alpha=0.8)  
ax[0,2].set(xlabel='SHMS Angle [deg]', ylabel=r'$x_{Bj}$')

# Pr vs. th_e
ax[1,0].plot( df['th_e'] , df['Pr'], linestyle='', marker='^', mec='b', mfc='w', alpha=0.8) 
ax[1,0].set(xlabel='SHMS Angle [deg]', ylabel=r'Recoil Momentum, $P_{r}$ [GeV/c]')

# th_rq vs. th_e
ax[1,1].plot( df['th_e'] , df['th_rq'], linestyle='', marker='^', mec='b', mfc='w', alpha=0.8) 
ax[1,1].set(xlabel='SHMS Angle [deg]', ylabel=r'$\theta_{rq}$ [deg]')

# Pf vs th_e
ax[1,2].plot( df['th_e'], df['Pf'], linestyle='', marker='^', mec='b', mfc='w', alpha=0.8)
ax[1,2].set(xlabel='SHMS Angle [deg]', ylabel='HMS Momentum [GeV/c]')

# Pf vs th_p
ax[2,0].plot( df['th_p'], df['Pf'], linestyle='', marker='^', mec='b', mfc='w', alpha=0.8)
ax[2,0].set(xlabel='HMS Angle [deg]', ylabel='HMS Momentum [GeV/c]')

# Q2 vs xbj
ax[2,1].plot( df['xbj'], df['Q2'], linestyle='', marker='^', mec='b', mfc='w', alpha=0.8)
ax[2,1].set(xlabel=r'$x_{bj}$', ylabel=r'$Q^{2}$ [GeV$^{2}$]')

# Pr vs th_rq
ax[2,2].plot( df['th_rq'], df['Pr'], linestyle='', marker='^', mec='b', mfc='w', alpha=0.8, label=r'E$_{b}$=14 GeV')  
ax[2,2].set(xlabel=r'$\theta_{rq}$ [deg]', ylabel=r'Recoil Momentum, $P_{r}$ [GeV/c]')
'''



E = [11,12,  14 , 16, 18, 20, 22]
clr  = ['grey', 'c',    'm',   'r',   'g',   'b', 'darkorange', 'violet', 'gold', 'lightcoral', 'olive', 'sandybrown'] #'darkgray']

#multi-file read
#------------------------
fig, ax  = plt.subplots(3, 3)

for i in np.arange(len(E)):

    print(E[i])
    df = pd.read_csv('idealistic_kinematics/Q2_4p5/kin_summary_Eb%.2f.txt'%(E[i]), comment='#')
    

    # kf vs. th_e
    if i==0:
        ax[0,0].axhline(y=2.0, c='red', linestyle='dashed', linewidth=2, zorder=0)
        ax[0,0].axhline(y=11.0, c='red', linestyle='dashed', linewidth=2, zorder=0)
        ax[0,0].axvline(x=5.5, c='red', linestyle='dashed', linewidth=2, zorder=0)
        #ax[0,0].axvline(x=40., c='red', linestyle='dashed', linewidth=2, zorder=0)

    ax[0,0].plot( df['th_e'] , df['kf'], linestyle='', marker='o', mec='k', mfc=clr[i], alpha=0.8)
    ax[0,0].set(xlabel='SHMS Angle [deg]', ylabel='SHMS Momentum [GeV/c]')

    
    # Q2 vs. th_e
    ax[0,1].plot( df['th_e'] , df['Q2'], linestyle='', marker='o', mec='k', mfc=clr[i], alpha=0.8)
    ax[0,1].set(xlabel='SHMS Angle [deg]', ylabel=r'$Q^{2}$ [GeV$^{2}$]')
    
    # xbj vs. th_e
    ax[0,2].plot( df['th_e'] , df['xbj'], linestyle='', marker='o', mec='k', mfc=clr[i], alpha=0.8)
    ax[0,2].set(xlabel='SHMS Angle [deg]', ylabel=r'$x_{Bj}$')
    
    # Pr vs. th_e
    ax[1,0].plot( df['th_e'] , df['Pr'], linestyle='', marker='o', mec='k', mfc=clr[i], alpha=0.8)
    ax[1,0].set(xlabel='SHMS Angle [deg]', ylabel=r'Recoil Momentum, $P_{r}$ [GeV/c]')
    
    # th_rq vs. th_e
    ax[1,1].plot( df['th_e'] , df['th_rq'], linestyle='', marker='o', mec='k', mfc=clr[i], alpha=0.8)
    ax[1,1].set(xlabel='SHMS Angle [deg]', ylabel=r'$\theta_{rq}$ [deg]')
    
    # Pf vs th_e
    ax[1,2].plot( df['th_e'], df['Pf'], linestyle='', marker='o', mec='k', mfc=clr[i], alpha=0.8)
    ax[1,2].set(xlabel='SHMS Angle [deg]', ylabel='HMS Momentum [GeV/c]')
    
    # Pf vs th_p
    ax[2,0].plot( df['th_p'], df['Pf'], linestyle='', marker='o', mec='k', mfc=clr[i], alpha=0.8)
    ax[2,0].set(xlabel='HMS Angle [deg]', ylabel='HMS Momentum [GeV/c]')
    
    # Q2 vs xbj
    ax[2,1].plot( df['xbj'], df['Q2'], linestyle='', marker='o', mec='k', mfc=clr[i], alpha=0.8)
    ax[2,1].set(xlabel=r'$x_{bj}$', ylabel=r'$Q^{2}$ [GeV$^{2}$]')
    
    # Pr vs th_rq
    ax[2,2].plot( df['th_rq'], df['Pr'], linestyle='', marker='o', mec='k', mfc=clr[i], alpha=0.8, label=r'E$_{b}$=%i GeV'%(E[i]))
    ax[2,2].set(xlabel=r'$\theta_{rq}$ [deg]', ylabel=r'Recoil Momentum, $P_{r}$ [GeV/c]')



plt.tight_layout()
plt.legend(bbox_to_anchor=(1.05, 1.0), loc='upper left')
plt.show()
