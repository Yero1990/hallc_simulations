import os
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# ipython make_plots_d2po.py 200,300 4.5 H_Q2   --> passes pmiss = 200 and 300 to pm_user and Q2 = 4.5 to q2_user arrays

# --- user arguments ---
if len(sys.argv) == 1:
    pm_user = [0]
    q2_user = [0]
    hist_user = [0]
    print('No arguments passed')

else:
    pm_user = sys.argv[1].split(',')  # central missing momentum setting
    q2_user = sys.argv[2].split(',')  # central Q2 setting
    hist_user = sys.argv[3]           # user histogram (must match filename hist prefix)
    
    pm_user = [int(s) for s in pm_user] # convert string to ints
    q2_user = [float(s) for s in q2_user] # convert string to ints

    
print('pm_user:', pm_user)
print('q2_user:', q2_user)

for i, ipm in enumerate(pm_user):
    
    for j, jq2 in enumerate(q2_user):

        fname = '%s_yield_d2pol_pm%s_Q2_%s.txt'%(hist_user, ipm, jq2)

        #-----------------------------
        if os.path.exists(fname):
            print('input file:',fname)
        else:
            print('The file does not exist, continue.')
            continue
        #-----------------------------
        
        df = pd.read_csv(fname, comment='#')

        #df.x0, df.ycont, np.sqrt(df.ycont)
        plt.hist(df.x0, bins=len(df.xb), weights=df.ycont, alpha=0.2, density=False, label=r'$Q^{2}=%s$, $P_{m}$=%s' % (jq2, ipm) )
        #plt.plot([1,2,3], [1,2,3], linestyle='None', mec='k', markersize=20, label=r'$Q^{2}=%s$, $P_{m}$=%s' % (jq2, ipm))
        plt.xlabel(r'x-label [units]')
        plt.ylabel(r'y-label [units]')
        plt.legend(frameon=False)
        
plt.show()
