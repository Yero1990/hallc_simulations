import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

df = pd.read_csv('d2_pm_bins_Q2_2p9_rad_thrq35.csv', comment='#')

pm_bin       = np.array(df['x0']) #[GeV]
y_cnt       = np.array(df['ycont']) #[GeV]
y_cnt_2weeks = y_cnt *2
y_ref = y_cnt * 0.
rel_err =  np.sqrt(y_cnt) / y_cnt
rel_err_2weeks =  np.sqrt(y_cnt_2weeks) / y_cnt_2weeks


plt.errorbar(pm_bin, y_ref, rel_err*100., marker='o', color='g', label='1 week (168 hrs)')
#plt.errorbar(pm_bin, y_ref, rel_err_2weeks*100., marker='o', color='r', label='2 weeks')
plt.xlabel('Missing Momentumm, $P_{m}$ [GeV/c]', fontsize=16)
plt.ylabel('Relative Error [%]', fontsize=16)
plt.grid(True)
plt.show()
