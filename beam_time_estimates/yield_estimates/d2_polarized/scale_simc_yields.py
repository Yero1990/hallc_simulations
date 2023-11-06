'''
Brief:
This script reads Misak's calculation
of pol, unpol d(e,e'p) theory cross sections
and attempts to scale the SIMC unpol yields by
these cross sections to get the polarized yield
at SIMC

To do the scaling correctly, I think one should do the following:
1) make a dense grid from Misak's code, i.e., calculate the cross sections
at very fine binning, maybe 10 MeV bins, and write to output: basic kinematics,
like:  Eb, Ef, th_e, Pf, th_p, Pm, ... cross sections . .. .,
2) from SIMC event loop, for every event of a particular (Eb, Ef, th_e, th_p, etc.) that
define the complete state of the particle, pick the corresponding cross section from the theory
grid

Now, there might be tricky parts, like: what if the simulation includes radiative effects + energy loss
which change the kinematics? In that case, the SIMC kinematics, affected by radiative effects might not
correspond to the true kinematics as calculated by the theory grid, and one might end up picking the
cross sections at the incrorect kinematics. But maybe these effects can be neglected, if one is doing
a rough estimation of rates ? 


'''

import numpy as np
import pandas as pd
from itertools import takewhile
import matplotlib.pyplot as plt
import sys
import uncertainties
from uncertainties import ufloat
from uncertainties import unumpy

  
# set filename
fname_theory='misak_calculations_by_Nathaly/Q2_3.50_x_1.30.dat'

fname_simc='d2_pm_bins_Q2_3p5_rad.csv'

# read .csv file
df_theory = pd.read_csv(fname_theory, comment='#')

df_simc = pd.read_csv(fname_simc, comment='#')

# set parameters
Pzz = 0.3  # tensor polarization (~30%)

pm_bin_th       = df_theory['pm_bin'] #[GeV]
sig_upol_paris  = df_theory['sigma_unpol_PARIS']  # [nb/GeV/str]
sig_pol_paris   = df_theory['sigma_ten_PARIS']  # [nb/GeV/str]

pm_bin_simc = df_simc['x0']


