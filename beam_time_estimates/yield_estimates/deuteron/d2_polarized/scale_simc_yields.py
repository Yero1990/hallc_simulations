'''
Brief:
This script reads Misak's calculation
of pol, unpol d(e,e'p) theory cross sections
and attempts to scale the SIMC unpol yields by
these cross sections to get the polarized yield
at SIMC
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

#fname_simc=''

# read .csv file
df_theory = pd.read_csv(fname_theory, comment='#')

#df_simc = pd.read_csv(fname_simc, comment='#')

# set parameters
Pzz = 0.3  # tensor polarization (~30%)

pm_bin_th       = df_theory['pm_bin'] #[GeV]
sig_upol_paris  = df_theory['sigma_unpol_PARIS']  # [nb/GeV/str]
sig_pol_paris   = df_theory['sigma_ten_PARIS']  # [nb/GeV/str]


