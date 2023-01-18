import matplotlib.pyplot as plt
import pandas as pd
import numpy as np


#Use Nice Fonts 
from matplotlib import rc
rc('text', usetex=True)
plt.rcParams["font.family"] = "Times New Roman"

# histogram utilitie for slicing 2d histogram in X 
def sliceX():

           
    # get generic file name
    file_path = ''
    
    # read .csv file
    df = pd.read_csv(summary_file_path, comment='#')

 
    
    #read data from file
    xb = f['xb']    # bin number
    x0 = f['x0']    # central x-bin value
    xlo = f['xlow'] # bin low edge value
    xup = f['xup']  # bin up edge value
    nbins = max(f['xb']) # total no. of bins
    xmin = min(xlo)
    xmax = max(xup)
    
    # the scaler factor is to scale the bin contents by beam time (hrs), if needed)
    bin_cnt = f['bin_cnt'] * scl_factor[idx]            # bin content
    bin_cnt_err = f['bin_cnt_err'] * scl_factor[idx]  # absolute error in bin content
