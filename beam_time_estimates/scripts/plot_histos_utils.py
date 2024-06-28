import sys
import os
import numpy as np
import pandas as pd
import numpy.ma as ma
import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from scipy.interpolate import interp1d

plt.rcParams["font.family"] = "Times New Roman"

# Brief: generic histos plotting utility functions script that loops over numerical histo files
# in a specified directory, and  makes and saves the plots
def get_label(label, ifile):
    with open(ifile, "r") as fp:
        for line in fp:
            if label in line:
                label = (line.split(':')[1]).strip()
                return label



def plot_rates():

    fig, axs = plt.subplots(1,3, figsize=(15,5))
    fig.suptitle('polarized d(e,e\'p) rates')


    #f1 = 'yield_estimates/d2_pol/smallFSI/output_rates_d2pol.txt'
    f2 = 'yield_estimates/d2_pol/smallFSI/optimized/output_rates_d2pol_optim.txt'  # optimized daata (centered hms delta)
    
    #df1 = pd.read_csv(f1, comment='#')
    df2 = pd.read_csv(f2, comment='#')
    
    Q2_set = list(set(df2.Q2_set))
    time = list(df2.time) # hrs
    charge = list(set(df2.charge)) # mC
    print(time)
    
    clr = ['b', 'r', 'g', 'm']
    
    for idx, iq2 in enumerate(Q2_set):
        print('idx:',idx)
        #axs[0].plot(df1.pm_set[df1.Q2_set==Q2_set[idx]], df1.deep_rates[df1.Q2_set==Q2_set[idx]]*3600., linestyle='None', marker='o', mfc=clr[idx], ms=6,  mec='k',  label='$Q^{2}$=%.1f'%(Q2_set[idx]))
        axs[0].plot(df2.pm_set[df2.Q2_set==Q2_set[idx]], df2.deep_rates[df2.Q2_set==Q2_set[idx]]*3600,  linestyle='None', marker='*', mfc=clr[idx], ms=12, mec='k',  label='$Q^{2}$=%.1f'%(Q2_set[idx]))
        axs[0].set_xlabel('central $p_{m}$ [MeV/c]')
        axs[0].set_ylabel('d(e,e\'p) rates [counts / hr]')
        axs[0].set_title('d(e,e\'p) yield rates')
        axs[0].legend()
        
        #axs[1].plot(df1.pm_set[df1.Q2_set==Q2_set[idx]], df1.daq_rates[df1.Q2_set==Q2_set[idx]]*3600, linestyle='None', marker='o', mfc=clr[idx], ms=6,  mec='k',  label='$Q^{2}$=%.1f'%(Q2_set[idx]))
        axs[1].plot(df2.pm_set[df2.Q2_set==Q2_set[idx]], df2.daq_rates[df2.Q2_set==Q2_set[idx]]*3600, linestyle='None', marker='*', mfc=clr[idx], ms=12, mec='k',  label='$Q^{2}$=%.1f'%(Q2_set[idx]))
        axs[1].set_xlabel('central $p_{m}$ [MeV/c]')
        axs[1].set_ylabel('DAQ rates [counts / hr]')
        axs[1].set_title('DAQ rates')
        axs[1].legend()
        
        #axs[2].plot(df1.pm_set[df1.Q2_set==Q2_set[idx]], df1.pm_counts[df1.Q2_set==Q2_set[idx]], linestyle='None', marker='o', mfc=clr[idx], ms=6, mec='k',  label='$Q^{2}$=%.1f'%(Q2_set[idx]))
        axs[2].plot(df2.pm_set[df2.Q2_set==Q2_set[idx]], df2.pm_counts[df2.Q2_set==Q2_set[idx]], linestyle='None', marker='*', mfc=clr[idx], ms=12, mec='k',  label='$Q^{2}$=%.1f, time=%.1f (hrs)'%(Q2_set[idx], time[idx] ))

        if(Q2_set[idx]==3.5):
            axs[2].plot(df2.pm_set[df2.Q2_set==Q2_set[idx]], 3.*df2.pm_counts[df2.Q2_set==Q2_set[idx]], linestyle='None', marker='*', mfc='g', ms=12, mec='k',  label='$Q^{2}$=%.1f, time=%.1f (hrs)'%(Q2_set[idx], 3.*time[idx] ))

        axs[2].set_xlabel('central $p_{m}$ [MeV/c]')
        axs[2].set_ylabel('d(e,e\'p) yield (counts)')
        axs[2].set_title('integrated yield')
        axs[2].legend()
        
    plt.show()
 

# original method
def combine_sets(kin_set=[], tgt='', hist_name='', model='', plot_flag=''):
    
    #Brief: function to overlay and combine multiple sets of the form kin_set[pm, Q2, scale] central values for a given target 'tgt'
    #e.g. combine_sets([ [200,3.7, 3], [300,4.0, 1]] ) which can loop over each set and overlay them
    
    rel_err_thrs = 0.5 # mask >30 % relative error

    clr = ['b', 'orange']  # colors depending on the number of elements in kin_set

    if plot_flag=='proj' or plot_flag=='tot_proj':
        # set figure subplots for the 1d projections
        fig, ax = plt.subplots(4, 3, tight_layout=True)
        
        fig.text(0.5, 0.002, 'missing momentum, p$_{m}$ [GeV/c]', ha='center', fontsize=12)
        #fig.text(0.005, 0.5, 'counts', va='center', rotation='vertical', fontsize=12)
        subplot_title =  r"Combined $(P_{m}, Q^{2})$ Kin. Settings"
        plt.suptitle(subplot_title, fontsize=15);
        fig.set_size_inches(12,10, forward=True)

    if plot_flag=='proj_err' or plot_flag=='tot_proj_err':
        # set figure subplots for the 1d projections (relative error)
        fig2, ax2 = plt.subplots(4, 3, tight_layout=True)
        
        fig2.text(0.5, 0.002, 'missing momentum, p$_{m}$ [GeV/c]', ha='center', fontsize=12)
        fig2.text(0.005, 0.5, 'relative stats. error', va='center', rotation='vertical', fontsize=12)
        subplot_title2 =  r"Combined $(P_{m}, Q^{2})$ Relative Stats. Error"
        plt.suptitle(subplot_title2, fontsize=15);
        fig2.set_size_inches(12,10, forward=True)

        
        
    offset=0 # offset for applying projections overlay 

    # sum elemtn by element of mutiple arrays
    # total = [sum(i) for i in zip(*pm_cnts)]
    
    print('len(kin_set):', len(kin_set))

    total_counts_per_xbin = 0
    
    # loop over each [pm, Q2, scale] set
    for i, ikin in enumerate(kin_set):
        
        pm = ikin[0]
        Q2 = ikin[1]
        scale = ikin[2]
    
        
        # increment offset for every new kin_set setting
        if i>0:
            offset = offset + 0.005  
        
        hist_basename    = 'H_%s_yield_d2pol_pm%d_Q2_%.1f.txt'%(hist_name, pm, Q2)   # generic histogram basename
        hist_file        = 'yield_estimates/d2_pol/smallFSI/optimized/histogram_data/tightEmiss_Cut/%s_pm%d_Q2_%.1f_%s/%s'%(tgt, pm, Q2, model, hist_basename)

        
        if not os.path.isfile(hist_file): continue
        
        df = pd.read_csv(hist_file, comment='#')

        xlabel = get_label('xlabel',     hist_file)
        ylabel = get_label('ylabel',     hist_file)
        title  = get_label('title',      hist_file)
        ybinw  = float(get_label('ybin_width', hist_file))
        xbinw  = float(get_label('xbin_width', hist_file))
        nxbins = len(df.xb[df.yb==df.yb[0]])
        nybins = len(df.yb[df.xb==df.xb[0]])
        ybc = (df.y0[df.x0==df.x0[0]]).to_numpy() # y-bin central value
        xbc = (df.x0[df.y0==df.y0[0]]).to_numpy() # x-bin central value

        df.zcont = df.zcont * scale    # scale by beam time (in multiples of the standard beamtime)
        counts = np.sum(df.zcont)

        
        if i==0:
            print('initializing counter')
            print('xbc, ybc', len(xbc), len(ybc))
            total_counts_per_xbin = np.zeros((len(xbc), len(ybc)))  # create a 2D zero array of length (xbins, ybins) to initialize counter

        
        jdx=0 #counter for subplots which have non-zero counts
        
        #loop over x-bins (for y-projections)
        for idx, xbin in enumerate( xbc ):
                           
            count_per_xbin     = df.zcont[df.x0==xbin]
            count_per_xbin_err = np.sqrt( count_per_xbin )
            
            count_per_xbin     = ma.masked_where(count_per_xbin==0, count_per_xbin)
            count_per_xbin_err = ma.masked_where(count_per_xbin==0, count_per_xbin_err)
            
            cnts = np.sum(count_per_xbin ) 


            # count sum  bin content for corresponding bins of a given idx
            total_counts_per_xbin[idx] = total_counts_per_xbin[idx] + np.array([count_per_xbin])

                
            # variables for plotting relative error
            count_per_xbin_rel_err = count_per_xbin_err / count_per_xbin
            y_const = np.zeros(len(count_per_xbin_rel_err))


            #------------------------------------------------------------------------
            # CALCULATE / PLOT COMBINE KINEMATICS (ONLY IF HAS REACHED END OF LOOP
            #------------------------------------------------------------------------
            # check if this is the last in both the kin_set and xbin loops
            if (i == len(kin_set) - 1) and (idx == (len(xbc) - 1)):

                if plot_flag=='tot_proj' or plot_flag=='tot_proj_err':
                    
                    jdx = 0

                    #loop over x-bins (for y-projections) -- again
                    for idx, xbin in enumerate( xbc ):
                    
                        # calculate the absolute and relative total errors
                        total_counts      = np.array(total_counts_per_xbin[idx])
                        total_counts_err  = np.sqrt(total_counts)
                        total_counts_rel_err  = 1. / np.sqrt(total_counts)
                        y_const               = np.zeros(len(total_counts_rel_err))
                        
                        # sum all counts of a particular xbin
                        cnts = np.sum(total_counts)
                        
                        if(np.ma.is_masked(cnts)): continue
                        
                        if(not np.ma.is_masked(cnts)):
                            
                            if plot_flag=='tot_proj' :
                                
                                #plot combined settings
                                ax = plt.subplot(5, 4, jdx+1)
                                ax.hist(ybc, bins=len(ybc), weights=total_counts, range=[min(df.ylow), max(df.yup)], alpha=0.5, ec='k', color="gray", density=False, label=r'%d counts'%(cnts))
                                plt.title(r'$\theta_{nq}$ = %d $\pm$ %d deg'%(xbin, xbinw/2.))
                                plt.legend(frameon=False, loc='upper right')
                                jdx = jdx+1

                            if plot_flag=='tot_proj_err' :

                                #plot combined settings relative error
                                ax2 = plt.subplot(5, 4, jdx+1)
                                total_counts_rel_err_m = ma.masked_where(total_counts_rel_err>rel_err_thrs, total_counts_rel_err)
                                y_const_m = ma.masked_where(total_counts_rel_err>rel_err_thrs, y_const)
                                ybc_m = ma.masked_where(total_counts_rel_err>rel_err_thrs, ybc)

                                
                                ax2.errorbar(ybc_m, y_const_m, total_counts_rel_err_m, marker='o', mfc='gray', mec='gray', markersize=4, ecolor='gray', linestyle='None', label=r'(%d MeV, %.1f GeV$^{2}$)'%(pm, Q2)) #//, label=r'%d counts'%(cnts))
                                plt.axhline(y = 0.20, color = 'r', linestyle = '--')
                                plt.axhline(y = -0.20, color = 'r', linestyle = '--')
                                plt.ylim(-0.6,0.6)
                                plt.xlim(0., 0.7)
                                plt.title(r'$\theta_{nq}$ = %d $\pm$ %d deg'%(xbin, xbinw/2.))
                                jdx = jdx+1
                                
                    plt.show()

                    
            if(np.ma.is_masked(cnts)): continue
            
            if(not np.ma.is_masked(cnts)):
            

                if plot_flag=='proj':
                    
                    #plot overlay of settings
                    ax = plt.subplot(4, 3, jdx+1)
                    ax.hist(ybc, bins=len(ybc), weights=count_per_xbin, range=[min(df.ylow), max(df.yup)], color=clr[i], alpha=0.5, ec='k', density=False, label=r'%d counts'%(cnts)+"\n"+r"setting:(%d MeV, %.1f GeV$^{2}$)"%(pm, Q2))
                    plt.title(r'$\theta_{nq}$ = %d $\pm$ %d deg'%(xbin, xbinw/2.))
                    plt.legend(frameon=False, loc='upper right')

                    jdx = jdx+1


                if plot_flag=='proj_err':
                    
                    #plot overlay of settings relative error
                    ax2 = plt.subplot(4, 3, jdx+1)
                    count_per_xbin_rel_err_m = ma.masked_where(count_per_xbin_rel_err>rel_err_thrs, count_per_xbin_rel_err)
                    y_const_m = ma.masked_where(count_per_xbin_rel_err>rel_err_thrs, y_const)
                    ybc_m = ma.masked_where(count_per_xbin_rel_err>rel_err_thrs, ybc)

                    ybc_m = ybc_m + offset
                    ax2.errorbar(ybc_m, y_const_m, count_per_xbin_rel_err_m, marker='o', markersize=4, linestyle='None', color=clr[i], elinewidth=2, label=r'(%d MeV, %.1f GeV$^{2}$)'%(pm, Q2)) #//, label=r'%d counts'%(cnts))
                    plt.axhline(y = 0.20, color = 'r', linestyle = '--')
                    plt.axhline(y = -0.20, color = 'r', linestyle = '--')
                    plt.ylim(-0.6,0.6)
                    plt.xlim(0.0, 0.7)
                    plt.title(r'$\theta_{nq}$ = %d $\pm$ %d deg'%(xbin, xbinw/2.))

                    #print('------> relative errors for (pm, Q2):', pm, Q2)
                    #print('rel_err:',  count_per_xbin_rel_err_m,)

                    jdx = jdx+1



              
                    
        plt.legend(frameon=False, loc='upper right')
                                
    plt.show()



# this alternative version of the method takes an additional input for combining multiple targets
def combine_sets_alt(kin_set=[], Q2 = 0, hist_name='', model='', plot_flag=''):
    
    #Brief: function to overlay and combine multiple sets of the form kin_set[tgt, pm, Q2, scale] central values
    #e.g. combine_sets([ [300,3.5, 2, 'n14'], [300, 3.5, 2, 'd2'], [ 300, 3.5, 2, 'he4'] ] ) which can loop over each set and overlay them
    
    rel_err_thrs = 0.5 # mask >30 % relative error

    #clr = ['orange', 'b', 'magenta']  # colors depending on the number of elements in kin_set
    clr = ['orange', 'b', 'magenta']

    if plot_flag=='proj' or plot_flag=='tot_proj':
        # set figure subplots for the 1d projections
        fig, ax = plt.subplots(4, 3, tight_layout=True)
        
        fig.text(0.5, 0.002, 'missing momentum, p$_{m}$ [GeV/c]', ha='center', fontsize=12)
        #fig.text(0.005, 0.5, 'counts', va='center', rotation='vertical', fontsize=12)
        subplot_title =  r"Combined $(P_{m}, Q^{2})$ Kin. Settings"
        plt.suptitle(subplot_title, fontsize=15);
        fig.set_size_inches(12,10, forward=True)

    if plot_flag=='proj_err' or plot_flag=='tot_proj_err':
        # set figure subplots for the 1d projections (relative error)
        fig2, ax2 = plt.subplots(4, 3, tight_layout=True)
        
        fig2.text(0.5, 0.002, 'missing momentum, p$_{m}$ [GeV/c]', ha='center', fontsize=12)
        fig2.text(0.005, 0.5, 'relative stats. error', va='center', rotation='vertical', fontsize=12)
        subplot_title2 =  r"Combined $(P_{m}, Q^{2})$ Relative Stats. Error"
        plt.suptitle(subplot_title2, fontsize=15);
        fig2.set_size_inches(12,10, forward=True)

        
        
    offset=0 # offset for applying projections overlay 

    # sum elemtn by element of mutiple arrays
    # total = [sum(i) for i in zip(*pm_cnts)]
    
    print('len(kin_set):', len(kin_set))

    total_counts_per_xbin = 0
    
    # loop over each [pm, Q2] set
    for i, ikin in enumerate(kin_set):

        
        pm    = ikin[0]
        scale = ikin[1]
        tgt   = ikin[2]
        
        # increment offset for every new kin_set setting
        if i>0:
            offset = offset + 0.005  
        
        hist_basename    = 'H_%s_yield_d2pol_pm%d_Q2_%.1f.txt'%(hist_name, pm, Q2)   # generic histogram basename
        hist_file        = 'yield_estimates/d2_pol/smallFSI/optimized/histogram_data/tightEmiss_Cut/%s_pm%d_Q2_%.1f_%s/%s'%(tgt, pm, Q2, model, hist_basename)
        #hist_file        = 'yield_estimates/d2_pol/smallFSI/optimized/histogram_data/tightEmiss_Cut/%s_pm%d_Q2_%.1f_%s/%s'%(tgt, pm, Q2, model, hist_basename)

        if not os.path.isfile(hist_file): continue
        
        df = pd.read_csv(hist_file, comment='#')

        xlabel = get_label('xlabel',     hist_file)
        ylabel = get_label('ylabel',     hist_file)
        title  = get_label('title',      hist_file)
        ybinw  = float(get_label('ybin_width', hist_file))
        xbinw  = float(get_label('xbin_width', hist_file))
        nxbins = len(df.xb[df.yb==df.yb[0]])
        nybins = len(df.yb[df.xb==df.xb[0]])
        ybc = (df.y0[df.x0==df.x0[0]]).to_numpy() # y-bin central value
        xbc = (df.x0[df.y0==df.y0[0]]).to_numpy() # x-bin central value

        df.zcont = df.zcont * scale    # scale by beam time (in multiples of the standard beamtime)
        counts = np.sum(df.zcont)

        
        if i==0:
            print('initializing counter')
            print('xbc, ybc', len(xbc), len(ybc))
            total_counts_per_xbin = np.zeros((len(xbc), len(ybc)))  # create a 2D zero array of length (xbins, ybins) to initialize counter

        
        jdx=0 #counter for subplots which have non-zero counts
        
        #loop over x-bins (for y-projections)
        for idx, xbin in enumerate( xbc ):
                           
            count_per_xbin     = df.zcont[df.x0==xbin]
            count_per_xbin_err = np.sqrt( count_per_xbin )
            
            count_per_xbin     = ma.masked_where(count_per_xbin==0, count_per_xbin)
            count_per_xbin_err = ma.masked_where(count_per_xbin==0, count_per_xbin_err)
            
            cnts = np.sum(count_per_xbin ) 


            # count sum  bin content for corresponding bins of a given idx
            total_counts_per_xbin[idx] = total_counts_per_xbin[idx] + np.array([count_per_xbin])

                
            # variables for plotting relative error
            count_per_xbin_rel_err = count_per_xbin_err / count_per_xbin
            y_const = np.zeros(len(count_per_xbin_rel_err))


            #------------------------------------------------------------------------
            # CALCULATE / PLOT COMBINE KINEMATICS (ONLY IF HAS REACHED END OF LOOP
            #------------------------------------------------------------------------
            # check if this is the last in both the kin_set and xbin loops
            if (i == len(kin_set) - 1) and (idx == (len(xbc) - 1)):

                if plot_flag=='tot_proj' or plot_flag=='tot_proj_err':
                    
                    jdx = 0

                    #loop over x-bins (for y-projections) -- again
                    for idx, xbin in enumerate( xbc ):
                    
                        # calculate the absolute and relative total errors
                        total_counts      = np.array(total_counts_per_xbin[idx])
                        total_counts_err  = np.sqrt(total_counts)
                        total_counts_rel_err  = 1. / np.sqrt(total_counts)
                        y_const               = np.zeros(len(total_counts_rel_err))
                        
                        # sum all counts of a particular xbin
                        cnts = np.sum(total_counts)
                        
                        if(np.ma.is_masked(cnts)): continue
                        
                        if(not np.ma.is_masked(cnts)):
                            
                            if plot_flag=='tot_proj' :
                                
                                #plot combined settings
                                ax = plt.subplot(5, 4, jdx+1)
                                ax.hist(ybc, bins=len(ybc), weights=total_counts, range=[min(df.ylow), max(df.yup)], alpha=0.5, ec='k', color="gray", density=False, label=r'%d counts'%(cnts))
                                plt.title(r'$\theta_{nq}$ = %d $\pm$ %d deg'%(xbin, xbinw/2.))
                                plt.legend(frameon=False, loc='upper right')
                                jdx = jdx+1

                            if plot_flag=='tot_proj_err' :

                                #plot combined settings relative error
                                ax2 = plt.subplot(5, 4, jdx+1)
                                total_counts_rel_err_m = ma.masked_where(total_counts_rel_err>rel_err_thrs, total_counts_rel_err)
                                y_const_m = ma.masked_where(total_counts_rel_err>rel_err_thrs, y_const)
                                ybc_m = ma.masked_where(total_counts_rel_err>rel_err_thrs, ybc)

                                
                                ax2.errorbar(ybc_m, y_const_m, total_counts_rel_err_m, marker='o', mfc='gray', mec='gray', markersize=4, ecolor='gray', linestyle='None', label=r'(%d MeV, %.1f GeV$^{2}$)'%(pm, Q2)) #//, label=r'%d counts'%(cnts))
                                plt.axhline(y = 0.20, color = 'r', linestyle = '--')
                                plt.axhline(y = -0.20, color = 'r', linestyle = '--')
                                plt.ylim(-0.6,0.6)
                                plt.xlim(0., 0.7)
                                plt.title(r'$\theta_{nq}$ = %d $\pm$ %d deg'%(xbin, xbinw/2.))
                                jdx = jdx+1
                                
                    plt.show()

                    
            if(np.ma.is_masked(cnts)): continue
            
            if(not np.ma.is_masked(cnts)):
            

                if plot_flag=='proj':
                    
                    #plot overlay of settings
                    ax = plt.subplot(4, 3, jdx+1)
                    ax.hist(ybc, bins=len(ybc), weights=count_per_xbin, range=[min(df.ylow), max(df.yup)], color=clr[i], alpha=0.5, ec='k', density=False, label=r'%d counts'%(cnts)+"\n"+r"setting:(%s, %d MeV, %.1f GeV$^{2}$)"%(tgt, pm, Q2))
                    plt.title(r'$\theta_{nq}$ = %d $\pm$ %d deg'%(xbin, xbinw/2.))
                    plt.legend(frameon=False, loc='upper right')

                    jdx = jdx+1


                if plot_flag=='proj_err':
                    
                    #plot overlay of settings relative error
                    ax2 = plt.subplot(4, 3, jdx+1)
                    count_per_xbin_rel_err_m = ma.masked_where(count_per_xbin_rel_err>rel_err_thrs, count_per_xbin_rel_err)
                    y_const_m = ma.masked_where(count_per_xbin_rel_err>rel_err_thrs, y_const)
                    ybc_m = ma.masked_where(count_per_xbin_rel_err>rel_err_thrs, ybc)

                    ybc_m = ybc_m + offset
                    ax2.errorbar(ybc_m, y_const_m, count_per_xbin_rel_err_m, marker='o', markersize=4, linestyle='None', color=clr[i], elinewidth=2, label=r'(%s, %d MeV, %.1f GeV$^{2}$)'%(tgt, pm, Q2)) #//, label=r'%d counts'%(cnts))
                    plt.axhline(y = 0.20, color = 'r', linestyle = '--')
                    plt.axhline(y = -0.20, color = 'r', linestyle = '--')
                    plt.ylim(-0.6,0.6)
                    plt.xlim(0.0, 0.7)
                    plt.title(r'$\theta_{nq}$ = %d $\pm$ %d deg'%(xbin, xbinw/2.))

                    #print('------> relative errors for (pm, Q2):', pm, Q2)
                    #print('rel_err:',  count_per_xbin_rel_err_m,)

                    jdx = jdx+1



              
                    
        plt.legend(frameon=False, loc='upper right')
                                
    plt.show()



    
    
def overlay_d2fsi(pm_set, thrq_set, hist_name, model, scale=[1,1,1]):
    '''
    Brief: generic function to overlay 1d histograms from multiple kin. files for deut fsi
    pm_set and thrq_set are lists of values representing the different kinematic settings
    scale: scale factor by beam time
    pm_set: central missing momentum setting, e.g. pm_set = [800]
    thrq_set : recoin neutron angle setting, e.g. thrq_set = [49, 60, 72]
    hist_name: histogram base name
    model: 'pwia', 'fsi'
    scale: default set to 1; scale counts by a multiple of beam_time (need to know what the beam time was during simulation)
    '''

    fig, axs = plt.subplots(1, figsize=(6,5))

    idx = 0  # counter for scale index 
    # loop over central missing momentum setting
    for ipm in pm_set:

    
        # loop over central q2 setting
        for ithrq in thrq_set:

            
            # set histogram file path
            #histos_file_path = 'path/to/histogram_data/pm%d_q2%d_%s/histo_name_pm_set_q2_set.txt'%(pm_set, q2_set, model, hist_name)

            hist_file = 'yield_estimates/d2_fsi/histogram_data/pm%d_thrq%d_%s/H_%s_yield_d2fsi_pm%d_thrq%d.txt'%(ipm, ithrq, model, hist_name, ipm, ithrq)

            df = pd.read_csv(hist_file, comment='#')

            # make plot of 1D histogram
            print("VALID 1D HIST")
            
            xlabel = get_label('xlabel', hist_file)
            ylabel = get_label('ylabel', hist_file) 
            title  = get_label('title', hist_file) 

            x = df.x0
            df.ycont = df.ycont * scale[idx]
            N = df.ycont
            Nerr = np.sqrt(N)
            
            x =   ma.masked_where(N==0, x)
            N =   ma.masked_where(N==0, N)
            Nerr = ma.masked_where(N==0, Nerr)

            
            plt.hist(df.x0, bins=len(df.x0), weights=df.ycont, alpha=0.5, ec='k', density=False, label=r"$\theta_{nq}=%d$ deg"%(ithrq)+"\n"+"$P_{m}$=%d MeV"%(ipm))


            #plt.errorbar(x, N, Nerr, linestyle='None', marker='o', mec='k', label=r"$\theta_{nq}=%d$ deg"%(ithrq)+"\n"+"$P_{m}$=%d MeV"%(ipm))

            #plt.legend()
            #plt.xlabel(xlabel, fontsize=18)
            #plt.ylabel(ylabel, fontsize=18)
            plt.xticks(fontsize = 22)
            plt.yticks(fontsize = 22)
            plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0), useMathText=True)
            axs.yaxis.offsetText.set_fontsize(18) # increase fontsie of y-axis sci not.
            #plt.title(title, fontsize=15)


            # limit the number of ticks
            max_xticks = 4
            max_yticks = 4
            xloc = plt.MaxNLocator(max_xticks)
            yloc = plt.MaxNLocator(max_yticks)
            axs.xaxis.set_major_locator(xloc)
            axs.yaxis.set_major_locator(yloc)

            
            # for drawing vertical lines of cuts

            # missing energy
            #plt.axvline(x = -0.02, color = 'r', linestyle = '--', linewidth=2)
            #plt.axvline(x = 0.04, color = 'r', linestyle = '--', linewidth=2)

            # hms delta
            #plt.axvline(x = -10, color = 'r', linestyle = '--', linewidth=2)
            #plt.axvline(x = 10, color = 'r', linestyle = '--', linewidth=2)

            # shms delta
            #plt.axvline(x = -10, color = 'r', linestyle = '--', linewidth=2)
            #plt.axvline(x = 22, color = 'r', linestyle = '--', linewidth=2)

            # Q2
            #plt.axvline(x = 4, color = 'r', linestyle = '--', linewidth=2)
            #plt.axvline(x = 5, color = 'r', linestyle = '--', linewidth=2)

            

            
            idx = idx + 1

            

    plt.show()


def overlay_d2pol(tgt_set, pm_set, Q2_set, hist_name, model, field, scale=1):
    '''
    Brief: generic function to overlay 1d histograms from multiple kin. files for deut pol. proposal

    -----------
    user input
    -----------
    pm_set: central missing momentum, e.g.  pm_set=[200], etc. (ONLY single-value can be used)
            since the code was structured so that for a single pm_set, loop over all Q2_set
            (and overlay all Q2 settings) 
    
    Q2_set: central Q2 setting, e.g. Q2_set=[3.5, 4.0, 4.5], (a single valued or multi-valued list can be provided)

    hist_name: histogram base name (if unfamiliar, will need to check the name in the histogram directory
               the generic histogram name has the form: H_[hist_name]_yield_d2pol_pm[pm_set]_Q2_[Q2_set].txt
               where the variable in brackets is replaced by user input

    model: can be "fsi" or "pwia" depending on the simulation done
    
    scale: time scale factor to scale yield by a multiple of the simulated beam-time, for example, if 168 hrs
    is the standard simulation beam time, and user sets scale = 3, then yield will be scaled to
    3 * yield, which is esentially triple the beam time for the yield (3 * 168 hrs)
    '''


    rel_err_thrs = 1000   #  relative stat. error threshold for masking

    #fig, axs = plt.subplots(2, sharex=True, figsize=(5,10))
    #fig, axs = plt.subplots(2, figsize=(6,7))
    fig, axs = plt.subplots(1, figsize=(6,5))

    #clr = ['orange', 'b', 'magenta']
    #clr = ['orange', 'b', 'magenta']
    clr = ['orange', 'deepskyblue', 'magenta']
 
    i=-1 # index for color counting
    
    offset=0 # for applying to overlayed data for easy visual
    # loop over the different target nuclei
    for itgt in tgt_set:
        i = i+1
        # loop over central missing momentum setting
        for ipm in pm_set:

            #offset=0 # for applying to overlayed data for easy visual
            # loop over central q2 setting
            for iq2 in Q2_set:
                
                # set histogram file path
                
                #hist_file = 'yield_estimates/d2_pol/smallFSI/optimized/histogram_data/%s_pm%d_Q2_%.1f_%s_%s/H_%s_yield_d2pol_pm%d_Q2_%.1f.txt'%(itgt, ipm, iq2, model, field, hist_name, ipm, iq2)

                # phi = 0 config
                hist_file = 'yield_estimates/d2_pol/smallFSI/optimized/histogram_data/phi0/%s_pm%d_Q2_%.1f_%s_%s/H_%s_yield_d2pol_pm%d_Q2_%.1f.txt'%(itgt, ipm, iq2, model, field, hist_name, ipm, iq2)




                
                print('tgt_set:', itgt, 'Pm:', ipm, 'Q2:', iq2)
                
                print('file %s exists?:'%(hist_file), os.path.isfile(hist_file))
                if not os.path.isfile(hist_file): continue
                
                df = pd.read_csv(hist_file, comment='#')
                
                # make plot of 1D histogram
                print("VALID 1D HIST")
                
                xlabel = get_label('xlabel', hist_file)
                ylabel = get_label('ylabel', hist_file) 
                title  = get_label('title', hist_file) 
                xbins  = get_label('xbins', hist_file) 
                xbinw  = float(get_label('xbin_width', hist_file))

                # increment offset for every new Q2 setting
                #if iq2 > Q2_set[0]:
                #    offset = offset + 0.1*xbinw

                print('tgt:', itgt, 'tgt_set[0]:',tgt_set[0], 'offset:', offset)
                # increment offset for every new tgt setting
                if itgt != tgt_set[0]:
                    offset = offset + 0.1*xbinw

                xbc = df.x0
                df.ycont = df.ycont * scale
                N = df.ycont
                Nerr = np.sqrt(N)
                
                xbc =   ma.masked_where(N==0, xbc)
                N =   ma.masked_where(N==0, N)
                Nerr = ma.masked_where(N==0, Nerr)
                Nerr_rel = Nerr/N
                y_const = np.zeros(len(N))
                
                counts = N.sum()
                # masked data to require rel. error > error threshold
                xbc        = ma.masked_where(Nerr_rel>rel_err_thrs, xbc)
                y_const  = ma.masked_where(Nerr_rel>rel_err_thrs, y_const)
                Nerr_rel = ma.masked_where(Nerr_rel>rel_err_thrs, Nerr_rel)
                N        = ma.masked_where(N>rel_err_thrs, N)
                Nerr     = ma.masked_where(N>rel_err_thrs, Nerr)
                
                
                # --- plot histogram ---
                print('min(df.xlow), max(df.xup) --> ', min(df.xlow), max(df.xup))
                print('xbin_center:', xbc)
                print('bins:', len(xbc))
                #axs.set_title(title)
                axs.hist(x=xbc, bins=len(xbc), range=[min(df.xlow), max(df.xup)], weights=N, alpha=0.15, ec='k',  color=clr[i], density=False, label="$P_{m}$=%d MeV"%(ipm)+"\n"+ "$Q^{2}$=%.1f GeV$^{2}$ (%d)"%(iq2, counts)+"\n"+"%s"%(itgt))
                #axs.set_ylabel(ylabel)
                #axs.set_xlabel(xlabel)

                

                #axs.legend(fontsize=12)
                axs.set_yscale('log')
                axs.xaxis.set_tick_params(labelbottom=True)
                #plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0), useMathText=True)  # comment out if using log
                axs.yaxis.offsetText.set_fontsize(18)
                plt.xticks(fontsize = 22)
                plt.yticks(fontsize = 22)

                #plt.ylim(1e-1, 5e5)

                # To specify the number of ticks on both or any single axes
                nxbins = len(xbc)

                #plt.locator_params(axis='y')
                #plt.locator_params(axis='x', tight=True)
                max_xticks = 7
                xloc = plt.MaxNLocator(max_xticks)
                axs.xaxis.set_major_locator(xloc)

                #plt.axvline(x = -0.01, color = 'r', linestyle = '--', linewidth=2)
                #plt.axvline(x = 0.04, color = 'r', linestyle = '--', linewidth=2)

                #plt.axvline(x = -10, color = 'r', linestyle = '--', linewidth=2)
                #plt.axvline(x = 10, color = 'r', linestyle = '--', linewidth=2)


                '''
                # apply offset and plot relative error
                xbc_off = xbc + offset
                axs[1].errorbar(xbc_off, y_const, Nerr_rel, linestyle='None', marker='o', mec='k', mfc=clr[i], ecolor=clr[i], label=r"$P_{m}$=%d MeV"%(ipm)+"\n"+"$Q^{2}$=%.1f GeV$^{2}$"%(iq2)+"\n"+"%s"%(itgt))
                axs[1].set_ylabel('Relative Error')
                axs[1].set_xlabel(xlabel)
                
                #plt.setp(axs, xlim=(np.min(xbc)-0.5*np.min(xbc), np.max(xbc)+0.5*np.min(xbc)))
                plt.axhline(y = 0.20, color = 'r', linestyle = '--')
                plt.axhline(y = -0.20, color = 'r', linestyle = '--')
                '''
                
    plt.show()    
            
def make_1d_Xprojections(h2_hist_name, pm_user, thrq_user, model):
    
    '''
    Brief: generic function makes 1D projections along x-axis (slicing ybins) for selected 2D histos
    '''

    histos_file_path = 'yield_estimates/d2_fsi/histogram_data/pm%d_thrq%d_%s/'%(pm_user, thrq_user, model)

    for fname in os.listdir (histos_file_path):
        
        # check if histo is 2D
        if ("_vs_" in fname):

            # check if specific yield histo is present
            if(h2_hist_name in fname):
            #if(1):
                hist_file = histos_file_path + fname

                print('Opening file -----> ', hist_file)
        
                df = pd.read_csv(hist_file, comment='#')

                xlabel = get_label('xlabel',     hist_file)
                ylabel = get_label('ylabel',     hist_file)
                title  = get_label('title',      hist_file)
                ybinw  = float(get_label('ybin_width', hist_file))
                

                ybc = (df.y0[df.x0==df.x0[0]]).to_numpy() # y-bin central value
               
                # count  to decide how many subplots to make
                X = round( np.sqrt( np.count_nonzero(df.y0[df.x0==df.x0[0]] )) )
                Y = X+1
                print('X:', X, 'Y:', Y)
                # set figure subplots
                fig, ax = plt.subplots(X, Y, sharex='col', sharey='row')
                fig.text(0.5, 0.007, xlabel, ha='center', fontsize=12)
                fig.text(0.01, 0.5, 'Counts', va='center', rotation='vertical', fontsize=12)
                subplot_title = title+': 1d x-projection (%s), setting: (%d MeV, %d deg)'%(model, pm_user, thrq_user)
                plt.suptitle(subplot_title, fontsize=15);
                fig.set_size_inches(14,10, forward=True)

    

                jdx = 0
                #loop over y-bins (for x-projections)
                for idx, ybin in enumerate( ybc ):
                    
                    xbins              = df.x0[df.y0==ybin]
                    count_per_ybin     = df.zcont[df.y0==ybin]
                    count_per_ybin_err = np.sqrt( count_per_ybin )
                
                    count_per_ybin     = ma.masked_where(count_per_ybin==0, count_per_ybin)
                    count_per_ybin_err = ma.masked_where(count_per_ybin==0, count_per_ybin_err)

                    cnts = np.sum(count_per_ybin )
                    
                    #---------------
                    if(np.ma.is_masked(cnts)): continue
                    #ax = plt.subplot(X, Y, idx+1)

                    # can think of additional conditons for plotting, for example, could put a constraint on the relatie error as well, and reduced
                    # the number of useless subplots
                    if(not np.ma.is_masked(cnts)):
                        jdx = jdx + 1
                        ax = plt.subplot(X, Y, jdx+1)
                        ax.errorbar(xbins, count_per_ybin, count_per_ybin_err, marker='o', markersize=4, linestyle='None', label=r'%d counts'%(cnts))

                    plt.title('$p_{m}$ = %d $\pm$ %d MeV'%(ybin*1000, ybinw*1000/2.))
                    plt.xlim([xbins.min(), xbins.max()])
                    plt.legend(frameon=False, loc='upper right')

                    
                plt.tight_layout()
                plt.show()


def make_all_plots(pm_user, thrq_user, model):

    '''
    Brief: loops thru all stored numerical histogram files (generated by the analysis code)
    and makes appropiate 1D or 2D histos
    '''

    # pm_user   : central missing momentum setting (e.g. 500)
    # thrq_user :  central recoil angle (e.g. 28)
    # model     : Laget 'pwia' or 'fsi' 

    #histos_file_path = 'yield_estimates/d2_fsi/histogram_data/pm%d_thrq%d_%s/'%(pm_user, thrq_user, model)
    histos_file_path = 'yield_estimates/d2_pol/smallFSI/optimized/histogram_data/pm%d_thrq%d_%s/'%(pm_user, thrq_user, model)

    for fname in os.listdir (histos_file_path):
        
        hist_file = histos_file_path + fname

        print('Opening file -----> ', hist_file)
        
        df = pd.read_csv(hist_file, comment='#')
        
        # check if histo is 2D
        if ("_vs_" in fname):
            print("VALID 2D HIST")
            xbins = len(df.xb[df.yb==df.yb[0]])
            ybins = len(df.yb[df.xb==df.xb[0]])
            zcont = np.array(df.zcont)
            counts =  np.sum(df.zcont)
            
            #if counts==0: break
            xlabel = get_label('xlabel', hist_file)
            ylabel = get_label('ylabel', hist_file) 
            title  = get_label('title', hist_file) 

            hist2d = plt.hist2d(x=df.x0 ,y=df.y0, bins=(xbins, ybins), range=[ [min(df.xlow), max(df.xup)], [min(df.ylow), max(df.yup)]], weights=zcont, cmap = 'viridis', norm=mcolors.LogNorm(vmin=0.1, vmax=np.sqrt(counts) ))
            plt.xlabel(xlabel, fontsize=12)
            plt.ylabel(ylabel, fontsize=12)
            plt.title(title,   fontsize=14)
             
            plt.text(0.6*(df.x0[df.y0==df.y0[0]]).max(), 0.7*(df.y0[df.x0==df.x0[0]]).max(), r"$\theta_{nq}=%d$ deg"%(thrq_user)+"\n"+"$P_{m}$=%d MeV"%(pm_user)+"\n"+"(counts = %d)"%(counts), fontsize=12)
            plt.colorbar(extend='max')

            plt.show()

        elif ("_2Davg" in fname):
            
            print("VALID 2D AVG")
            xbins = len(df.xb[df.yb==df.yb[0]])
            ybins = len(df.yb[df.xb==df.xb[0]])
            zcont = np.array(df.zcont)

            xlabel = get_label('xlabel', hist_file)
            ylabel = get_label('ylabel', hist_file) 
            title  = get_label('title', hist_file) 

            if(zcont.min()>0):
                hist2d = plt.hist2d(df.x0 ,df.y0, weights=zcont, bins=(xbins, ybins), range=[ [min(df.xlow), max(df.xup)], [min(df.ylow), max(df.yup)]], cmap = 'viridis', norm=mcolors.LogNorm())
            else:
                hist2d = plt.hist2d(df.x0 ,df.y0, weights=zcont, bins=(xbins, ybins), range=[ [min(df.xlow), max(df.xup)], [min(df.ylow), max(df.yup)]], cmap = 'viridis')
                #plt.scatter(df.x0, df.y0, c=zcont, s = 1, cmap = 'viridis', vmin = zcont.min(), vmax = zcont.max())
            plt.text(0.6*(df.x0[df.y0==df.y0[0]]).max(), 0.7*(df.y0[df.x0==df.x0[0]]).max(), r"$\theta_{nq}=%d$ deg"%(thrq_user)+"\n"+"$P_{m}$=%d MeV"%(pm_user), fontsize=12)
            plt.colorbar(extend='max')
            plt.xlabel(xlabel, fontsize=12)
            plt.ylabel(ylabel, fontsize=12)
            plt.title(title,   fontsize=14)

            plt.show()


            
        elif ("_2Davg_" and "_vs_") not in fname:
            # make plot of 1D histogram
            print("VALID 1D HIST")
            counts =  np.sum(df.ycont)

            xlabel = get_label('xlabel', hist_file)
            ylabel = get_label('ylabel', hist_file) 
            title  = get_label('title', hist_file) 

      
            plt.hist(df.x0, bins=len(df.xb), weights=df.ycont, range=[min(df.xlow), max(df.xup)], alpha=0.2, density=False, label=r"$\theta_{nq}=%d$ deg"%(thrq_user)+"\n"+"$P_{m}$=%d MeV"%(pm_user)+ "\n"+"(counts = %d)"%(counts))
            plt.xlabel(xlabel, fontsize=12)
            plt.ylabel(ylabel, fontsize=12)
            plt.title(title,   fontsize=14)
            plt.legend(frameon=False, fontsize=12)
            plt.show()


            
'''
Brief: Plotting histos utilities specialized for
d(e,e'p) fsi studies proposal
'''

def make_ratios_d2fsi(pm_set, thrq_set, scale, plot_flag=''):


    
                
    if plot_flag=='ratio' or plot_flag=='ratio_err':
        
        # set figure subplots for ratio
        #fig, ax = plt.subplots(5, 8, sharex='col', sharey='row')
        # fig, ax = plt.subplots(5, 8) #original

        fig, ax = plt.subplots(4, 3)  # only pm ~520 - 960 (12 plots)
        #fig, ax = plt.subplots(3, 3)  # only pm = 500 setting (9 plots)
        
        #fig.text(0.5, 0.01, r'Recoil Angle $\theta_{nq}$ [deg]', ha='center', fontsize=12)
        #fig.text(0.01, 0.5, r'R = FSI / PWIA', va='center', rotation='vertical', fontsize=12)
        #subplot_title ='angular distributions FSI/PWIA ratio'   #setting: (%d MeV, %d deg)'%(pm_set, thrq_set)
        #plt.suptitle(subplot_title, fontsize=15);
        fig.set_size_inches(8,10, forward=True)

    offset=0
    # loop over central missing momentum kin. setting 
    for ipm in pm_set:

             
        scl_idx = 0  # scale index
        
        # loop over central recoil angle kin. setting for a given central momentum
        for ithrq in thrq_set:

            offset = offset + 1
            
            print('ithrq: ', ithrq)
            hist_file                 = 'H_Pm_vs_thrq_yield_d2fsi_pm%d_thrq%d.txt'%(ipm, ithrq)  # histogram file with numerical info
            histos_file_path_pwia = 'yield_estimates/d2_fsi/histogram_data/pm%d_thrq%d_pwia/%s'%(ipm, ithrq, hist_file)
            histos_file_path_fsi  = 'yield_estimates/d2_fsi/histogram_data/pm%d_thrq%d_fsi/%s'%(ipm, ithrq, hist_file)

            # read histogram param
            pm_binw   = float(get_label('ybin_width', histos_file_path_pwia))
            thrq_binw = float(get_label('xbin_width', histos_file_path_pwia))

                        
            rel_err_thrs = 0.4   #  relative stat. error threshold for masking

            # read dataframe
            df_fsi  = pd.read_csv(histos_file_path_fsi,  comment='#')
            df_pwia = pd.read_csv(histos_file_path_pwia, comment='#')

            # get central bin values arrays
            thrq_bins = (df_fsi.x0).to_numpy()
            pm_bins   = (df_fsi.y0[thrq_bins==thrq_bins[0]]).to_numpy()

            
            # get bin content / bin content error
            fsi_N       = df_fsi.zcont * scale[scl_idx]
            fsi_Nerr    = np.sqrt(fsi_N) 
            fsi_rel_err = fsi_Nerr / fsi_N

            fsi_N    = ma.masked_where(fsi_N==0, fsi_N)
            fsi_Nerr = ma.masked_where(fsi_N==0, fsi_Nerr)

            pwia_N       = df_pwia.zcont * scale[scl_idx]
            pwia_Nerr    = np.sqrt(pwia_N)
            pwia_rel_err = pwia_Nerr / pwia_N
            
            pwia_N    = ma.masked_where(pwia_N==0, pwia_N)
            pwia_Nerr = ma.masked_where(pwia_N==0, pwia_Nerr)

            
            ratio = fsi_N / pwia_N
            ratio_err = ratio * np.sqrt((fsi_Nerr/fsi_N)**2 + (pwia_Nerr/pwia_N)**2)
            
            ratio_rel_err = ratio_err / ratio

            
            ratio = ma.masked_where(ratio_rel_err>rel_err_thrs,     ratio)
            ratio_err = ma.masked_where(ratio_rel_err>rel_err_thrs, ratio_err)

            rel_error = ratio_err / ratio
            y_const = np.zeros(len(rel_error))

      
            scl_idx = scl_idx + 1  # increment index for every increment in ipm

          
            #print('ratio:', ratio)
            #print('thrq_bins:', thrq_bins)
            #print('ratio[thrq=37.5]:', ratio[thrq_bins==37.5])

            idx_plot = 0 
            for idx, pm_bin in enumerate(pm_bins):

      

                cnts = np.sum(count_per_xbin )
                
                print('idx_plot', idx_plot)
                if pm_bin<0.520 or pm_bin>0.960: continue  # only for 12 plots (for pm=800 setting)

                #if pm_bin<=0.160 or pm_bin>0.520:  continue  # only for pm=500 setting

                
                if plot_flag=='ratio':

                
                  
                    
                    # do interpolation    (need to do this, but for combined thrq, not individually)    
                    #f = interp1d(thrq_bins[df_fsi.y0==pm_bin], ratio[df_fsi.y0==pm_bin], kind='linear', bounds_error=True)
                    #x_interp = np.linspace(thrq_min, thrq_max, num=100)
                    #y_interp = f(x_interp)
                    #print('x_interp;', x_interp)
                    #print('y_interp;', y_interp)


                    
                    # ---- plot ratio fsi/pwia -----
                    # ax = plt.subplot(5, 8, idx+1) original
                    ax = plt.subplot(4, 3, idx_plot+1)  # pm=800
                    #ax = plt.subplot(3, 3, idx_plot+1)  # pm=500

               
                    ax.errorbar(thrq_bins[df_fsi.y0==pm_bin], ratio[df_fsi.y0==pm_bin], ratio_err[df_fsi.y0==pm_bin], marker='o', linestyle='None', ms=5, label=r'$\theta_{nq}=%.1f$ deg'%ithrq)


                    # alternatively plot a yield (try)
                    #ax.hist(fsi_N[df_fsi.y0==pm_bin], bins=len(ybc), weights=count_per_xbin, range=[min(df.ylow), max(df.yup)], ec='k', density=False, label=r'%d counts (%.1f GeV$^{2}$)'%(cnts, jq2))


                    
                    ax.set_title('$p_{m}$ = %d $\pm$ %d MeV'%(pm_bin*1000, pm_binw*1000/2.), fontsize=16)
                    plt.axhline(1, linestyle='--', color='gray')
                    
                    plt.vlines(x = 70, ymin=1, ymax=4.5, color = 'r', linestyle = '--', linewidth=1.5) # reference line at 70 deg

                    ax.set_xlim(20,90)
                    ax.set_ylim(0,4.)
                    plt.xticks(fontsize = 16)
                    plt.yticks(fontsize = 16)

                    if pm_bin==0.520:
                        plt.legend(loc='upper left')
 
                    # limit the number of ticks
                    max_xticks = 4
                    max_yticks = 4
                    xloc = plt.MaxNLocator(max_xticks)
                    yloc = plt.MaxNLocator(max_yticks)
                    ax.xaxis.set_major_locator(xloc)
                    ax.yaxis.set_major_locator(yloc)

             
                if plot_flag=='ratio_err':

                    ax = plt.subplot(4, 3, idx_plot+1)
                    
                    # mask if exceed rel error thresh
                    rel_error_m = ma.masked_where(rel_error>rel_err_thrs, rel_error)
                    y_const_m = ma.masked_where(rel_error>rel_err_thrs, y_const)
                    thrq_bins_m = ma.masked_where(rel_error>rel_err_thrs,thrq_bins)
                    
                    
                    thrq_bins_m = thrq_bins_m + offset
                    ax.errorbar(thrq_bins_m[df_fsi.y0==pm_bin], y_const_m[df_fsi.y0==pm_bin], rel_error_m[df_fsi.y0==pm_bin], marker='o', linestyle='None', ms=5, label=r'$\theta_{nq}=%.1f$ deg'%ithrq)
                    
                    
                    ax.set_title('$p_{m}$ = %d $\pm$ %d MeV'%(pm_bin*1000, pm_binw*1000/2.), fontsize=20)
                    plt.axhline(0, linestyle='--', color='gray')
                
                    plt.axhline(y = 0.20, color = 'r', linestyle = '--')
                    plt.axhline(y = -0.20, color = 'r', linestyle = '--')
                        
                    ax.set_xlim(20,90)
                    ax.set_ylim(-0.5,0.5)
                    plt.xticks(fontsize = 20)
                    plt.yticks(fontsize = 20)
                    
                    if pm_bin==0.520:
                        plt.legend(loc='lower right')
                        
                    # limit the number of ticks
                    max_xticks = 4
                    max_yticks = 4
                    xloc = plt.MaxNLocator(max_xticks)
                    yloc = plt.MaxNLocator(max_yticks)
                    ax.xaxis.set_major_locator(xloc)
                    ax.yaxis.set_major_locator(yloc)

                    
                         
                idx_plot = idx_plot+1
   
    plt.tight_layout()
   
    plt.show()
    plt.savefig('test.png')

def make_projY_d2fsi(h2_hist_name, pm_user, thrq_user, model, plot_flag, scale=[1,1,1]):


    ifig = 1 # counter for 2d histogram figures

    # define collimator polygon shape (for plotting contour lines) 
    shms_hsize = 8.5
    shms_vsize = 12.5
    
    hms_hsize = 4.575
    hms_vsize = 11.646
    
    
    coord_shms = [[  shms_hsize,     shms_vsize/2.],
                  [  shms_hsize/2.,  shms_vsize   ],
                  [ -shms_hsize/2.,  shms_vsize   ],
                  [ -shms_hsize,     shms_vsize/2.],
                  [ -shms_hsize,    -shms_vsize/2.],
                  [ -shms_hsize/2., -shms_vsize   ],
                  [  shms_hsize/2., -shms_vsize   ],
                  [  shms_hsize,    -shms_vsize/2.],
                  [  shms_hsize,     shms_vsize/2.]]
    
    
    
    coord_hms = [[  hms_hsize,     hms_vsize/2.],
                 [  hms_hsize/2.,  hms_vsize   ],
                 [ -hms_hsize/2.,  hms_vsize   ],
                 [ -hms_hsize,     hms_vsize/2.],
                 [ -hms_hsize,    -hms_vsize/2.],
                 [ -hms_hsize/2., -hms_vsize   ],
                 [  hms_hsize/2., -hms_vsize   ],
                 [  hms_hsize,    -hms_vsize/2.],
                 [  hms_hsize,     hms_vsize/2.]]


    coord_shms.append(coord_shms[0])
    coord_hms.append(coord_hms[0])
    x_shms, y_shms = zip(*coord_shms)
    x_hms, y_hms = zip(*coord_hms)
 
    # loop over central pm setting
    for ipm in pm_user:

        idx=0
        for ithrq in thrq_user:
            
            # set histo base name and file path
            h2_hist_basename = 'H_%s_yield_d2fsi_pm%d_thrq_%d.txt'%(h2_hist_name, ipm, ithrq)   # generic histogram name
            hist_file_path = 'yield_estimates/d2_fsi/histogram_data/pm%d_thrq%d_%s/%s'%(ipm, ithrq, model, h2_hist_basename)  #optimized kinematics

            print('hist_file_path:', hist_file_path)
            # check if file exists, else continue reading next file
            if not os.path.isfile(hist_file_path):

                # try again 
                h2_hist_basename = 'H_%s_yield_d2fsi_pm%d_thrq%d.txt'%(h2_hist_name, ipm, ithrq)   # generic histogram name
                hist_file_path = 'yield_estimates/d2_fsi/histogram_data/pm%d_thrq%d_%s/%s'%(ipm, ithrq, model, h2_hist_basename)  #optimized kinematics

                if not os.path.isfile(hist_file_path):
                    print('FILE DOES NOT EXIST !')
                    continue

            df = pd.read_csv(hist_file_path, comment='#')

            # get histo parameters from .txt file
            xlabel = get_label('xlabel',     hist_file_path)
            ylabel = get_label('ylabel',     hist_file_path)
            title  = get_label('title',      hist_file_path)
            ybinw  = float(get_label('ybin_width', hist_file_path))
            xbinw  = float(get_label('xbin_width', hist_file_path))
            nxbins = len(df.xb[df.yb==df.yb[0]])
            nybins = len(df.yb[df.xb==df.xb[0]])
            ybc = (df.y0[df.x0==df.x0[0]]).to_numpy() # y-bin central value
            xbc = (df.x0[df.y0==df.y0[0]]).to_numpy() # x-bin central value
            if "_2Davg" in h2_hist_basename:
                print('not setting up counts')
            else:
                df.zcont = df.zcont * scale[idx]
                counts = np.sum(df.zcont)

            # index counter for thrq
            idx = idx+1
            
            if plot_flag=='2d':

                print('PASS1')
                # plotting the 2d histo
                fig2d = plt.figure(ifig, figsize=(9,6))

                if "_2Davg" in h2_hist_basename:
                    plt.hist2d(df.x0 ,df.y0, weights=df.zcont, bins=(nxbins, nybins), range=[ [min(df.xlow), max(df.xup)], [min(df.ylow), max(df.yup)]], cmap = 'viridis')
                else:                        
                    plt.hist2d(df.x0 ,df.y0, weights=df.zcont, bins=(nxbins, nybins), range=[ [min(df.xlow), max(df.xup)], [min(df.ylow), max(df.yup)]], cmap = 'viridis', norm=mcolors.LogNorm())
                    print('TRIED TO PLOT')
                    # plot the countour lines of collimator when appropiate
                    if "hXColl_vs_hYColl" in h2_hist_basename:
                        plt.plot(x_hms,y_hms, color='r', linewidth=2)  # plot hms collimator geometry contour line
                    if "eXColl_vs_eYColl" in h2_hist_basename:
                        plt.plot(x_shms,y_shms, color='r', linewidth=2)  # plot shms collimator geometry contour line
                        
                # limit the number of ticks
                #max_xticks = 4
                #max_yticks = 4
                #xloc = plt.MaxNLocator(max_xticks)
                #yloc = plt.MaxNLocator(max_yticks)
                #axs.xaxis.set_major_locator(xloc)
                #axs.yaxis.set_major_locator(yloc)

                plt.xlabel(xlabel, fontsize=12)
                plt.ylabel(ylabel, fontsize=12)
                plt.title(title,   fontsize=14)
                plt.xticks(fontsize = 30)
                plt.yticks(fontsize = 30)
                cbar = plt.colorbar()
                cbar.ax.tick_params(labelsize=30)

         

            
                
                ifig = ifig+1
            
                plt.grid()
                plt.show()
                
# original method (without multiple targets option)
def make_projY_d2pol(h2_hist_name, pm_user, Q2_user, model, plot_flag, scale=1):
    # NEED TO FIX THIS FUNCTION, AS IT CURRENTLY DISPLAYS MULTIPLE SUBPLOTS, WHERE ONLY ONE IS NEEDED
    '''
    Brief: generic function makes 1D projections along y-axis (slicing xbins) for selected 2D histos,

    -----------
    user input
    -----------
    h2_hist_name: 2D histogram base name (if unfamiliar, will need to check the name in the histogram directory
               the generic histogram name has the form: H_[h2_hist_name]_yield_d2pol_pm[pm_set]_Q2_[Q2_set].txt  
               where the variable in brackets is replaced by user input
    
    pm_set: central missing momentum, e.g.  pm_set=[200], etc. (ONLY single-value can be used)
            since the code was structured so that for a single pm_set, loop over all Q2_set
            (and overlay all Q2 settings) 
    
    Q2_set: central Q2 setting, e.g. Q2_set=[3.5, 4.0, 4.5], (a single valued or multi-valued list can be provided)

    model: can be "fsi" or "pwia" depending on the simulation done

    plot_flag:  plot_flag="2d"       -> produces a standard 2D histogram (with log or linear z scale, depending on which type of 2D found)
                plot_flag="proj"     -> produces a sliced set of subplots, where each subplot is the Y-projection of the 2D histogram
                plot_flag="proj_err" -> produces a sliced set of subplots, where each subplot is the Y-projection relative error (dN/N) of the 2D histogram
    
    scale: time scale factor to scale yield by a multiple of the simulated beam-time, for example, if 168 hrs
    is the standard simulation beam time, and user sets scale = 3, then yield will be scaled to
    3 * yield, which is esentially triple the beam time for the yield (3 * 168 hrs)
    
    '''
    
    #histos_file_path = 'yield_estimates/d2_fsi/histogram_data/pm%d_thrq%d_%s/'%(pm_user, thrq_user, model)

    rel_err_thrs = 0.5 # mask >30 % relative error
   
    ifig = 1 # counter for 2d histogram figures

    
    # loop over central pm setting
    for ipm in pm_user:

        #ipm = ikin_set[0]
        #jq2 = ikin_set[1]
        
        if plot_flag=='proj' or plot_flag=='proj_err':
            # set figure subplots for the 1d projections
            #fig, ax = plt.subplots(6, 3, sharex='col', sharey='row', tight_layout=True)
            fig, ax = plt.subplots(3, 2, tight_layout=True)

            fig.text(0.5, 0.002, 'missing momentum, p$_{m}$ [GeV/c]', ha='center', fontsize=12)
            fig.text(0.005, 0.5, 'counts', va='center', rotation='vertical', fontsize=12)
            subplot_title =  r"p$_{m}$ vs. $\theta_{nq}$ 1d projection (%s), central p$_{m}$ setting: %d MeV"%(model, ipm)
            plt.suptitle(subplot_title, fontsize=15);
            fig.set_size_inches(12,10, forward=True)
        

        offset=0 # offset for applying projections overlay 
        # loop over central q2 setting
        for jq2 in Q2_user:

            # increment offset for every new Q2 setting
            if jq2 > Q2_user[0]:
                offset = offset + 0.003                                
                
            
            # set histo base name and file path
            h2_hist_basename = 'H_%s_yield_d2pol_pm%d_Q2_%.1f.txt'%(h2_hist_name, ipm, jq2)   # generic histogram name
            #hist_file_path = 'yield_estimates/d2_pol/smallFSI/phi_0deg/histogram_data/pm%d_Q2_%.1f_%s/%s'%(ipm, jq2, model, h2_hist_basename)             #original kinematics
            hist_file_path = 'yield_estimates/d2_pol/smallFSI/optimized/histogram_data/pm%d_Q2_%.1f_%s/%s'%(ipm, jq2, model, h2_hist_basename)  #optimized kinematics

            print('hist_file_path:', hist_file_path)
            # check if file exists, else continue reading next file
            if not os.path.isfile(hist_file_path): continue

                    
            df = pd.read_csv(hist_file_path, comment='#')

            # get histo parameters from .txt file
            xlabel = get_label('xlabel',     hist_file_path)
            ylabel = get_label('ylabel',     hist_file_path)
            title  = get_label('title',      hist_file_path)
            ybinw  = float(get_label('ybin_width', hist_file_path))
            xbinw  = float(get_label('xbin_width', hist_file_path))
            nxbins = len(df.xb[df.yb==df.yb[0]])
            nybins = len(df.yb[df.xb==df.xb[0]])
            ybc = (df.y0[df.x0==df.x0[0]]).to_numpy() # y-bin central value
            xbc = (df.x0[df.y0==df.y0[0]]).to_numpy() # x-bin central value
            if "_2Davg" in h2_hist_basename:
                print('not setting up counts')
            else:
                df.zcont = df.zcont * scale
                counts = np.sum(df.zcont)

      

            if plot_flag=='2d':
                
                # plotting the 2d histo
                fig2d = plt.figure(ifig, figsize=(9,6))

                if "_2Davg" in h2_hist_basename:
                    plt.hist2d(df.x0 ,df.y0, weights=df.zcont, bins=(nxbins, nybins), range=[ [min(df.xlow), max(df.xup)], [min(df.ylow), max(df.yup)]], cmap = 'viridis')
                else:
                    plt.hist2d(df.x0 ,df.y0, weights=df.zcont, bins=(nxbins, nybins), range=[ [min(df.xlow), max(df.xup)], [min(df.ylow), max(df.yup)]], cmap = 'viridis', norm=mcolors.LogNorm())
                    plt.text(0.6*(df.x0[df.y0==df.y0[0]]).max(), 0.7*(df.y0[df.x0==df.x0[0]]).max(), r"Q$^{2}$=%.1f GeV$^{2}$"%(jq2)+"\n"+"$P_{m}$=%d MeV"%(ipm)+"\n"+"(counts = %d)"%(counts), fontsize=12)


                plt.xlabel(xlabel, fontsize=12)
                plt.ylabel(ylabel, fontsize=12)
                plt.title(title,   fontsize=14)
                plt.xticks(fontsize = 22)
                plt.yticks(fontsize = 22)

                
                plt.colorbar()
                ifig = ifig+1
            
                plt.grid()
                plt.show() 

            jdx=0 #counter for subplots which have non-zero counts
            
            #loop over x-bins (for y-projections)
            for idx, xbin in enumerate( xbc ):
                           
                count_per_xbin     = df.zcont[df.x0==xbin]
                count_per_xbin_err = np.sqrt( count_per_xbin )

            
                count_per_xbin     = ma.masked_where(count_per_xbin==0, count_per_xbin)
                count_per_xbin_err = ma.masked_where(count_per_xbin==0, count_per_xbin_err)


                cnts = np.sum(count_per_xbin )

                    
                # variables for plotting relative error
                count_per_xbin_rel_err = count_per_xbin_err / count_per_xbin
                y_const = np.zeros(len(count_per_xbin_rel_err))
                

                #---------------

                if(np.ma.is_masked(cnts)): continue
            
                if(not np.ma.is_masked(cnts)):

                    if plot_flag=='proj':
                        ax = plt.subplot(3, 2, jdx+1)

                        if "_2Davg" in h2_hist_basename:
                            ax.errorbar(ybc, count_per_xbin, count_per_xbin_err, marker='o', markersize=4, linestyle='None', label=r'%.1f GeV$^{2}$'%(jq2)) #//, label=r'%d counts'%(cnts))
                            #ax.hist(ybc, bins=len(ybc), weights=count_per_xbin, range=[min(df.ylow), max(df.yup)], alpha=0.5, ec='k', density=False, label=r'(%.1f GeV$^{2}$)'%(jq2))
                            ax.set_ylim(-20,20)
                        else:
                            ax.hist(ybc, bins=len(ybc), weights=count_per_xbin, range=[min(df.ylow), max(df.yup)], ec='k', density=False, label=r'%d counts (%.1f GeV$^{2}$)'%(cnts, jq2))

                        plt.title(r'$\theta_{nq}$ = %d $\pm$ %d deg'%(xbin, xbinw/2.))
                        plt.legend(frameon=False, loc='upper right')
                        jdx = jdx+1
                        
                    # ---- relative error plots ----
                    if plot_flag=='proj_err':

                        ax = plt.subplot(3, 2, jdx+1)
                        count_per_xbin_rel_err_m = ma.masked_where(count_per_xbin_rel_err>rel_err_thrs, count_per_xbin_rel_err)

                        y_const_m = ma.masked_where(count_per_xbin_rel_err>rel_err_thrs, y_const)
                        ybc_m = ma.masked_where(count_per_xbin_rel_err>rel_err_thrs, ybc)

                        ybc_m = ybc_m + offset
                        ax.errorbar(ybc_m, y_const_m, count_per_xbin_rel_err_m, marker='o', markersize=4, linestyle='None', label=r'%.1f GeV$^{2}$'%(jq2)) #//, label=r'%d counts'%(cnts))
                        plt.axhline(y = 0.20, color = 'r', linestyle = '--')
                        plt.axhline(y = -0.20, color = 'r', linestyle = '--')
                        plt.xlim(0.15, 0.72)
                        plt.title(r'$\theta_{nq}$ = %d $\pm$ %d deg'%(xbin, xbinw/2.))
                        jdx = jdx+1
            plt.legend(frameon=False, loc='upper right')

                  
       
    if plot_flag=='proj' or plot_flag=='proj_err':
        plt.show()




def calc_dilution(pm_user, Q2_user, model, field, scale=1):

    '''
    Brief: function to calculate the d2pol dilution factor: dilution = d2 / (d2 + he4 + n14)
    '''

    overlay_flag = True
    
    # create file to write dilution factors
    ofname = 'd2pol_dilution_factors_pm%d_Q2_%.1f_phi0_%s.csv' %(pm_user, Q2_user, field)
    ofile = open(ofname, 'w+')
    ofile.write('# d2pol dilution factors (%s) \n'%(field))
    ofile.write('# \n'
                '# Header Definitions: \n'
                '#        \n'
                '# pm_bin      : missing momentum bin (GeV/c)\n'
                '# thrq_bin    : neutron recoil angle bin (deg)'
                '# d2          : d(e ep) counts  \n'
                '# d2_err       : d(e ep) counts absolute error  \n'
                '# he4          : he4(e ep) counts  \n'
                '# he4_err       : he4(e ep) counts absolute error  \n'
                '# n14          : n14(e ep) counts  \n'
                '# n14_err       : n14(e ep) counts absolute error  \n'
                '# dilution       : dilution factor defined as N_D / (N_D + N_14 + N_He)\n'
                '# dilution_err   : absolute error on dilution factor\n'
        
                )
    ofile.write('thrq_bin,pm_bin,d2,d2_err,he4,he4_err,n14,n14_err,dilution,dilution_err\n') 



    # if overlay_flag == False:
    
    # uncomment comment for plotting dilution factor into subplots binned in thrq 
    
    fig, ax = plt.subplots(3, 3)
    fig.set_size_inches(10,10, forward=True)
    fig.text(0.5, 0.01, 'missing momentum, p$_{m}$ [GeV/c]', ha='center', fontsize=14)
    fig.text(0.01, 0.5, 'dilution', va='center', rotation='vertical', fontsize=14)
    subplot_title =  r"dilution factor, central p$_{m,cent}$ = %d MeV, $Q^{2}$ = %.1f GeV$^{2}$"%(pm_user, Q2_user)
    plt.suptitle(subplot_title, fontsize=15);
    plt.tight_layout()
    
    #------------------------------------------------------
    
    
    # h2 -> 2d histo (not related to hydrogrn)
    #h2_hist_basename = 'H_Pm_vs_thrq_dil_yield_d2pol_pm%d_Q2_%.1f.txt'%(pm_user, Q2_user)   # generic histogram name  (finer bins)
    h2_hist_basename = 'H_Pm_vs_thrq_yield_d2pol_pm%d_Q2_%.1f.txt'%(pm_user, Q2_user)   # generic histogram name (coarser bins)
    
    d2_file_path = 'yield_estimates/d2_pol/smallFSI/optimized/histogram_data/phi0/d2_pm%d_Q2_%.1f_%s_%s/%s'%(pm_user, Q2_user, model, field, h2_hist_basename)  #optimized kinematics
    he4_file_path = 'yield_estimates/d2_pol/smallFSI/optimized/histogram_data/phi0/he4_pm%d_Q2_%.1f_%s_%s/%s'%(pm_user, Q2_user, model, field, h2_hist_basename)  #optimized kinematics
    n14_file_path = 'yield_estimates/d2_pol/smallFSI/optimized/histogram_data/phi0/n14_pm%d_Q2_%.1f_%s_%s/%s'%(pm_user, Q2_user, model, field, h2_hist_basename)  #optimized kinematics
    
    df_d2  = pd.read_csv(d2_file_path, comment='#')
    df_he4 = pd.read_csv(he4_file_path, comment='#')
    df_n14 = pd.read_csv(n14_file_path, comment='#')
    
    # get histo parameters from .txt file (general params)
    xlabel = get_label('xlabel',     d2_file_path)
    ylabel = get_label('ylabel',     d2_file_path)
    title  = get_label('title',      d2_file_path)
    ybinw  = float(get_label('ybin_width', d2_file_path))
    xbinw  = float(get_label('xbin_width', d2_file_path))
    nxbins = len(df_d2.xb[df_d2.yb==df_d2.yb[0]])
    nybins = len(df_d2.yb[df_d2.xb==df_d2.xb[0]])
    ybc = (df_d2.y0[df_d2.x0==df_d2.x0[0]]).to_numpy() # y-bin central value (pm)
    xbc = (df_d2.x0[df_d2.y0==df_d2.y0[0]]).to_numpy() # x-bin central value (thrq)
    
    # scale the counts per bin (default scale = 1)
    df_d2.zcont  = df_d2.zcont * scale
    df_he4.zcont = df_he4.zcont * scale
    df_n14.zcont = df_n14.zcont * scale
    

    jdx=0 #counter for subplots which have non-zero counts

    #loop over x-bins i.e., thrq bins (for y-projections)
    for idx, xbin in enumerate( xbc ):

        # ignore thetarq >110
        if xbin>=110: continue

        # uncomment for plotting dilution factor into subplots binned in thrq 
        ax = plt.subplot(3, 3, jdx+1)

        
        d2_count_per_xbin     = df_d2.zcont[df_d2.x0==xbin]
        d2_count_per_xbin_err = np.sqrt( d2_count_per_xbin )
        
        
        d2_count_per_xbin     = ma.masked_where(d2_count_per_xbin==0, d2_count_per_xbin)
        d2_count_per_xbin_err = ma.masked_where(d2_count_per_xbin==0, d2_count_per_xbin_err)


        he4_count_per_xbin     = df_he4.zcont[df_he4.x0==xbin]
        he4_count_per_xbin_err = np.sqrt( he4_count_per_xbin )
        
        # mask values
        he4_count_per_xbin     = ma.masked_where(he4_count_per_xbin==0, he4_count_per_xbin)
        he4_count_per_xbin_err = ma.masked_where(he4_count_per_xbin==0, he4_count_per_xbin_err)

        n14_count_per_xbin     = df_n14.zcont[df_n14.x0==xbin]
        n14_count_per_xbin_err = np.sqrt( n14_count_per_xbin )
        
        
        n14_count_per_xbin     = ma.masked_where(n14_count_per_xbin==0, n14_count_per_xbin)
        n14_count_per_xbin_err = ma.masked_where(n14_count_per_xbin==0, n14_count_per_xbin_err)

        
        
        # rename variables for simplicity
        x     = d2_count_per_xbin
        sigx  = d2_count_per_xbin_err
        y     = he4_count_per_xbin
        sigy  = he4_count_per_xbin_err
        z     = n14_count_per_xbin
        sigz  = n14_count_per_xbin_err
        
        dilution     =  x / ( x + y + z)
        dilution_err =  np.sqrt( ((sigx**2 * (y+z)**2 ) + (x**2 * (sigy**2 + sigz**2))) / (x + y + z)**4 )

        dsum = np.sum( dilution )
            
        if(np.ma.is_masked(dsum)): continue

        if(not np.ma.is_masked(dsum)):


            # get min/max values (ignoring pm values correponding to masked dilution) -- for interpolation purpose
            pm_min = ma.min(ma.masked_where(ma.getmask(dilution), ybc) )
            pm_max = ma.max(ma.masked_where(ma.getmask(dilution), ybc) )
            
            
            # do interpolation        
            f = interp1d(ybc, dilution, kind='linear', bounds_error=True)
            x_interp = np.linspace(pm_min, pm_max, num=100)
            y_interp = f(x_interp)
            
            

        
            #if overlay_flag == True:

            
            # -- plotting option: overlay dilution factors for all thrq_bins ---
            '''
            # uncomment if plotting overlay for all thrq_bins 
            fig= plt.subplot()
            plt.plot(x_interp,  y_interp, marker='None', linestyle='--', label=r'$\theta_{nq}$ = %d $\pm$ %d deg'%(xbin, xbinw/2.) )
            plt.xticks([0.0, 0.2, 0.4, 0.6], fontsize = 32)
            plt.yticks(fontsize = 32)
            plt.ylim(0, 1.0)
            plt.xlim(0, 0.65)
            plt.legend(fontsize=16, loc='lower right')
            '''
            #--------------------------------------------------------------------
            
            
            # loop over y-bins (pm_bins) to write to file
            for ipm, ybin in enumerate( ybc ):
            
                ofile.write("%.1f,%.3f,%.1f,%.3f,%.1f,%.3f,%.1f,%.3f,%.3f,%.3f\n" % (xbin, ybin, x[ipm], sigx[ipm], y[ipm], sigy[ipm], z[ipm], sigz[ipm], dilution[ipm], dilution_err[ipm] ))


            
            #uncomment if plotting individual subplots
            # plot interpolated function
            ax.plot(x_interp,  y_interp, marker='None', alpha=0.9, linestyle='--', color='r')

            # plot data
            ax.errorbar(ybc, dilution, dilution_err, marker='o', markersize=8, alpha=0.4, linestyle='None', color='r', label=r'%.1f GeV$^{2}$'%(Q2_user))

            
            plt.title(r'$\theta_{nq}$ = %d $\pm$ %d deg'%(xbin, xbinw/2.), fontsize=15)
            plt.xticks(fontsize = 18)
            plt.yticks(fontsize = 18)
            plt.ylim(0.,1.0)
            plt.xlim(0, 0.7)
            
            
            jdx = jdx+1
       
    ofile.close()
    plt.show()

def overlay_dilution():

    print('overlay_dilution')
    fig, ax = plt.subplots(3, 3)
    fig.set_size_inches(10,9, forward=True)
    #fig.text(0.5, 0.01, 'missing momentum, p$_{m}$ [GeV/c]', ha='center', fontsize=14)   
    #fig.text(0.01, 0.5, 'dilution factor', va='center', rotation='vertical', fontsize=14)
    #subplot_title =  r"dilution factor"
    #plt.suptitle(subplot_title, fontsize=15);
    plt.tight_layout()
            
    df_fieldON = pd.read_csv('d2pol_dilution_factors_pm350_Q2_2.5_fieldON.csv', comment='#')
    df_fieldOFF = pd.read_csv('d2pol_dilution_factors_pm350_Q2_2.5_fieldOFF.csv', comment='#')

    
    pm_bin   = df_fieldON.pm_bin[df_fieldON.thrq_bin==df_fieldON.thrq_bin[0]]
    thrq_bin = df_fieldON.thrq_bin[df_fieldON.pm_bin==df_fieldON.pm_bin[0]]

    jdx = 0
    #loop over x-bins (thrq) 
    for idx, xbin in enumerate( thrq_bin ):

        
        ax = plt.subplot(3, 3, jdx+1)
        
        dilution_per_xbin_fieldON     = df_fieldON.dilution[df_fieldON.thrq_bin==xbin]
        dilution_per_xbin_err_fieldON = df_fieldON.dilution_err[df_fieldON.thrq_bin==xbin]
        
        dilution_per_xbin_fieldOFF     = df_fieldOFF.dilution[df_fieldOFF.thrq_bin==xbin]
        dilution_per_xbin_err_fieldOFF = df_fieldOFF.dilution_err[df_fieldOFF.thrq_bin==xbin]
       
        
        print('pm_bin:',pm_bin, 'dilution (field ON):', dilution_per_xbin_fieldON, 'dilution_err (field ON):', dilution_per_xbin_err_fieldON, )
        ax.errorbar(pm_bin, dilution_per_xbin_fieldOFF , dilution_per_xbin_err_fieldOFF, alpha=0.6, marker='s', markersize=8, linestyle='None', color='k', label=r'field OFF')
        ax.errorbar(pm_bin, dilution_per_xbin_fieldON , dilution_per_xbin_err_fieldON, alpha=0.6, marker='o', markersize=8, linestyle='None', color='r', label=r'field ON')                   
       
        plt.ylim(0.3,1.2)
        plt.xlim(0.0,1.08)
        plt.xticks(fontsize = 18)
        plt.yticks(fontsize = 18)
        plt.title(r'$\theta_{nq}$ = %d $\pm$ %d deg'%(xbin, 10), fontsize=14)
        jdx = jdx+1
    plt.legend()
    plt.show()
        
    
def make_projY_d2pol(h2_hist_name, tgt_set, pm_user, Q2_user, model, field, plot_flag, scale=1):
    '''
    Brief: generic function makes 1D projections along y-axis (slicing xbins) for selected 2D histos,
    (this version of the method takes multiple targets and a single value of Q2)
    -----------
    user input
    -----------
    h2_hist_name: 2D histogram base name (if unfamiliar, will need to check the name in the histogram directory
               the generic histogram name has the form: H_[h2_hist_name]_yield_d2pol_pm[pm_set]_Q2_[Q2_set].txt  
               where the variable in brackets is replaced by user input
    
    tgt_set: target nuclei, e.g. tgt_set=['n14', 'd2', 'he4']  provide multiple targets will overlay

    pm_set: central missing momentum, e.g.  pm_set=[200], etc. (ONLY single-value can be used)
            since the code was structured so that for a single pm_set, loop over all Q2_set
            (and overlay all Q2 settings) 
    
    Q2_set: central Q2 setting, e.g. Q2_set=3.5, (a single valued MUST be provided)

    model: can be "fsi" or "pwia" depending on the simulation done
    field: "fieldON" or "fieldOFF" (target magnetic field)

    plot_flag:  plot_flag="2d"       -> produces a standard 2D histogram (with log or linear z scale, depending on which type of 2D found)
                plot_flag="proj"     -> produces a sliced set of subplots, where each subplot is the Y-projection of the 2D histogram
                plot_flag="proj_err" -> produces a sliced set of subplots, where each subplot is the Y-projection relative error (dN/N) of the 2D histogram
    
    scale: time scale factor to scale yield by a multiple of the simulated beam-time, for example, if 168 hrs
    is the standard simulation beam time, and user sets scale = 3, then yield will be scaled to
    3 * yield, which is esentially triple the beam time for the yield (3 * 168 hrs)
    
    '''
    

    rel_err_thrs = 0.5 # mask >30 % relative error

    clr = ['orange', 'deepskyblue', 'magenta']
    
    ifig = 1 # counter for 2d histogram figures


    # define collimator polygon shape (for plotting contour lines) 
    shms_hsize = 8.5
    shms_vsize = 12.5
    
    hms_hsize = 4.575
    hms_vsize = 11.646
    
    
    coord_shms = [[  shms_hsize,     shms_vsize/2.],
                  [  shms_hsize/2.,  shms_vsize   ],
                  [ -shms_hsize/2.,  shms_vsize   ],
                  [ -shms_hsize,     shms_vsize/2.],
                  [ -shms_hsize,    -shms_vsize/2.],
                  [ -shms_hsize/2., -shms_vsize   ],
                  [  shms_hsize/2., -shms_vsize   ],
                  [  shms_hsize,    -shms_vsize/2.],
                  [  shms_hsize,     shms_vsize/2.]]
    
    
    
    coord_hms = [[  hms_hsize,     hms_vsize/2.],
                 [  hms_hsize/2.,  hms_vsize   ],
                 [ -hms_hsize/2.,  hms_vsize   ],
                 [ -hms_hsize,     hms_vsize/2.],
                 [ -hms_hsize,    -hms_vsize/2.],
                 [ -hms_hsize/2., -hms_vsize   ],
                 [  hms_hsize/2., -hms_vsize   ],
                 [  hms_hsize,    -hms_vsize/2.],
                 [  hms_hsize,     hms_vsize/2.]]


    coord_shms.append(coord_shms[0])
    coord_hms.append(coord_hms[0])
    x_shms, y_shms = zip(*coord_shms)
    x_hms, y_hms = zip(*coord_hms)
    
    # loop over central pm setting
    for ipm in pm_user:

        #ipm = ikin_set[0]
        #jq2 = ikin_set[1]
        
        if plot_flag=='proj' or plot_flag=='proj_err':
            # set figure subplots for the 1d projections
            #fig, ax = plt.subplots(6, 3, sharex='col', sharey='row', tight_layout=True)
            fig, ax = plt.subplots(3, 3)
            fig.set_size_inches(12,10, forward=True)
           
            #fig.text(0.5, 0.01, 'missing momentum, p$_{m}$ [GeV/c]', ha='center', fontsize=14)
            #if plot_flag=='proj':
            #    fig.text(0.01, 0.5, 'counts', va='center', rotation='vertical', fontsize=14)
            #elif plot_flag=='proj_err':
            #    fig.text(0.01, 0.5, 'relative error', va='center', rotation='vertical', fontsize=14)

            subplot_title =  r"p$_{m}$ vs. $\theta_{nq}$ 1d projection (%s), central p$_{m}$ setting: %d MeV"%(model, ipm)
            #plt.suptitle(subplot_title, fontsize=15);
            plt.tight_layout()
            
        offset=0 # offset for applying projections overlay
        i = -1 # color index counter
        # loop over central targets
        for itgt in tgt_set:

            i = i+1
            
            # increment offset for every new Q2 setting
            if itgt != tgt_set[0]:
                offset = offset + 0.008                                
                
            
            # set histo base name and file path
            h2_hist_basename = 'H_%s_yield_d2pol_pm%d_Q2_%.1f.txt'%(h2_hist_name, ipm, Q2_user)   # generic histogram name
            #hist_file_path = 'yield_estimates/d2_pol/smallFSI/optimized/histogram_data/%s_pm%d_Q2_%.1f_%s_%s/%s'%(itgt, ipm, Q2_user, model, field, h2_hist_basename)  #optimized kinematics

            # phi = 0 config
            hist_file_path = 'yield_estimates/d2_pol/smallFSI/optimized/histogram_data/phi0/%s_pm%d_Q2_%.1f_%s_%s/%s'%(itgt, ipm, Q2_user, model, field, h2_hist_basename)  #optimized kinematics

            
            print('hist_file_path:', hist_file_path)
            # check if file exists, else continue reading next file
            if not os.path.isfile(hist_file_path): continue

                    
            df = pd.read_csv(hist_file_path, comment='#')

            # get histo parameters from .txt file
            xlabel = get_label('xlabel',     hist_file_path)
            ylabel = get_label('ylabel',     hist_file_path)
            title  = get_label('title',      hist_file_path)
            ybinw  = float(get_label('ybin_width', hist_file_path))
            xbinw  = float(get_label('xbin_width', hist_file_path))
            nxbins = len(df.xb[df.yb==df.yb[0]])
            nybins = len(df.yb[df.xb==df.xb[0]])
            ybc = (df.y0[df.x0==df.x0[0]]).to_numpy() # y-bin central value
            xbc = (df.x0[df.y0==df.y0[0]]).to_numpy() # x-bin central value
            if "_2Davg" in h2_hist_basename:
                print('not setting up counts')
            else:
                df.zcont = df.zcont * scale
                counts = np.sum(df.zcont)

      

            if plot_flag=='2d':
                
                # plotting the 2d histo
                fig2d = plt.figure(ifig, figsize=(9,6))

                if "_2Davg" in h2_hist_basename:
                    plt.hist2d(df.x0 ,df.y0, weights=df.zcont, bins=(nxbins, nybins), range=[ [min(df.xlow), max(df.xup)], [min(df.ylow), max(df.yup)]], cmap = 'viridis')
                else:
                    plt.hist2d(df.x0 ,df.y0, weights=df.zcont, bins=(nxbins, nybins), range=[ [min(df.xlow), max(df.xup)], [min(df.ylow), max(df.yup)]], cmap = 'viridis', norm=mcolors.LogNorm())

                    #plt.plot(x_shms,y_shms, color='r', linewidth=2)  # plot shms collimator geometry contour line
                    #plt.plot(x_hms,y_hms, color='r', linewidth=2)  # plot hms collimator geometry contour line
                    #plt.text(0.6*(df.x0[df.y0==df.y0[0]]).max(), 0.7*(df.y0[df.x0==df.x0[0]]).max(), r"Q$^{2}$=%.1f GeV$^{2}$"%(Q2_user)+"\n"+"$P_{m}$=%d MeV"%(ipm)+"\n"+"(counts = %d)"%(counts)+"\n"+"target: %s"%(itgt), fontsize=12)



                
                plt.xlabel(xlabel, fontsize=12)
                plt.ylabel(ylabel, fontsize=12)
                plt.title(title,   fontsize=14)
                plt.xticks(fontsize = 22)
                plt.yticks(fontsize = 22)
                cbar = plt.colorbar()
                cbar.ax.tick_params(labelsize=22)
                
                ifig = ifig+1
            
                plt.grid()
                plt.show() 

            jdx=0 #counter for subplots which have non-zero counts
            
            #loop over x-bins (for y-projections)
            for idx, xbin in enumerate( xbc ):
                           
                count_per_xbin     = df.zcont[df.x0==xbin]
                count_per_xbin_err = np.sqrt( count_per_xbin )

            
                count_per_xbin     = ma.masked_where(count_per_xbin==0, count_per_xbin)
                count_per_xbin_err = ma.masked_where(count_per_xbin==0, count_per_xbin_err)


                cnts = np.sum(count_per_xbin )

                    
                # variables for plotting relative error
                count_per_xbin_rel_err = count_per_xbin_err / count_per_xbin
                y_const = np.zeros(len(count_per_xbin_rel_err))
                

                #---------------

                if(np.ma.is_masked(cnts)): continue
            
                if(not np.ma.is_masked(cnts)):

                    if plot_flag=='proj':
                        ax = plt.subplot(3, 3, jdx+1)
                        
                        if "_2Davg" in h2_hist_basename:
                            ax.errorbar(ybc, count_per_xbin, count_per_xbin_err, marker='o', markersize=4, linestyle='None', label=r'%.1f GeV$^{2}$'%(Q2_user)) #//, label=r'%d counts'%(cnts))
                            #ax.hist(ybc, bins=len(ybc), weights=count_per_xbin, range=[min(df.ylow), max(df.yup)], alpha=0.5, ec='k', density=False, label=r'(%.1f GeV$^{2}$)'%(Q2_user))
                            ax.set_ylim(-20,20)
                        else:
                            ax.hist(ybc, bins=len(ybc), weights=count_per_xbin, range=[min(df.ylow), max(df.yup)], alpha=0.2,  color=clr[i], ec='k', density=False, label=r'%d counts (%s)'%(cnts, itgt))
                            plt.axhline(y = 500, color = 'r', linestyle = '--')
                            ax.set_yscale('log')
                            plt.ylim(0.1, 2e5)
                            plt.xticks(fontsize = 14)
                            plt.yticks(fontsize = 14)
                        plt.title(r'$\theta_{nq}$ = %d $\pm$ %d deg'%(xbin, xbinw/2.), fontsize=15)
                        plt.legend(frameon=False, loc='upper right')
                        jdx = jdx+1
                        
                    # ---- relative error plots ----
                    if plot_flag=='proj_err':

                        ax = plt.subplot(3, 3, jdx+1)
                        count_per_xbin_rel_err_m = ma.masked_where(count_per_xbin_rel_err>rel_err_thrs, count_per_xbin_rel_err)

                        y_const_m = ma.masked_where(count_per_xbin_rel_err>rel_err_thrs, y_const)
                        ybc_m = ma.masked_where(count_per_xbin_rel_err>rel_err_thrs, ybc)

                        ybc_m = ybc_m + offset
                        ax.errorbar(ybc_m, y_const_m, count_per_xbin_rel_err_m, marker='o', markersize=8,  alpha=0.4, elinewidth=2, color=clr[i], linestyle='None', label=r'%s'%(itgt)) #//, label=r'%d counts'%(cnts))
                        plt.axhline(y = 0.20, color = 'r', linestyle = '--')
                        plt.axhline(y = -0.20, color = 'r', linestyle = '--')
                        plt.xlim(0.15, 0.65)
                        plt.xticks(fontsize = 14)
                        plt.yticks(fontsize = 14)
                        plt.title(r'$\theta_{nq}$ = %d $\pm$ %d deg'%(xbin, xbinw/2.), fontsize=15)
                        jdx = jdx+1
            plt.legend(frameon=False, loc='upper right')

                  
       
    if plot_flag=='proj' or plot_flag=='proj_err':
        plt.show()






        
# call functions here (can later be passed thru steering code)
# make_plots(800, 79, 'pwia')
# make_1d_Xprojections('H_Pm_vs_thrq_yield', 800, 79, 'pwia')
# make_ratios_d2fsi([800], [28, 49, 55, 60, 66, 72, 79], 'ratio')
# overlay_d2fsi([800], [28, 49, 55], 'thrq', 'fsi')


#*************************************
#
# D2 FSI STUDIES PROPOSAL PLOTS  *
#
#*************************************

# pm=800 setting scale is in hours
#           blue, orange, green 
#thrq_set =  [72,  60,  49]
#scale_set = [160, 144, 200]  # hrs

# pm=500 setting scale is in hours
#thrq_set =  [70]
#scale_set = [12]  # hrs


#scale_set = [1, 1, 1]

# spectrometer kinematics
#overlay_d2fsi([800], thrq_set, 'Pf',  'fsi', scale_set)
#overlay_d2fsi([800], thrq_set, 'thp', 'fsi', scale_set)
#overlay_d2fsi([800], thrq_set, 'kf',  'fsi', scale_set)
#overlay_d2fsi([800], thrq_set, 'the', 'fsi', scale_set)


# electron kinematics
#overlay_d2fsi([800], thrq_set, 'nu',      'fsi', scale_set)
#overlay_d2fsi([800], thrq_set, 'xbj',     'fsi', scale_set)
#overlay_d2fsi([800], thrq_set, 'Q2_nsc',  'fsi', scale_set)
#overlay_d2fsi([800], thrq_set, 'q',       'fsi', scale_set)


# missing variables
#overlay_d2fsi([800], thrq_set, 'Pm',             'fsi', scale_set)
#overlay_d2fsi([800], thrq_set, 'Em_nuc_nsc',     'fsi', scale_set)


# angle distributions
#overlay_d2fsi([800], thrq_set, 'thq',    'fsi', scale_set)
#overlay_d2fsi([800], thrq_set, 'thpq',   'fsi', scale_set)
#overlay_d2fsi([800], thrq_set, 'thrq',   'fsi', scale_set)
#overlay_d2fsi([800], thrq_set, 'phi_pq', 'fsi', scale_set)

#===================
# plot yield ratio
#===================
#make_projY_d2fsi('Pm_vs_thrq',[800], [72], 'pwia', '2d', [160])
#make_projY_d2fsi('Pm_vs_thrq',[800], [60], 'pwia', '2d', [144])
#make_projY_d2fsi('Pm_vs_thrq',[800], [49], 'pwia', '2d', [200])

#make_ratios_d2fsi([800], [49, 60, 72], scale=[200,144,160], plot_flag='ratio_err')  # scale represent hrs
#make_ratios_d2fsi([500], [70], scale=[24], plot_flag='ratio')  # scale represent hrs



# plot 2d correlations
#make_projY_d2fsi('hXColl_vs_hYColl_nsc',[800], [72], 'fsi', '2d', [160])
#make_projY_d2fsi('eXColl_vs_eYColl_nsc',[800], [72], 'fsi', '2d', [160])


#make_projY_d2fsi('hxptar_vs_exptar',[800], [72], 'fsi', '2d', [160])
#make_projY_d2fsi('hyptar_vs_eyptar',[800], [72], 'fsi', '2d', [160])
#make_projY_d2fsi('hdelta_vs_edelta',[800], [72], 'fsi', '2d', [160])
#make_projY_d2fsi('Q2_vs_xbj_nsc',       [800], [72], 'fsi', '2d', [160])
#make_projY_d2fsi('Em_nuc_vs_Pm_nsc',[800], [72], 'fsi', '2d', [160])




#overlay_d2fsi([800], [72, 60, 49], 'Q2_nsc', 'fsi', scale=[1,1,1])
#overlay_d2fsi([800], [72, 60, 49], 'Em_nuc_nsc', 'fsi', scale=[1,1,1])

#overlay_d2fsi([500], [70], 'Pm', 'fsi', scale=[12])
#overlay_d2fsi([800], [72, 60, 49], 'Pm', 'fsi', scale=[120,132,240])

# acceptance
#overlay_d2fsi([800], thrq_set, 'edelta_nsc', 'fsi',  scale_set)
#overlay_d2fsi([800], thrq_set, 'hdelta_nsc', 'fsi', scale_set)
#make_projY_d2fsi('hXColl_vs_hYColl_nsc',[800], [72], 'fsi', '2d', [160])
#make_projY_d2fsi('eXColl_vs_eYColl_nsc',[800], [72], 'fsi', '2d', [160])






#*************************************
#
# D2 POLARIZED PROPOSAL PLOTS  *
#
#*************************************

# for overlay_2dpol() and make_projY_d2pol(), select the single-valued central momentum setting and multi-value Q2 setting for plotting
pm_set = [400]
q2_set = [2.0]
tgt_set = ['d2', 'n14', 'he4' ]
#tgt_set = ['n14']

field = 'fieldON'


scale = 1 # in multiple of weeks ( defaults to scale=1 - 2 week, if scale = 2 -> 4 weeks, . . . )

# ----- plot the kinematics variables in which a cut is used (without the self cut, ie nsc or no self cut) -----





#make_projY_d2pol('Em_nuc_vs_Pm_nsc',     ['d2'], pm_set, 2.5, 'fsi', field, '2d')
#make_projY_d2pol('hXColl_vs_hYColl_nsc',     ['d2'], pm_set, 2.0, 'fsi', field, '2d')
#make_projY_d2pol('eXColl_vs_eYColl',     ['d2'], pm_set, 2.0, 'fsi', field, '2d')


# ------ plot kinematic variables -----
# no self-cut 
#overlay_d2pol(tgt_set, pm_set, q2_set, 'Q2_nsc',     'fsi', field, scale)
#overlay_d2pol(tgt_set, pm_set, q2_set, 'Em_nuc_nsc', 'fsi', field, scale)
#overlay_d2pol(tgt_set, pm_set, q2_set, 'edelta_nsc', 'fsi', field, scale)
#overlay_d2pol(tgt_set, pm_set, q2_set, 'hdelta_nsc', 'fsi', field, scale)



# spectrometer kinematics
#overlay_d2pol(tgt_set, pm_set, q2_set, 'Pf',     'fsi', field,  scale)   # proton momentum
#overlay_d2pol(tgt_set, pm_set, q2_set, 'thp',    'fsi', field,  scale)  # proton angle
#overlay_d2pol(tgt_set, pm_set, q2_set, 'kf',     'fsi', field,  scale)   # e- momentum
#overlay_d2pol(tgt_set, pm_set, q2_set, 'the',    'fsi', field,  scale)  # e- angle


# electron kinematics
#overlay_d2pol(tgt_set, pm_set, q2_set, 'nu',     'fsi', field,  scale)   # energy transfer
#overlay_d2pol(tgt_set, pm_set, q2_set, 'xbj',    'fsi', field,  scale)  # x-bjorken
#overlay_d2pol(tgt_set, pm_set, q2_set, 'Q2_nsc',     'fsi', field, scale)
#overlay_d2pol(tgt_set, pm_set, q2_set, 'q',      'fsi', field,  scale)    # 3-momentum (q) transfer



# missing variables
#overlay_d2pol(tgt_set, pm_set, q2_set, 'Pm',     'fsi', field,  scale)   # missing momentum
#overlay_d2pol(tgt_set, pm_set, q2_set, 'Em_nuc_nsc', 'fsi', field, scale)

# angle distributions
#overlay_d2pol(tgt_set, pm_set, q2_set, 'thq',    'fsi', field,  scale)  # 3-momentum (q) angle
#overlay_d2pol(tgt_set, pm_set, q2_set, 'thpq',   'fsi', field,  scale)    # in-plane angle between (proton,q)
#overlay_d2pol(tgt_set, pm_set, q2_set, 'thrq',   'fsi', field,  scale)    # in-plane angle between (recoil,q)
#overlay_d2pol(tgt_set, pm_set, q2_set, 'phi_pq', 'fsi', field,  scale)  # out-of-plane angle between (proton, q)
#overlay_d2pol(tgt_set, pm_set, q2_set, 'cphi_pq' 'fsi', field,  scale)


# ----- plot acceptance variables ----

# reconstructed variables

#overlay_d2pol(tgt_set, pm_set, q2_set, 'exptar', 'fsi', field, scale)
#overlay_d2pol(tgt_set, pm_set, q2_set, 'eyptar', 'fsi', field, scale)
#overlay_d2pol(tgt_set, pm_set, q2_set, 'eytar',  'fsi', field, scale)
#overlay_d2pol(tgt_set, pm_set, q2_set, 'edelta', 'fsi', field, scale)

#overlay_d2pol(tgt_set, pm_set, q2_set, 'hxptar', 'fsi', field, scale)
#overlay_d2pol(tgt_set, pm_set, q2_set, 'hyptar', 'fsi', field, scale)
#overlay_d2pol(tgt_set, pm_set, q2_set, 'hytar',  'fsi', field, scale)
#overlay_d2pol(tgt_set, pm_set, q2_set, 'hdelta_nsc', 'fsi', field, scale)

# acceptance correlations

#make_projY_d2pol('hXColl_vs_hYColl_nsc', ['d2'], pm_set, 2.5, 'fsi', field, '2d')
#make_projY_d2pol('eXColl_vs_eYColl_nsc', ['d2'], pm_set, 2.5, 'fsi', field, '2d')

#make_projY_d2pol('hxptar_vs_exptar',     ['d2'], pm_set, 2.5, 'fsi', field, '2d')
#make_projY_d2pol('hyptar_vs_eyptar',     ['d2'], pm_set, 2.5, 'fsi', field, '2d')
#make_projY_d2pol('hdelta_vs_edelta',     ['d2'], pm_set, 2.5, 'fsi', field, '2d')
#make_projY_d2pol('Q2_vs_xbj',     ['d2'], pm_set, 2.0, 'fsi', field, '2d')
#make_projY_d2pol('Em_nuc_vs_Pm_nsc',  ['d2'], pm_set, 2.0, 'fsi', field, '2d')



# focal plane
#make_projY_d2pol('hxfp_vs_hyfp', tgt_set, pm_set, 2.5, 'fsi', field, '2d')
#make_projY_d2pol('exfp_vs_eyfp', tgt_set, pm_set, 2.5, 'fsi', field, '2d')


# ------ Pm vs theta_rq yield projections and errors -----

#make_projY_d2pol('Pm_vs_thrq', ['d2'], pm_set, 2.0, 'fsi', 'fieldON', '2d', 1)
#make_projY_d2pol('Pm_vs_thrq', tgt_set, pm_set, 2.0, 'fsi', 'fieldON', 'proj', 1)
#make_projY_d2pol('Pm_vs_thrq', tgt_set, pm_set, 2.0, 'fsi', 'fieldON', 'proj_err', 1)

#calc_dilution(400, 2.0, 'fsi', 'fieldON', scale=4)

#calc_dilution(350, 2.5, 'fsi', 'fieldOFF', scale=1)

#overlay_dilution()




# ----- plot 2D average kinematics --- 
#make_projY_d2pol('kf_2Davg', pm_set, q2_set, 'fsi', '2d')
#make_projY_d2pol('kf_2Davg', pm_set, q2_set, 'fsi', 'proj')

#make_projY_d2pol('the_2Davg', pm_set, q2_set, 'fsi', '2d')
#make_projY_d2pol('the_2Davg', pm_set, q2_set, 'fsi', 'proj')

#make_projY_d2pol('Q2_2Davg', pm_set, q2_set, 'fsi', '2d')
#make_projY_d2pol('Q2_2Davg', pm_set, q2_set, 'fsi', 'proj')

#make_projY_d2pol('phi_pq_2Davg', pm_set, q2_set, 'fsi', '2d')
#make_projY_d2pol('phi_pq_2Davg', pm_set, q2_set, 'fsi', 'proj')

#make_projY_d2pol('cphi_pq_2Davg', pm_set, q2_set, 'fsi', '2d')
#make_projY_d2pol('cphi_pq_2Davg', pm_set, q2_set, 'fsi', 'proj')


# plot rates
#plot_rates()
