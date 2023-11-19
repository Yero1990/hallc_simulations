import sys
import os
import numpy as np
import pandas as pd
import numpy.ma as ma
import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

# Brief: generic histos plotting utility functions script that loops over numerical histo files
# in a specified directory, and  makes and saves the plots

def get_label(label, ifile):
    with open(ifile, "r") as fp:
        for line in fp:
            if label in line:
                label = (line.split(':')[1]).strip()
                return label
    
def overlay_d2fsi(pm_set, thrq_set, hist_name, model):
    '''
    Brief: generic function to overlay 1d histograms from multiple kin. files for deut fsi
    pm_set and thrq_set are lists of values representing the different kinematic settings
    '''

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
            N = df.ycont
            Nerr = np.sqrt(N)
            
            x =   ma.masked_where(N==0, x)
            N =   ma.masked_where(N==0, N)
            Nerr = ma.masked_where(N==0, Nerr)
         
            plt.hist(df.x0, bins=len(df.x0), weights=df.ycont, alpha=0.5, ec='k', density=False, label=r"$\theta_{rq}=%d$ deg"%(ithrq)+"\n"+"$P_{m}$=%d MeV"%(ipm))
            #plt.errorbar(x, N, Nerr, linestyle='None', marker='o', mec='k', label=r"$\theta_{rq}=%d$ deg"%(ithrq)+"\n"+"$P_{m}$=%d MeV"%(ipm))

            plt.legend()
            plt.xlabel(xlabel, fontsize=15)
            plt.ylabel(ylabel, fontsize=15)
            plt.title(title, fontsize=15)

    plt.show()


def overlay_d2pol(pm_set, Q2_set, hist_name, model):
    '''
    Brief: generic function to overlay 1d histograms from multiple kin. files for deut pol. proposal
    pm_set and Q2_set are lists of values representing the different kinematic settings (needs to be checked if it works)
    '''

    rel_err_thrs = 0.3   #  relative stat. error threshold for masking

    fig, axs = plt.subplots(2, sharex=True, figsize=(7,9))
    
    # loop over central missing momentum setting
    for ipm in pm_set:

        offset=0 # for applying to overlayed data for easy visual
        # loop over central q2 setting
        for iq2 in Q2_set:

     
                
            # set histogram file path
            #histos_file_path = 'path/to/histogram_data/pm%d_q2%d_%s/histo_name_pm_set_q2_set.txt'%(pm_set, q2_set, model, hist_name)

            hist_file = 'yield_estimates/d2_pol/smallFSI/phi_180deg/histogram_data/pm%d_Q2_%.1f_%s/H_%s_yield_d2pol_pm%d_Q2_%.1f.txt'%(ipm, iq2, model, hist_name, ipm, iq2)

            print('Pm:', ipm, 'Q2:', iq2)

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
            if iq2 > Q2_set[0]:
                offset = offset + 0.1*xbinw

            
            xbc = df.x0
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
            axs[0].set_title(title)
            axs[0].hist(x=xbc, bins=len(xbc), range=[min(df.xlow), max(df.xup)], weights=N, alpha=0.5, ec='k', density=False, label="$P_{m}$=%d MeV, $Q^{2}$=%.1f GeV$^{2}$ (%d)"%(ipm, iq2, counts))
            axs[0].set_ylabel(ylabel)
            #axs[0].set_ylim(0, )  # modify ylim in case legend needs to fit better
            
            axs[0].legend()
            axs[0].xaxis.set_tick_params(labelbottom=True)
            #plt.xlabel(xlabel, fontsize=15)
            #plt.ylabel(ylabel, fontsize=15)
            #plt.title(title, fontsize=15)

            # apply offset 
            xbc_off = xbc + offset
            axs[1].errorbar(xbc_off, y_const, Nerr_rel, linestyle='None', marker='o', mec='k', label=r"$P_{m}$=%d MeV"%(ipm)+"\n"+"$Q^{2}$=%.1f GeV$^{2}$"%(iq2))
            axs[1].set_ylabel('Relative Error')
            axs[1].set_xlabel(xlabel)

           
            plt.setp(axs, xlim=(np.min(xbc)-0.5*np.min(xbc), np.max(xbc)+0.5*np.min(xbc)))
            plt.axhline(y = 0.20, color = 'r', linestyle = '--')
            plt.axhline(y = -0.20, color = 'r', linestyle = '--')
            #plt.legend()
            #plt.xlabel(xlabel, fontsize=15)
            #plt.ylabel(ylabel, fontsize=15)
            #plt.title(title, fontsize=15)
                
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

    histos_file_path = 'yield_estimates/d2_fsi/histogram_data/pm%d_thrq%d_%s/'%(pm_user, thrq_user, model)

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
             
            plt.text(0.6*(df.x0[df.y0==df.y0[0]]).max(), 0.7*(df.y0[df.x0==df.x0[0]]).max(), r"$\theta_{rq}=%d$ deg"%(thrq_user)+"\n"+"$P_{m}$=%d MeV"%(pm_user)+"\n"+"(counts = %d)"%(counts), fontsize=12)
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
            plt.text(0.6*(df.x0[df.y0==df.y0[0]]).max(), 0.7*(df.y0[df.x0==df.x0[0]]).max(), r"$\theta_{rq}=%d$ deg"%(thrq_user)+"\n"+"$P_{m}$=%d MeV"%(pm_user), fontsize=12)
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

      
            plt.hist(df.x0, bins=len(df.xb), weights=df.ycont, range=[min(df.xlow), max(df.xup)], alpha=0.2, density=False, label=r"$\theta_{rq}=%d$ deg"%(thrq_user)+"\n"+"$P_{m}$=%d MeV"%(pm_user)+ "\n"+"(counts = %d)"%(counts))
            plt.xlabel(xlabel, fontsize=12)
            plt.ylabel(ylabel, fontsize=12)
            plt.title(title,   fontsize=14)
            plt.legend(frameon=False, fontsize=12)
            plt.show()


            
'''
Brief: Plotting histos utilities specialized for
d(e,e'p) fsi studies proposal
'''

def make_ratios_d2fsi(pm_set, thrq_set, plot_flag=''):


    
                
    if plot_flag=='ratio':
        
        # set figure subplots for ratio
        fig, ax = plt.subplots(5, 8, sharex='col', sharey='row')
        fig.text(0.5, 0.01, r'Recoil Angle $\theta_{rq}$ [deg]', ha='center', fontsize=12)
        fig.text(0.01, 0.5, r'R = FSI / PWIA', va='center', rotation='vertical', fontsize=12)
        subplot_title ='angular distributions FSI/PWIA ratio'   #setting: (%d MeV, %d deg)'%(pm_set, thrq_set)
        plt.suptitle(subplot_title, fontsize=15);
        fig.set_size_inches(14,10, forward=True)
        
    # loop over central missing momentum kin. setting 
    for ipm in pm_set:

        # loop over central recoil angle kin. setting for a given central momentum
        for ithrq in thrq_set: 
            print('ithrq: ', ithrq)
            hist_file                 = 'H_Pm_vs_thrq_yield_d2fsi_pm%d_thrq%d.txt'%(ipm, ithrq)  # histogram file with numerical info
            histos_file_path_pwia = 'yield_estimates/d2_fsi/histogram_data/pm%d_thrq%d_pwia/%s'%(ipm, ithrq, hist_file)
            histos_file_path_fsi  = 'yield_estimates/d2_fsi/histogram_data/pm%d_thrq%d_fsi/%s'%(ipm, ithrq, hist_file)

            # read histogram param
            pm_binw   = float(get_label('ybin_width', histos_file_path_pwia))
            thrq_binw = float(get_label('xbin_width', histos_file_path_pwia))

                        
            rel_err_thrs = 0.3   #  relative stat. error threshold for masking

            # read dataframe
            df_fsi  = pd.read_csv(histos_file_path_fsi,  comment='#')
            df_pwia = pd.read_csv(histos_file_path_pwia, comment='#')

            # get central bin values arrays
            thrq_bins = (df_fsi.x0).to_numpy()
            pm_bins   = (df_fsi.y0[thrq_bins==thrq_bins[0]]).to_numpy()

            
            # get bin content / bin content error
            fsi_N       = df_fsi.zcont 
            fsi_Nerr    = np.sqrt(fsi_N) 
            fsi_rel_err = fsi_Nerr / fsi_N

            fsi_N    = ma.masked_where(fsi_N==0, fsi_N)
            fsi_Nerr = ma.masked_where(fsi_N==0, fsi_Nerr)

            pwia_N       = df_pwia.zcont
            pwia_Nerr    = np.sqrt(pwia_N)
            pwia_rel_err = pwia_Nerr / pwia_N
            
            pwia_N    = ma.masked_where(pwia_N==0, pwia_N)
            pwia_Nerr = ma.masked_where(pwia_N==0, pwia_Nerr)

            
            ratio = fsi_N / pwia_N
            ratio_err = ratio * np.sqrt((fsi_Nerr/fsi_N)**2 + (pwia_Nerr/pwia_N)**2)
            
            ratio_rel_err = ratio_err / ratio

            
            ratio = ma.masked_where(ratio_rel_err>rel_err_thrs,     ratio)
            ratio_err = ma.masked_where(ratio_rel_err>rel_err_thrs, ratio_err)

            #print('ratio:', ratio)
            #print('thrq_bins:', thrq_bins)
            #print('ratio[thrq=37.5]:', ratio[thrq_bins==37.5])

            for idx, pm_bin in enumerate(pm_bins):
                
                if plot_flag=='ratio':
                    
                    # ---- plot ratio fsi/pwia -----
                    ax = plt.subplot(5, 8, idx+1)
                    ax.errorbar(thrq_bins[df_fsi.y0==pm_bin], ratio[df_fsi.y0==pm_bin], ratio_err[df_fsi.y0==pm_bin], marker='o', linestyle='None', ms=5, label=r'$\theta_{rq}=%.1f$ deg'%ithrq)
                    ax.set_title('$p_{m}$ = %d $\pm$ %d MeV'%(pm_bin*1000, pm_binw*1000/2.), fontsize=10)
                    plt.axhline(1, linestyle='--', color='gray')
                                  
   
    plt.tight_layout()
    plt.legend()
    plt.show()
    #plt.savefig('test.png')

def make_projY_d2pol(h2_hist_name, pm_user, Q2_user, model, plot_flag):
    # NEED TO FIX THIS FUNCTION, AS IT CURRENTLY DISPLAYS MULTIPLE SUBPLOTS, WHERE ONLY ONE IS NEEDED
    '''
    Brief: generic function makes 1D projections along y-axis (slicing xbins) for selected 2D histos,
    
    h2_hist_name: histogram base name
    pm_user: central missing momentu setting
    Q2_user: central Q2 setting
    model: "pwia" or "fsi"
    plot_flag: "2d" or "proj" or "proj_err"
    '''

    #histos_file_path = 'yield_estimates/d2_fsi/histogram_data/pm%d_thrq%d_%s/'%(pm_user, thrq_user, model)

    rel_err_thrs = 0.3 # mask >30 % relative error
   
    ifig = 1 # counter for 2d histogram figures

    
    # loop over central pm setting
    for ipm in pm_user:

        
        if plot_flag=='proj' or plot_flag=='proj_err':
            # set figure subplots for the 1d projections
            fig, ax = plt.subplots(6, 3, sharex='col', sharey='row', tight_layout=True)
            fig.text(0.5, 0.002, 'missing momentum, p$_{m}$ [GeV/c]', ha='center', fontsize=12)
            fig.text(0.005, 0.5, 'counts', va='center', rotation='vertical', fontsize=12)
            subplot_title =  r"p$_{m}$ vs. $\theta_{rq}$ 1d projection (%s), central p$_{m}$ setting: %d MeV"%(model, ipm)
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
            hist_file_path = 'yield_estimates/d2_pol/smallFSI/phi_180deg/histogram_data/pm%d_Q2_%.1f_%s/%s'%(ipm, jq2, model, h2_hist_basename)


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
            counts = np.sum(df.zcont)

      

            if plot_flag=='2d':
                
                # plotting the 2d histo
                fig2d = plt.figure(ifig, figsize=(9,6))
                plt.hist2d(df.x0 ,df.y0, weights=df.zcont, bins=(nxbins, nybins), range=[ [min(df.xlow), max(df.xup)], [min(df.ylow), max(df.yup)]], cmap = 'viridis', norm=mcolors.LogNorm())
                plt.xlabel(xlabel, fontsize=12)
                plt.ylabel(ylabel, fontsize=12)
                plt.title(title,   fontsize=14)
                
                plt.text(0.6*(df.x0[df.y0==df.y0[0]]).max(), 0.7*(df.y0[df.x0==df.x0[0]]).max(), r"Q$^{2}$=%.1f GeV$^{2}$"%(jq2)+"\n"+"$P_{m}$=%d MeV"%(ipm)+"\n"+"(counts = %d)"%(counts), fontsize=12)
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
                        ax = plt.subplot(6, 3, jdx+1)
                        ax.hist(ybc, bins=len(ybc), weights=count_per_xbin, range=[min(df.ylow), max(df.yup)], alpha=0.5, ec='k', density=False, label=r'%d counts (%.1f GeV$^{2}$)'%(cnts, jq2))

                        plt.title(r'$\theta_{rq}$ = %d $\pm$ %d deg'%(xbin, xbinw/2.))
                        plt.legend(frameon=False, loc='upper right')
                        jdx = jdx+1
                        
                    # ---- relative error plots ----
                    if plot_flag=='proj_err':

                        ax = plt.subplot(6, 3, jdx+1)
                        count_per_xbin_rel_err_m = ma.masked_where(count_per_xbin_rel_err>rel_err_thrs, count_per_xbin_rel_err)

                        y_const_m = ma.masked_where(count_per_xbin_rel_err>rel_err_thrs, y_const)
                        ybc_m = ma.masked_where(count_per_xbin_rel_err>rel_err_thrs, ybc)

                        ybc_m = ybc_m + offset
                        ax.errorbar(ybc_m, y_const_m, count_per_xbin_rel_err_m, marker='o', markersize=4, linestyle='None', label=r'%.1f GeV$^{2}$'%(jq2)) #//, label=r'%d counts'%(cnts))
                       
                        plt.title(r'$\theta_{rq}$ = %d $\pm$ %d deg'%(xbin, xbinw/2.))
                        jdx = jdx+1
            plt.legend(frameon=False, loc='upper right')

                  
       
    if plot_flag=='proj' or plot_flag=='proj_err':
        plt.show()
            
# call functions here (can later be passed thru steering code)
# make_plots(800, 79, 'pwia')
# make_1d_Xprojections('H_Pm_vs_thrq_yield', 800, 79, 'pwia')
# make_ratios_d2fsi([800], [28, 49, 55, 60, 66, 72, 79], 'ratio')
# overlay_d2fsi([800], [28, 49, 55], 'thrq', 'fsi')



# ---- example of d2 polarized plotting scripts -----
#overlay_d2pol([300], [3.5, 4.0, 4.5], 'Pm', 'fsi')
#make_projY_d2pol('Pm_vs_thrq', [300], [3.5, 4.0, 4.5], 'fsi', '2d')
#make_projY_d2pol('Pm_vs_thrq', [300], [3.5], 'fsi', 'proj')
#make_projY_d2pol('Pm_vs_thrq', [300], [3.5], 'fsi', 'proj_err')
