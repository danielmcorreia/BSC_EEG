#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr  1 19:40:48 2022

@author: Daniel
"""

import mne
from scipy.io import loadmat
import numpy as np
import matplotlib.pyplot as plt

data_folder_path = '/Volumes/LaCie Armaz/DPrataLab/ISC/FINAL_PREPROCESSING/ISC_CALC/'

xlabels = ['Neutral', 'Neutral', 'High +', 'High +', 'Low -', 'Low -', 'Low +', 'Low +', 'High -', 'High -', 'Neutral', 'Neutral', 'Low -', 'Low -', 'High +', 'High +', 'High -', 'High -', 'Low +', 'Low +']

NComps = 3

scalp_drug1 = {}
scalp_drug2 = {}

# loads data
for video_num in range(1,21):
    
    scalp_drug1["scalp{}".format(video_num)] = \
        (loadmat(data_folder_path + 'DRUG1/VIDEO' + str(video_num) + '_scalp_projection')['A_drug1'].T)[0:3]
    
    scalp_drug2["scalp{}".format(video_num)] = \
        (loadmat(data_folder_path + 'DRUG2/VIDEO' + str(video_num) + '_scalp_projection')['A_drug2'].T)[0:3]
        
'-----------------------------------------------------------------------------'

# multiplies vectores by +1 or -1 for visualization purposes
for video_num in range(1,21):
    for c in range(3):
        if scalp_drug1['scalp{}'.format(video_num)][c][0:6].mean() < 0:
            scalp_drug1['scalp{}'.format(video_num)][c] = scalp_drug1['scalp{}'.format(video_num)][c]*-1
        if scalp_drug2['scalp{}'.format(video_num)][c][0:6].mean() < 0:
            scalp_drug2['scalp{}'.format(video_num)][c] = scalp_drug2['scalp{}'.format(video_num)][c]*-1
    
# throubleshoots specific component for specific movie
scalp_drug1['scalp2'][1] = scalp_drug1['scalp2'][1]*-1

'-----------------------------------------------------------------------------'

# calculates averages/movies of topoplots

topo_neutral1 = np.mean(np.stack((scalp_drug1['scalp1'], scalp_drug1['scalp2'], scalp_drug1['scalp11'], scalp_drug1['scalp12']), axis=2), axis=2)
topo_HighPos1 = np.mean(np.stack((scalp_drug1['scalp3'], scalp_drug1['scalp4'], scalp_drug1['scalp15'], scalp_drug1['scalp16']), axis=2), axis=2)
topo_HighNeg1 = np.mean(np.stack((scalp_drug1['scalp9'], scalp_drug1['scalp10'], scalp_drug1['scalp17'], scalp_drug1['scalp18']), axis=2), axis=2)
topo_LowPos1 = np.mean(np.stack((scalp_drug1['scalp7'], scalp_drug1['scalp8'], scalp_drug1['scalp19'], scalp_drug1['scalp20']), axis=2), axis=2)
topo_LowNeg1 = np.mean(np.stack((scalp_drug1['scalp5'], scalp_drug1['scalp6'], scalp_drug1['scalp13'], scalp_drug1['scalp14']), axis=2), axis=2)
    
topo_neutral2 = np.mean(np.stack((scalp_drug2['scalp1'], scalp_drug2['scalp2'], scalp_drug2['scalp11'], scalp_drug2['scalp12']), axis=2), axis=2)
topo_HighPos2 = np.mean(np.stack((scalp_drug2['scalp3'], scalp_drug2['scalp4'], scalp_drug2['scalp15'], scalp_drug2['scalp16']), axis=2), axis=2)
topo_HighNeg2 = np.mean(np.stack((scalp_drug2['scalp9'], scalp_drug2['scalp10'], scalp_drug2['scalp17'], scalp_drug2['scalp18']), axis=2), axis=2)
topo_LowPos2 = np.mean(np.stack((scalp_drug2['scalp7'], scalp_drug2['scalp8'], scalp_drug2['scalp19'], scalp_drug2['scalp20']), axis=2), axis=2)
topo_LowNeg2 = np.mean(np.stack((scalp_drug2['scalp5'], scalp_drug2['scalp6'], scalp_drug2['scalp13'], scalp_drug2['scalp14']), axis=2), axis=2)

'-----------------------------------------------------------------------------'

def get_info():
    # extract channels from raw data and define info/montage
    path_import = '/Volumes/LaCie Armaz/DPrataLab/ISC/Raw_data/EEG/bsc049isc.vhdr'; # any random subject will do
    #path_import = '/Volumes/disk/Work/Raw_Data/EEG/bsc049isc.vhdr';
    raw = mne.io.read_raw_brainvision(path_import, preload = True, verbose = False);
    ch_names = raw.ch_names[:63]
    mne.add_reference_channels(raw, ref_channels=["FCz"], copy=False)
    mne.set_eeg_reference(raw, ref_channels='average', copy=False)
    
    info = mne.create_info(ch_names, 250, ['eeg']*63) #info structure; removes EMG channel here
    info.set_montage('standard_1005')
    
    return info

'-----------------------------------------------------------------------------'

def get_correlations():
    
    labels = ["Neutral", "High Positive", "High Negative", "Low Positive", "Low Negative"]
    
    # correlations btw drugs1/2    
    print("\n\nCorrelations between Drugs")
    
    for c in range(3):
        
        print("\nComponent{}".format(c+1) + '\n')
        
        print('Neutral:')
        print(np.corrcoef(topo_neutral1[c], topo_neutral2[c])[0][1])
        print('HighPos:')
        print(np.corrcoef(topo_HighPos1[c], topo_HighPos2[c])[0][1])
        print('HighNeg:')
        print(np.corrcoef(topo_HighNeg1[c], topo_HighNeg2[c])[0][1])
        print('LowPos:')
        print(np.corrcoef(topo_LowPos1[c], topo_LowPos2[c])[0][1])
        print('LowNeg:')
        print(np.corrcoef(topo_LowNeg1[c], topo_LowNeg2[c])[0][1])
    
    print("\n\nCorrelations between Movies, per component")
    
    # pairwise correlations between movies, per component
    for c in range(3):
        print("\nComponent{}".format(c+1) + '\n')
        
        print("\n Drug1 \n")
        
        for video11, i in zip([topo_neutral1,topo_HighPos1,topo_HighNeg1,topo_LowPos1,topo_LowNeg1], range(5)):
            for video12, j in zip([topo_neutral1,topo_HighPos1,topo_HighNeg1,topo_LowPos1,topo_LowNeg1], range(5)):
                print(labels[i] + ' --- ' + labels[j])
                print(np.corrcoef(video11[c],video12[c])[0][1])
                
        print("\n Drug2 \n")
                
        for video21, i in zip([topo_neutral2,topo_HighPos2,topo_HighNeg2,topo_LowPos2,topo_LowNeg2], range(5)):
            for video22, j in zip([topo_neutral2,topo_HighPos2,topo_HighNeg2,topo_LowPos2,topo_LowNeg2], range(5)):
                print(labels[i] + ' --- ' + labels[j])
                print(np.corrcoef(video21[c],video22[c])[0][1])

'-----------------------------------------------------------------------------'
# TOPOPLOTS

def topoplot(scalp_drug1, scalp_drug2):
    
    info = get_info()

    for video_num in range(1,21):
    
        fig, axs = plt.subplots(nrows = 2, ncols = 3, figsize=(9,4))

        fig.suptitle('Video{}'.format(video_num) + ': ' + str(xlabels[video_num-1]))
    
        for c in range(3):
            
            # fix topoplots
            if scalp_drug1['scalp{}'.format(video_num)][c][0:6].mean() < 0:
                scalp_drug1['scalp{}'.format(video_num)][c] = scalp_drug1['scalp{}'.format(video_num)][c]*-1
            if scalp_drug2['scalp{}'.format(video_num)][c][0:6].mean() < 0:
                scalp_drug2['scalp{}'.format(video_num)][c] = scalp_drug2['scalp{}'.format(video_num)][c]*-1
            
            scalp_drug1['scalp2'][1] = scalp_drug1['scalp2'][1]*-1
            
            # plot topo maps
            im,cm = mne.viz.plot_topomap(scalp_drug1['scalp{}'.format(video_num)][c], info, vmin=-0.5, vmax=0.5, axes = axs[0][c])
            im,cm = mne.viz.plot_topomap(scalp_drug2['scalp{}'.format(video_num)][c], info, vmin=-0.5, vmax=0.5, axes = axs[1][c])
        
            # component info
            axs[0][c].set_title('Component {}'.format(c+1))

            # manually fiddle the position of colorbar
            ax_x_start = 0.93
            ax_x_width = 0.02
            ax_y_start = 0.25
            ax_y_height = 0.5
            cbar_ax = fig.add_axes([ax_x_start, ax_y_start, ax_x_width, ax_y_height])
            clb = fig.colorbar(im, cax=cbar_ax)
            clb.ax.set_title("Weight (A.U.)", fontsize = 8) # title on top of colorbar

            # Drug info
            axs[0][0].set_ylabel('OT', fontsize = 12.5)
            axs[1][0].set_ylabel('PLC', fontsize = 12.5)
            
        #plt.savefig("/Volumes/LaCie Armaz/DPrataLab/ISC/FINAL_PREPROCESSING/RESULTS/TOPOPLOTS/Figure{}.pdf".format(video_num), format="pdf", bbox_inches="tight")
            
        fig.tight_layout()

'-----------------------------------------------------------------------------'
# AVERAGED TOPOPLOTS

def averaged_topoplot():
    
    info = get_info()
    
    labels = ["Neutral", "High +", "High -", "Low +", "Low -"]
    i = 0
    for video_cat1, video_cat2, cat_num in zip([topo_neutral1,topo_HighPos1,topo_HighNeg1,topo_LowPos1,topo_LowNeg1], [topo_neutral2,topo_HighPos2,topo_HighNeg2,topo_LowPos2,topo_LowNeg2], np.arange(1,6)):
    
        fig, axs = plt.subplots(nrows = 2, ncols = 3, figsize=(9,4))

        fig.suptitle(str(labels[cat_num-1]), weight='bold')
    
        for c in range(3):
            
            # plot topo maps
            im,cm = mne.viz.plot_topomap(video_cat1[c], info, vmin=-0.5, vmax=0.5, axes = axs[0][c])
            im,cm = mne.viz.plot_topomap(video_cat2[c], info, vmin=-0.5, vmax=0.5, axes = axs[1][c])
        
            # component info
            axs[0][c].set_title('Component {}'.format(c+1))

            # manually fiddle the position of colorbar
            ax_x_start = 0.93
            ax_x_width = 0.02
            ax_y_start = 0.25
            ax_y_height = 0.5
            cbar_ax = fig.add_axes([ax_x_start, ax_y_start, ax_x_width, ax_y_height])
            clb = fig.colorbar(im, cax=cbar_ax)
            clb.ax.set_title("Weight (A.U.)", fontsize = 8) # title on top of colorbar

            # Drug info
            axs[0][0].set_ylabel('OT', fontsize = 12.5)
            axs[1][0].set_ylabel('PLC', fontsize = 12.5)
            
        
            
        fig.tight_layout()
        
        plt.savefig("/Volumes/LaCie Armaz/DPrataLab/ISC/FINAL_PREPROCESSING/RESULTS_out/TOPOPLOTS_AVG/{}.pdf".format(labels[i]), format="pdf", bbox_inches="tight")
        i +=1
            
#topoplot(scalp_drug1, scalp_drug2)
averaged_topoplot()
get_correlations()
