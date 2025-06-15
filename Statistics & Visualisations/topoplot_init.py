#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 21 16:36:21 2022

@author: Daniel
"""
#from mne import channels
import mne
import matplotlib.pyplot as plt
from scipy.io import loadmat

''' to identify where data used in components is moslty coming from in the scalp'''

# extract channels from raw data and define info/montage
path_import = '/Volumes/LaCie Armaz/DPrataLab/ISC/Raw_data/EEG/bsc049isc.vhdr'; # any random subject will do
#path_import = '/Volumes/disk/Work/Raw_Data/EEG/bsc049isc.vhdr'; # any random subject will do
raw = mne.io.read_raw_brainvision(path_import, preload = True, verbose = False);
ch_names = raw.ch_names[:63]
mne.add_reference_channels(raw, ref_channels=["FCz"], copy=False)
mne.set_eeg_reference(raw, ref_channels='average', copy=False)



info = mne.create_info(ch_names, 250, ['eeg']*63) #info structure; removes EMG channel here
info.set_montage('standard_1005')


# load data from all videos

path_import_project = '/Volumes/LaCie Armaz/DPrataLab/ISC/FINAL_PREPROCESSING/ISC_CALC/';
#path_import_project = '/Volumes/disk/Work/FINAL_PREPROCESSING/ISC_CALC/';



scalp_drug1 = {}
scalp_drug2 = {}

for video_num in range(1,21):
    
    #drug1
    scalp_drug1["scalp{}".format(video_num)] = loadmat(path_import_project + 'DRUG1/VIDEO' + str(video_num) + '_scalp_projection')['A_drug1'].T
    
    #drug2
    scalp_drug2["scalp{}".format(video_num)] = loadmat(path_import_project + 'DRUG2/VIDEO' + str(video_num) + '_scalp_projection')['A_drug2'].T

xlabels = ['Neutral', 'Neutral', 'Erotic', 'Erotic', 'Social-', 'Social-', 'Social+', 'Social+', 'Horror', 'Horror', 'Neutral', 'Neutral', 'Social-', 'Social-', 'Erotic', 'Erotic', 'Horror', 'Horror', 'Social+', 'Social+']

for video_num in range(1,21):
    
    fig, axs = plt.subplots(nrows = 2, ncols = 3, figsize=(9,4))

    fig.suptitle('Video{}'.format(video_num) + ': ' + str(xlabels[video_num-1]))
    
    for c in range(3):
        
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
        #clb.ax.set_title(unit_label,fontsize=fontsize) # title on top of colorbar

        # Drug info
        axs[0][0].set_ylabel('Drug1', fontsize = 12.5)
        axs[1][0].set_ylabel('Drug2', fontsize = 12.5)
            
    fig.tight_layout()

'----------------'

#ten_twenty_montage = mne.channels.make_standard_montage('standard_1005')

#fig = ten_twenty_montage.plot(kind='3d')
#fig.gca().view_init(azim=70, elev=15)  # set view angle
#ten_twenty_montage.plot(kind='topomap', show_names=False)    



'''ten_twenty_montage = mne.channels.make_standard_montage('standard_1020')
ten_twenty_montage.plot(show_names=False)
fig = ten_twenty_montage.plot(kind='3d')
fig.gca().view_init(azim=70, elev=15)


fig = plt.figure()
ax2d = fig.add_subplot(121)
ax3d = fig.add_subplot(122, projection='3d')
raw.plot_sensors(ch_type='eeg', axes=ax2d)
raw.plot_sensors(ch_type='eeg', axes=ax3d, kind='3d')
ax3d.view_init(azim=70, elev=15)'''
