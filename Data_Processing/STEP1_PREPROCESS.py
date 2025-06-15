#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov  5 15:20:18 2021

@author: Daniel
"""

import os
import numpy as np
import mne
import matplotlib.pyplot as plt
from DEF_PREPROCESSSING import run_ICA, remove_bad_channels, remove_outliers
from SEGMENTATION import segment
mne.set_log_level("WARNING")

path_import = '/Volumes/LaCie Armaz/DPrataLab/ISC/Raw_data/EEG/'

path_save_ICA = '/Volumes/LaCie Armaz/DPrataLab/ISC/FINAL_PREPROCESSING/CHECKPOINTS/ICA_checkpoint/'
path_save_segment = '/Volumes/LaCie Armaz/DPrataLab/ISC/FINAL_PREPROCESSING/CHECKPOINTS/segment_checkpoint/'
path_save_bach_ch = '/Volumes/LaCie Armaz/DPrataLab/ISC/FINAL_PREPROCESSING/CHECKPOINTS/bad_ch_checkpoint/'
path_save_final_preprocessing = '/Volumes/LaCie Armaz/DPrataLab/ISC/FINAL_PREPROCESSING/CHECKPOINTS/final_preprocessing_checkpoint/'
path_save_final_preprocessing_out = '/Volumes/LaCie Armaz/DPrataLab/ISC/FINAL_PREPROCESSING/CHECKPOINTS/final_preprocessing_checkpoint_out/'


checkpoint_fail = []

#subjects_list = [file for file in os.listdir(path_import) if (file.endswith('.vhdr') and not file.startswith("."))]
#subjects_list = [file for file in os.listdir(path_save_ICA) if (file.endswith('.fif') and not file.startswith("."))]
subjects_list = [file for file in os.listdir(path_save_ICA) if (file.endswith('.fif') and not file.startswith("."))]


rejected_datapoints = []

for subject in range(43,len(subjects_list)): #len(range(subjects_list))
    
    
    print('\n\nSUBJECT#{}_READING_DATA--------------------------------------------'.format(subject))

    # reads & plots data for current subject
    #current_subject = mne.io.read_raw_brainvision(os.path.join(path_import, subjects_list[subject]), preload = True, verbose = False);
    #current_subject = mne.io.read_raw_fif(path_save_ICA + subjects_list[subject], preload = True, verbose = None)
    
    rejected_datapoints.append([subjects_list[subject]])
    
    '''print('\n\nSUBJECT#{}_HIGH-PASS_FILTERING-------------------------------------'.format(subject))

    # high-pass data at 1Hz
    current_subject = current_subject.filter(l_freq=1, h_freq=None, picks = 'eeg', verbose=False);
    
    
    print('\n\nSUBJECT#{}_NOTCH_FILTERING-----------------------------------------'.format(subject))
    
    # Notch Filter 50Hz harmonic to remove powerline interference
    current_subject = current_subject.notch_filter(np.arange(50, 499, 50), picks = 'eeg', verbose=False);
    
    
    print('\n\nSUBJECT#{}_RE-REFERENCING-------------------------------------------'.format(subject))
    
    # add "fake" referencing channel to mantain rank, since ours is missing
    mne.add_reference_channels(current_subject, ref_channels=["FCz"], copy=False);
    
    # re-reference to average
    mne.set_eeg_reference(current_subject, ref_channels='average', copy=False);
    
    
    print('\n\nSUBJECT#{}_APPLYING_ICA--------------------------------------------'.format(subject))
    
    current_subject = run_ICA(current_subject)

    # saves preprocessed data ICA (1) to allow further modification in subsequent steps
    current_subject.save(path_save_ICA+subjects_list[subject][:-5]+".fif", overwrite=True, fmt='single', picks='all');
    
    # ATTENTION: Manual modification of bad channels is injected here; run code from here if there's a need to modify them manually
    
    current_subject.plot()
    plt.pause(3)
    
    print('\n\nSUBJECT#{}_SEGMENTING----------------------------------------------'.format(subject))
    
    signal_triggers, timepoints_video_begin, check_fail = segment(current_subject)
    
    if check_fail:
        checkpoint_fail.append(subject) #these subjects' segmentations should be dealth wit manually
        continue'''
    
    for video_num in range(20):
        

        print('\nVIDEO{}_SEGMENTING'.format(video_num+1))
        
        # segments video; len of videos is 40 sec so i'm taking sfreq*40sec
        '''video_trace = signal_triggers[timepoints_video_begin[video_num]:timepoints_video_begin[video_num]+(1000*40), 0:63]
        
        #creates indiv structure
        ch_names = current_subject.ch_names[:63]
        info = mne.create_info(ch_names, 1000, ['eeg']*63) #info structure; removes EMG channel here
        info.set_montage('standard_1005')
        data = mne.io.RawArray(video_trace.T, info)
        
        
        print('\nVIDEO{}_DOWNSAMPLING'.format(video_num+1))
        
        data = data.resample(250)
        
        # saves preprocessed data segmented (2) to allow further modification in subsequent steps
        #data.save(path_save_segment+"VIDEO"+str(video_num+1)+"/"+subjects_list[subject][:-4]+".fif", overwrite=True);

        
        print('\nVIDEO{}_ID_BAD_CHANNELS'.format(video_num+1))
        
        remove_channels_trace = data.get_data() # removes unecessary EMG channel
        info = mne.create_info(ch_names, 250, ['eeg']*63) # info needs updating
        data.info['bads'] = current_subject.info['bads'] #includes manually rejected channels, if any

        find_bad_channels = remove_bad_channels(remove_channels_trace, kIQDp = 4)
        
        bad_channels = find_bad_channels[0]
        
        # mark bad channels
        if bad_channels.any(): #if not empty
            bad_channels = bad_channels.flatten()
            data.info['bads'] = list(np.array(ch_names).take(bad_channels)) #just to plot
            
            
        for bad_channel in data.info['bads']:
            bad_channel = ch_names.index(bad_channel)
            remove_channels_trace[bad_channel] = [0]*len(remove_channels_trace[bad_channel]) # replaces channel by 0'''
        
        #data = mne.io.RawArray(remove_channels_trace, info);
        #data.save(path_save_bach_ch+"VIDEO"+str(video_num+1)+"/"+subjects_list[subject][:-4]+".fif", overwrite=True);
        
        
        current_subject = mne.io.read_raw_fif(path_save_bach_ch + 'VIDEO{}/'.format(video_num+1) + subjects_list[subject], preload = True, verbose = None)
        remove_channels_trace = current_subject.get_data()
        
        #creates indiv structure
        ch_names = current_subject.ch_names[:63]
        info = mne.create_info(ch_names, 250, ['eeg']*63) #info structure; removes EMG channel here
        info.set_montage('standard_1005')
        
        print('\nVIDEO{}_REMOVE_OUTLIERS'.format(video_num+1))
        
        find_outliers = remove_outliers(remove_channels_trace)
        
        correct_outliers_trace = find_outliers[0]
        
        
        print('\nVIDEO{}_PLOTTING & EXPORTING'.format(video_num+1))
        
        # saves final preprocessed (4) to allow further modification in subsequent steps
        data = mne.io.RawArray(correct_outliers_trace*10e5, info);
        data.save(path_save_final_preprocessing_out+"VIDEO"+str(video_num+1)+"/"+subjects_list[subject][:-4]+".fif", overwrite=True, fmt='double');
        

        print('\n\nSUBJECT#{}_PREPROCESS_COMPLETED--------------------------------'.format(subject))
    
        #data.plot();
        
        #plt.pause(15)'''
        
        #rejected_datapoints[-1].append((video_num, find_bad_channels[1], find_outliers[1]))
        rejected_datapoints[-1].append((video_num, find_outliers[1]))

    
    
    plt.close('all')
    
print('\n\nThe following subjects have insonsistencies in triggers and should be dealt with manually: ', checkpoint_fail )