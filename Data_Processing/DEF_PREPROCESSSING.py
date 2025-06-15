#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov  5 14:54:40 2021

@author: Daniel
"""

import numpy as np
import matplotlib.pyplot as plt
from mne.preprocessing import ICA


def run_ICA(subject):
    
    '''Fits ICA to subject; waits for manual input of components to remove (should befed manually);
    removes chosen components from data; returns subject's signal, eye-blink & eye-movement free
    
    Requires:
        mne.io.brainvision.brainvision.RawBrainVision data format of a subject's data toanalyse
    Ensures:
        mne.io.brainvision.brainvision.RawBrainVision data format of subject's signal, eye-blink & eye-movement free'''
    
    ica = ICA(max_iter='auto', random_state=23, method = 'picard') #allow_ref_meg False
    ica.fit(subject, decim = 3) #decim to speed up things #reject to reject based on amplitude
    ica.plot_sources(subject)
    ica.plot_components(inst = subject)
    
    # some time to interact with components and decide which components to remove
    plt.pause(30)
    
    # asks for input to exclude eye related components
    components_exclude = []
    component_num = input('Insert component index to exclude (0-63): ')
    
    # continues asking until no input is given
    while component_num:
        components_exclude.append(int(component_num))
        component_num = input('Insert component index to exclude (0-63): ')
    
    ica.plot_overlay(subject, exclude=components_exclude, picks='eeg')
    
    # ica.apply() changes the Raw object in-place, so let's make a copy first:
    reconst_subject = subject.copy()
    ica.apply(reconst_subject, exclude=components_exclude)
    
    # plots original and ICA-removed
    subject.plot(title = 'original timeseries')
    #plt.pause(5)
    reconst_subject.plot(title = 'ICA-removed timeseries')
    plt.pause(5)
    #plt.close()
    
    return reconst_subject


def remove_bad_channels(signal, kIQDp = 4):
            
    '''Channels whose average power exceeded mean channel power by kIQDp stds were excluded, according to Dmochowski et al. (2014, 2018); unless manually selected
    
    Requires
        vector signal to correct; kIQDp (predefined 3) multiple of interquartile differences to mark as outliers channels
    Ensures:
        signal w/ found bad channels replaced by 0'''
    
    # Find bad channels based on power outliers; plots data for additional hand selection
    # according to Dmochoski et al., 2012, 2014, 2018
    logpower = np.log(np.std(signal, axis=1))
    Q=np.percentile(np.log(np.std(signal, axis=1)),[25, 50, 75])
    badchannels = np.argwhere(logpower-Q[1]>kIQDp*(Q[2]-Q[0]))

    return badchannels


def remove_outliers(signal, thresh_stds = 4, iters = 4):
    
    #!# Very inneficient. Make sure to re-write this func. #!#
    
    '''Removes outliers above defined amplitude threshold of thresh_stds * std signal and 100ms before/after, according to Dmochowski et al. (2012)
    Requires:
        vector signal to correct; threshold_stst (predefined 4) and num of iterations (predefined 4)
    Ensures:
        cleaned signal, outlier-free'''
        
    #done by Zhu et al.
    
    for j in range(iters): 
        for i, channel in enumerate(signal):
            
            mean = np.mean(channel)
            std = np.std(channel)
    
            final_list = np.array([x if abs(x - mean) < thresh_stds * std else np.nan for x in channel])
            
            try: nans = np.stack(np.argwhere(np.isnan(final_list)), axis=1)[0] # replaces nans by 0s
            except: continue
        
        
            #removes 100ms before and after events = 25 sampling points #according to Dmochoski et al., 2014
            for ind in nans:
                final_list[ind-25:ind]=np.nan
                try: final_list[ind:ind+25]=np.nan #cases of end vector
                except: continue
                
            final_list = np.nan_to_num(final_list) #converts nans to 0
        
            signal[i] = final_list
        
            #print('channel {}, round {} check'.format(i+1,j+1))
            
    return signal



def run_rPCA(signal):
    
    '''implemented in MATLAB; see '''
    
    return None