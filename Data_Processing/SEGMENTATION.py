#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov  5 16:36:33 2021

@author: Daniel
"""

import mne
import numpy as np


def segment(subject):
    
    '''following function is designed to find the indexes where the desired TTL_code
    if found in the stimulus column of the array (see below)
    Requires: x - TTL code of interest; z - vector of stimulus id, with same length
    as eeg trace (either 0 or specific TTL code)
    Ensures: list of indexes of TTL code of interest in vector z'''
    
    get_indexes = lambda x, z: [i for (y, i) in zip(z, range(len(z))) if x == y]

    signal = subject.get_data().T

    
    print('\nExtracting trigger info...')
    
    #identifies events in data; gets array in format ( [ [timepoint, 0, event_id], [...], ... ], { TTL Code info } )
    events = mne.events_from_annotations(subject)

    event_TTL_code = []
    event_time = []
    #loops through array and gets lists of TTL codes and timestamps
    for event in events[0]:
        event_TTL_code.append(event[2])
        event_time.append(event[0])
    
    #print(event_time)
        
    print('\nAligning triggers with recordings...')
    
    #creates array of 0's to store trigger codes in respective timepoint
    data_triggers = [[0]] * signal.shape[0]

    for i, timepoint in enumerate(event_time):
        
        if not event_TTL_code[i] == 99999: #ignores 1st TTL code, not relevant
            current_code = event_TTL_code[i] #gets current timepoint associated TTL code
            data_triggers[timepoint-1] = [current_code] #replaces value 0 by trigger code in respective timepoint
    
    #appends column to data signal
    signal_triggers = np.append(signal, data_triggers, axis=1)
    
    #print(signal_triggers.shape)
    
    
    print('\nTrigger checkpoint...')
    
    if np.count_nonzero(signal_triggers[:,-1] == 31) == 20:
        check_fail = 0
        pass
        
    else:
        '''if current subject doesn't have the trigger scheme as expected (20 30's, 31'2 and 32's - corresponding to
        a normal video segmentation scheme), these subjects are passed and are not included in this initial analysis
        (for now! this should be dealth with manually, before moving to next step; ask Vasco what happened in these cases)'''
        
        check_fail = 1
    
    
    print('\nSegmenting timeseries...')
    
    #gets timepoints of video_begin and video_end
    timepoints_video_begin = get_indexes(31,signal_triggers[:,-1])   
    
    return signal_triggers, timepoints_video_begin, check_fail
