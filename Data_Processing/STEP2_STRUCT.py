#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Nov  6 01:38:52 2021

@author: Daniel
"""


'''to run this script, data in import folder path should be divided by videos, as outputted from previous step;
you should have a new folder named 'PREPROCESSED', with subfolders 'DRUG1' and 'DRUG2', each of them containing
20 folders, entitled VIDEO#, one for each of the 20 movies '''

import os
import mne
import scipy.io
import numpy as np
import shutil
import pandas as pd

path_import = '/Volumes/LaCie Armaz/DPrataLab/ISC/FINAL_PREPROCESSING/CHECKPOINTS/final_preprocessing_checkpoint_out/'

path_export_drug1 = '/Volumes/LaCie Armaz/DPrataLab/ISC/FINAL_PREPROCESSING/PREPROCESSED_out/DRUG1/'
path_export_drug2 = '/Volumes/LaCie Armaz/DPrataLab/ISC/FINAL_PREPROCESSING/PREPROCESSED_out/DRUG2/'

rejected_trials = pd.read_excel ('/Volumes/LaCie Armaz/DPrataLab/ISC/FINAL_PREPROCESSING/rejected_trials.xlsx')[['participant','per_x','Trial','Rejected_by_Anottation']]
rejected_trials = rejected_trials.loc[(rejected_trials['Rejected_by_Anottation'] == 1) | (rejected_trials['per_x'] >= 50)] # 11.2% of the data was rejected!

#subjects_list = [file for file in os.listdir(path_import) if (file.endswith('.fif'))]


# divide subjects by drug, manually

subjects_DRUG1 = ['BSC_002.fif', 'BSC_003.fif', 'BSC_006.fif', 'BSC_009.fif', 'BSC_011.fif', 'BSC_012.fif', 'BSC_017.fif',
                  'BSC_020.fif', 'BSC_021.fif', 'BSC_022.fif', 'BSC_027.fif', 'BSC_028.fif', 'BSC_031.fif', 'BSC_032.fif',
                  'BSC_037.fif', 'BSC_038.fif', 'BSC_039.fif', 'BSC_040.fif', 'BSC_043.fif', 'BSC_044.fif', 'BSC_047.fif',
                  'BSC_048.fif', 'BSC_057.fif', 'BSC_058.fif', 'BSC_059.fif', 'BSC_060.fif', 'BSC_061.fif', 'BSC_062.fif']

subjects_DRUG2 = ['BSC_007.fif', 'BSC_008.fif', 'BSC_013.fif', 'BSC_014.fif', 'BSC_015.fif', 'BSC_016.fif', 'BSC_023.fif',
                  'BSC_024.fif', 'BSC_025.fif', 'BSC_026.fif', 'BSC_029.fif', 'BSC_030.fif', 'BSC_033.fif', 'BSC_034.fif',
                  'BSC_035.fif', 'BSC_036.fif', 'BSC_041.fif', 'BSC_042.fif', 'BSC_045.fif', 'BSC_046.fif', 'BSC_049.fif',
                  'BSC_050.fif', 'BSC_055.fif', 'BSC_056.fif']


# copy data from individuals for each of the corresponding folders
    
for video_num in range(1,21):
    for subdir, dirs, files in os.walk(path_import+'VIDEO{}'.format(video_num)):
        for subject in files:
          
            reject_trials = rejected_trials[rejected_trials['participant'].isin([subject[:-4]])].get('Trial').to_list()
 
            if subject in subjects_DRUG1 and video_num not in reject_trials:
                shutil.copyfile(path_import + 'VIDEO{video}/'.format(video=video_num) + subject, path_export_drug1+'VIDEO{}/'.format(video_num) + subject)
            if subject in subjects_DRUG2 and video_num not in reject_trials:
                shutil.copyfile(path_import + 'VIDEO{video}/'.format(video=video_num) + subject, path_export_drug2+'VIDEO{}/'.format(video_num) + subject)
                
 
keep_track_rej_drug1 = []
keep_track_rej_drug2 = []     

# organizes data in format ready to be read by ISC scripts (either for Python or MATLAB scripts)
for drug in range(1,3):
    
    if drug == 1:
        subjects_list = subjects_DRUG1
        keep_track_rej = keep_track_rej_drug1
        
    if drug == 2:
        subjects_list = subjects_DRUG2
        keep_track_rej = keep_track_rej_drug2
    
    path_import = '/Volumes/LaCie Armaz/DPrataLab/ISC/FINAL_PREPROCESSING/PREPROCESSED_out/DRUG{}/'.format(drug)

    path_save_mat = '/Volumes/LaCie Armaz/DPrataLab/ISC/FINAL_PREPROCESSING/ISC_STRUCT_out/DRUG{}/'.format(drug)
    
    for video_num in range(20):
        
        keep_track_rej.append([])
        
        video_array_mat = []
        video_array_py = []
    
        for subject in subjects_list:
            
            # troubleshoots for cases when videos were removed;
            try:
                current_subject = mne.io.read_raw_fif(path_import + 'VIDEO'+str(video_num+1) + "/" + subject, preload = True, verbose = None)
                
                current_subject_array = current_subject.get_data()
                
                #for MATLAB code
                video_array_mat.append(current_subject_array.T)
        
                #for Python code
                #video_array_py.append(current_subject_array)
                
                keep_track_rej[-1].append(1)

            except:
                
                keep_track_rej[-1].append(0)
                pass

    
        #for MATLAB code:
        video_data_mat = np.array(np.stack(video_array_mat, axis = 2)) 
        dict_mat = {'X': video_data_mat, 'fs': 250}
        scipy.io.savemat(path_save_mat+"VIDEO"+str(video_num+1)+'.mat', mdict=dict_mat, appendmat=True, format='5', long_field_names=False)
    
        #for Python code
        #video_data_py = np.array(np.stack(video_array_py, axis = 0))
        #np.save(path_save_py+"VIDEO"+str(video_num+1)+".npy", video_data_py)'''

np.save("/Volumes/LaCie Armaz/DPrataLab/ISC/FINAL_PREPROCESSING/PREPROCESSED_out/"+"keep_track_rej_drug1.npy", keep_track_rej_drug1)
np.save("/Volumes/LaCie Armaz/DPrataLab/ISC/FINAL_PREPROCESSING/PREPROCESSED_out/"+"keep_track_rej_drug2.npy", keep_track_rej_drug2)

    