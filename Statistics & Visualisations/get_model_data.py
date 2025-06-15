#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 18 11:08:10 2022

@author: Daniel
"""

from scipy.io import loadmat
import pandas as pd
import numpy as np

subjects_DRUG1 = np.array(['BSC_002.fif', 'BSC_003.fif', 'BSC_006.fif', 'BSC_009.fif', 'BSC_011.fif', 'BSC_012.fif', 'BSC_017.fif',
                  'BSC_020.fif', 'BSC_021.fif', 'BSC_022.fif', 'BSC_027.fif', 'BSC_028.fif', 'BSC_031.fif', 'BSC_032.fif',
                  'BSC_037.fif', 'BSC_038.fif', 'BSC_039.fif', 'BSC_040.fif', 'BSC_043.fif', 'BSC_044.fif', 'BSC_047.fif',
                  'BSC_048.fif', 'BSC_057.fif', 'BSC_058.fif', 'BSC_059.fif', 'BSC_060.fif', 'BSC_061.fif', 'BSC_062.fif'])

subjects_DRUG2 = np.array(['BSC_007.fif', 'BSC_008.fif', 'BSC_013.fif', 'BSC_014.fif', 'BSC_015.fif', 'BSC_016.fif', 'BSC_023.fif',
                  'BSC_024.fif', 'BSC_025.fif', 'BSC_026.fif', 'BSC_029.fif', 'BSC_030.fif', 'BSC_033.fif', 'BSC_034.fif',
                  'BSC_035.fif', 'BSC_036.fif', 'BSC_041.fif', 'BSC_042.fif', 'BSC_045.fif', 'BSC_046.fif', 'BSC_049.fif',
                  'BSC_050.fif', 'BSC_055.fif', 'BSC_056.fif'])


path_import_drug1 = '/Volumes/LaCie Armaz/PAPER_REVISIONS_1902/ISC_CALC_out/keep_track_rej_drug1_afterISC.npy'
path_import_drug2 = '/Volumes/LaCie Armaz/PAPER_REVISIONS_1902/ISC_CALC_out/keep_track_rej_drug2_afterISC.npy'

rj_drug1 = np.load(path_import_drug1).astype(np.float32)
rj_drug2 = np.load(path_import_drug2).astype(np.float32)

xlabels = ['Neutral', 'Neutral', 'High +', 'High +', 'Low -', 'Low -', 'Low +', 'Low +', 'High -', 'High -', 'Neutral', 'Neutral', 'Low -', 'Low -', 'High +', 'High +', 'High -', 'High -', 'Low +', 'Low +']

data_folder_path = '/Volumes/LaCie Armaz/PAPER_REVISIONS_1902/ISC_CALC_out/'

ISC_persubject_sum_drug1 = {}
ISC_persubject_sum_drug2 = {}

for video_num in range(1,21):
    
    ISC_persubject_sum_drug1["ISC_persubject_sum{}".format(video_num)] = \
            loadmat(data_folder_path + 'DRUG1/VIDEO' + str(video_num) + '_ISC_persubject_sum')['ISC_persubject_sum_drug1'][0]
        
    ISC_persubject_sum_drug2["ISC_persubject_sum{}".format(video_num)] = \
            loadmat(data_folder_path + 'DRUG2/VIDEO' + str(video_num) + '_ISC_persubject_sum')['ISC_persubject_sum_drug2'][0]

dict_ISC = {1: ISC_persubject_sum_drug1, 2: ISC_persubject_sum_drug2}

#print(ISC_persubject_sum_drug1)

#print(len(ISC_persubject_sum_drug1))

#print(rj_drug1)

###LONG FORMAT DATA###

df_columns = ['Subject', 'ISC', 'Video', 'Drug']

# get list of subjects per video
#dict_subjects = {}
#for video in range(0,20):
    #dict_subjects[video] = {1: subjects_DRUG1[rj_drug1[video] == 2], 2: subjects_DRUG2[rj_drug2[video] == 2]}
    #dict_subjects['drug1'] = subjects_DRUG1[rj_drug1[video] == 2]
    #dict_subjects['drug2'] = subjects_DRUG2[rj_drug2[video] == 2]
    #list_subjects_DRUG1.append(subjects_DRUG1[rj_drug1[video] == 2])
    #list_subjects_DRUG2.append(subjects_DRUG2[rj_drug2[video] == 2])

#print(dict_subjects)

# get list of subjects per video
dict_ISC_DRUG1 = {}
dict_ISC_DRUG2 = {}
for video in range(0,20):
    #dict_subjects = {1: subjects_DRUG1[rj_drug1[video] == 2], 2: subjects_DRUG2[rj_drug2[video] == 2]}

    dict_ISC_DRUG1[video] = {'subjects': subjects_DRUG1[rj_drug1[video] == 2], 'drug': 1, 'video': video+1, 'ISC': dict_ISC[1]["ISC_persubject_sum{}".format(video+1)]}
    dict_ISC_DRUG2[video] = {'subjects': subjects_DRUG2[rj_drug2[video] == 2], 'drug': 2, 'video': video+1, 'ISC': dict_ISC[2]["ISC_persubject_sum{}".format(video+1)]}


df = (pd.DataFrame.from_dict(dict_ISC_DRUG1).T).append(pd.DataFrame.from_dict(dict_ISC_DRUG2).T)
df = df.explode(column=['subjects', 'ISC'])

df.to_excel("/Volumes/LaCie Armaz/PAPER_REVISIONS_1902/ISC_CALC_out/long_format_data.xlsx")



#print(dict_ISC_DRUG2)

#print(list_subjects_DRUG1)
#print(list_subjects_DRUG2)
#df_data['Video'] = pd.Series([1]*len(subjects_DRUG1) + [2]*len(subjects_DRUG2))



### SHORT FORMAT DATA CODE ###            

'''df_final = pd.DataFrame(subjects_DRUG1+subjects_DRUG2, columns = ['subjects'])

df_final['drug'] = pd.Series([1]*len(subjects_DRUG1) + [2]*len(subjects_DRUG2))

for video_num in range(1,21):
    i1 = 0
    i2 = 0
    for sub1 in range(len(subjects_DRUG1)):
        if rj_drug1[video_num-1][sub1] == 1:
            rj_drug1[video_num-1][sub1] = ISC_persubject_sum_drug1["ISC_persubject_sum{}".format(video_num)][sub1+i1]
        else:
            i1-=1
            
    for sub2 in range(len(subjects_DRUG2)):
        if rj_drug2[video_num-1][sub2] == 1:
            rj_drug2[video_num-1][sub2] = ISC_persubject_sum_drug2["ISC_persubject_sum{}".format(video_num)][sub2+i2]
        else:
            i2-=1
    


df=pd.DataFrame()

#replace 0 with NaN

for label, id_movie in zip(xlabels, range(1,21)):
    
    df[label + '_{}'.format(id_movie)] = pd.Series(np.concatenate((rj_drug1[id_movie-1], rj_drug2[id_movie-1])))
    
    #df[label + '_{}'.format(id_movie)] = pd.Series([1]*len(subjects_DRUG1) + [2]*len(subjects_DRUG2))

#df['Neutral'] = pd.DataFrame.from_dict(ISC_persubject_sum_drug1['ISC_persubject_sum1']+ISC_persubject_sum_drug2['ISC_persubject_sum1'])

df.replace(0.0, np.nan, inplace=True)

df_final["Neutral"] = df[["Neutral_1","Neutral_2","Neutral_11", "Neutral_12"]].mean(axis=1)
df_final["HighPos"] = df[["High +_3","High +_4","High +_15", "High +_16"]].mean(axis=1)
df_final["LowNeg"] = df[["Low -_5","Low -_6","Low -_13", "Low -_14"]].mean(axis=1)
df_final["LowPos"] = df[["Low +_7","Low +_8","Low +_19", "Low +_20"]].mean(axis=1)
df_final["HighNeg"] = df[["High -_9","High -_10","High -_17", "High -_18"]].mean(axis=1)

df_final.replace(np.nan, '', inplace=True)
#.fillna('')


#df_final.to_excel("/Volumes/LaCie Armaz/DPrataLab/ISC/FINAL_PREPROCESSING/MODELS_out/ANOVA_data.xlsx")'''