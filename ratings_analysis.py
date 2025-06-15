#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec  6 12:24:43 2021

@author: danielcorreia
"""

import os
import numpy as np
from scipy.io import loadmat
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib import rcParams

# LOADS PERSUBJECT_SUM and PERSECOND DATA

data_folder_ISC = '/Volumes/LaCie Armaz/DPrataLab/ISC/ISC_CALC/MATLAB/'

ISC_persecond_drug1 = []
ISC_persubject_sum_drug1 = []

ISC_persecond_drug2 = []
ISC_persubject_sum_drug2 = []


for video_num in range(1,21):
    
    #drug1
    #ISC_persecond_drug1["ISC_persecond{}".format(video_num)] = loadmat(data_folder_ISC + 'DRUG1/video' + str(video_num) + '_ISC_persecond')
    ISC_persubject_sum_drug1.append(loadmat(data_folder_ISC + 'DRUG1/video' + str(video_num) + '_ISC_persubject_sum')['ISC_persubject_sum'][0].tolist())
    
    #drug2
    #ISC_persecond_drug2["ISC_persecond{}".format(video_num)] = loadmat(data_folder_ISC + 'DRUG2/video' + str(video_num) + '_ISC_persecond')
    ISC_persubject_sum_drug2.append(loadmat(data_folder_ISC + 'DRUG2/video' + str(video_num) + '_ISC_persubject_sum')['ISC_persubject_sum'][0].tolist())


# order from subjects here is mantained from ISC calculation, so it's secure to associate ratings to ISC from prev exported arrays

subjects_DRUG1 = ['BSC_002_IA', 'BSC_003_IA', 'BSC_006_IA', 'BSC_009_IA', 'BSC_011_IA', 'BSC_012_IA',
                  'BSC_017_IA', 'BSC020ISC', 'BSC021ISC', 'BSC_022_ISC', 'BSC_027_ISC', 'bsc028isc', 'bsc031isc',
                  'BSC_032_ISC', 'BSC_037_ISC', 'bsc038isc', 'BSC_039_ISC', 'BSC_040_ISC', 'bsc043isc', 'BSC_044_ISC',
                  'bsc047isc', 'bsc048isc', 'bsc057isc', 'bsc058isc', 'bsc059isc', 'bsc060isc', 'bsc061isc', 'bsc062isc']

subjects_DRUG2 = ['BSC_007_IA', 'BSC_008_IA', 'BSC_013_IA', 'BSC_014_IA', 'BSC_015_IA', 'BSC_016_IA', 'BSC_023_ISC',
                  'BSC_024_ISC', 'BSC025ISC', 'BSC026ISC', 'bsc029isc', 'bsc030isc', 'BSC_033_ISC', 'bsc034isc',
                  'bsc035isc', 'BSC_036_ISC', 'bsc041isc', 'bsc042isc', 'BSC_045_ISC', 'bsc046isc', 'bsc049isc',
                  'bsc050isc', 'bsc055isc', 'BSC_056_ISC']

# rejected according to annotations taken during data acquisition
rejected_ratings = ['BSC_015_IA', 'bsc038isc', 'bsc043isc', 'bsc058isc']

# data does not exist
missing_ratings = ['bsc062isc']


# Movie IDs

movie_categories = ['Neutral 5001', 'Neutral 5002', 'Neutral 5000', 'Neutral 5008',
                    'Erotic 2002', 'Erotic 2000', 'Erotic 2009', 'Erotic 2003',
                    'Horror 1002', 'Horror 1001', 'Horror 1007', 'Horror 1006',
                    'Social- 3004', 'Social- 3003', 'Social- 3002', 'Social- 3007',
                    'Social+ 4004', 'Social+ 4006', 'Social+ 4001', 'Social+ 4009']

movie_categories_idx = [[0,1,10,11], [2,3,14,15], [8,9,16,17], [5,6,13,14], [6,7,18,19]]

'-----------------------------------------------------------------------------'

# SHORT RATINGS ANALYSIS (Arousal and Valence)

# loads short ratings data

#path_import = '/Users/danielcorreia/Desktop/short_ratings/'

path_import = '/Volumes/LaCie Armaz/DPrataLab/ISC/short_ratings/'

# rejected movies for specific subjects --> this was done mannualy, in the raw data (by replacing by 0s)

# prepares data

arousal_ratings_drug1 = []
arousal_ratings_drug2 = []

valence_ratings_drug1 = []
valence_ratings_drug2 = []

for drug in range(1,3):
    
    if drug == 1:
        subjects_list = subjects_DRUG1
        arousal_ratings = arousal_ratings_drug1
        valence_ratings = valence_ratings_drug1
        
    if drug == 2:
        subjects_list = subjects_DRUG2
        arousal_ratings = arousal_ratings_drug2
        valence_ratings = valence_ratings_drug2
        
    for subject in subjects_list:
        
        if subject in subjects_list and subject not in rejected_ratings and subject not in missing_ratings:
            
            # reads lines from current txt
            arousal_lines = [x.split('\t') for x in open(path_import + subject + '/RESULTS_AROUSAL_VA.txt').readlines()]
            valence_lines = [x.split('\t') for x in open(path_import + subject + '/RESULTS_VALENCE_VA.txt').readlines()]

            
            # gets rating for each line/movie in current subject
            # for some reason, there are 10-marked ratings; must have been misuderstanding by participants; those are replaced by 9s here
            arousal_ratings_subject = [int(rating[5][:-1]) if int(rating[5][:-1]) != 10 else 9 for rating in arousal_lines[1::]]
            valence_ratings_subject = [int(rating[5][:-1]) if int(rating[5][:-1]) != 10 else 9 for rating in valence_lines[1::]]

        if subject in missing_ratings or subject in rejected_ratings:
            arousal_ratings_subject = [np.nan]*20
            valence_ratings_subject = [np.nan]*20
            pass
    
        arousal_ratings.append(arousal_ratings_subject)
        valence_ratings.append(valence_ratings_subject)
    
arousal_ratings_drug1 = (np.array(arousal_ratings_drug1).T).tolist()
arousal_ratings_drug2 = (np.array(arousal_ratings_drug2).T).tolist()

valence_ratings_drug1 = (np.array(valence_ratings_drug1).T).tolist()
valence_ratings_drug2 = (np.array(valence_ratings_drug2).T).tolist()


# i manupulated data manually for cases when secific movies where rejected in specific subjects, replacing rating by 0. To match ISC data, let's remove 0s
#can't use np.array here coz of != lenghts that come from ISC matrices; gotta do this inneficiently w/ lists
for movie_arousal_drug1, movie_arousal_drug2, movie_valence_drug1, movie_valence_drug2 in zip(arousal_ratings_drug1, arousal_ratings_drug2, valence_ratings_drug1, valence_ratings_drug2):

    while 0 in movie_arousal_drug1: movie_arousal_drug1.remove(0)
    while 0 in movie_arousal_drug2: movie_arousal_drug2.remove(0)
    
    while 0 in movie_valence_drug1: movie_valence_drug1.remove(0)
    while 0 in movie_valence_drug2: movie_valence_drug2.remove(0)
    


# - BOXPLOT _ AROUSAL RATINGS

# prepare data for boxplot - Arousal
df_boxplot_ratings = pd.DataFrame(data = {"Drug1": arousal_ratings_drug1, "Drug2": arousal_ratings_drug2})
df_boxplot_ratings = df_boxplot_ratings.stack().reset_index()
df_boxplot_ratings.columns = ["Video", "Drug", "Rating"]
df_boxplot_ratings = df_boxplot_ratings.explode("Rating")

plt.figure(1)
plt.suptitle('Arousal Ratings', fontsize=16)
rcParams['figure.figsize'] = 16,8

df_boxplot_ratings.index = np.rint((df_boxplot_ratings.index/2)-0.1)
ax = sns.boxplot(x=df_boxplot_ratings.index, y="Rating", hue="Drug", data=df_boxplot_ratings) #, palette = "ch:s=.25,rot=.5"
ax = sns.swarmplot(x=df_boxplot_ratings.index, y="Rating", hue="Drug", data=df_boxplot_ratings, alpha = 0.6, dodge = True, color = ".3", size = 3.5)
xlabels = ['Neutral', 'Neutral', 'Erotic', 'Erotic', 'Social-', 'Social-', 'Social+', 'Social+', 'Horror', 'Horror', 'Neutral', 'Neutral', 'Social-', 'Social-', 'Erotic', 'Erotic', 'Horror', 'Horror', 'Social+', 'Social+']
ax.set_xticklabels(xlabels, rotation=70)
ax.set_xlabel('Movie')
ax.set(ylim=(0, 10))
plt.show()


# - BOXPLOT _ VALENCE RATINGS

# prepare data for boxplot - Valence
df_boxplot_ratings = pd.DataFrame(data = {"Drug1": valence_ratings_drug1, "Drug2": valence_ratings_drug2})
df_boxplot_ratings = df_boxplot_ratings.stack().reset_index()
df_boxplot_ratings.columns = ["Video", "Drug", "Rating"]
df_boxplot_ratings = df_boxplot_ratings.explode("Rating")

plt.figure(2)
plt.suptitle('Valence Ratings', fontsize=16)
rcParams['figure.figsize'] = 16,8

df_boxplot_ratings.index = np.rint((df_boxplot_ratings.index/2)-0.1)
ax = sns.boxplot(x=df_boxplot_ratings.index, y="Rating", hue="Drug", data=df_boxplot_ratings) #, palette = "ch:s=.25,rot=.5"
ax = sns.swarmplot(x=df_boxplot_ratings.index, y="Rating", hue="Drug", data=df_boxplot_ratings, alpha = 0.6, dodge = True, color = ".3", size = 3.5)
xlabels = ['Neutral', 'Neutral', 'Erotic', 'Erotic', 'Social-', 'Social-', 'Social+', 'Social+', 'Horror', 'Horror', 'Neutral', 'Neutral', 'Social-', 'Social-', 'Erotic', 'Erotic', 'Horror', 'Horror', 'Social+', 'Social+']
ax.set_xticklabels(xlabels, rotation=70)
ax.set_xlabel('Movie')
ax.set(ylim=(0, 10))
plt.show()



# - LINEAR REG. PLOTS

# PLOTS AROUSAL AND VALENCE AGAINST PERSUBJECT_SUM DATA (to see possible trend between both)

'''i = 0

for movie_cat in movie_categories_idx:
    
    # For Arousal
    
    fig, axes = plt.subplots(1, 4, figsize=(16, 4))
    
    fig.suptitle('Short Ratings', fontsize=14)
    
    xlim = [0.02,0.31]
    ylim = [0,10]
    for ax in axes:
        ax.set_xlim(xlim)
        ax.set_ylim(ylim)

    sns.regplot(x=ISC_persubject_sum_drug1[movie_cat[0]], y=arousal_ratings_drug1[movie_cat[0]], ci=None, ax=axes[0], truncate=False)
    sns.regplot(x=ISC_persubject_sum_drug2[movie_cat[0]], y=arousal_ratings_drug2[movie_cat[0]], ci=None, ax=axes[0], truncate=False)
    axes[0].set_title(movie_categories[i])
    
    sns.regplot(x=ISC_persubject_sum_drug1[movie_cat[1]], y=arousal_ratings_drug1[movie_cat[1]], ci=None, ax=axes[1], truncate=False)
    sns.regplot(x=ISC_persubject_sum_drug2[movie_cat[1]], y=arousal_ratings_drug2[movie_cat[1]], ci=None, ax=axes[1], truncate=False)
    axes[1].set_title(movie_categories[i+1])
    
    sns.regplot(x=ISC_persubject_sum_drug1[movie_cat[2]], y=arousal_ratings_drug1[movie_cat[2]], ci=None, ax=axes[2], truncate=False)
    sns.regplot(x=ISC_persubject_sum_drug2[movie_cat[2]], y=arousal_ratings_drug2[movie_cat[2]], ci=None, ax=axes[2], truncate=False)
    axes[2].set_title(movie_categories[i+2])
    
    sns.regplot(x=ISC_persubject_sum_drug1[movie_cat[3]], y=arousal_ratings_drug1[movie_cat[3]], ci=None, ax=axes[3], truncate=False)
    sns.regplot(x=ISC_persubject_sum_drug2[movie_cat[3]], y=arousal_ratings_drug2[movie_cat[3]], ci=None, ax=axes[3], truncate=False)
    axes[3].set_title(movie_categories[i+3])
    
    for ax in axes:
        ax.legend([],[], frameon=False)
        ax.set_yticks([i for i in range(1,10,1)])
        ax.set_ylabel('Arousal Rating')
        ax.set_xlabel('ISC')
        ax.set_xticks([-0.01,.04,.09,.14,.19,.24,.29,.34])
         
    fig.tight_layout()
    
    # For Valence
    
    fig, axes = plt.subplots(1, 4, figsize=(16, 4))
    
    fig.suptitle('Short Ratings', fontsize=14)
    
    xlim = [0.02,0.31]
    ylim = [0,10]
    for ax in axes:
        ax.set_xlim(xlim)
        ax.set_ylim(ylim)
    
    sns.regplot(x=ISC_persubject_sum_drug1[movie_cat[0]], y=valence_ratings_drug1[movie_cat[0]], ci=None, ax=axes[0], truncate=False)
    sns.regplot(x=ISC_persubject_sum_drug2[movie_cat[0]], y=valence_ratings_drug2[movie_cat[0]], ci=None, ax=axes[0], truncate=False)
    axes[0].set_title(movie_categories[i])
    
    sns.regplot(x=ISC_persubject_sum_drug1[movie_cat[1]], y=valence_ratings_drug1[movie_cat[1]], ci=None, ax=axes[1], truncate=False)
    sns.regplot(x=ISC_persubject_sum_drug2[movie_cat[1]], y=valence_ratings_drug2[movie_cat[1]], ci=None, ax=axes[1], truncate=False)
    axes[1].set_title(movie_categories[i+1])
    
    sns.regplot(x=ISC_persubject_sum_drug1[movie_cat[2]], y=valence_ratings_drug1[movie_cat[2]], ci=None, ax=axes[2], truncate=False)
    sns.regplot(x=ISC_persubject_sum_drug2[movie_cat[2]], y=valence_ratings_drug2[movie_cat[2]], ci=None, ax=axes[2], truncate=False)
    axes[2].set_title(movie_categories[i+2])
    
    sns.regplot(x=ISC_persubject_sum_drug1[movie_cat[3]], y=valence_ratings_drug1[movie_cat[3]], ci=None, ax=axes[3], truncate=False)
    sns.regplot(x=ISC_persubject_sum_drug2[movie_cat[3]], y=valence_ratings_drug2[movie_cat[3]], ci=None, ax=axes[3], truncate=False) #, truncate=False
    axes[3].set_title(movie_categories[i+3])
    
    for ax in axes:
        ax.legend([],[], frameon=False)
        ax.set_yticks([i for i in range(1,10,1)])
        ax.set_ylabel('Valence Rating')
        ax.set_xlabel('ISC')
        ax.set_xticks([-0.01,.04,.09,.14,.19,.24,.29,.34])
         
    fig.tight_layout()
    
    i+=4'''
    
'-----------------------------------------------------------------------------'

### SHOULD I DO BARPLOT OF EACH VALENCE/AROUSAL SCALE (1-9) AGAINST AVERAGE ISC (PER MOVIE, DRUG)? ###

'-----------------------------------------------------------------------------'

# LIVE RATINGS ANALYSIS (Arousal and Valence)

'''path_import = '/Users/danielcorreia/Desktop/live_ratings/'

# objetivo é ver se há maior corr em drug1 vs. drug2 com ISC_persecond of drug1/2

movie_id = ['5001.mp4','5002.mp4','2002.mp4','2000.mp4','3004.mp4','3003.mp4','4004.mp4','4006.mp4',
            '1002.mp4','1001.mp4','5000.mp4','5008.mp4','3002.mp4','3007.mp4','2009.mp4','2003.mp4',
            '1007.mp4','1006.mp4','4001.mp4','4009.mp4',]

live_ratings_drug1 = []
live_ratings_drug2 = []

subjects_w_live_ratings = []

for (dirpath, dirnames, filenames) in os.walk(path_import):
    
    subjects_w_live_ratings.extend(dirnames)


for subject in subjects_w_live_ratings:
    
    if subject not in rejected_ratings:
    
        lines = [x.split('\t') for x in open(path_import + subject + '/RESULTS_VIDEO_VB.txt').readlines()]

        i=0
        subject_ratings = []
        movie_ratings = np.zeros(195)
    
        for line in lines[1::]:

            if movie_id[i] == line[3]:

                timepoint = int(((float(line[8])/1000)*195)/40) #gets timepoint (ms to sec) adjusted to a 195 long vector to match ISC_persecond
                rating = int(line[7]) #gets rating
                movie_ratings[timepoint-1::] = rating 
                
            else:
                #saves ratings of last movie
                subject_ratings.append(movie_ratings.tolist())
                    
                #moves on to next movie category
                i+=1
                movie_ratings = np.zeros(195)
                    
                timepoint = int(((float(line[8])/1000)*195)/40)
                rating = int(line[7])
                    
                movie_ratings[timepoint-1::] = rating 
        
        subject_ratings.append(movie_ratings.tolist()) #saves last movie

        if len(subject_ratings) == 19:
            # troubleshooting cases where last rating is missing due to no change (only happens in last, otherwise is resgisted as a single 0)
        
            movie_ratings = np.zeros(195)
            subject_ratings.append(movie_ratings.tolist())
            
        if subject in subjects_DRUG1:
            live_ratings_drug1.append(subject_ratings)
    
        if subject in subjects_DRUG2:
            live_ratings_drug2.append(subject_ratings)


mean_live_ratings1 = np.mean(live_ratings_drug1, axis=0)
mean_live_ratings2 = np.mean(live_ratings_drug2, axis=0)


i = 0

for movie_cat in movie_categories_idx:
    
    fig, axes = plt.subplots(1, 4, figsize=(16, 4))
    
    fig.suptitle('Live Ratings', fontsize=14)
    
    xlim = [0,40]
    ylim = [-5,5]
    for ax in axes:
        ax.set_xlim(xlim)
        ax.set_ylim(ylim)
    
    x = np.arange(0, 40, 40/195)
    
    sns.lineplot(x=x, y=mean_live_ratings1[movie_cat[0]], ax=axes[0])
    sns.lineplot(x=x, y=mean_live_ratings2[movie_cat[0]], ax=axes[0])
    axes[0].set_title(movie_categories[i])
    
    sns.lineplot(x=x, y=mean_live_ratings1[movie_cat[1]], ax=axes[1])
    sns.lineplot(x=x, y=mean_live_ratings2[movie_cat[1]], ax=axes[1])
    axes[1].set_title(movie_categories[i+1])
    
    sns.lineplot(x=x, y=mean_live_ratings1[movie_cat[2]], ax=axes[2])
    sns.lineplot(x=x, y=mean_live_ratings2[movie_cat[2]], ax=axes[2])
    axes[2].set_title(movie_categories[i+2])
    
    sns.lineplot(x=x, y=mean_live_ratings1[movie_cat[3]], ax=axes[3])
    sns.lineplot(x=x, y=mean_live_ratings2[movie_cat[3]], ax=axes[3])
    axes[3].set_title(movie_categories[i+3])
    
    
    for ax in axes:
        ax.legend([],[], frameon=False)
        ax.set_yticks([i for i in range(-5,6,1)])
        ax.set_ylabel('Average Rating')
        ax.set_xlabel('Time (sec)')
        #ax.set_xticks([-0.01,.04,.09,.14,.19,.24,.29,.34])
         
    fig.tight_layout()
    
    i += 4'''
