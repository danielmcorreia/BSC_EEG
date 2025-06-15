#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 11 15:01:39 2021

@author: Daniel
"""

#import mne
import pandas as pd
import seaborn as sns
from scipy.io import loadmat
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams
import math
import sys
#from statannot import add_stat_annotation
from prop_test import prop_ztest, eval_pval, wilcoxon_paired

data_folder_path = '/Volumes/LaCie Armaz/PAPER_REVISIONS_1902/ISC_CALC_out/'
#data_folder_path = '/Volumes/LaCie Armaz/FINAL_PREPROCESSING_PAPER/ISC_CALC_out/'
#data_folder_path = '/Volumes/disk/Work/FINAL_PREPROCESSING/ISC_CALC/'

# define number of components to use
NComps = 3

'-----------------------------------------------------------------------------'


# 1 pd.dataframe would be so much more efficient than a thousand dicts; change this

# Import DRUG1

ISC_persecond_drug1 = {}
pvalues_persecond_drug1 = {}
chance_val_perwindow_drug1 = {}

#ISC_persubject_drug1 = {}
ISC_persubject_sum_drug1 = {}
pvalues_persubject_drug1 = {}
chance_val_persubject_drug1 = {}

ISC_spatial_filter_eig_drug1 = {}
pvalues_spatial_filter_eigval_drug1 = {}
chance_spatial_filter_eigval_drug1 = {}

# Import DRUG2

ISC_persecond_drug2 = {}
pvalues_persecond_drug2 = {}
chance_val_perwindow_drug2 = {}

#ISC_persubject_drug2 = {}
ISC_persubject_sum_drug2 = {}
pvalues_persubject_drug2 = {}
chance_val_persubject_drug2 = {}

ISC_spatial_filter_eig_drug2 = {}
pvalues_spatial_filter_eigval_drug2 = {}
chance_spatial_filter_eigval_drug2 = {}

# import pvalues differences
pvalues_dif = {}

# Loading Data

for video_num in range(1,21):
    
    # DRUG1
    
    ISC_persecond_drug1["ISC_persecond{}".format(video_num)] = \
        loadmat(data_folder_path + 'DRUG1/VIDEO' + str(video_num) + '_ISC_persecond')
    pvalues_persecond_drug1["pvalues{}".format(video_num)] = \
        (loadmat(data_folder_path + 'DRUG1/VIDEO' + str(video_num) + '_pvals_signific_persec')['pvals'].T)[0:3]
    chance_val_perwindow_drug1["chance_val_perwindow{}".format(video_num)] = \
        np.percentile(np.transpose(loadmat(data_folder_path + 'DRUG1/VIDEO' + str(video_num) + '_chance_val_perwindow')['chance_val_perwindow1'][:, 0:3]),95,0) #95th percentile (point at which 5% of a population set will exceed the referenced value) of the null distribution served as the threshold for statistical significance, so we are taking that here from the permut data for visualizations purposes (p-values already reflect that)
    
    ###ISC_persubject_drug1["ISC_persubject{}".format(video_num)] = \
        ###loadmat(data_folder_path + 'DRUG1/VIDEO' + str(video_num) + '_ISC_persubject')['ISC_persubject1']
    ISC_persubject_sum_drug1["ISC_persubject_sum{}".format(video_num)] = \
        loadmat(data_folder_path + 'DRUG1/VIDEO' + str(video_num) + '_ISC_persubject_sum')['ISC_persubject_sum_drug1'][0]
    pvalues_persubject_drug1["pvalues{}".format(video_num)] = \
        (loadmat(data_folder_path + 'DRUG1/VIDEO' + str(video_num) + '_pvals_signific_persubject_sum')['pvals'].T).mean()
    chance_val_persubject_drug1["chance_persubject_sum{}".format(video_num)] = \
        np.percentile(loadmat(data_folder_path + 'DRUG1/VIDEO' + str(video_num) + '_chance_val_persubject_sum')['chance_val_persubject_sum1'],95,2).flatten()
    
    ISC_spatial_filter_eig_drug1["ISC_spatial_filter_eig{}".format(video_num)] = \
        loadmat(data_folder_path + 'DRUG1/VIDEO' + str(video_num) + '_spatial_filter_eigval')['ISC_drug1'][0:NComps]
    pvalues_spatial_filter_eigval_drug1["pvalues{}".format(video_num)] = \
        (loadmat(data_folder_path + 'DRUG1/VIDEO' + str(video_num) + '_pvals_adjust_signific_eigen')['pval_adjust'].T)[0:3]
    chance_spatial_filter_eigval_drug1["chance_spatial_filter_eigval{}".format(video_num)] = \
        np.percentile(loadmat(data_folder_path + 'DRUG1/VIDEO' + str(video_num) + '_chance_spatial_filter_eigval')['chance_val_eig1'][0][0:NComps],95,1)
    

    # DRUG2
    
    ISC_persecond_drug2["ISC_persecond{}".format(video_num)] = \
        loadmat(data_folder_path + 'DRUG2/VIDEO' + str(video_num) + '_ISC_persecond')
    pvalues_persecond_drug2["pvalues{}".format(video_num)] = \
        (loadmat(data_folder_path + 'DRUG2/VIDEO' + str(video_num) + '_pvals_signific_persec')['pvals'].T)[0:3]
    chance_val_perwindow_drug2["chance_val_perwindow{}".format(video_num)] = \
        np.percentile(np.transpose(loadmat(data_folder_path + 'DRUG2/VIDEO' + str(video_num) + '_chance_val_perwindow')['chance_val_perwindow2'][:, 0:3]),95,0)
    
    ###ISC_persubject_drug2["ISC_persubject{}".format(video_num)] = \
        ###loadmat(data_folder_path + 'DRUG2/VIDEO' + str(video_num) + '_ISC_persubject')['ISC_persubject2']
    ISC_persubject_sum_drug2["ISC_persubject_sum{}".format(video_num)] = \
        loadmat(data_folder_path + 'DRUG2/VIDEO' + str(video_num) + '_ISC_persubject_sum')['ISC_persubject_sum_drug2'][0]
    pvalues_persubject_drug2["pvalues{}".format(video_num)] = \
        (loadmat(data_folder_path + 'DRUG2/VIDEO' + str(video_num) + '_pvals_signific_persubject_sum')['pvals'].T).mean()
    chance_val_persubject_drug2["chance_persubject_sum{}".format(video_num)] = \
        np.percentile(loadmat(data_folder_path + 'DRUG2/VIDEO' + str(video_num) + '_chance_val_persubject_sum')['chance_val_persubject_sum2'],95,2).flatten()

    ISC_spatial_filter_eig_drug2["ISC_spatial_filter_eig{}".format(video_num)] = \
        loadmat(data_folder_path + 'DRUG2/VIDEO' + str(video_num) + '_spatial_filter_eigval')['ISC_drug2'][0:NComps]
    pvalues_spatial_filter_eigval_drug2["pvalues{}".format(video_num)] = \
        (loadmat(data_folder_path + 'DRUG2/VIDEO' + str(video_num) + '_pvals_adjust_signific_eigen')['pval_adjust'].T)[0:3]
    chance_spatial_filter_eigval_drug2["chance_spatial_filter_eigval{}".format(video_num)] = \
        np.percentile(loadmat(data_folder_path + 'DRUG2/VIDEO' + str(video_num) + '_chance_spatial_filter_eigval')['chance_val_eig2'][0][0:NComps],95,1)
    

    # pvalues differences
    
    pvalues_dif["pvalues{}".format(video_num)] = \
        (loadmat(data_folder_path + 'stats/pvals_difference/VIDEO' + str(video_num) + '_pvals_difference')['pvals'].T)#[0:3]


# Loading External Statistics

#pvals_wsrt = pd.read_csv(data_folder_path+'stats/pvals_WSRT.txt', sep='\t', header=None)
#corr_vals = pd.read_csv(data_folder_path+'stats/corr_vals.txt', sep='\t', header=None)

# movie labels
xlabels = ['Neutral', 'Neutral', 'High +', 'High +', 'Low -', 'Low -', 'Low +', 'Low +', 'High -', 'High -', 'Neutral', 'Neutral', 'Low -', 'Low -', 'High +', 'High +', 'High -', 'High -', 'Low +', 'Low +']
#xlabels = ['Neutral', 'Neutral', 'Erotic', 'Erotic', 'Social-', 'Social-', 'Social+', 'Social+', 'Horror', 'Horror', 'Neutral', 'Neutral', 'Social-', 'Social-', 'Erotic', 'Erotic', 'Horror', 'Horror', 'Social+', 'Social+']


'-----------------------------------------------------------------------------'  
# BOXPLOT PERSUBJECT

def boxplot_persubject(ISC1, ISC2, chance_visible=False):
    
    '''Plots boxplot of...
    Requires:
    Ensures: '''

    # prepare data for boxplot
    df_boxplot_persub = pd.DataFrame(data = {"Drug1": ISC_persubject_sum_drug1, "Drug2": ISC_persubject_sum_drug2})
    df_boxplot_persub = df_boxplot_persub.stack().reset_index()
    df_boxplot_persub.columns = ["Video", "Drug", "ISC"]
    df_boxplot_persub = df_boxplot_persub.explode("ISC")

    plt.figure(1)
    rcParams['figure.figsize'] = 27,9

    #df_boxplot_persub.index = np.rint((df_boxplot_persub.index/2)-0.01)
    new_list=[]
    original_list = df_boxplot_persub.index
    
    print(len(original_list))

    df_boxplot_persub['x_pos'] = df_boxplot_persub.index

    for index, row in df_boxplot_persub.iterrows():
        
        if row['Drug'] == 1:
            row['x_pos'] = row['Video']
        if row['Drug'] == 2:
            row['x_pos'] = row['Video'] + 1

    print(df_boxplot_persub)



    #df_boxplot_persub['x_pos'] = df_boxplot_persub.index
    #df_boxplot_persub = df_boxplot_persub.reset_index(drop=True)

    ax = sns.boxplot(x="x_pos", y="ISC", hue="Drug", width=0.95, data=df_boxplot_persub, showfliers = False, color='gold', hue_order = ['Drug2', 'Drug1'], palette=sns.color_palette(['darkgray', 'gold']))
    ax = sns.swarmplot(x="x_pos", y="ISC", hue="Drug", data=df_boxplot_persub, alpha = 0.7, dodge = True, color = "0", size = 3.5, hue_order = ['Drug2', 'Drug1'], palette=sns.color_palette(['black', 'black']))
    
    # chance boxplots
    if chance_visible:


        df_boxplot_chance_persub = pd.DataFrame(data = {"Drug1": chance_val_persubject_drug1, "Drug2": chance_val_persubject_drug2})
        df_boxplot_chance_persub = df_boxplot_chance_persub.stack().reset_index()
        df_boxplot_chance_persub.columns = ["Video", "Drug", "ISC"]
        df_boxplot_chance_persub = df_boxplot_chance_persub.explode("ISC")

        df_boxplot_chance_persub['x_pos'] = df_boxplot_chance_persub.index

        for index, row in df_boxplot_chance_persub.iterrows():
        
            if row['Drug'] == 1:
                row['x_pos'] = row['Video']
            if row['Drug'] == 2:
                row['x_pos'] = row['Video'] + 1
    
        #df_boxplot_chance_persub.index = np.rint((df_boxplot_chance_persub.index/2)-0.01)
        #df_boxplot_chance_persub['x_pos'] = df_boxplot_chance_persub.index
        df_boxplot_chance_persub = df_boxplot_chance_persub.reset_index(drop=True)
        
        ax = sns.boxplot(x="x_pos", y="ISC", hue="Drug", width=0.95, data=df_boxplot_chance_persub, showfliers = False, color='snow', hue_order = ['Drug2', 'Drug1'], boxprops=dict(facecolor=(0,0,0,0)))
        ax = sns.swarmplot(x="x_pos", y="ISC", hue="Drug", data=df_boxplot_chance_persub, alpha = 0.2, dodge = True, color = "0", size = 3.5, hue_order = ['Drug2', 'Drug1'], palette=sns.color_palette(['black', 'black']))
        
    # add significance annotations
    x_pos_1 = -1.8
    x_pos_2 = -1.8
    for movie in range(1,21):
        x_pos_1 += 2
        x_pos_2 += 2
        plt.text(x_pos_1, 0.003,  eval_pval(pvalues_persubject_drug1['pvalues{}'.format(movie)])[1], fontsize = 7, rotation=90)    
        plt.text(x_pos_2 +0.4, 0.003,  eval_pval(pvalues_persubject_drug2['pvalues{}'.format(movie)])[1], fontsize = 7, rotation=90)  
    
    # make it pretty
    xtick_loc = [0.5,2.5,4.5,6.5,8.5,10.5,12.5,14.5,16.5,18.5,20.5,22.5,24.5,26.5,28.5,30.5,32.5,34.5,36.5,38.5]
    ax.set_xticks(xtick_loc)
    ax.set_xticklabels(xlabels, rotation=70)
    ax.set_xlabel('Movie', weight='bold')
    ax.set_ylabel('ISC', weight='bold')
    ax.set(ylim=(0, 0.35))
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    plt.rcParams['axes.linewidth'] = 2
    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles[:2], ['PLC', 'OT'], frameon=False)
    plt.show()
        
'-----------------------------------------------------------------------------'
# COMPONENT COEFFICIENTS

def spatial_filter_eigval(ISC1, ISC2):

    # TOO INNEFICIENT; FIX THIS
    
    # prepare real data for plot
    df_eig = pd.DataFrame(data = {"Drug1": ISC_spatial_filter_eig_drug1, "Drug2": ISC_spatial_filter_eig_drug2})
    df_eig = df_eig.stack().reset_index()
    df_eig.columns = ["Video", "Drug", "eig"]
    df_eig = df_eig.explode("eig").explode("eig")
    df_eig = df_eig.astype({"eig":float})

    # prepare randomized data for plot
    df_chance_eig = pd.DataFrame(data = {"Drug1": chance_spatial_filter_eigval_drug1, "Drug2": chance_spatial_filter_eigval_drug2})
    df_chance_eig = df_chance_eig.stack().reset_index()
    df_chance_eig.columns = ["Video", "Drug", "eig"]
    df_chance_eig = df_chance_eig.explode("eig").explode("eig")

    
    for video_num in range(1,21):
    
        fig, axs = plt.subplots(1,2, figsize=(6,4), sharey=True)
        
        # make it pretty
        fig.suptitle("Movie{}".format(video_num), weight='bold')
        axs[0].set(ylim=((0, 0.25)))
        
        # filters dataframes
        df_eig_selected = df_eig.query("Video == 'ISC_spatial_filter_eig{}'".format(video_num))
        df_eig_selected['comp_nums'] = list(range(1, NComps*2+1))
        df_chance_eig_selected = df_chance_eig.query("Video == 'chance_spatial_filter_eigval{}'".format(video_num))
        df_chance_eig_selected['comp_nums'] = list(range(1, NComps+1)) + list(range(1, NComps+1))
        
        
        for drug in range(1,3):
            
            # plots line and bar plots
            sns.pointplot(data=df_eig_selected.groupby(['Drug']).get_group('Drug{}'.format(drug)), x="comp_nums", y="eig", color='black', ax = axs[drug-1],)
            sns.barplot(x="comp_nums", y="eig", data=df_chance_eig_selected.groupby(['Drug']).get_group('Drug{}'.format(drug)), color = 'grey', ax = axs[drug-1])
            
            # plots significance
            for line in range(0,NComps):
                if df_eig_selected.groupby(['Drug']).get_group('Drug{}'.format(drug)).iloc[line]['eig'] > df_chance_eig_selected.groupby(['Drug']).get_group('Drug{}'.format(drug)).iloc[line]['eig']:
                    axs[drug-1].text(line, 0.23, '*')
                else: axs[drug-1].text(line, 0.23, 'n.s.')
         
                
        # make it pretty
        axs[0].set_xlabel('Component #')
        axs[1].set_xlabel('Component #')
        axs[0].set_title('OT')
        axs[1].set_title('PLC')
        axs[1].set_yticks([])
        axs[1].yaxis.label.set_visible(False)
        axs[0].set_ylabel('Coefficient')
        axs[0].spines['top'].set_visible(False)
        axs[1].spines['top'].set_visible(False)
        axs[0].spines['right'].set_visible(False)
        axs[1].spines['right'].set_visible(False)
        axs[0].set_yticks(np.arange(0,0.26,0.05))
        
        plt.savefig("/Volumes/LaCie Armaz/DPrataLab/ISC/FINAL_PREPROCESSING/RESULTS_out/EIGENVALS/Figure{}.pdf".format(video_num), format="pdf", bbox_inches="tight")
        
        #plt.show()
    
        fig.tight_layout()

'-----------------------------------------------------------------------------'
# TIMESERIES PERSECOND AND RESPECTIVE PROP. PLOT(S) ON THE SIDE

def timeseries_persecond(ISC_persecond_drug1, ISC_persecond_drug2, pvalues_drug1, pvalues_drug2, pvalues_dif, pvals_wsrt, corr_values, plot_time_windows=True, plot_chance_levels=True, plot_type=None, thresh=.01):
    
    '''
    plot_type = "merged, separated, None" '''
    
    
    global p_vals_signif
    global p_vals_diff
    global p_vals_merged
    
    p_vals_signif = []
    p_vals_diff = []
    p_vals_merged = []

    video_num = 1
    for (k1,v1), (k2,v2) in zip(ISC_persecond_drug1.items(), ISC_persecond_drug2.items()):
    
        if plot_type == "separated":
            fig, axs = plt.subplots(3, 3, figsize=(12,7), gridspec_kw={'width_ratios': [5, 1, 1]}, constrained_layout=True)
            fig.suptitle('Video{}'.format(video_num) + ': ' + str(xlabels[video_num-1]), weight='bold', fontsize='large')
        
            p_vals_comp_signif = []
            p_vals_comp_diff = []
            
        if plot_type == "merged":
            fig, axs = plt.subplots(3, 2, figsize=(12,7), gridspec_kw={'width_ratios': [5, 1]}, constrained_layout=True)
            fig.suptitle('Video{}'.format(video_num) + ': ' + str(xlabels[video_num-1]), weight='bold', fontsize='large')
        
            p_vals_comp_merged = []
            
        if plot_type == "":
            fig, axs = plt.subplots(3, 1, figsize=(18,4.5), constrained_layout=True)
            fig.suptitle('Video{}'.format(video_num) + ': ' + str(xlabels[video_num-1]), weight='bold', fontsize='large')
        
            p_vals_comp_merged = []
            
        
        for c in range(0,3):
            
            # PLOT TIMESERIES PER SECOND
            
            persecond_drug1 = v1['ISC_persecond_drug1'][c]
            persecond_drug2 = v2['ISC_persecond_drug2'][c]
        
            #x = np.arange(0, 40, 40/195) #adjust x axis coz ISC is calculated over a time window of 5 secs
            x = np.arange(1, 40, 39/190) #NEW
            
            if plot_type == "":
                
                # line plot of persecond
                axs[c].set_title('Component {}:'.format(c+1), loc='left', fontsize=11)
                axs[c].plot(x, persecond_drug1, label='OT', color='gold')
                axs[c].plot(x, persecond_drug2, label='PLC', color='black')
            
                # plot legend 
                leg = axs[0].legend(bbox_to_anchor=(0.9, 1.45), loc='upper left', ncol=2, mode = "expand", frameon=False)
                for line in leg.get_lines():
                    line.set_linewidth(5)
                
                # plot Pearson corr and wsrt stats info
                #axs[c].text(33, 0.35, 'corr. = ' + str(round(corr_vals[c][video_num-1],2)), fontsize=9)
                #axs[c].text(33, 0.30, 'w.s.r.t = ' + eval_pval(pvals_wsrt[c][video_num-1])[1] + ' ({})'.format(eval_pval(pvals_wsrt[c][video_num-1])[0]), fontsize=9)
            
    
                # to plot chance levels of significance Drug1 and Drug2
                if plot_chance_levels:
                    axs[c].fill_between(x, -0.03, chance_val_perwindow_drug1['chance_val_perwindow{}'.format(video_num)][c], facecolor='gold', alpha=0.4)
                    axs[c].fill_between(x, -0.03, chance_val_perwindow_drug2['chance_val_perwindow{}'.format(video_num)][c], facecolor='black', alpha=0.3)       
            
                # to plot significant time-windows based on permutation tests on the difference of ISC drug1-drug2
                if plot_time_windows:
                    axs[c].fill_between(x, -0.03, 0.4, facecolor='grey', alpha=0.2, where = pvalues_dif['pvalues{}'.format(video_num)][c]<thresh)
                    axs[c].fill_between(x, -0.03, 0.4, facecolor='grey', alpha=0.2, where = pvalues_dif['pvalues{}'.format(video_num)][c]>1-thresh)
            
                
                # makes plots pretty
            
                axs[c].set_yticks([0.05,0.2,0.35])
                axs[c].set_ylim(-0.03, 0.4)
            
                axs[c].spines['top'].set_visible(False)
                axs[c].spines['right'].set_visible(False)
            
                axs[c].set_ylabel('ISC')
                axs[c].set_xlabel('Time (sec)')
            
            
            
            else:
                # makes plots pretty
            
                axs[c][0].set_yticks([0,0.1,0.2,0.3, 0.4])
                axs[c][0].set_ylim(-0.03, 0.4)
            
                axs[c][0].spines['top'].set_visible(False)
                axs[c][0].spines['right'].set_visible(False)
            
                axs[c][0].set_ylabel('ISC')
                axs[c][0].set_xlabel('Time (sec)')
            
                # line plot of persecond
                axs[c][0].set_title('Component {}:'.format(c+1), loc='left', fontsize=11)
                axs[c][0].plot(x, persecond_drug1, label='OT', color='gold')
                axs[c][0].plot(x, persecond_drug2, label='PLC', color='black')
            
                # plot legend 
                leg = axs[0][0].legend(bbox_to_anchor=(0.135, 1.45), loc='upper left', ncol=2, mode = "expand", frameon=False)
                for line in leg.get_lines():
                    line.set_linewidth(5)
                
                # plot Pearson corr and wsrt stats info
                axs[c][0].text(33, 0.35, 'corr. = ' + str(round(corr_vals[c][video_num-1],2)), fontsize=9)
                axs[c][0].text(33, 0.30, 'w.s.r.t = ' + eval_pval(pvals_wsrt[c][video_num-1])[1] + ' ({})'.format(eval_pval(pvals_wsrt[c][video_num-1])[0]), fontsize=9)
            
    
                # to plot chance levels of significance Drug1 and Drug2
                if plot_chance_levels:
                    axs[c][0].fill_between(x, -0.03, chance_val_perwindow_drug1['chance_val_perwindow{}'.format(video_num)][c], facecolor='gold', alpha=0.4)
                    axs[c][0].fill_between(x, -0.03, chance_val_perwindow_drug2['chance_val_perwindow{}'.format(video_num)][c], facecolor='black', alpha=0.3)       
            
                # to plot significant time-windows based on permutation tests on the difference of ISC drug1-drug2
                if plot_time_windows:
                    axs[c][0].fill_between(x, -0.03, 0.4, facecolor='grey', alpha=0.2, where = pvalues_dif['pvalues{}'.format(video_num)][c]<thresh)
                    axs[c][0].fill_between(x, -0.03, 0.4, facecolor='grey', alpha=0.2, where = pvalues_dif['pvalues{}'.format(video_num)][c]>1-thresh)
            
            
                # makes plots pretty
            
                axs[c][0].set_yticks([0,0.1,0.2,0.3, 0.4])
                axs[c][0].set_ylim(-0.03, 0.4)
            
                axs[c][0].spines['top'].set_visible(False)
                axs[c][0].spines['right'].set_visible(False)
            
                axs[c][0].set_ylabel('ISC')
                axs[c][0].set_xlabel('Time (sec)')
            
            if plot_type == "separated":
            
                # PLOT % SIGNIFICANT WINDOWS
            
                # get data
                perc_drug1 = ((sum(pvalues_drug1['pvalues{}'.format(video_num)][c]<thresh)) / len(pvalues_dif['pvalues{}'.format(video_num)][c]))*100
                perc_drug2 = ((sum(pvalues_drug2['pvalues{}'.format(video_num)][c]<thresh)) / len(pvalues_dif['pvalues{}'.format(video_num)][c]))*100
            
                axs[c][1].bar([0,1], [perc_drug1, perc_drug2], color='grey', alpha=0.5, tick_label=['OT', 'PLC'], align='center', width=0.9)
            
                axs[c][1].spines['top'].set_visible(False)
                axs[c][1].spines['right'].set_visible(False)
            
                axs[c][1].set_ylim(0, 80)

                axs[c][1].set_ylabel('% significant windows\n(p<{}, uncorrected)'.format(thresh))
            
                # plot significance
                pval = float(prop_ztest(sum(pvalues_drug1['pvalues{}'.format(video_num)][c]<thresh), len(pvalues_dif['pvalues{}'.format(video_num)][c]), perc_drug2/100)) #gets pval
                p_vals_comp_signif.append(pval)
                if not math.isnan(pval):
                    axs[c][1].plot([0, 0, 1, 1], [78, 80, 80, 78], lw=1, color='black')
                    axs[c][1].text((0+1)*.5, 80, eval_pval(pval)[0], ha='center', va='bottom')

            
                # PLOT % SIGNIFICANT DIFFERENCES
            
                # get data
                perc_drug1 = ((sum(pvalues_dif['pvalues{}'.format(video_num)][c]<thresh)) / len(pvalues_dif['pvalues{}'.format(video_num)][c]))*100
                perc_drug2 = ((sum(pvalues_dif['pvalues{}'.format(video_num)][c]>1-thresh)) / len(pvalues_dif['pvalues{}'.format(video_num)][c]))*100
            
                axs[c][2].bar([0,1], [perc_drug1, perc_drug2], color='grey', alpha=0.5, tick_label=['OT>PLC', 'PLC>OT'], align='center', width=0.9)
            
                axs[c][2].spines['top'].set_visible(False)
                axs[c][2].spines['right'].set_visible(False)
            
                axs[c][2].set_ylim(0, 80)
            
                # plot significance
                pval = float(prop_ztest(sum(pvalues_dif['pvalues{}'.format(video_num)][c]<thresh), len(pvalues_dif['pvalues{}'.format(video_num)][c]), perc_drug2/100)) #gets pval
                p_vals_comp_diff.append(pval)
                if not math.isnan(pval):
                    axs[c][2].plot([0, 0, 1, 1], [78, 80, 80, 78], lw=1, color='black')
                    axs[c][2].text((0+1)*.5, 80, eval_pval(pval)[0], ha='center', va='bottom')
            
            if plot_type == "merged":
                
                #get data
                perc_drug1 = ((sum(np.bitwise_and(pvalues_drug1['pvalues{}'.format(video_num)][c]<thresh, pvalues_dif['pvalues{}'.format(video_num)][c]<thresh))) / sum(pvalues_drug1['pvalues{}'.format(video_num)][c]<thresh))*100
                perc_drug2 = ((sum(np.bitwise_and(pvalues_drug2['pvalues{}'.format(video_num)][c]<thresh, pvalues_dif['pvalues{}'.format(video_num)][c]>1-thresh))) / sum(pvalues_drug2['pvalues{}'.format(video_num)][c]<thresh))*100
                
                #perc_drug1 = ((sum(pvalues_dif['pvalues{}'.format(video_num)][c]<thresh)) / len(pvalues_dif['pvalues{}'.format(video_num)][c]))*100
                #perc_drug2 = ((sum(pvalues_dif['pvalues{}'.format(video_num)][c]>1-thresh)) / len(pvalues_dif['pvalues{}'.format(video_num)][c]))*100
                
                axs[c][1].bar([0,1], [perc_drug1, perc_drug2], color='grey', alpha=0.5, tick_label=['OT>PLC', 'PLC>OT'], align='center', width=0.9)
            
                axs[c][1].spines['top'].set_visible(False)
                axs[c][1].spines['right'].set_visible(False)
            
                axs[c][1].set_ylim(0, 80)
                
                axs[c][1].set_ylabel('% significant windows\n(p<{}, uncorrected)'.format(thresh))
            
                # plot significance
                pval = float(prop_ztest(sum(np.bitwise_and(pvalues_drug1['pvalues{}'.format(video_num)][c]<thresh, pvalues_dif['pvalues{}'.format(video_num)][c]<thresh)), sum(pvalues_drug1['pvalues{}'.format(video_num)][c]<thresh), perc_drug2/100)) #gets pval
                p_vals_comp_merged.append(pval)
                if not math.isnan(pval):
                    axs[c][1].plot([0, 0, 1, 1], [79, 80, 80, 79], lw=1, color='black')
                    axs[c][1].text((0+1)*.5, 80, eval_pval(pval)[0], ha='center', va='bottom')
                
                
        if plot_type == "separated":
            p_vals_signif.append(p_vals_comp_signif)    
            p_vals_diff.append(p_vals_comp_diff)   
            #plt.savefig("/Volumes/LaCie Armaz/DPrataLab/ISC/FINAL_PREPROCESSING/RESULTS/PERSECOND_SEP/Figure{}_thresh{}.pdf".format(video_num, str(int(thresh*100))), format="pdf", bbox_inches="tight")

        if plot_type == "merged":
            p_vals_merged.append(p_vals_comp_merged) 
            #plt.savefig("/Volumes/LaCie Armaz/DPrataLab/ISC/FINAL_PREPROCESSING/RESULTS/PERSECOND_MERGED/Figure{}_thresh{}.pdf".format(video_num, str(int(thresh*100))), format="pdf", bbox_inches="tight")
            
        if plot_type == "":
            p_vals_merged.append(p_vals_comp_merged) 
            plt.savefig("/Volumes/LaCie Armaz/DPrataLab/ISC/FINAL_PREPROCESSING/RESULTS_out/PERSECOND/Figure{}_thresh{}.pdf".format(video_num, str(int(thresh*100))), format="pdf", bbox_inches="tight")


        video_num+=1
        
'-----------------------------------------------------------------------------'
# % significant time windows

def perc_windows(ISC1, ISC2, pvalues_drug1, pvalues_drug2, thresh=0.05):
    
    ''' barplot of % significant time windows, per drug, per MOVIE CATEGORY'''
    
    xlabels = ['Neutral', 'Low -', 'Low +', 'High -', 'High +']
    
    # lazy inneficient way, change this whenever possible...
    
    # prepares data
    store_perc_drug1_avg = []; store_perc_drug2_avg = []
    store_perc_drug1 = []; store_perc_drug2 = []
    
    nobs = 195
    
    for video_num in range(1,21):
        perc_drug1 = 0; perc_drug2 = 0
        
        for c in range(0,3): 
            # get data
            perc_drug1 += ((sum(pvalues_drug1['pvalues{}'.format(video_num)][c]<thresh)) / nobs)*100
            perc_drug2 += ((sum(pvalues_drug2['pvalues{}'.format(video_num)][c]<thresh)) / nobs)*100
            
        store_perc_drug1.append(perc_drug1/3); store_perc_drug2.append(perc_drug2/3)
    
    for videos in [[0,1,10,11],[4,5,12,13],[6,7,18,19],[8,9,16,17],[2,3,14,15]]: # neutral, low-, low+, high-, high+
        store_perc_drug1_avg.append( ((store_perc_drug1[videos[0]]) + (store_perc_drug1[videos[1]]) + (store_perc_drug1[videos[2]]) + (store_perc_drug1[videos[3]]))/4 )
        store_perc_drug2_avg.append( ((store_perc_drug2[videos[0]]) + (store_perc_drug2[videos[1]]) + (store_perc_drug2[videos[2]]) + (store_perc_drug2[videos[3]]))/4 )
        #print(wilcoxon_paired([store_perc_drug1[videos[0]], store_perc_drug1[videos[1]], store_perc_drug1[videos[2]], store_perc_drug1[videos[3]]], [store_perc_drug2[videos[0]], store_perc_drug2[videos[1]], store_perc_drug2[videos[2]], store_perc_drug2[videos[3]]]))
    
    # prepare data for barplot
    df = pd.DataFrame(data = {"Drug1": store_perc_drug1_avg, "Drug2": store_perc_drug2_avg}); df = df.stack().reset_index(); df.columns = ["Video", "Drug", "perc"]; df = df.explode("perc")
    
    # plots data
    plt.figure(1); rcParams['figure.figsize'] = 16,8
    df.index = np.rint((df.index/2)-0.1)
    ax = sns.barplot(x=df.index, y="perc", hue="Drug", data=df, color='gold', hue_order = ['Drug2', 'Drug1'], palette=sns.color_palette(['darkgray', 'gold']))
    
    # add significance annotations
    for video in range(1,6):
        pval = prop_ztest((store_perc_drug1_avg[video-1]/100) * nobs, nobs, store_perc_drug2_avg[video-1]/100)
        plt.text(video-1-0.06, 65,  str(eval_pval(pval)[0]) + str(eval_pval(pval)[1]), fontsize = 12)
        plt.plot([video-1.1, video-1.1, video-1+0.1, video-1+0.1], [62, 62.5, 62.5, 62], lw=1, color='black')
    
    # make comparisons btw movies to report
    for video1 in range(1,6):
        for video2 in range(1,6):
            print("{} x {}".format(xlabels[video1-1],xlabels[video2-1]))
            print(eval_pval(prop_ztest((store_perc_drug1_avg[video1-1]/100) * nobs, nobs, store_perc_drug1_avg[video2-1]/100))[0])
            print(eval_pval(prop_ztest((store_perc_drug2_avg[video1-1]/100) * nobs, nobs, store_perc_drug2_avg[video2-1]/100))[0])
    
    # make it pretty
    ax.set_xticklabels(xlabels, rotation=70)
    ax.set_xlabel('Movie', weight='bold')
    ax.set_ylabel('% significant timepoints \n (averaged across components and movie categories)', weight='bold', fontsize=12)
    ax.set(ylim=(0, 70))
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    plt.rcParams['axes.linewidth'] = 2
    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles[:2], ['PLC', 'OT'], frameon=False)
    plt.show()
    
'-----------------------------------------------------------------------------'
# % significant different time windows

def perc_windows_differences(ISC1, ISC2, pvalues_drug1, pvalues_drug2, thresh=0.05):
    
    ''' barplot of % significant time windows, per drug, per MOVIE CATEGORY'''
    
    xlabels = ['Neutral', 'Low -', 'Low +', 'High -', 'High +']
    
    # lazy inneficient way, change this whenever possible...
    
    # prepares data
    store_perc_drug1_avg = []; store_perc_drug2_avg = []
    store_count_drug1 = []; store_count_drug2 = []
    
    for video_num in range(1,21):
        count_drug1 = 0; count_drug2 = 0
        
        for c in range(0,3): 
            
            # get data ###from the significant timepoints, what's the proportion of them where there's a signfiicant difference between drugs, and in which direction (ie OT>Plc or Plc>OT)
            count_drug1 += sum(np.bitwise_and(np.bitwise_or(pvalues_persecond_drug1['pvalues{}'.format(video_num)][c]<thresh, pvalues_persecond_drug2['pvalues{}'.format(video_num)][c]<thresh), pvalues_dif['pvalues{}'.format(video_num)][c]<thresh))
            count_drug2 += sum(np.bitwise_and(np.bitwise_or(pvalues_persecond_drug1['pvalues{}'.format(video_num)][c]<thresh, pvalues_persecond_drug2['pvalues{}'.format(video_num)][c]<thresh), pvalues_dif['pvalues{}'.format(video_num)][c]>1-thresh))
        
        store_count_drug1.append(count_drug1/3); store_count_drug2.append(count_drug2/3)
    
        
    nobs = 195
        
    for videos in [[0,1,10,11],[4,5,12,13],[6,7,18,19],[8,9,16,17],[2,3,14,15]]: # neutral, low-, low+, high-, high+
        
        avg_counts = ((store_count_drug1[videos[0]]) + (store_count_drug1[videos[1]]) + (store_count_drug1[videos[2]]) + (store_count_drug1[videos[3]]))/4
        store_perc_drug1_avg.append( (avg_counts / nobs)*100 )

        avg_counts = ((store_count_drug2[videos[0]]) + (store_count_drug2[videos[1]]) + (store_count_drug2[videos[2]]) + (store_count_drug2[videos[3]]))/4
        store_perc_drug2_avg.append( (avg_counts / nobs)*100 )
    
    # prepare data for barplot
    df = pd.DataFrame(data = {"Drug1": store_perc_drug1_avg, "Drug2": store_perc_drug2_avg}); df = df.stack().reset_index(); df.columns = ["Video", "Drug", "perc"]; df = df.explode("perc")
    
    # plots data
    plt.figure(1); rcParams['figure.figsize'] = 7,5
    df.index = np.rint((df.index/2)-0.1)
    ax = sns.barplot(x=df.index, y="perc", hue="Drug", data=df, color='gold', hue_order = ['Drug2', 'Drug1'], palette=sns.color_palette(['darkgray', 'gold']))
    
    # add significance annotations
    for video in range(1,6):
        pval = prop_ztest((store_perc_drug1_avg[video-1]/100) * nobs, nobs, store_perc_drug2_avg[video-1]/100)
        plt.text(video-1-0.2, 19,  str(eval_pval(pval)[0]) + str(eval_pval(pval)[1]), fontsize = 14)
        plt.plot([video-1.2, video-1.2, video-1+0.2, video-1+0.2], [17.5, 17.8, 17.8, 17.5], lw=2, color='black')
    
    # make comparisons btw movies to report
    for video1 in range(1,6):
        for video2 in range(1,6):
            print("{} x {}".format(xlabels[video1-1],xlabels[video2-1]))
            print(eval_pval(prop_ztest((store_perc_drug1_avg[video1-1]/100) * nobs, nobs, store_perc_drug1_avg[video2-1]/100))[0])
            print(eval_pval(prop_ztest((store_perc_drug2_avg[video1-1]/100) * nobs, nobs, store_perc_drug2_avg[video2-1]/100))[0])
    
    # make it pretty
    ax.set_xticklabels(xlabels, rotation=70)
    ax.set_xlabel('Movie', weight='bold', fontsize=14)
    ax.set_ylabel('% significant differences \n (averaged across components and movie categories)', weight='bold', fontsize=14)
    ax.set(ylim=(0, 20))
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    plt.rcParams['axes.linewidth'] = 3
    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles[:2], ['PLC>OT', 'OT>PLC'], frameon=False)
    plt.show()
    

'-----------------------------------------------------------------------------'

# mean of significant different time windows

def mean_perc_windows_differences(ISC1, ISC2, pvalues_drug1, pvalues_drug2, pvalues_dif, thresh=0.05):
    
    ''' barplot of % significant time windows, per drug, per MOVIE CATEGORY'''

    #xlabels = ['Neutral', 'Low -', 'Low +', 'High -', 'High +']
    
    # lazy inneficient way, change this whenever possible...
    
    # prepares data
    store_perc_drug1_avg = []; store_perc_drug2_avg = []
    store_perc_drug1_avg_ = []; store_perc_drug2_avg_ = []
    
    store_perc_drug1 = []; store_perc_drug2 = []
    store_perc_drug1_ = []; store_perc_drug2_ = []
    
    for video_num in range(1,21):
        #perc_drug1 = 0; perc_drug2 = 0
        compare1 = []; compare2 =[]
        compare1_ = []; compare2_ =[]
        
        for c in range(0,3): 
            # get data
            
            # for statistically significant ISC for at least 1 of the drugs and for where's a signfiicant difference between them
            compare1.append( np.where( np.bitwise_and(np.bitwise_or(pvalues_persecond_drug1['pvalues{}'.format(video_num)][c]<thresh, pvalues_persecond_drug2['pvalues{}'.format(video_num)][c]<thresh), np.bitwise_or(pvalues_dif['pvalues{}'.format(video_num)][c]<thresh, pvalues_dif['pvalues{}'.format(video_num)][c]>1-thresh)), ISC_persecond_drug1['ISC_persecond{}'.format(video_num)]['ISC_persecond_drug1'][c], [np.nan]*190))
            compare2.append( np.where( np.bitwise_and(np.bitwise_or(pvalues_persecond_drug2['pvalues{}'.format(video_num)][c]<thresh, pvalues_persecond_drug1['pvalues{}'.format(video_num)][c]<thresh), np.bitwise_or(pvalues_dif['pvalues{}'.format(video_num)][c]>1-thresh, pvalues_dif['pvalues{}'.format(video_num)][c]<thresh)), ISC_persecond_drug2['ISC_persecond{}'.format(video_num)]['ISC_persecond_drug2'][c], [np.nan]*190))
            
            # for statistically significant ISC for BOTH drugs (same time)
            compare1_.append( np.where( np.bitwise_and(pvalues_persecond_drug1['pvalues{}'.format(video_num)][c]<thresh, pvalues_persecond_drug2['pvalues{}'.format(video_num)][c]<thresh), ISC_persecond_drug1['ISC_persecond{}'.format(video_num)]['ISC_persecond_drug1'][c], [np.nan]*190))
            compare2_.append( np.where( np.bitwise_and(pvalues_persecond_drug1['pvalues{}'.format(video_num)][c]<thresh, pvalues_persecond_drug2['pvalues{}'.format(video_num)][c]<thresh), ISC_persecond_drug2['ISC_persecond{}'.format(video_num)]['ISC_persecond_drug2'][c], [np.nan]*190))


            
            #compare1p += np.nanmean(np.where(np.bitwise_and(np.bitwise_or(pvalues_persecond_drug1['pvalues{}'.format(video_num)][c]<thresh, pvalues_persecond_drug2['pvalues{}'.format(video_num)][c]<thresh), pvalues_dif['pvalues{}'.format(video_num)][c]<thresh), ISC_persecond_drug1['ISC_persecond{}'.format(video_num)]['ISC_persecond_drug1'][c], [np.nan]*195))
            #compare2p += np.nanmean(np.where(np.bitwise_and(np.bitwise_or(pvalues_persecond_drug2['pvalues{}'.format(video_num)][c]<thresh, pvalues_persecond_drug1['pvalues{}'.format(video_num)][c]<thresh), pvalues_dif['pvalues{}'.format(video_num)][c]>1-thresh), ISC_persecond_drug2['ISC_persecond{}'.format(video_num)]['ISC_persecond_drug2'][c], [np.nan]*195))

            # antigo
            #compare1 += np.nanmean(np.where(np.bitwise_and(np.bitwise_or(pvalues_persecond_drug1['pvalues{}'.format(video_num)][c]<thresh, pvalues_persecond_drug2['pvalues{}'.format(video_num)][c]<thresh), pvalues_dif['pvalues{}'.format(video_num)][c]<thresh), ISC_persecond_drug1['ISC_persecond{}'.format(video_num)]['ISC_persecond_drug1'][c], [np.nan]*195))
            #compare2 += np.nanmean(np.where(np.bitwise_and(np.bitwise_or(pvalues_persecond_drug2['pvalues{}'.format(video_num)][c]<thresh, pvalues_persecond_drug1['pvalues{}'.format(video_num)][c]<thresh), pvalues_dif['pvalues{}'.format(video_num)][c]>1-thresh), ISC_persecond_drug2['ISC_persecond{}'.format(video_num)]['ISC_persecond_drug2'][c], [np.nan]*195))

            
            
            
            #compare1 = [value for idx, value in enumerate(ISC1['ISC_persecond{}'.format(video_num)][c]) if np.bitwise_and(pvalues_drug1['pvalues{}'.format(video_num)][c]<thresh, pvalues_dif['pvalues{}'.format(video_num)][c]<thresh)[idx] is True]
            #compare2 = [value for idx, value in enumerate(ISC2['ISC_persecond{}'.format(video_num)][c]) if np.bitwise_and(pvalues_drug2['pvalues{}'.format(video_num)][c]<thresh, pvalues_dif['pvalues{}'.format(video_num)][c]>1-thresh)[idx] is True]
            #print(compare1)
            #if len(compare1) != 0:
                #perc_drug1 += np.mean(compare1) #sum(pvalues_drug1['pvalues{}'.format(video_num)][c]<thresh)
            #if len(compare2) != 0:
                #perc_drug2 += np.mean(compare2) #sum(pvalues_drug2['pvalues{}'.format(video_num)][c]<thresh)
            #print([value for idx, value in enumerate(ISC1) if np.bitwise_and(pvalues_drug1['pvalues{}'.format(video_num)][c]<thresh, pvalues_dif['pvalues{}'.format(video_num)][c]<thresh)[idx] is True])
            
        store_perc_drug1.append(compare1); store_perc_drug2.append(compare2)
        store_perc_drug1_.append(compare1_); store_perc_drug2_.append(compare2_)
        
    print("on differences")
    
    for videos in [[0,1,10,11],[4,5,12,13],[6,7,18,19],[8,9,16,17],[2,3,14,15]]: # neutral, low-, low+, high-, high+
        
        #concat data for videos
        movie1 = np.array(store_perc_drug1[videos[0]] + store_perc_drug1[videos[1]] + store_perc_drug1[videos[2]] + store_perc_drug1[videos[3]]).flatten()
        movie1 = movie1[~np.isnan(movie1)]
        store_perc_drug1_avg.append(movie1)
        
        movie2 = np.array(store_perc_drug2[videos[0]] + store_perc_drug2[videos[1]] + store_perc_drug2[videos[2]] + store_perc_drug2[videos[3]]).flatten()
        movie2 = movie2[~np.isnan(movie2)]
        store_perc_drug2_avg.append(movie2)
    
    print(wilcoxon_paired(store_perc_drug1_avg[0], store_perc_drug2_avg[0]))
    print(wilcoxon_paired(store_perc_drug1_avg[1], store_perc_drug2_avg[1]))
    print(wilcoxon_paired(store_perc_drug1_avg[2], store_perc_drug2_avg[2]))
    print(wilcoxon_paired(store_perc_drug1_avg[3], store_perc_drug2_avg[3]))
    print(wilcoxon_paired(store_perc_drug1_avg[4], store_perc_drug2_avg[4]))
    
    
    print("au naturale")
    
    for videos in [[0,1,10,11],[4,5,12,13],[6,7,18,19],[8,9,16,17],[2,3,14,15]]: # neutral, low-, low+, high-, high+
        
        #concat data for videos
        movie1 = np.array(store_perc_drug1_[videos[0]] + store_perc_drug1_[videos[1]] + store_perc_drug1_[videos[2]] + store_perc_drug1_[videos[3]]).flatten()
        movie1 = movie1[~np.isnan(movie1)]
        store_perc_drug1_avg_.append(movie1)
        
        movie2 = np.array(store_perc_drug2_[videos[0]] + store_perc_drug2_[videos[1]] + store_perc_drug2_[videos[2]] + store_perc_drug2_[videos[3]]).flatten()
        movie2 = movie2[~np.isnan(movie2)]
        store_perc_drug2_avg_.append(movie2)
    
    print(wilcoxon_paired(store_perc_drug1_avg_[0], store_perc_drug2_avg_[0]))
    print(wilcoxon_paired(store_perc_drug1_avg_[1], store_perc_drug2_avg_[1]))
    print(wilcoxon_paired(store_perc_drug1_avg_[2], store_perc_drug2_avg_[2]))
    print(wilcoxon_paired(store_perc_drug1_avg_[3], store_perc_drug2_avg_[3]))
    print(wilcoxon_paired(store_perc_drug1_avg_[4], store_perc_drug2_avg_[4]))
    
    
    '''global df
    # prepare data for barplot
    df = pd.DataFrame(data = {"Drug1": [np.nanmean(x) for x in store_perc_drug1_avg], "Drug2": [np.nanmean(x) for x in store_perc_drug2_avg]}); df = df.stack().reset_index(); df.columns = ["Video", "Drug", "perc"]; df = df.explode("perc")
    
    # plots data
    plt.figure(1); rcParams['figure.figsize'] = 16,8
    df.index = np.rint((df.index/2)-0.1)
    ax = sns.barplot(x=df.index, y="perc", hue="Drug", data=df, color='gold', hue_order = ['Drug2', 'Drug1'], palette=sns.color_palette(['darkgray', 'gold']))
    
    # add significance annotations
    for video in range(1,6):
        pval = wilcoxon_paired(store_perc_drug1_avg[video-1], store_perc_drug2_avg[video-1])
        plt.text(video-1-0.17, 0.33,  str(eval_pval(pval)[0]) + str(eval_pval(pval)[1]), fontsize = 10)
        #sum Trues, 195*4
        plt.plot([video-1.1, video-1.1, video-1+0.1, video-1+0.1], [0.326, 0.328, 0.328, 0.326], lw=1, color='black')
    
    # make comparisons btw movies to report
    #for video1 in range(1,6):
        #for video2 in range(1,6):
            #print("{} x {}".format(xlabels[video1-1],xlabels[video2-1]))
            #print(eval_pval(prop_ztest(store_perc_drug1_avg[video1-1], 100, store_perc_drug1_avg[video2-1]/100))[0])
            #print(eval_pval(prop_ztest(store_perc_drug2_avg[video1-1], 100, store_perc_drug2_avg[video2-1]/100))[0])
    
    # make it pretty
    ax.set_xticklabels(xlabels, rotation=70)
    ax.set_xlabel('Movie', weight='bold')
    ax.set_ylabel('ISC \n (summed across components)', weight='bold', fontsize=12)
    ax.set(ylim=(0, 0.35))
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    plt.rcParams['axes.linewidth'] = 2
    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles[:2], ['PLC>OT', 'OT>PLC'], frameon=False)
    plt.show()'''

'-----------------------------------------------------------------------------'

# mean of significant time windows

'''def mean_perc_windows(ISC1, ISC2, pvalues_drug1, pvalues_drug2, thresh=0.05):
    
    #barplot of mean ISC in significantly different time windows of OT > OLC or OLC>OT,
    #for significant time points, per drug, per MOVIE CATEGORY
    
    xlabels = ['Neutral', 'Low -', 'Low +', 'High -', 'High +']
    
    # lazy inneficient way, change this whenever possible...
    
    # prepares data
    store_perc_drug1_avg = []; store_perc_drug2_avg = []
    store_perc_drug1 = []; store_perc_drug2 = []
    
    for video_num in range(1,21):
        #perc_drug1 = 0; perc_drug2 = 0
        compare1 = 0; compare2 = 0
        
        for c in range(0,3): 
            # get data
            
            
            
            compare1 += np.nanmean(np.where(pvalues_persecond_drug1['pvalues{}'.format(video_num)][c]<thresh, ISC_persecond_drug1['ISC_persecond{}'.format(video_num)]['ISC_persecond_drug1'][c], [np.nan]*195))
            compare2 += np.nanmean(np.where(pvalues_persecond_drug2['pvalues{}'.format(video_num)][c]<thresh, ISC_persecond_drug2['ISC_persecond{}'.format(video_num)]['ISC_persecond_drug2'][c], [np.nan]*195))

            
            
            
            #compare1 = [value for idx, value in enumerate(ISC1['ISC_persecond{}'.format(video_num)][c]) if np.bitwise_and(pvalues_drug1['pvalues{}'.format(video_num)][c]<thresh, pvalues_dif['pvalues{}'.format(video_num)][c]<thresh)[idx] is True]
            #compare2 = [value for idx, value in enumerate(ISC2['ISC_persecond{}'.format(video_num)][c]) if np.bitwise_and(pvalues_drug2['pvalues{}'.format(video_num)][c]<thresh, pvalues_dif['pvalues{}'.format(video_num)][c]>1-thresh)[idx] is True]
            #print(compare1)
            #if len(compare1) != 0:
                #perc_drug1 += np.mean(compare1) #sum(pvalues_drug1['pvalues{}'.format(video_num)][c]<thresh)
            #if len(compare2) != 0:
                #perc_drug2 += np.mean(compare2) #sum(pvalues_drug2['pvalues{}'.format(video_num)][c]<thresh)
            #print([value for idx, value in enumerate(ISC1) if np.bitwise_and(pvalues_drug1['pvalues{}'.format(video_num)][c]<thresh, pvalues_dif['pvalues{}'.format(video_num)][c]<thresh)[idx] is True])
            
        store_perc_drug1.append(compare1); store_perc_drug2.append(compare2)
    
    for videos in [[0,1,10,11],[4,5,12,13],[6,7,18,19],[8,9,16,17],[2,3,14,15]]: # neutral, low-, low+, high-, high+
        store_perc_drug1_avg.append( ((store_perc_drug1[videos[0]]) + (store_perc_drug1[videos[1]]) + (store_perc_drug1[videos[2]]) + (store_perc_drug1[videos[3]]))/4 )
        store_perc_drug2_avg.append( ((store_perc_drug2[videos[0]]) + (store_perc_drug2[videos[1]]) + (store_perc_drug2[videos[2]]) + (store_perc_drug2[videos[3]]))/4 )
    
    # prepare data for barplot
    df = pd.DataFrame(data = {"Drug1": store_perc_drug1_avg, "Drug2": store_perc_drug2_avg}); df = df.stack().reset_index(); df.columns = ["Video", "Drug", "perc"]; df = df.explode("perc")
    
    # plots data
    plt.figure(1); rcParams['figure.figsize'] = 16,8
    df.index = np.rint((df.index/2)-0.1)
    ax = sns.barplot(x=df.index, y="perc", hue="Drug", data=df, color='gold', hue_order = ['Drug2', 'Drug1'], palette=sns.color_palette(['darkgray', 'gold']))
    
    # add significance annotations
    for video in range(1,6):
        plt.text(video-1-0.06, 0.33,  eval_pval(prop_ztest(store_perc_drug1_avg[video-1], 100, store_perc_drug2_avg[video-1]/100))[0], fontsize = 10)
        #sum Trues, 195*4
        plt.plot([video-1.1, video-1.1, video-1+0.1, video-1+0.1], [0.326, 0.328, 0.328, 0.326], lw=1, color='black')
    
    # make comparisons btw movies to report
    for video1 in range(1,6):
        for video2 in range(1,6):
            print("{} x {}".format(xlabels[video1-1],xlabels[video2-1]))
            print(eval_pval(prop_ztest(store_perc_drug1_avg[video1-1], 100, store_perc_drug1_avg[video2-1]/100))[0])
            print(eval_pval(prop_ztest(store_perc_drug2_avg[video1-1], 100, store_perc_drug2_avg[video2-1]/100))[0])
    
    # make it pretty
    ax.set_xticklabels(xlabels, rotation=70)
    ax.set_xlabel('Movie', weight='bold')
    ax.set_ylabel('ISC \n (summed across components)', weight='bold', fontsize=12)
    ax.set(ylim=(0, 0.35))
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    plt.rcParams['axes.linewidth'] = 2
    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles[:2], ['PLC>OT', 'OT>PLC'], frameon=False)
    plt.show()'''




pvals_wsrt = []
corr_vals = []

# call vizs here

boxplot_persubject(ISC_persubject_sum_drug1, ISC_persubject_sum_drug2, chance_visible=True)  
#spatial_filter_eigval(ISC_spatial_filter_eig_drug1, ISC_spatial_filter_eig_drug2)
#timeseries_persecond(ISC_persecond_drug1, ISC_persecond_drug2, pvalues_persecond_drug1, pvalues_persecond_drug2, pvalues_dif, pvals_wsrt, corr_vals, plot_type="", thresh=0.05)

#perc_windows(ISC_persecond_drug1, ISC_persecond_drug2, pvalues_persecond_drug1, pvalues_persecond_drug2, thresh=0.05)
###mean_perc_windows(ISC_persecond_drug1, ISC_persecond_drug2, pvalues_persecond_drug1, pvalues_persecond_drug2, thresh=0.05)

#perc_windows_differences(ISC_persecond_drug1, ISC_persecond_drug2, pvalues_persecond_drug1, pvalues_persecond_drug2, thresh=0.05)
#mean_perc_windows_differences(ISC_persecond_drug1, ISC_persecond_drug2, pvalues_persecond_drug1, pvalues_persecond_drug2, pvalues_dif, thresh=0.05)

#Topoplots:
#Faz paiwise correlations between ALL MOVIES; produce figure matrix movies x components x drug for supplementary data 


'-----------------------------------------------------------------------------'
#pearson corr scalp
#harry = np.array([scalp_drug2['scalp{}'.format(video)][1] for video in range(1,21)])
#np.corrcoef(harry)

'-----------------------------------------------------------------------------'

# % significant time windows for drug1 and drug2, compare both; show % significant
# time windows overlapping btw drugs, just to show that it's not too small to calculate average'
# average ISC during significant time windows for drug1 AND drug2, simultenously'''

'''jui = pd.DataFrame(list(ISC_persubject_sum_drug1.items()))
jui = jui.transpose()
jui.columns = jui.iloc[0]
jui = jui.iloc[1: , :]'''
#jui = jui.apply(lambda x: pd.Series(x['samples']),axis=1).stack().reset_index(level=1, drop=True)
