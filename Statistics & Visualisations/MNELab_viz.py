#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 17 10:38:19 2021

@author: Daniel
"""

#from copy import deepcopy
#import mne
from mnelab.io import read_raw
import os
from matplotlib import pyplot as plt


drug = 1
video = 9

path_import = '/Volumes/LaCie Armaz/DPrataLab/ISC/PREPROCESSED/DRUG{}/VIDEO{}/'.format(drug, video)

subjects_list = [file for file in os.listdir(path_import) if (file.endswith('.fif') and not file.startswith("."))]

datasets = []

for subject in range(10): 
    data = read_raw(path_import+'/'+subjects_list[subject], preload=True)
    datasets.insert(subject, data)

for subject in datasets:
    subject.plot(n_channels=62)
    figManager = plt.get_current_fig_manager()
    figManager.window.showMaximized()

