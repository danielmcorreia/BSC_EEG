# EEG-ISC_analysis
 
**INSTRUCTIONS**

To run this analysis pipeline:

**1. RUN STEP1_PREPROCESS**

   STEP1_PREPROCESS performs all the necessary steps for the preprocessing of EEG data, namely:
   - Applies high-pass filter to the data (1Hz);
   - Applies notch filter to the data (50Hz, harmonic);
   - Fits ICA to the data, in order to remove eye-blink and eye-movement related artifacts (requires manual selection of corresponding componentes);
   - Calls segmenting function (see below) and applies it to the data. At this step, the data is divided into 20 folders (1 for each analysed video);
   - Downsamples the data to 250Hz;
   - Calls bad_channel identification function (see below) and applies it to the data;
   - Calls remove_outliers function (see below) and applies it to the data.
   
   In between each major step (ICA, segmentation, bad_channel removal and outlier removal), checkpoints are performed and data is saved for further pick up without the need of running all the analysis from the beginning again.

--------------------------------------------------------------------------------------------------------------------------------------------------------------------
 
**2. RUN STEP2_STRUCT**

   This script defines the distributions of subjects by drugs and manipulates the data so that it is ready for calculating ISC (either with MATLAB or PYTHON code, see below).

--------------------------------------------------------------------------------------------------------------------------------------------------------------------

**3. RUN ISC_calculation** (only applicable if you want to calculate ISC using MATLAB code, otherwise proceed to step 4)
 
--------------------------------------------------------------------------------------------------------------------------------------------------------------------
   
**4. RUN ISC_visualization**

   Reads data outputted from ISC calculation (if calcualted in MATLAB) and creates visualization to observe results; if ISC is meant to be calculated using the Python code, it calls ISC_calculation.py, calculates it and outputs results;

--------------------------------------------------------------------------------------------------------------------------------------------------------------------

**Note:**
Make sure to follow the guidelines (found at the top of each of the scripts, if applicable), regarding the organization of import and export paths required by each script.

--------------------------------------------------------------------------------------------------------------------------------------------------------------------

Other necessary scripts for the analysis (though no interaction with them is necessary):

- DEF_PREPROCESSSING: Defines major functions necessary for data preprocessing, namely:
   - run_ICA: fits ICA to data; waits for manual input of eye-artifact related components to remove; removes them from data and returns it;
   - remove_bad_channels: identifies bad channels based on average power exceeding mean channel power by a multiple of interquartile differences;
   - remove_outliers: removes outliers above defined amplitude threshold of thresh_stds * std signal and 100ms before/after;
- SEGMENTATION: Segments the data by identifying TTL code related to movie beginning and takes out 40 sec after that time code; returns segmented data;
- ISC_calculation.py: Python implementation of ISC (inspired by code from the Parra Lab)
- ISC_calculation.mat: MATLAB implementation of ISC (inspired by code from the Parra Lab)

