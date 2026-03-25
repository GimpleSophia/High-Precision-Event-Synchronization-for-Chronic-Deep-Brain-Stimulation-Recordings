## Code for paper: High Precision Event Synchronization for Chronic Deep Brain Stimulation Local Field Potential Recordings
In the paper "High-Precision Event Synchronization for Chronic Deep Brain Stimulation Local Field Potential Recordings" we show and evaluate the synchronization of Medtronic Percept DBS local field potential (LFP) recordings to EEG recordings using Transcutaneous subcutaneous stimulation (TES) to induce matching artifacts in both recordings.

This repository contains 3 folders:
- Python_synchronization_functions: functions for loading, preprocessing and synchronizing LFP recordings to EEG recordings based on TES artefacts in Python
- Matlab_synchronization_funcitons: functions for synchronizing LFP and EEG recordings based on TES artefacts in Matlab, for loading and preprocessing of the LFP data from   ".json" files we recommend using an already avaliable toolbox e.g. https://github.com/jgvhabets/PyPerceive
- Paper_analysis_code: Code for replicating our analysis from the paper

One example recording for testing the synchronization is provided via the Open Science Framework (OSF) repository (https://osf.io/dnvgb)


