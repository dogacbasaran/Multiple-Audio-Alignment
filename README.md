# Multiresolution Alignment for Multiple Unsynchronized Audio Sequences using Sequential Monte Carlo Samplers

This software is an implementation of 

## Folder Structure

There are 4 main folders;

### audio_data

Contains the audio files to be synchronized.

### SMC_based_alignment

Contains the Matlab codes for the multiresolution audio alignment software. 

### Fingerprinting_based_alignment

Contains the Matlab codes for the baseline-fingerprint based audio alignment software.

### Evaluation 

Contains the python codes for the evaluation of accuracy for both multiresolution alignment and baseline methods.

This folder contains three subfolders:  

-SMC_offset_estimation_results: Contains the resulting alignment estimates of the multiresolution alignment system. 

-fingerprinting_offset_estimation_results: Contains the results alignment estimates fo the baseline system. There are several results for different thresholds that determines matching sequences.

-ground_truth: Contains provided ground-truth file "Jiku_GT_090912.xml", the "groundtruth.txt" that provides the ground-truth matching sequences and their offsets at each row.



