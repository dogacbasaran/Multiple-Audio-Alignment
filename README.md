# Multiresolution Alignment for Multiple Unsynchronized Audio Sequences using Sequential Monte Carlo Samplers

## How to use the software

There are 3 separate parts of the software; Multiresolution Alignment, Fingerprinting based alignment and Evaluation.

### Multiresolution Alignment Software

The software consists of 5 Matlab files and is located under the folder "SMC_based_alignment".

To run the software first provide the path of the input dataset to feature_extract_module;

	dataset_features = feature_extract_module(<Path to the dataset>)

Then the resulting struct 'dataset_features' is fed into the SMC_main_module that finds the best alignment estimates using the multiresolution alignment scheme proposed in DSP Manuscript,
	
	[Clusters, r_clusters, time_elapsed] = SMC_main_module(dataset_features)

'Clusters' is a cell array and a non-empty cell with index 'i' holds the list of connected sequences of the cluster 'i'. 

'r_clusters' is a cell array and a non-empty cell with index 'i' holds the relative offsets of the connected sequences in cluster 'i'

'time_elapsed' variable holds the elapsed time in seconds for the procedure to be completed i.e., all the sequences are aligned.
 
For the demonstration of the software, please run "SMC_demonstration.m" file. 
Here, the software takes the audio files for GT_090912 event of the Jiku dataset as input (available under folder "audio_data") and returns the list of connected sequences (Clusters), their relative offset information (r_clusters) and the elapsed time information (time_elapsed). The results are written to a text file for each connected pair separately in a format 

	<Sequence 1> <Sequence 2> <Relative offset>

The results text file is available under the folder "Evaluation/SMC_offset_estimation_results".

The detailed explanations and block diagrams of the main modules 'feature_extract_module', 'SMC_main_module', 'sequential_alignment_module' and 'SMC_core_module' can be found in the SoftwareX manuscript. 

### Fingerprinting based alignment

The fingerprinting software is obtained from 

D. Ellis, Robust Landmark based Audio Fingerprinting, web source, available (2009).
URL http://labrosa.ee.columbia.edu/matlab/fingerprint/

The codes are available under the folder "Fingerprinting_based_alignment".

This software is not directly applicable to the alignment setting since it is a query-by-example based audio fingerprinting software. Here, we simply try to match all pairs of sequences from the number of exactly matching ngerprints and compute the relative offset from the time
information of  fingerprints. A more detailed description can be found in Sec.4.1 in DSP manuscript. 

For demonstration of the fingerprinting based alignment software, please run "fingerprint_based_audio_alignment.m" file. Here, the software takes the audio files for GT_090912 event of the Jiku dataset as input and the connected pair of sequences are estimated with their relative offsets using a threshold for the exact matching hashes. The threshold value is searched with a grid search method. The results are written to a separate text file for each threshold consisting each connected pair in a format 

	<Sequence 1> <Sequence 2> <Relative offset>

The results text file is available under the folder 
"Evaluation/fingerprinting_offset_estimation_results".

### Evaluation 

Contains the python codes for the evaluation of accuracy for both multiresolution alignment and baseline methods.

This folder contains three subfolders:  

-SMC_offset_estimation_results: Contains the resulting alignment estimates of the multiresolution alignment system. 

-fingerprinting_offset_estimation_results: Contains the results alignment estimates fo the baseline system. There are several results for different thresholds that determines matching sequences.

-ground_truth: Contains provided ground-truth file "Jiku_GT_090912.xml", the "groundtruth.txt" that provides the ground-truth matching sequences and their offsets at each row.



