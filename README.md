# Multiresolution Alignment for Multiple Unsynchronized Audio Sequences using Sequential Monte Carlo Samplers

This is the distributed software for the article "Multiresolution Alignment for Multiple Unsynchronized Audio Sequences using Sequential Monte Carlo Samplers" submitted to Digital Signal Processing - SoftwareX joint special issue on Reproducible Research in signal processing.

## How to use the software

There are 3 separate parts of the software; Multiresolution Alignment, Baseline - Fingerprinting based alignment and Evaluation.

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

The results text file is available under the folder "Evaluation/SMC_offset_estimation_results". The name convention is;

	offset_estimation_SMC_result_<dd>_<mm>_<yyyy>_<hh>h_<mm>m.txt

As an example,

	offset_estimation_SMC_result_16_12_2016_13h_26m.txt

The detailed explanations and block diagrams of the main modules 'feature_extract_module', 'SMC_main_module', 'sequential_alignment_module' and 'SMC_core_module' can be found in the SoftwareX manuscript. 

### Baseline - Fingerprinting based alignment

The fingerprinting software is obtained from 

D. Ellis, Robust Landmark based Audio Fingerprinting, web source, available (2009).
URL http://labrosa.ee.columbia.edu/matlab/fingerprint/

The codes are available under the folder "Fingerprinting_based_alignment".

This software is not directly applicable to the alignment setting since it is a query-by-example based audio fingerprinting software. Here, we simply try to match all pairs of sequences from the number of exactly matching ngerprints and compute the relative offset from the time
information of  fingerprints. A more detailed description can be found in Sec.4.1 in DSP manuscript. 

For demonstration of the fingerprinting based alignment software, please run "fingerprint_based_audio_alignment.m" file. Here, the software takes the audio files for GT_090912 event of the Jiku dataset as input and the connected pair of sequences are estimated with their relative offsets using a threshold for the exact matching hashes. The threshold value is searched with a grid search method. The results are written to a separate text file for each threshold consisting each connected pair in a format 

	<Sequence 1> <Sequence 2> <Relative offset>

The results text file is available under the folder 
"Evaluation/fingerprinting_offset_estimation_results". The name convention is;

	offset_estimation_fingerprinting_thr_<value>_result.txt

As an example,

	offset_estimation_fingerprinting_thr_20_result.txt


### Evaluation 

The evaluation software is written in python 2.7 and contains 3 files "groundtruth.py", "compute_accuracy.py" and "fingerprinting_evaluation.py" available under the folder "Evaluation".

1- groundtruth.py

The ground-truth for the Jiku dataset is given in the "Jiku_GT_090912.xml" file available under the folder "Evaluation/ground_truth". The ground-truth is obtained from

M. Guggenberger, M. Lux, L. Boszormenyi, A Synchronization Ground
Truth for the Jiku Mobile Video Dataset 

"Jiku_GT_090912.xml" file contains the offset of each sequence on the universal time line that is not compatible with the evaluation procedure as in Sec. 4.2 in DSP Manuscript. We reformat the ground-truth information with the "groundtruth.py" file into a text file "ground_truth.txt" available under "Evaluation/ground_truth". 

Note that "groundtruth.py" file have to be run for one time. Once the "ground_truth.txt" file is created, there is no need to run it again. 

2- compute_accuracy.py

When an experiment is conducted with the multiresolution alignment software, the results are written into a text file and saved under the folder "Evaluation/SMC_offset_estimation_results". Similarly when an experiment is conducted with the baseline method, the results are written into a text file and saved under the folder "Evaluation/fingerprinting_offset_estimation_results".  

"compute_accuracy.py" file computes the accuracy of an alignment estimate result. It requires the ground-truth information from "ground_truth.txt" file and the text file that contains the resulting alignment estimates.  

The variable 'offset_estimation_result_filename' has to be set to the name of the results text file, as an example ;

	offset_estimation_result_filename = 'offset_estimation_SMC_result_16_11_2016.txt'

3- fingerprinting_evaluation.py

The baseline method requires a threshold to decide a matching/not matching decision between two sequences. A grid search is applied to tune the threshold for best accuracy result. The alignment estimation results are computed for thresholds {10 , 20 , ... , 150} with the baseline. "fingerprinting_evaluation.py" file computes accuracy for each estimation result using the "compute_accuracy.py", and prints the best accuracy and the threshold value. It also plots a figure with two subplots; the accuracy for each threshold, FP, FN_1 and FN_2 values for each threshold. 

## For reviewers

To test all the aspects provided by software, one should follow the steps below; 

1- Run "/Evaluation/groundtruth.py" 

2- Run "/SMC_based_alignment/SMC_demonstration.m"

3- Run "/Evaluation/compute_accuracy.py" by first setting the filename as the most recent files name,

i.e., offset_estimation_SMC_result_16_11_2016_13h_26m.txt

4- Run "/Fingerprinting_based_alignment/fingerprint_based_audio_alignment.m"

5- Run "/Evaluation/fingerprinting_evaluation.py"

Note that "ground_truth.txt" and baseline results for each threshold are already available in their respective paths. One can also skip the steps 1 and 4, and apply the rest of the steps in the same order for a shorter test. 


