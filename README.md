# Multiresolution Alignment for Multiple Unsynchronized Audio Sequences using Sequential Monte Carlo Samplers

This is the distributed software for the article "Multiresolution Alignment for Multiple Unsynchronized Audio Sequences using Sequential Monte Carlo Samplers" submitted to Digital Signal Processing - SoftwareX joint special issue on Reproducible Research in signal processing.

## How to use the software

There are 3 separate parts of the software; 

1- Multiresolution multiple audio alignment software. 

2- Fingerprinting based multiple audio alignment software (Baseline).

3- Evaluation.

Here, we demostrate the software with using the audio dataset from GT_090912 event of the Jiku dataset (available [here](https://www.dropbox.com/sh/ktirf3t6b8lhs7d/AADlCm24Xw2A_qru5sUP71jFa?dl=0)). Note that the audio files have a name convention for microphones with multiple recordings. In the provided link, the audio filenames are modified accordingly.

### Multiresolution Multiple Audio Alignment Software

This is the main software that computes the alignment estimates of unsynchronized audio files using the SMC based multiresolution multiple audio alignment method. The software consists of 5 Matlab files and is located under the folder "SMC_based_alignment". The documentation of the software is available [here](www.dogacbasaran.com/Software\_documentation/SMC\_based\_alignment\_documentation/index.html)
 
For the demonstration of the software, please apply the following steps,

1- Download and decompress Jiku dataset with the provided download link.

2- Download and decompress the project from the github repo
    
3- To run SMC based multiresolution multiple audio alignment software, simply run 

	/path/to/project/SMC\_based\_alignment/SMC\_demonstration.m 

and choose /path/to/audio_data as input in the browse menu. The resulting alignment estimates will be written in a text file under 

	/path/to/project/Evaluation/SMC_offset_estimation_results 

with the name convention,

	offset_estimation_SMC_result_<dd>_<mm>_<yyyy>_<hh>h_<mm>m.txt

4-  To run the fingerprinting based multiple audio alignment system, simply run
 
	/path/to/project/Fingerprinting\_based\_alignment/fingerprinting\_based\_audio\_alignment.m

and choose /path/to/audio_data as input in the browse menu. The resulting alignment estimates will be written in a text file under 

	/path/to/project/Evaluation/fingerprinting_offset_estimation_results 

with the name convention,

	offset_estimation_fingerprinting_thr_<value>_result.txt

5- To evaluate any result, simply run

	/path/to/project/Evaluation/compute\_accuracy.py

and choose /path/to/alignment_result to see the accuracy, precision, recall and F-measure scores.

For the demonstration of the software, please run "SMC_demonstration.m" file and choose the path to the downloaded and decompressed audio dataset as input.

 Once the user obtained the dataset, the user can hence any input audio dataset should be  and returns the list of connected sequences (Clusters), their relative offset information (r_clusters) and the elapsed time information (time_elapsed). The results are written to a text file for each connected pair separately in a format 

	<Sequence 1> <Sequence 2> <Relative offset>

The results text file is available under the folder "Evaluation/SMC_offset_estimation_results". The name convention is;

	offset_estimation_SMC_result_<dd>_<mm>_<yyyy>_<hh>h_<mm>m.txt

As an example,

	offset_estimation_SMC_result_16_12_2016_13h_26m.txt

The detailed explanations and block diagrams of the main modules 'feature_extract_module', 'SMC_main_module', 'sequential_alignment_module' and 'SMC_core_module' can be found in the SoftwareX manuscript. 

### Fingerprinting Based Multiple Audio Alignment Software (Baseline)

The fingerprinting software is obtained from 

D. Ellis, Robust Landmark based Audio Fingerprinting, web source, available (2009).
URL http://labrosa.ee.columbia.edu/matlab/fingerprint/

The codes for fingerprinting software are available under the subfolder "Fingerprinting_based_alignment/fingerprint_labrosa".

This software is not directly applicable to the alignment setting since it is a query-by-example based audio fingerprinting software. Here, we simply try to match all pairs of sequences from the number of exactly matching fingerprints and compute the relative offset from the time
information of fingerprints. A more detailed description can be found in Sec.4.1 in DSP manuscript. 

For demonstration of the fingerprinting based alignment software, please run "Fingerprinting_based_alignment/fingerprint_based_audio_alignment.m" file. Here, the software takes the audio files for GT_090912 event of the Jiku dataset as input and the connected pair of sequences are estimated with their relative offsets using a threshold for the exact matching hashes. The threshold value is searched with a grid search method. The results are written to a separate text file for each threshold consisting each connected pair in a format 

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

The baseline method requires a threshold to decide a matching/not matching decision between two sequences. A grid search is applied to tune the threshold for best accuracy result. The alignment estimation results are computed for thresholds {10 , 20 , ... , 150} with the baseline. "fingerprinting_evaluation.py" file computes accuracy for each estimation result using the "compute_accuracy.py", and prints the best accuracy and the threshold value. It also plots a figure with two subplots; the accuracy for each threshold and precision, recall and F-measure values for each threshold. 

## For reviewers

To test all the aspects provided by software, one should follow the steps below; 

1- Run "/Evaluation/groundtruth.py" 

2- Run "/SMC_based_alignment/SMC_demonstration.m"

3- Run "/Evaluation/compute_accuracy.py" by first setting the filename as the most recent files name,

i.e., offset_estimation_SMC_result_16_11_2016_13h_26m.txt

4- Run "/Fingerprinting_based_alignment/fingerprint_based_audio_alignment.m"

5- Run "/Evaluation/fingerprinting_evaluation.py"

Note that "ground_truth.txt" and baseline results for each threshold are already available in their respective paths. One can also skip the steps 1 and 4, and apply the rest of the steps in the same order for a shorter test. 


