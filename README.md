# Multiresolution Alignment for Multiple Unsynchronized Audio Sequences using Sequential Monte Carlo Samplers

This is the distributed software for the article "Multiresolution Alignment for Multiple Unsynchronized Audio Sequences using Sequential Monte Carlo Samplers" submitted to Digital Signal Processing - SoftwareX joint special issue on Reproducible Research in signal processing.

There are 3 separate parts of the software; 

1- Multiresolution multiple audio alignment software. 

2- Fingerprinting based multiple audio alignment software (Baseline).

3- Evaluation.

Here, we demostrate the software with using the audio dataset from GT_090912 event of the Jiku dataset (available [here](https://www.dropbox.com/sh/ktirf3t6b8lhs7d/AADlCm24Xw2A_qru5sUP71jFa?dl=0)). Note that the audio files have a name convention for microphones with multiple recordings. In the provided link, the audio filenames are modified accordingly.

## Multiresolution Multiple Audio Alignment Software

This is the main software that computes the alignment estimates of unsynchronized audio files using the SMC based multiresolution multiple audio alignment method. The software is written in Matlab and is located under the folder "/SMC_based_alignment". The documentation of the software is available [here](http://www.dogacbasaran.com/Software_documentation/SMC_based_alignment_documentation/index.html).
 
## Fingerprinting Based Multiple Audio Alignment Software (Baseline)

As a baseline, we use a fingerprinting based alignment approach. The fingerprinting software is obtained from 

D. Ellis, Robust Landmark based Audio Fingerprinting, web source, available (2009).
URL http://labrosa.ee.columbia.edu/matlab/fingerprint/

The codes for fingerprinting software are available under the subfolder "/Fingerprinting_based_alignment/fingerprint_labrosa".

Note that this software is not directly applicable to the alignment setting since it is a query-by-example based audio fingerprinting software. For alignment purposes, we simply count the number of exact hash(fingerprint) matches between each pair of sequences. Then by thresholding according to the number of hash matches, we decide if the sequences are matching. The time information of the matching hashes are then used to compute the relative offset between sequences. A more detailed description can be found in Sec.4.1 in the DSP manuscript.  

## Evaluation 

The evaluation software is written in python 2.7 and is located under the folder "/Evaluation". The documentation of the software is available [here](http://www.dogacbasaran.com/Software_documentation/Evaluation_documentation/index.html).

The ground-truth for the Jiku dataset is given in the "Jiku_GT_090912.xml" file available under the subfolder "/Evaluation/ground_truth". The ground-truth is obtained from

M. Guggenberger, M. Lux, L. Boszormenyi, A Synchronization Ground
Truth for the Jiku Mobile Video Dataset 

"Jiku_GT_090912.xml" file contains the offset of each sequence on the universal time line that is not compatible with the evaluation procedure where the relative offsets are considered. Hence, we reformat the ground-truth information with the "groundtruth.py" file into a text file "ground_truth.txt" available under "Evaluation/ground_truth". 

Note that "groundtruth.py" file have to be run for one time. Once the "ground_truth.txt" file is created, there is no need to run it again. 

The baseline method requires a threshold to decide a matching/not matching decision between two sequences. A grid search is applied to tune the threshold for best accuracy result using the "fingerprinting_evaluation.py" module. The threshold that results in highest accuracy is chosen for comparison. The software prints the best accuracy with the respective threshold value and plots a figure with two subplots; the accuracy for each threshold in the first plot and precision, recall and F-measure values for each threshold in the second plot. 

## How to use the software

For the demonstration of the software, please apply the following steps,

1- Download and decompress Jiku dataset with the provided download link.

2- Download and decompress the project from the github repo
    
3- To run SMC based multiresolution multiple audio alignment software, simply run 

	/path/to/project/SMC_based_alignment/SMC_demonstration.m 

and choose /path/to/audio_data as input in the browse menu. The resulting alignment estimates will be written in a text file under 

	/path/to/project/Evaluation/SMC_offset_estimation_results 

with the name convention,

	offset_estimation_SMC_result_<dd>_<mm>_<yyyy>_<hh>h_<mm>m.txt

4- To evaluate the estimation results, simply run

	/path/to/project/Evaluation/compute_accuracy.py

and choose /path/to/project/SMC_offset_estimation_results/result_file to compute the accuracy, precision, recall and F-measure scores.

4-  To run the fingerprinting based multiple audio alignment system, simply run
 
	/path/to/project/Fingerprinting_based_alignment/fingerprinting_based_audio_alignment.m

and choose /path/to/audio_data as input in the browse menu. The estimation results written to a separate text file for each threshold under

	/path/to/project/Evaluation/fingerprinting_offset_estimation_results 

with the name convention,

	offset_estimation_fingerprinting_thr_<value>_result.txt

5- To evaluate the estimation results, simply run

	/path/to/project/Evaluation/fingerprinting_evaluation.py
	
the best threshold (highest accuracy) will be computed and the respective evalutaion metrics will be printed on the screen.



