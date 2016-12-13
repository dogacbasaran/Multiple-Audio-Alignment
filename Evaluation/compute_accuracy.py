""" 
compute_accuracy.py: This is the main file that computes the accuracy of an alignment
estimate. 

How to use:
Manually write the name of the result text file from the result folder to the variable
"offset_estimation_result_filename"
Example:
    offset_estimation_result_filename = 'offset_estimation_SMC_result_16_11_2016.txt'

The code can work under Windows or Linus based OS. 

Copyright (C) 2016  Dogac Basaran

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>
"""


import numpy as np
import json
import os
from os import listdir

# Sets the key name for ground_truth dictionary entry, order the sequences in 
# ascending order of microphone number and record number
# Example: If sequence = 'mic4_rec2.wav', sequence_ = 'mic1_rec3.wav'
#           the key_name = 'mic1_rec3.wav mic4_rec2.wav'
# Input: sequence - First filename
#        sequence_ - Second filename
# Output: key_name - Either 'filename filename_' or 'filename_ filename'
def set_key_name(sequence, sequence_):
    mic_number = int(sequence[sequence.find('mic')+3:sequence.find('_')])
    mic_number_ = int(sequence_[sequence_.find('mic')+3:sequence_.find('_')])
    rec_number = int(sequence[sequence.find('rec')+3:sequence.find('.wav')])
    rec_number_ = int(sequence_[sequence_.find('rec')+3:sequence_.find('.wav')])
    
    if mic_number > mic_number_:
        key_name = sequence_ + ' ' + sequence
    elif mic_number == mic_number_:
        if rec_number > rec_number_:
            key_name = sequence_ + ' ' + sequence
        else:
            key_name = sequence + ' ' + sequence_
    else:
        key_name = sequence + ' ' + sequence_
    return key_name
    
def set_relative_offset(key_name, sequence, relative_offset_):
    if key_name.find(sequence)==0:
        relative_offset_distance = relative_offset_
    else:
        relative_offset_distance = -relative_offset_
    return relative_offset_distance
    
def extract_estimated_pairs(path, offset_estimation_result_filename, coeff=1.):
    estimations = {}
    # Start reading the estimation results 
    f = open(path + '\\' + offset_estimation_result_filename,'r') 
    for line in f.readlines():
        estimated_sequences = line.split(" ")
        sequence = estimated_sequences[0] # The first sequence is current sequence
        if len(estimated_sequences) > 2: # This means no aligned sequences with the current sequence
            # For each estimated sequence
            for sequence_relative_offset in estimated_sequences[1:-1]:
                # Estimated sequence and its relative offset to current sequence
                sequence_ = sequence_relative_offset.split(',')[0]
                relative_offset_ = coeff * float(sequence_relative_offset.split(',')[1])
                key_name = set_key_name(sequence, sequence_)
                relative_offset_ = set_relative_offset(key_name, sequence, relative_offset_)
                if estimations.has_key(key_name)==False:
                    estimations[key_name] = relative_offset_
    f.close()
    return estimations

def compute_accuracy(path, offset_estimation_result_filename, print_results = False):
    path1 = path[0]
    path2 = path[1]
    path3 = path[2]
    path4 = path[3]

    FN_1 = 0. # False Negative Type 1
    FN_2 = 0. # False Negative Type 2
    FP = 0. # False Positive
    
    # Read precomputed ground truth dictionary
    if path1.find('/')!=-1: # Windows based
        ground_truth = json.load(file(path1 + '\\' + 'ground_truth.txt')) 
    else: # Linux based
        ground_truth = json.load(file(path1 + '/' + 'ground_truth.txt')) 
    
    # Extract estimated pairs dictionary
    if offset_estimation_result_filename.find('fingerprint') == -1: 
        estimations = extract_estimated_pairs(path3, offset_estimation_result_filename)
    else:
        estimations = extract_estimated_pairs(path4, offset_estimation_result_filename)
        
    # Scan estimations to find FN_2 and FP types of errors
    tolerance = 10 # The alignment is acceptable in neighborhood of ground truth with amount of tolerance
    for key, value in estimations.iteritems():
        if ground_truth.has_key(key):
            true_relative_offset = ground_truth[key]
            estimated_relative_offset = value
            #print(key + ' ' + np.str(true_relative_offset) + ' ' + np.str(estimated_relative_offset) + '   ' + np.str(np.abs(true_relative_offset-estimated_relative_offset)))
            if np.abs(true_relative_offset - estimated_relative_offset) > tolerance:
                FN_2+=1.            
        else:
            FP+=1.
    
    # Scan ground_truth to find FN_1 type of error
    for key in ground_truth.iterkeys():
        if estimations.has_key(key) == False:
            FN_1+=1.
    
    # Find number of recordings for each microphone
    number_of_microphones = 4
    number_of_recordings = np.zeros(number_of_microphones,dtype='float64')
    for filename in listdir(path2):
        mic_number = int(filename[filename.find('mic')+3:filename.find('_')])
        rec_number = int(filename[filename.find('rec')+3:filename.find('.wav')])
        if number_of_recordings[mic_number-1]<rec_number:
            number_of_recordings[mic_number-1] = rec_number
    
    # Compute the number of included pairs as total number - excluded number
    # Note that we assume that the recordings from the same microphone do not
    # overlap in the algorithm, hence remove those pairs from the total number
    K = np.sum(number_of_recordings)
    total_number_of_pairs = K * (K-1)/2
    excluded_number_of_pairs = 0                 
    for i in range(len(number_of_recordings)):
        K_tmp = number_of_recordings[i]
        excluded_number_of_pairs += K_tmp * (K_tmp-1)/2            
    
    # Compute the final evaluation score OMEGA
    Omega = 1. - (FN_1 + FN_2 + FP)/(total_number_of_pairs-excluded_number_of_pairs)
    if print_results == True:
        print("\nThe evaluation results for " + offset_estimation_result_filename)
        print(('False Negative Type 1 - FN_1 = {0}').format(FN_1))
        print(('False Negative Type 2 - FN_2 = {0}').format(FN_2))
        print(('False Positive - FP = {0}').format(FP))
        print(('Omega = {0}').format(Omega))
    
    return (Omega, FN_1, FN_2, FP)  

if __name__ == '__main__':
    cw_path = os.getcwd();
    if cw_path.find('/')==-1:
        cw_path_parent = cw_path[:cw_path.find('\\Evaluation')]
        path1 = cw_path + '\\ground_truth'
        path2 = cw_path_parent + '\\audio_data'
        path3 = cw_path + '\\SMC_offset_estimation_results'
        path4 = cw_path + '\\fingerprinting_offset_estimation_results'
    else:
        cw_path_parent = cw_path[:cw_path.find('/Evaluation')]
        path1 = cw_path + 'ground_truth'
        path2 = cw_path_parent + '/audio_data'
        path3 = cw_path + 'SMC_offset_estimation_results'
        path4 = cw_path + 'fingerprinting_offset_estimation_results'

    path = [path1, path2, path3, path4]

    #offset_estimation_result_filename = 'offset_estimation_SMC_result_16_11_2016.txt'
    
    offset_estimation_result_filename = 'offset_estimation_fingerprinting_thr_20_result_14_11_2016.txt'
    
    Omega, FN_1, FN_2, FP = compute_accuracy(path, offset_estimation_result_filename, True) 
    