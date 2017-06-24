""" 
compute_accuracy.py
~~~~~~~~~~~~~~~~~~~ 
.. topic:: Content

    This is the main file that computes the accuracy of an alignment
    estimate. 

    How to use:
        Manually write the name of the result text file from the result folder to the variable
        "offset_estimation_result_filename"
    Example:
        offset_estimation_result_filename = 'offset_estimation_SMC_result_16_11_2016_13h_26m.txt'

    Copyright (C) 2016  Author: Dogac Basaran

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>"""


import numpy as np
import json
import os
from os import listdir
import Tkinter
import tkFileDialog

def set_key_name(sequence1, sequence2):
    
    """ Sets the key name for ground_truth dictionary entry, order the sequences in 
        ascending order of microphone number and record number
                 
    **Parameters**
    
    sequence1: String
        The name of the first file
    sequence2: String
        The name of the second file
        
    **Returns**
    
    key_name: String
        Key name for the ground_truth dictionary."""
        
    mic_number = int(sequence1[sequence1.find('mic')+3:sequence1.find('_')])
    mic_number_ = int(sequence2[sequence2.find('mic')+3:sequence2.find('_')])
    rec_number = int(sequence1[sequence1.find('rec')+3:sequence1.find('.wav')])
    rec_number_ = int(sequence2[sequence2.find('rec')+3:sequence2.find('.wav')])
    
    if mic_number > mic_number_:
        key_name = sequence2 + ' ' + sequence1
    elif mic_number == mic_number_:
        if rec_number > rec_number_:
            key_name = sequence2 + ' ' + sequence1
        else:
            key_name = sequence1 + ' ' + sequence2
    else:
        key_name = sequence1 + ' ' + sequence2
    return key_name
    
def set_relative_offset(key_name, sequence, relative_offset):
    
    """Sets the estimated relative offset between two sequences according to the
       ordering of the sequences
       
   **Parameters**
   
   key_name: String
       The key name that is consistent with the ground_truth dictionary
   sequence: String
       The sequence name
   relative_offset: float 
       Estimated relative offset
   
   **Returns**
   
   relative_offset_distance: float
       Relative offset according to the ordering of the aligned sequences"""
           
    if key_name.find(sequence)==0:
        relative_offset_distance = relative_offset
    else:
        relative_offset_distance = -relative_offset
    return relative_offset_distance
    
def extract_estimated_pairs(offset_estimation_result_filename, coeff=1.):
    
    """ Gets the estimated alignments from a txt file and extracts all the pair of alignment in format convenient with ground_truth
    dictionary
        
    **Parameters**
    
    offset_estimation_result_filename: String
        Name of the estimation file with full path
    coeff: Float
            
    
    **Returns**
    
    estimations: Dictionary
            The dictionary containing the estimated pair of alignments"""
            
    estimations = {}
    
    # Start reading the estimation results 
    f = open(offset_estimation_result_filename,'r')

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

def compute_accuracy(path, offset_estimation_result_filename, verbose = False):
    
    """Computes the main evaluation metrics namely Accuracy, Precision, Recall and F-measure
      
    **Parameters**
    
    path: List
        List contains the path to the ground truth and path to the audio dataset
    offset_estimation_result_filename: String
        Name of the estimation file
    verbose: Boolean (default False)
        Prints the resulting metrics as well as the number of TP, TN, FP and FN 
        if it is True
        
    **Returns**
    
    Accuracy: Float
        The accuracy of the estimation
    Precision: Float
        The precision of the estimation
    Recall: Float
        The recall of the estimation
    F_measure: Float 
        The F-measure of the estimation
    TP: Integer 
        The number of the true positives
    TN: Integer 
        The number of the true negatives
    FP: Integer 
        The number of the false positives
    FN: Integer
        The number of the false negatives"""
    
    path_ground_truth = path[0]
    path_audio_data = path[1]

    FN = 0. # False Negative
    FP = 0. # False Positive
    TP = 0. # True Positive
    TN = 0. # True Negative
    
    # Read precomputed ground truth dictionary
    if path_ground_truth.find('/')==-1: # Windows based
        ground_truth = json.load(file(path_ground_truth + '\\' + 'ground_truth.txt')) 
    else: # Linux based
        ground_truth = json.load(file(path_ground_truth + '/' + 'ground_truth.txt')) 
    
    # Extract estimated pairs dictionary
    if offset_estimation_result_filename.find('fingerprint') == -1: 
        estimations = extract_estimated_pairs(offset_estimation_result_filename)
    else:
        estimations = extract_estimated_pairs(offset_estimation_result_filename)
        
    # Scan estimations to find FN, FP types of errors and TP
    tolerance = 10 # The alignment is acceptable in neighborhood of ground truth with amount of tolerance
    for key, value in estimations.iteritems():
        if ground_truth.has_key(key):
            true_relative_offset = ground_truth[key]
            estimated_relative_offset = value
            if np.abs(true_relative_offset - estimated_relative_offset) > tolerance:
                FN+=1.            
            else:
                TP+=1
        else:
            FP+=1.
    
    # Scan ground_truth to find FN type of error
    for key in ground_truth.iterkeys():
        if estimations.has_key(key) == False:
            FN+=1.
    
    # Find number of recordings for each microphone
    number_of_microphones = 4
    number_of_recordings = np.zeros(number_of_microphones,dtype='float64')
    for filename in listdir(path_audio_data):
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
    
    # True Negatives are computed
    TN = total_number_of_pairs - (FN + FP + TP) - excluded_number_of_pairs
    
    # Compute the evaluation metrics: Accuracy, Precision, Recall, F-measure
    Accuracy = (TP + TN)/(TP + TN + FP + FN)
    Precision = TP / (TP + FP)
    Recall = TP / (TP + FN)
    F_measure = 2 * (Precision * Recall)/(Precision + Recall)
    
    if verbose == True:
        print("\nThe evaluation results for " + offset_estimation_result_filename)
        print(('\nFalse Negative - FN = {0}').format(FN))
        print(('False Positive - FP = {0}').format(FP))
        print(('True Positive - TP = {0}').format(TP))
        print(('True Negative - TN = {0}').format(TN))
        print(('\nAccuracy = {0}').format(Accuracy))
        print(('Precision = {0}').format(Precision))
        print(('Recall = {0}').format(Recall))
        print(('F-measure = {0}').format(F_measure))
    
    return (Accuracy, Precision, Recall, F_measure, TP, TN, FP, FN)  

if __name__ == '__main__':
    cw_path = os.getcwd();
    if cw_path.find('/')==-1:
        cw_path_parent = cw_path[:cw_path.find('\\Evaluation')]
        path_ground_truth = cw_path + '\\ground_truth'
        #path_audio_data = cw_path_parent + '\\audio_data'
    else:
        cw_path_parent = cw_path[:cw_path.find('/Evaluation')]
        path_ground_truth = cw_path + '/ground_truth'
        #path_audio_data = cw_path_parent + '/audio_data'
        
    root = Tkinter.Tk()
    root.withdraw() #use to hide tkinter window
    
    currdir = os.getcwd()
    path_audio_data = tkFileDialog.askdirectory(parent=root, initialdir=currdir, title='Please select a directory')

    offset_estimation_result_filename = tkFileDialog.askopenfilename(parent=root, initialdir=currdir, title='Please select a result file')
   
    path = [path_ground_truth, path_audio_data]

    Accuracy, Precision, Recall, F_measure, TP, TN, FP, FN = compute_accuracy(path, offset_estimation_result_filename, True) 
    