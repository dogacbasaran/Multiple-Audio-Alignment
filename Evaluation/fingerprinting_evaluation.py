""" 
fingerprinting_evaluation.py  
~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.. topic:: Contents

    This function searches for the best threshold value for the fingerprinting
    based alignment system by evaluating the Accuracy, Precision, Recall and 
    F-measure metrics. The best threshold is found using the aacuracy measure.
    
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
import compute_accuracy
import matplotlib.pyplot as plt
import os
import Tkinter
import tkFileDialog

if __name__ == '__main__':
    
    cw_path = os.getcwd() 
    if cw_path.find('/')==-1:
        cw_path_parent = cw_path[:cw_path.find('\\Evaluation')]
        path_ground_truth = cw_path + '\\ground_truth'
        #path_audio_data = cw_path_parent + '\\audio_data'
        path_fingerprinting_results = cw_path + '\\fingerprinting_offset_estimation_results'
    else:
        cw_path_parent = cw_path[:cw_path.find('/Evaluation')]
        path_ground_truth = cw_path + '/ground_truth'
        #path_audio_data = cw_path_parent + '/audio_data'
        path_fingerprinting_results = cw_path + '/fingerprinting_offset_estimation_results'
    
    root = Tkinter.Tk()
    root.withdraw() #use to hide tkinter window
    
    currdir = os.getcwd()
    path_audio_data = tkFileDialog.askdirectory(parent=root, initialdir=currdir, title='Please select a directory')
 
    path = [path_ground_truth, path_audio_data]
    
    accuracy_list = []
    precision_list = []
    recall_list = []
    F_measure_list = []
    
    TP_list = []
    TN_list = []
    FP_list = []
    FN_list = []
    
    # Illustration of the results of fingerprinting based alignment for 
    # threholds: 10, 20, ..., 150
    thresholds = range(10,151,10)
    for thr in thresholds:        
        offset_estimation_result_filename = '{0}/offset_estimation_fingerprinting_thr_'.format(path_fingerprinting_results) + np.str(thr) + '_result.txt'
        Accuracy, Precision, Recall, F_measure, TP, TN, FP, FN = compute_accuracy.compute_accuracy(path, offset_estimation_result_filename)
        accuracy_list.append(100*Accuracy)
        precision_list.append(Precision)
        recall_list.append(Recall)
        F_measure_list.append(F_measure)
        TP_list.append(TP)
        TN_list.append(TN)
        FN_list.append(FN)
        FP_list.append(FP)
        
    best_accuracy_fingerprinting_based = max(accuracy_list)    
    index_of_best = accuracy_list.index(best_accuracy_fingerprinting_based)
        
    print('The best accuracy is obtained for threshold = {0} with accuracy = {1}'.format(thresholds[index_of_best], best_accuracy_fingerprinting_based))
    
    print('\nOther metrics:')
    print(('\nFalse Negative - FN = {0}').format(FN_list[index_of_best]))
    print(('False Positive - FP = {0}').format(FP_list[index_of_best]))
    print(('True Positive - TP = {0}').format(TP_list[index_of_best]))
    print(('True Negative - TN = {0}').format(TN_list[index_of_best]))
    print(('\nAccuracy = {0}').format(accuracy_list[index_of_best]))
    print(('Precision = {0}').format(precision_list[index_of_best]))
    print(('Recall = {0}').format(recall_list[index_of_best]))
    print(('F-measure = {0}').format(F_measure_list[index_of_best]))
    
    fig, axes = plt.subplots(1,2)
    axes[0].plot(thresholds, accuracy_list)
    axes[0].set_xlabel('Thresholds')
    axes[0].set_ylabel('Accuracy (%)')
    axes[0].set_ylim([80,100])
        
    axes[1].plot(thresholds, precision_list, label = '$Precision$')
    axes[1].plot(thresholds, recall_list, '+',label = '$Recall$')
    axes[1].plot(thresholds, F_measure_list, '--',label = '$F-measure$')
    axes[1].set_xlabel('Thresholds')
    axes[1].legend()    
    
    fig.tight_layout()
    if cw_path.find('/')==-1:
        fig.savefig(path_fingerprinting_results + '\\' + 'Fingerprinting_based_estimation_results.png')
    else:
        fig.savefig(path_fingerprinting_results + '/' + 'Fingerprinting_based_estimation_results.png')
    
    
    
    
