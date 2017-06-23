""" 
fingerprinting_evaluation.py  

The code can work under Windows or Linus OS. 

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
    along with this program.  If not, see <http://www.gnu.org/licenses/>
"""

import numpy as np
import compute_accuracy
import matplotlib.pyplot as plt
import os

if __name__ == '__main__':
    cw_path = os.getcwd()
    
    if cw_path.find('/')==-1:
        cw_path_parent = cw_path[:cw_path.find('\\Evaluation')]
        path1 = cw_path + '\\ground_truth'
        path2 = cw_path_parent + '\\audio_data'
        path3 = cw_path + '\\SMC_offset_estimation_results'
        path4 = cw_path + '\\fingerprinting_offset_estimation_results'
    else:
        cw_path_parent = cw_path[:cw_path.find('/Evaluation')]
        path1 = cw_path + '/ground_truth'
        path2 = cw_path_parent + '/audio_data'
        path3 = cw_path + '/SMC_offset_estimation_results'
        path4 = cw_path + '/fingerprinting_offset_estimation_results'
    
    path = [path1, path2, path3, path4]
    
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
        offset_estimation_result_filename = 'offset_estimation_fingerprinting_thr_' + np.str(thr) + '_result.txt'
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
        
    print('The best accuracy for fingerprinting based method \n   Threshold = {0}  Accuracy = %{1}'.format(thresholds[index_of_best], best_accuracy_fingerprinting_based))
        
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
##    
#    axis_font = {'fontname':'Arial', 'size':'15'}
#
#    plt.Figure()
#    plt.plot(thresholds, accuracy_list)
#    plt.xlabel('Thresholds', **axis_font)
#    plt.ylabel('Accuracy (%)', **axis_font)
#    plt.ylim([80,100])
#    
#    plt.Figure()
#    plt.rc('legend', fontsize=14)
#    plt.plot(thresholds, precision_list, label = '$Precision$')
#    plt.plot(thresholds, recall_list, '+',label = '$Recall$')
#    plt.plot(thresholds, F_measure_list, '--',label = '$F-measure$')
#    plt.xlabel('Thresholds', **axis_font)
#    plt.legend()
    
    
    fig.tight_layout()
    if cw_path.find('/')==-1:
        fig.savefig(path4 + '\\' + 'Fingerprinting_based_estimation_results.png')
    else:
        fig.savefig(path4 + '/' + 'Fingerprinting_based_estimation_results.png')
    
    
    
    
