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
    
    accuracy_omega = [] # Accuracy results list
    accuracy_FN_1 = [] # False Negative Type_1 results list
    accuracy_FN_2 = [] # False Negative Type_2 results list
    accuracy_FP = [] # False Positive results list

    # Illustration of the results of fingerprinting based alignment for 
    # threholds: 10, 20, ..., 150
    thresholds = range(10,151,10)
    for thr in thresholds:        
        offset_estimation_result_filename = 'offset_estimation_fingerprinting_thr_' + np.str(thr) + '_result.txt'
        Omega, FN_1, FN_2, FP = compute_accuracy.compute_accuracy(path, offset_estimation_result_filename)
        accuracy_omega.append(100*Omega)
        accuracy_FN_1.append(FN_1)
        accuracy_FN_2.append(FN_2)
        accuracy_FP.append(FP)
    
    best_accuracy_fingerprinting_based = max(accuracy_omega)    
    index_of_best = accuracy_omega.index(best_accuracy_fingerprinting_based)
        
    print('The best accuracy for fingerprinting based method \n   Threshold = {0}  Accuracy = %{1}'.format(thresholds[index_of_best], best_accuracy_fingerprinting_based))
    
    fig, axes = plt.subplots(1,2)
    axes[0].plot(thresholds, accuracy_omega)
    axes[0].set_xlabel('Thresholds')
    axes[0].set_ylabel('Accuracy (%)')
    axes[0].set_ylim([80,100])
    
    axes[1].plot(thresholds, accuracy_FN_1, label = '$FN_1$')
    axes[1].plot(thresholds, accuracy_FN_2, label = '$FN_2$')
    axes[1].plot(thresholds, accuracy_FP, label = '$FP$')
    axes[1].set_xlabel('Thresholds')
    axes[1].set_ylabel('Counts')
    axes[1].set_ylim([-5, 65])
    axes[1].legend()
    
    fig.tight_layout()
    fig.savefig(path4 + '\\' + 'Fingerprinting_based_estimation_results.png')
    
    
    
    
    
