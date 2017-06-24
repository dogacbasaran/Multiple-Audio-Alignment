""" 
groundtruth.py
~~~~~~~~~~~~~~~

.. topic:: Contents:
    
    The groundtruth module pulls the offset values of the sequences from the Jiku_GT_090912.xml 
    file and forms the pairwise connected sequence list. The list is then written in
    a text file "grount_truth.txt". 

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
import scipy.io.wavfile
import json
import os
from os import listdir

def find_filelengths(path):        
    """ Find the length of each sequence, in addition finds the silence parts in the
        beginning of each sequence for calibration of the offsets accordingly 
        
    **Parameters**
    
    path: String 
        The path to the audio dataset

    **Returns**
    
    filelengths: Dictionary
        The dictionary that contains the file lengths in seconds. Keys are the filenames
    silence_at_the_beginning: Dictionary
        The dictionary that contains the silence part in seconds for each file. Keys are the filenames"""
        
    filelengths = {}
    silence_at_the_beginning = {}
    for filename in listdir(path):
        Fs, x = scipy.io.wavfile.read(path + '\\' + filename)
        if np.where(x!=0)[0][0]-1 >=0:
            silence_at_the_beginning[filename] = np.float64(np.where(x!=0)[0][0]-1)/np.float64(Fs)
        else:
            silence_at_the_beginning[filename] = 0.
        filelengths[filename] = np.float64(x.shape[0])/np.float64(Fs)    
    return filelengths, silence_at_the_beginning

 
def extract_offsets(path, filelengths, silence_in_the_beginning, hopsize = 0.02):
    
    """Extract offset information of each file from ground truth file GT_090912.xml 
    and convert length and offset information from seconds to frames (STFT based).
    Note that this piece of code is unique for the GT_090912.xml file. For other
    Jiku datasets it has to be modified!!
 
    **Parameters**
    
    path: String
        path of the GT_090912.xml file
    filelengths: Dictionary
        The dictionary that contains the file lengths in seconds. Keys are the filenames
    silence_at_the_beginning: Dictionary
        The dictionary that contains the silence part in seconds for each file. Keys are the filenames
    hopsize: Float (default 0.02)
        Hop size in the STFT
        
    **Returns**
    
    offset_length_in_frames: List
        The list of offset lengths in frames    
    """
    
    # cnti is the number of recordings for microphone i
    cnt1 = 1
    cnt2 = 1
    cnt3 = 1
    cnt4 = 1
    offset_length_in_frames = []
    for line in open(path1 + 'Jiku_GT_090912.xml','r').readlines():
        if line.find('offset') != -1:                        
            # Extract the length of the sequence
            mic_number = np.int(line.split('_')[2])
            if mic_number == 0: # represents 1st mic
                filename = 'mic1_rec' + np.str(cnt1) + '.wav'
                N = np.floor(filelengths[filename]/hopsize) # Length in frames
                cnt1+=1
            elif mic_number == 2: # represents 2nd mic
                filename = 'mic2_rec' + np.str(cnt2) + '.wav'
                N = np.floor(filelengths[filename]/hopsize) # Length in frames
                cnt2+=1
            elif mic_number == 3: # represents 3rd mic
                filename = 'mic3_rec' + np.str(cnt3) + '.wav'
                N = np.floor(filelengths[filename]/hopsize) # Length in frames
                cnt3+=1
            elif mic_number == 4: # represents 4th mic
                filename = 'mic4_rec' + np.str(cnt4) + '.wav'
                N = np.floor(filelengths[filename]/hopsize) # Length in frames
                cnt4+=1               
            
            # Extract the offset of the sequence
            tmp_offset = line.split(':') 
            offset = np.float(tmp_offset[1]) * 3600. + np.float(tmp_offset[2]) * 60. + np.float(tmp_offset[3][:6])
            offset = offset + silence_in_the_beginning[filename] # Calibration for silence parts
            offset = np.floor(offset/hopsize) # Convert offset from seconds to frames
            
            # Write Offset-Length tuple to the list    
            offset_length_in_frames.append((filename, offset, N))            
    return offset_length_in_frames

def set_key_name(filename1, filename2):
    
    """Sets the key name for ground_truth dictionary entry, order the sequences in 
    ascending order of microphone number and record number. 
        
    **Parameters**
    
    filename1: String
        The name of the first file
    filename2: String
        The name of the second file
        
    **Returns**
    
    key_name: String
        Key name for the ground_truth dictionary."""
    
    mic_number = int(filename1[filename1.find('mic')+3:filename1.find('_')])
    mic_number_ = int(filename2[filename2.find('mic')+3:filename2.find('_')])
    rec_number = int(filename1[filename1.find('rec')+3:filename1.find('.wav')])
    rec_number_ = int(filename2[filename2.find('rec')+3:filename2.find('.wav')])
    
    if mic_number > mic_number_:
        key_name = filename2 + ' ' + filename1
    elif mic_number == mic_number_:
        if rec_number > rec_number_:
            key_name = filename2 + ' ' + filename1
        else:
            key_name = filename1 + ' ' + filename2
    else:
        key_name = filename1 + ' ' + filename2
    return key_name
    
def set_relative_offset(key_name, filename, offset1, offset2):
    
    """Sets the true relative distance between two sequences according to the
       ordering of the sequences
       
    **Parameters**
    
    key_name: String
        The key name for the ground_truth dictionary entry
    filename: String
        The sequence name
    offset1: Integer
        Starting point of the first sequence in frames
    offset2: Integer
        Starting point of the second sequence in frames
   
    **Returns**
    
    relative_offset_distance: float
        Relative distance according to the ordering of the aligned sequences"""
    
    if key_name.find(filename)==0:
        relative_offset_distance = offset2-offset1
    else:
        relative_offset_distance = offset1-offset2
    return relative_offset_distance

def set_paths():
    
    """Set the paths to the ground_truth (Path 1) and audio data (Path 2). Detects
    the operating system and set the paths accordingly.
    
    **Parameters**
    
    **Returns**
    
    path1: String
        Path to the ground truth
    path2: String
        Path to the audio dataset
    wl: Integer
        1 represents Windows based and 2 represents Linux based systems."""
    
    cw_path = os.getcwd();
    if cw_path.find('/')==-1:
        cw_path_parent = cw_path[:cw_path.find('\\Evaluation')]
        path1 = cw_path + '\\ground_truth\\'
        path2 = cw_path_parent + '\\audio_data\\'
        wl = 1; # Defines windows based OS
    else:
        cw_path_parent = cw_path[:cw_path.find('/Evaluation/')]
        path1 = cw_path + 'ground_truth/'
        path2 = cw_path_parent + '/audio_data/'
        wl = 2; # Defines windows based OS
        
    return (path1,path2,wl)

if __name__ == '__main__':
    
    path1, path2, wl = set_paths()
    
    filelengths, silence_in_the_beginning = find_filelengths(path2)        
    
    offset_length_in_frames = extract_offsets(path2, filelengths, silence_in_the_beginning, hopsize = 0.032)        
    
    # Finding the overlapping sequences and relative distances between offsets   
    ground_truth = {}     
    res = 4 # The amount in frames that overlapping won't be counted
    f = open(path1 + 'ground_truth_sequence_ordering.txt','w') # write sequence ordering to a txt file
    for i in range(len(offset_length_in_frames)):
        f.write(offset_length_in_frames[i][0] + "\n")    
        
        filename = offset_length_in_frames[i][0]
        offset = offset_length_in_frames[i][1]
        length = offset_length_in_frames[i][2]
        for j in range(i+1,len(offset_length_in_frames)):
            filename_ = offset_length_in_frames[j][0]
            offset_ = offset_length_in_frames[j][1]
            length_ = offset_length_in_frames[j][2]
            
            if offset + length > offset_ + res:
                key_name = set_key_name(filename, filename_)
                relative_offset_distance = set_relative_offset(key_name, filename, offset, offset_)
                ground_truth[key_name] = relative_offset_distance
                print(key_name)       
                print(ground_truth[key_name])
            else:
                break
               
    json.dump(ground_truth,file(path1 + 'ground_truth.txt','w'))
    f.close()
