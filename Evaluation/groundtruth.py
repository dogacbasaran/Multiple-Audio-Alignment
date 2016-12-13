""" 
grount_truth.py reads the offset values of the sequences from the GT_090912.xml 
file and forms the pairwise connected sequence list. The list is then written in
a text file "grount_truth.txt".
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
import scipy.io.wavfile
import json
import os
from os import listdir


# Find the length of each sequence, in addition finds the silence parts in the
# beginning of each sequence for calibration of the offsets accordingly 
# Input: path - path of the files
# Output: filelengths - Type: Dictionary keys: filenames values: lengths of files
#         silence_in_the_beginning - Type Dictionary keys: filenames values: amount of silence in seconds
def find_filelengths(path):        
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

# Extract offset information of each file from ground truth file GT_090912.xml 
# and convert length and offset information from seconds to frames (STFT based)
# Note that this piece of code is unique for the GT_090912.xml file. For other
# Jiku datasets it has to be modified!!
# Input: path - path of the GT_090912.xml file
#        filelengths - dictionary keys: filenames values: lengths of files
# Output: offset_length_in_frames - Type:List in format "filename-offset-length"
def extract_offsets(path, filelengths, silence_in_the_beginning, hopsize = 0.02):
    # cnti is the number of recordings for microphone i
    cnt1 = 1
    cnt2 = 1
    cnt3 = 1
    cnt4 = 1
    offset_length_in_frames = []
    for line in open(path1 + '\\GT_090912.xml','r').readlines():
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

# Sets the key name for ground_truth dictionary entry, order the sequences in 
# ascending order of microphone number and record number
# Example: If filename = 'mic4_rec2.wav', filename_ = 'mic1_rec3.wav'
#           the key_name = 'mic1_rec3.wav mic4_rec2.wav'
# Input: filename - First filename
#        filename_ - Second filename
# Output: key_name - Either 'filename filename_' or 'filename_ filename'
def set_key_name(filename, filename_):
    mic_number = int(filename[filename.find('mic')+3:filename.find('_')])
    mic_number_ = int(filename_[filename_.find('mic')+3:filename_.find('_')])
    rec_number = int(filename[filename.find('rec')+3:filename.find('.wav')])
    rec_number_ = int(filename_[filename_.find('rec')+3:filename_.find('.wav')])
    
    if mic_number > mic_number_:
        key_name = filename_ + ' ' + filename
    elif mic_number == mic_number_:
        if rec_number > rec_number_:
            key_name = filename_ + ' ' + filename
        else:
            key_name = filename + ' ' + filename_
    else:
        key_name = filename + ' ' + filename_
    return key_name
    
def set_relative_offset(key_name, filename, offset, offset_):
    if key_name.find(filename)==0:
        relative_offset_distance = offset_-offset
    else:
        relative_offset_distance = offset-offset_
    return relative_offset_distance

def set_paths():
    cw_path = os.getcwd();
    if cw_path.find('/')!=-1:
        cw_path = cw_path[:cw_path.find('\\Evaluation\\')]
        path1 = cw_path + '\\ground_truth\\'
        path2 = cw_path + '\\audio_data\\'
        wl = 1; # Defines windows based OS
    else:
        cw_path = cw_path[:cw_path.find('/Evaluation/')]
        path1 = cw_path + '/ground_truth/'
        path2 = cw_path + '/audio_data/'
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
