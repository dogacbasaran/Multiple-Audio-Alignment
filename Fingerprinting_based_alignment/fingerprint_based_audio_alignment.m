% FINGERPRINTING_BASED_AUDIO_ALIGNMENT: Implements fingerprinting based
% audio alignment 
%   This function uses the fingerprinting scheme from 
%   D. Ellis (2009), "Robust Landmark-Based Audio Fingerprinting", 
%   web resource, available: http://labrosa.ee.columbia.edu/matlab/fingerprint/ .
%   
%   The function reads the audio files from the dataset, applies a grid
%   search to find the best threshold to detect matching sequences. The
%   estimated alignments with the highest accuracy are written in a text
%   file.

% Copyright (C) 2016  Author: Dogac Basaran
%
%    This program is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    This program is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with this program.  If not, see <http://www.gnu.org/licenses/>

tic
current_folder = pwd;
if isempty(strfind(current_folder,'/'))==1 % OS is Windows
    parent_folder = current_folder(1:strfind(current_folder,'\Fingerprinting_based_alignment'));
    load_path = [parent_folder 'audio_data\'];    
    addpath([current_folder '\fingerprint_labrosa']);
else % OS is Linux
    parent_folder = current_folder(1:strfind(current_folder,'/Fingerprinting_based_alignment'));
    load_path = [parent_folder 'audio_data/'];
    addpath([current_folder '/fingerprint_labrosa']);
end

% Read files
tks= myls([load_path '*.wav']);
 
filenames = cell(1, length(tks));
for k = 1:length(tks)   
    str = tks{k};
    filename = str(strfind(str,load_path)+length(load_path):end);
    filenames{k} = filename;
end 
 
% Initialize the hash table database array 
clear_hashtable
% Calculate the landmark hashes for each reference track and store
% it in the array 
add_tracks(tks);
 
fprintf('\nFeature extraction: %f secs\n',toc)

% A grid search for the threshold to find the most similar sequences
thresholds = 10:10:150;

estimated_alignments = cell(1, length(tks));
for n = 1:length(thresholds)
    tic
    thr = thresholds(n);
    fileID = fopen(['offset_estimation_fingerprinting_thr_' num2str(thr) '_result.txt'],'w');

    for k = 1:length(tks)
        [x,Fs] = audioread(tks{k});
        current_sequence = filenames{k};

        % Run the query
        R = match_query(x,Fs);
        % Find the aligned sequences with the current_sequence
        R_aligned = R(R(:,2)>thr,:);
        estimated_alignments{k} = R_aligned;
    end

    % Generate the output file for the evaluation part    
    for k = 1:length(tks)
        current_sequence = filenames{k};
        str_towrite = current_sequence;

        R_aligned = estimated_alignments{k};
        if size(R_aligned,1) ~= 1
            for i = 2:length(R_aligned(:,1))
                number_of_landmarks = R_aligned(i,2);
                R_aligned_ = estimated_alignments{R_aligned(i,1)};
                
                number_of_landmarks_ = R_aligned_(R_aligned_(:,1)==k,2);
                if isempty(number_of_landmarks_)
                    number_of_landmarks_ = 0;
                end
                if number_of_landmarks > number_of_landmarks_
                    relative_offset = -R_aligned(i,3);
                else
                    relative_offset = R_aligned_(R_aligned_(:,1)==k,3);
                end
                str_towrite = [str_towrite ' ' filenames{R_aligned(i,1)} ',' num2str(relative_offset)];          
            end       
        end
        str_towrite = [str_towrite ' \n'];
        fprintf(fileID,str_towrite);
    end
    
    fclose(fileID);
    fprintf('\nThreshold = %d, the elapsed time for procedure: %f secs\n\n',thr, toc)
end
