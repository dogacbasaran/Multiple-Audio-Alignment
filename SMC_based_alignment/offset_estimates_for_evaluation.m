%% Creates the text file where the results of the alignment are listed
% Copyright (C) 2016  Dogac Basaran

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

function offset_estimates_for_evaluation(Clusters, rs_fixed, dataset_features)

MicRec_sorted = dataset_features.MicRec_sorted;
Nsorted = dataset_features.Nsorted;

cw_path = pwd; % Get current working directory
% Decide whether it is windows or linux based system based on the slash orientation
if ~isempty(strfind(cw_path,'/')) % Linux based system
    cw_path = cw_path(1:strfind(cw_path, 'SMC_based_alignment')-1);
    gt_path = [cw_path 'Evaluation/ground_truth/'];
    results_path = [cw_path 'Evaluation/SMC_offset_estimation_results/'];
else % Windows based system
    cw_path = cw_path(1:strfind(cw_path, 'SMC_based_alignment')-1);
    gt_path = [cw_path 'Evaluation\\ground_truth\\'];
    results_path = [cw_path 'Evaluation\\SMC_offset_estimation_results\\'];
end

% Read the ground truth ordering of the sequences
fid = fopen([gt_path 'ground_truth_sequence_ordering.txt'],'r');
GT_ordering = textscan(fid, '%s', 'Delimiter', '\n'); % Ground truth ordering list
fclose(fid);

clk = clock;
date_of_experiment = [num2str(clk(3)) '_' num2str(clk(2)) '_' num2str(clk(1)) '_' num2str(clk(4)) 'h_' num2str(clk(5)) 'm'];
fileID = fopen([results_path 'offset_estimation_SMC_result_' date_of_experiment '.txt'],'w');

estimated_offsets = cell(length(GT_ordering{1}),1);
for k = 1:length(GT_ordering{1}) 
    filename = GT_ordering{1}(k);
    filename = filename{1};
    % Extract the mic number and rec number of the sequence
    mic_number = str2double(filename(strfind(filename,'mic')+3:strfind(filename,'rec')-2));
    rec_number = str2double(filename(strfind(filename,'rec')+3:strfind(filename,'.wav')-1));
    
    str_towrite = ['mic' num2str(mic_number) '_rec' num2str(rec_number) '.wav '];
    
    ind = find(MicRec_sorted(:,1)==mic_number & MicRec_sorted(:,2)==rec_number);
    res = 0; % The amount of overlap in frames for which the overlap doesn't count
    for c = 1:length(Clusters)
        if sum(Clusters{c}==ind) % Checks if the cluster has that sequence
            sequence_offsets = rs_fixed{c};
            sequence_list = Clusters{c};
            sequence_offsets(sequence_list == ind) = [];              
            sequence_list(sequence_list == ind) = [];                        
            
            % Check each sequence if they are overlapping with the current
            % sequence
            current_sequence_offset = rs_fixed{c}(Clusters{c}==ind);
            N_current_sequence = Nsorted(end,ind);
            %N_current_sequence = N(end,ind);
            
            for j = 1:length(sequence_list)
                test_sequence_offset = sequence_offsets(j);
                N_test_sequence = Nsorted(end,sequence_list(j));
                %N_test_sequence = N(end,sequence_list(j));
                if test_sequence_offset <=current_sequence_offset
                    if test_sequence_offset + N_test_sequence > current_sequence_offset + res
                        str_towrite = [str_towrite 'mic' num2str(MicRec_sorted(sequence_list(j),1)) '_rec' ...
                                       num2str(MicRec_sorted(sequence_list(j),2)) '.wav,'...
                                       num2str(test_sequence_offset - current_sequence_offset) ' '];
                    end
                else
                    if current_sequence_offset + N_current_sequence > test_sequence_offset + res                        
                        str_towrite = [str_towrite 'mic' num2str(MicRec_sorted(sequence_list(j),1)) '_rec' ...
                                       num2str(MicRec_sorted(sequence_list(j),2)) '.wav,'...
                                       num2str(test_sequence_offset - current_sequence_offset) ' '];
                    end
                end
            end
            
            str_towrite = [str_towrite '\n'];
            fprintf(fileID,str_towrite);
            break;
        end
    end
end
fclose(fileID);

