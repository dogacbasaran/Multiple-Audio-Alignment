%% SMC based alignment module
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


function [Clusters, rs_fixed, time_elapsed] = SMC_main_module(dataset_features, print_proc, display_results)

if nargin<2
    print_proc = true;
    display_results = false;
elseif nargin<3
    display_results = false;
end

%% SMC Sampler based Sequential Alignment
disp(' ')
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp('% SMC procedure started.. %')
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%')
tic

% Initialize the parameters of SMC procedure
NoC = 1; % Number of estimated clusters
NAS = 1; % Number of aligned sources for current cluster
K = dataset_features.K; % Number of audio sequences
Clusters = cell(K,1);  % All estimated clusters
Clusters{1} = 1;
rs_fixed = cell(K,1); % All estimated alignments
rs_fixed{1} = 1;
p = 1:K; % The list of sequence indices that are going to be aligned
r_fixed = 1; % The starting points of aligned sequences for current cluster
run_number = 1; % Number of runs of the sequential matching algorithm
terminate = false; % terminate = true ==> Terminate the process

% Unclustered list of sequence indices. At start, all the sequences are unclustered!
% This variable holds the unclustered list of sequence indices through all runs
Unclustered = 1:K; 

while ~terminate
    % Display the Run number and sequence indices to be aligned
    disp(' ')
    disp(['-Run #' num2str(run_number) '-']) 
    if print_proc
        fprintf('\nSequence Indices to be aligned: ');
        pu = p(Unclustered);
        for k=1:length(pu)
            fprintf('%d ',pu(k));
        end
        fprintf('\n');
    end
    
    % There are two types of runs of procedure:
    %     Type 1: There are no pre-aligned sequences in the cluster
    %     Type 2: There are pre-aligned sequences in the cluster
    % Go through all unclustered sequences once to update the aligned sequence list
    [r_fixed, Cluster, Unclustered] = sequential_alignment_module(p(Unclustered), rs_fixed{NoC}, Clusters, NoC, run_number, dataset_features);
    
    run_number = run_number +1; % Update the run number

    % CLusters{NoC} holds the sequence indices aligned before the procedure with number run_number
    % Cluster holds the sequence indices aligned after the procedure with number run_number
    % Comparing these values, one can decide whether to decide the next procedure run
    % would be Type 1 or Type 2 .
    if length(Cluster) ~= length(Clusters{NoC}) % Type 2
        % Update the current cluster information and appropriate sequence indices        
        [Clusters, rs_fixed] = update_clusters_offsets(Clusters, rs_fixed, NoC, Cluster, r_fixed);

        terminate = check_no_unclustered(Unclustered);
        if terminate == false 
            % When a few number of consecutive sequences starting from the 
            % first one are aligned at the end of a run. The remaining 
            % unclustered sequences don't align with this cluster. 
            % This case happens the minimum index of unclustered sequences 
            % and is larger than maximum index of clustered sequences
            mn_unclustered = min(Unclustered);
            mx_clustered = max(Cluster);
            if  mx_clustered < mn_unclustered
                % Switch to next cluster
                [Clusters, rs_fixed, NoC] = switch_to_next_cluster(Clusters, rs_fixed, NoC, Unclustered);
                terminate = check_no_unclustered(Unclustered);                
            end            
        end
    else % Type 1
        % Switch to next cluster
        [Clusters, rs_fixed, NoC] = switch_to_next_cluster(Clusters, rs_fixed, NoC, Unclustered);
        terminate = check_no_unclustered(Unclustered);        
    end
end

time_elapsed = toc;
offset_estimates_for_evaluation(Clusters, rs_fixed, dataset_features);


function [Clusters, rs_fixed] = update_clusters_offsets(Clusters, rs_fixed, NoC, Cluster, r_fixed)
    % Update the current cluster information and appropriate sequence indices        
    Clusters{NoC} = Cluster;
    rs_fixed{NoC} = r_fixed;

function [Clusters, rs_fixed, NoC] = switch_to_next_cluster(Clusters,  rs_fixed, NoC, Unclustered)
    NoC = NoC + 1;
    Clusters{NoC} = Unclustered(1);
    rs_fixed{NoC} = 1;

function terminate = check_no_unclustered(Unclustered)
    if isempty(Unclustered)
        terminate = true;
    elseif length(Unclustered)==1
        terminate = true;
    else
        terminate = false;
    end

function display_results(Clusters, rs_fixed, NoC, dataset_features)
    ss = dataset_features.ss;
    N = dataset_features.N;
    Fs = dataset_features.Fs;
    MicRec_sorted = dataset_features.MicRec_sorted;
    hopsize = dataset_features.Ws / 2;
     % Displaying the alignment results
    for i=1:NoC-1
        len = length(Clusters{i});
        figure,
        mx = max(rs_fixed{i}+N(end,Clusters{i}));
        for j = 1:len
            subplot(len,1,j),plot((rs_fixed{i}(j)*hopsize*Fs:rs_fixed{i}(j)*hopsize*Fs+length(ss{Clusters{i}(j)})-1)/Fs,ss{Clusters{i}(j)})
            title(['Seq. ' num2str(MicRec_sorted(Clusters{i}(j),1))])
            title(['Mic' num2str(MicRec_sorted(Clusters{i}(j),1)) 'Rec' num2str(MicRec_sorted(Clusters{i}(j),2))])
            set(gca, 'xlim', [0 (mx*hopsize*Fs+1)/Fs])
        end
    end
 
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
        