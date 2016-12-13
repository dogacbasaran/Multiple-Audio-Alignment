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
    