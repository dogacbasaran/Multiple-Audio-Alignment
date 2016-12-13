% function [r_fixed, Cluster, Uncluster, internal_struct] = SM_M5_F5(p, r_fixed, Cluster, data_struct, internal_struct)
% 
% Sequential Matching Algorithm for M5: Be-Be F5: Feature: Threshold first Difference Through Time and Frequency
% 
% The inputs and outputs of the function are as follows;
% 
% Inputs:
%   p                      : The sequence indices that are going to be aligned
%   r_fixed                : The starting points of aligned sequences
%   Cluster                : The sequence indices that are overlapping (that are in the same cluster)
%   S                      : The sequences at each resolution (1st row is the lowest resolution)
%   data_struct            : Structure contains the features of sequences of all resolution, length of each sequence,
%                                   number of subbands, sampling frequency
%                                   and minimum number of frame overlap, to be accepted as aligned
%   internal_struct        : Structure contains internal variables S1 and S2
% 
% Outputs:
%   r_fixed                : Updated starting points of aligned sequences
%   Cluster                : Updated sequence indices that are overlapping (that are in the same cluster)
%   Uncluster              : The sequence indices that are not in the same cluster (that can not be clustered with current cluster)
%   internal_struct        : Updated structure contains internal variables S1 and S2
% 
% The generative model is,
% 
% lambda_tau ~ Be(lambda_tau; alpha_L)  True time line index
% r_k ~ U[1,T-N_k+1] True time index for source k
% x_(k,n) | r_k, lambda_tau ~ PROD(1:T) P(x_(k,n) | r_k, lambda_tau)^[n=tau-r_k]
% 
% P(x_(k,n) | r_k, lambda_tau) = exp(Sum_i Sum_j [x_(k,n)=i][lambda_tau=j] log(w))
%
% The aim is to find the best alignment of sources applying the following
% methodology:
% - Fix the first source alignment variable assume that it is known.
% - Match the next source the first source as a pairwise manner. 
% - Then match the next source according to the fixed source alignments
% 
% - When aligning sources sequentially, if the current source do not
% overlap with the previous sources, the ordering of the sources changes 
% in a way that the non-overlapping source is moved to the end of the 
% ordering and the rest of the sources are circulary shift to left by 1
% Example:
% p = [4 2 3 5 1 6]
% Assume that source 4 and 2 overlaps the current search for source 3 and
% it does not overlap with the previous sources then the ordering is
% changed to p = [4 2 5 1 6 3]
%% Sequential Alignment module: Applies the sequential alignment method in Manuscript Sec. 3.4
% Copyright (C) 2016 Dogac Basaran

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

function [r_fixed, Cluster, Unclustered] = sequential_alignment_module(p, r_fixed, Clusters, NoC, run_number, dataset_features)

% startResolution = 1;
% N = data_struct.Nsteps(startResolution,:); % The last row of N -> The lengths of sequences in the lowest resolution
global data_struct
data_struct = dataset_features;

F = data_struct.F;
Nsteps = data_struct.Nsteps;
Cluster = Clusters{NoC};

wsteps = data_struct.wsteps2; % We choose exponential increase of w parameter through different resolutions as annealing

% Parameters are set for fresh Cluster with one sequence or a Cluster with more than one aligned sequences
%   NAS: Number of aligned sources
%   start: Starting sequence index from unclustered sequence list "p"
%   initialClusterLength: Initial length of the current cluster
if length(Cluster) == 1 % The unclustered sequences are aligned against each other (fresh cluster)
    NAS = 1; 
    start = 2; 
    initialClusterLength = 1;
else                    % The unclustered sequences are aligned agains the current cluster with pre-aligned sequences
    NAS = length(Cluster); 
    start = 1; 
    initialClusterLength = NAS;
end

% The list of unclustered sequence indices after one pass over all sequences is completed
Unclustered = [];

for s = start:length(p) % index of the source to be aligned
    
    currentSequence = p(s); % Current sequence index
    
    if initialClusterLength ~=1 && initialClusterLength == NAS
        if currentSequence > Cluster(end)
            Unclustered = [Unclustered p(s:end)];
            for k=s:length(p)
                fprintf('%d is NOT assigned to cluster!!\n',p(k));
            end
            break;
        end
    end
        
    startResolution = set_resolution(Cluster, p, s);

    [alignment_possible, r_start, r_end, startResolution, closest] = set_search_space(startResolution, Cluster, r_fixed, p, s);    
        
    if alignment_possible
        
        logL = initialize_samples(startResolution, r_start, r_end, currentSequence, wsteps, Cluster, r_fixed, p, s);
        
        
        closestRecLeft = closest{1};
        closestRecRight = closest{2};
        closestRecLeft_index = closest{3};
        closestRecRight_index = closest{4}; 
                      
        [r_max, Num_overlapping_frames] = SMC_core_module(currentSequence, data_struct.S, F, Nsteps, wsteps, logL, Cluster, r_fixed,...
                                                               startResolution, r_start, r_end, closestRecLeft, closestRecRight,...
                                                               closestRecLeft_index, closestRecRight_index);

        [Cluster, NAS, r_fixed, Unclustered] = check_amount_of_overlap(Num_overlapping_frames, currentSequence, Unclustered, Cluster, r_fixed, r_max, NAS );
    else % The case where the current sequence cant be aligned due to other sequences from the same microphone
        Unclustered = update_Unclustered(Unclustered, currentSequence);        
    end
end           
print_results(Cluster, Unclustered, run_number);

function print_results(Cluster, Unclustered, run_number)
    % Print the results of the current Run
    fprintf('\nResults of the Run: #%d', run_number);
    fprintf('\nClustered sequences: ');
    for k=1:length(Cluster)
        fprintf('%d ',Cluster(k));
    end
    fprintf('\nUnclustered sequences: ');
    if isempty(Unclustered)
        fprintf('ALL CLUSTERED\n')
    else    
        for k=1:length(Unclustered)
            fprintf('%d ',Unclustered(k));
        end
        fprintf('\n');
    end

function [Cluster, r_fixed, NAS] = update_Cluster(Cluster, currentSequence, r_max, r_fixed, NAS, print_proc)
    global data_struct
    Nsteps = data_struct.Nsteps;

    if nargin < 6
        print_proc = true;
    end
    
    Cluster = [Cluster currentSequence];
    NAS = NAS + 1;
    if r_max < Nsteps(end,currentSequence)+1
        r_fixed = r_fixed + Nsteps(end,currentSequence)+1-r_max;
        r_fixed = [r_fixed 1];                                
    else                
        r_fixed = [r_fixed r_max-Nsteps(end,currentSequence)];                
    end            

    if print_proc
        % Printing the alignment result
        fprintf('\n%d is assigned to cluster!!\n',currentSequence);
        fprintf('\nThe cluster: ');
        for i=1:length(Cluster)
            fprintf(' %d',Cluster(i));
        end
        fprintf('\n');
    end
    
function Unclustered = update_Unclustered(Unclustered, currentSequence, print_proc)

    if nargin < 3
        print_proc = true;
    end

    Unclustered = [Unclustered currentSequence];
    
    if print_proc
        % Printing the alignment result
        fprintf('\n%d is NOT assigned to cluster!!\n',currentSequence);            
    end
        
function [Cluster, NAS, r_fixed, Unclustered] = check_amount_of_overlap(Num_overlapping_frames, currentSequence, Unclustered, Cluster, r_fixed, r_max, NAS)
    global data_struct
    Nsteps = data_struct.Nsteps;
    
    if Num_overlapping_frames > data_struct.MinNum_overlapping_frames; % DECIDE that the current source is aligned with previous sources
        [Cluster, r_fixed, NAS] = update_Cluster(Cluster, currentSequence, r_max, r_fixed, NAS);
    else
        Unclustered = update_Unclustered(Unclustered, currentSequence);
    end

function logL = initialize_samples(startResolution, r_start, r_end, currentSequence, wsteps, Cluster, r_fixed, p, s)
    global data_struct
    
    N = data_struct.Nsteps(startResolution,:); % The last row of N -> The lengths of sequences in the lowest resolution
    L = data_struct.Lsteps(startResolution);  % Lowest resolution L
    w = wsteps(startResolution); % w parameter for the lowest resolution
    x = data_struct.S(startResolution,:); % Data in lowest resolution
    F = data_struct.F; % The number of frequency bins in the feature

    T = sum(N(Cluster)) + 2*N(p(s)) - 1;

    r_fixed_Low = ceil(r_fixed/L); % The fixed offset values in the current resolution
        
    % Intermediate values for computing likelihood function in the lowest
    % resolution
    S1 = zeros(1,T);
    S2 = zeros(F,T);

    for k = 1:length(Cluster)
        S1(r_fixed_Low(k):r_fixed_Low(k)+N(Cluster(k))-1) = S1(r_fixed_Low(k):r_fixed_Low(k)+N(Cluster(k))-1) + 1;        
        S2(:,r_fixed_Low(k):r_fixed_Low(k)+N(Cluster(k))-1) = S2(:,r_fixed_Low(k):r_fixed_Low(k)+N(Cluster(k))-1) + squeeze(x{Cluster(k)});
    end

    len_phi = r_end-r_start+1;
    logL = zeros(1, len_phi); % Initialization of the likelihood function        

    temp_S1 = circshift(S1,[0,N(currentSequence)]); 
    temp_S2 = circshift(S2,[0,N(currentSequence)]);

    xs = squeeze(x{currentSequence});

    for r_est = r_start:r_end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%% Computation of likelihood function %%%%%%%%
        %%%%%%%% for the current alignment r_est    %%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        temp1 = temp_S1;
        temp2 = temp_S2;
        temp1(r_est:r_est+N(currentSequence)-1) = temp1(r_est:r_est+N(currentSequence)-1) + 1;
        temp2(:,r_est:r_est+N(currentSequence)-1) = temp2(:,r_est:r_est+N(currentSequence)-1) + xs;

        ind = find(temp1>0);
        F1 = temp2(:,ind); % Number of ones
        F2 = (temp1(ind)'*ones(1,F)*L)' - temp2(:,ind); % Number of zeros        

        logL(r_est) = sum(sum(log(0.5 * w.^(F1) .* (1-w).^(F2) + 0.5 * w.^(F2) .* (1-w).^(F1))));        
    end
    

function [r_closestRecLeft, N_sequencesLeft, closestRecLeft, closestRecLeft_index] = find_closest_Rec_left(currentSequence, startResolution, Cluster, r_fixed, ind)
    global data_struct;
    
    L = data_struct.Lsteps(startResolution);
    r_fixed_Low = ceil(r_fixed/L);
    
    currentMic = data_struct.MicRec_sorted(currentSequence,1); % Current sequence: index of Microphone
    currentRec = data_struct.MicRec_sorted(currentSequence,2); % Current sequence: index of Recording by that Microphone
    clusterMics = data_struct.MicRec_sorted(Cluster,1); % The indices of Microphones in the current Cluster
    closestRecLeft = []; 
    closestRecLeft_index = [];

    N = data_struct.Nsteps(startResolution,:); % Gather the length of each sequence in current resolution level
    
    clusterRecs = data_struct.MicRec_sorted(Cluster,2);
    % First check the sequences from the same microphone
    %   1- Sequences recorded "before" the current sequence
    tempInd = find(clusterRecs(ind)<currentRec);
    if ~isempty(tempInd) % There are some recordings in the Cluster that are recorded before the current index
        % The constraint on the left
        [mx, index] = max(clusterRecs(ind(tempInd)));
        closestRecLeft_index = ind(tempInd(index));

        closestRecLeft = Cluster(closestRecLeft_index);                
        r_closestRecLeft = r_fixed_Low(closestRecLeft_index);
        ind1 = find(data_struct.MicRec_sorted(:,2)>=data_struct.MicRec_sorted(closestRecLeft,2) ...
                  & data_struct.MicRec_sorted(:,2)<data_struct.MicRec_sorted(currentSequence,2));
        ind2 = find(data_struct.MicRec_sorted(ind1,1) == data_struct.MicRec_sorted(closestRecLeft,1));
        N_sequencesLeft = N(ind1(ind2));    
    else % There are no recordings in the Cluster that are recorded before the current index
        r_closestRecLeft = 1;
        N_sequencesLeft = 0;
    end

function [r_closestRecRight, N_sequencesRight, closestRecRight, closestRecRight_index] = find_closest_Rec_right(currentSequence, startResolution, Cluster, r_fixed, ind)
    global data_struct;
    
    L = data_struct.Lsteps(startResolution);
    r_fixed_Low = ceil(r_fixed/L);
    
    currentMic = data_struct.MicRec_sorted(currentSequence,1); % Current sequence: index of Microphone
    currentRec = data_struct.MicRec_sorted(currentSequence,2); % Current sequence: index of Recording by that Microphone
    clusterMics = data_struct.MicRec_sorted(Cluster,1); % The indices of Microphones in the current Cluster
    closestRecRight = []; 
    closestRecRight_index = [];

    N = data_struct.Nsteps(startResolution,:); % Gather the length of each sequence in current resolution level
    
    clusterRecs = data_struct.MicRec_sorted(Cluster,2);
    % First check the sequences from the same microphone
    %   1- Sequences recorded "after" the current sequence
    tempInd = find(clusterRecs(ind)>currentRec);
    if ~isempty(tempInd) % There are some recordings in the Cluster that are recorded after the current index
        % The constraint on the right
        [mx, index] = max(clusterRecs(ind(tempInd)));
        closestRecRight_index = ind(tempInd(index));

        closestRecRight = Cluster(closestRecRight_index);                
        r_closestRecRight = r_fixed_Low(closestRecRight_index);
        ind1 = find(data_struct.MicRec_sorted(:,2)>=data_struct.MicRec_sorted(closestRecRight,2) ...
                  & data_struct.MicRec_sorted(:,2)<data_struct.MicRec_sorted(currentSequence,2));
        ind2 = find(data_struct.MicRec_sorted(ind1,1) == data_struct.MicRec_sorted(closestRecRight,1));
        N_sequencesRight = N(ind1(ind2));    
    else % There are no recordings in the Cluster that are recorded after the current index
        r_closestRecRight = 1;
        N_sequencesRight = 0;
    end
    
function [alignment_possible, r_start, r_end, closest] = apply_constraints(startResolution, Cluster, r_fixed, p, s)
    global data_struct;
    N = data_struct.Nsteps(startResolution,:); % Gather the length of each sequence in current resolution level
    L = data_struct.Lsteps(startResolution);
    r_fixed_Low = ceil(r_fixed/L); % The fixed offset values in the current resolution
        
    closest = cell(1,4); 

    % alignment_possible = 1 -> The sequence can't be aligned with current cluster 
    % alignment_possible = 0 -> The sequence might be aligned with current cluster
    % Start with assuming it can be aligned
    alignment_possible = true; 
    
    r_start = 1;
    length_of_cluster = max(r_fixed_Low + N(Cluster))-min(r_fixed_Low); % Length of current cluster
    r_end = N(p(s)) + length_of_cluster; % The end of search space is length of the cluster + length of sequence
    
    currentSequence = p(s); % Current sequence index
    currentMic = data_struct.MicRec_sorted(currentSequence,1); % Current sequence: index of Microphone
    currentRec = data_struct.MicRec_sorted(currentSequence,2); % Current sequence: index of Recording by that Microphone
    clusterMics = data_struct.MicRec_sorted(Cluster,1); % The indices of Microphones in the current Cluster
    
    if length(Cluster) == 1 % The case: Fresh cluster (trying to align two sequences)
        % The Microphone indices of the recordings are the same -> Constraint 1
        % See Manuscript: Sec. Experimental Results
        if clusterMics == currentMic  
            alignment_possible = false;            
        end 
        %closest = [0 0 0 0]; % This vector will be useless if alignment is not possible           
    else  % The case: Mature cluster (trying to align a sequence to a group of pre-aligned sequences)
        ind = find(clusterMics == currentMic);
        if ~isempty(ind)        
            [r_closestRecLeft, N_sequencesLeft, closestRecLeft, closestRecLeft_index] = find_closest_Rec_left(currentSequence, startResolution, Cluster, r_fixed, ind);
            r_start = r_closestRecLeft + sum(N_sequencesLeft) + 1;
            
            [r_closestRecRight, N_sequencesRight, closestRecRight, closestRecRight_index] = find_closest_Rec_right(currentSequence, startResolution, Cluster, r_fixed, ind);
            r_end = r_closestRecRight - N(currentSequence) ;
            
            if r_start >= r_end  % If the first possible alignment is out of search space
                alignment_possible = false;
                %closest = [0 0 0 0]; % This vector will be useless if alignment is not possible
            else
                closest{1} = closestRecLeft; 
                closest{2} = closestRecRight; 
                closest{3} = closestRecLeft_index;  
                closest{4} = closestRecRight_index;                 
            end          
        %else
        %    closest = [0 0 0 0]; % This vector will be useless if no recording with the same microphone index exist in the cluster            
        end
        % An extra check if r_start is smaller than r_end
        if r_end <= r_start
            alignment_possible = false;
            %closest = [0 0 0 0]; % This vector will be useless if alignment is not possible
        end        
    end
    
function enough_samples = check_number_of_samples(r_start, r_end, startResolution)
    global data_struct
    % If the number of samples is below the threshold (minNumberOfSamples)
    % then we start with a higher resolution to have enough samples
    numberOfSamples = r_end - r_start + 1;
    if numberOfSamples < data_struct.minNumberOfSamples
        enough_samples = false;    
    else
        enough_samples = true;
    end        
    
function new_startResolution = re_set_resolution(r_start, r_end, startResolution)
    global data_struct
    old_startResolution = startResolution;
    numberOfSamples = r_end - r_start + 1;
    new_startResolution = min([old_startResolution + ceil(log2(data_struct.minNumberOfSamples/numberOfSamples)) data_struct.numSteps]);
            
function [alignment_possible, r_start, r_end, startResolution, closest] = set_search_space(startResolution, Cluster, r_fixed, p, s)    
    
    [alignment_possible, r_start, r_end, closest] = apply_constraints(startResolution, Cluster, r_fixed, p, s);
    % After applying the constraints and decide that the Cluster and the current sequence can be aligned, 
    % one still needs to check if the phi() function has enough samples at the current resolution. 
    if alignment_possible
        enough_samples = check_number_of_samples(r_start, r_end);
        if ~enough_samples % If there is not enough samples re-set the resolution and r_start, r_end variables
            startResolution = re_set_resolution(r_start, r_end, startResolution);        
            % Note that by re-calling apply_constraints function, the start and end points of search
            % space are automatically set,
            [alignment_possible, r_start, r_end, closest] = apply_constraints(startResolution, Cluster, r_fixed, p, s);
        end
    end
    
function all_sequences_exist = check_all_sequences_exist(startResolution, Cluster, p, s)
    global data_struct;
    
    if sum([data_struct.Nsteps(startResolution,Cluster) data_struct.Nsteps(startResolution,p(s))] == 0) == 0
        all_sequences_exist = true;
    else
        all_sequences_exist = false;
    end

function startResolution = set_resolution(Cluster, p, s)
    global data_struct;
    
    startResolution = 1; % For each sequence, start with the lowest possible resolution
        
    % Note that some sequence are very long i.e., mic4_rec5.wav -> 17.57
    % and some are very short i.e., mic1_rec1.wav -> 00.29. Hence some lower level
    % resolutions for short sequences don't exist. We identify such cases from the 
    % length information of the k'th sequence N(k) by the following way; If none 
    % of the lengths at the current resolution are zero, it means all the sequences
    % to be aligned exist at that resolution. 
    all_sequences_exist = check_all_sequences_exist(startResolution, Cluster, p, s);
    if ~all_sequences_exist 
        % Increase the resolution level upto where all the sequences to be aligned exist
        all_sequences_exist = false;
        while(~all_sequences_exist)
            startResolution = startResolution + 1;
            all_sequences_exist = check_all_sequences_exist(startResolution, Cluster, p, s);
        end
    end
    