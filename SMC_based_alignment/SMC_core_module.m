%% SMC core module: The main implementation of the SMC procedure
%
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

function [r_max, Num_overlapping_frames] = SMC_core_(s, S, F, Nsteps,  wsteps, phi_L_low, Cluster, r_fixed, startResolution, r_start, r_end, closestRecLeft, closestRecRight , closestRecLeft_index, closestRecRight_index)

% Note: The first value of wsteps is for the lowest resolution where each offset is computed.
% That's why there are length(wsteps)-1 low resolution steps.
numSteps = length(wsteps)-1; % Number of low resolution steps 
endResolution = numSteps + 1; % Highest resolution

% The values of Phi(r) are the samples at the lowest resolution. The procedure starts from 1 resolution higher
if startResolution ~= endResolution
    startResolution = startResolution + 1; 
end

% The sparse smoothing kernel parameters. For details see Manuscript Section 3.2 Forward Markov Transition Kernel Design
winlen = [16 12 8 6 4 3 2 1]; % The choices of window length  for the smoothing kernel (Averaging)
hoplen = [8  4  4 2 2 1 1 1]; % The hop length between two averaging windows
% The index of window and hop lengths. As an example, if there are 10 resolution levels i.e., numSteps = 10
% Lprec_ind = [5 5 5 5 5 5 4 3 2 1]
% For Lprec_ind(1) = 5  ->  windowlength = winlen(5) = 4, hoplength = hoplen(5) = 2
Lprec_ind = [5*ones(1,numSteps-4) 4 3 2 1]; 
 
% r = 1 represents the non-overlapping alignment. That is just computed at the highest resolution
% Therefore in the intermediate resolutions the samples don't go to r=1.
if r_start == 1
    r_start =2;
end
NumParticles = r_end - r_start + 1; % Number of particles 
R = r_start:r_end; % Initial values of particles ( 2 3 ..) for the lowest resolution
fprintf('\nThe number of particles: %d \n',NumParticles);

% The initial distribution is shown as Pi_1
Pi_1 = phi_L_low(R)*0.001; % Annealing the initial stationary distribution
Pi_1 = Pi_1 - max(Pi_1);
Pi_1 = exp(Pi_1)/sum(exp(Pi_1));

weight = Pi_1; % Initial weights

for step = startResolution:endResolution  % Low resolution steps 1/256, 1/128, 1/64...
    
    P = winlen(Lprec_ind(step-1):end); % Window lengths at each intermediate distribution between resolution changes
    H = hoplen(Lprec_ind(step-1):end); % Hop lengths at each intermediate distribution between resolution changes

    w = wsteps(step); % The hyperparameter w is set according to annealing strategy
    N = Nsteps(step,:); % Lengths of sources at current resolution level
    
    L = 2^(endResolution -step); % Current resolution L
    
    % Alignments(offsets) of the fixed sequences for resolution L    
    r_fixed_L = r_fixed/(L);
    [mx,ind] = min(r_fixed_L);
    r_fixed_L(ind) = 1; % Note: Minimum value of r_fixed is always 1, hence it is fixed to 1.
    r_fixed_L = ceil(r_fixed_L);        
    
    T = max(r_fixed_L + N(Cluster)) + 2*N(s) - 1; % Length of the hidden feature vector lambda

    % S1 and S2 are precomputed values for effective computation of likelihood
    [S1, S2, temp_S1, temp_S2, length_of_likelihood] = compute_intermediate_values_for_likelihood(T, F, Cluster, r_fixed_L, N, s, S, step);
    
    % Set the boundaries of the offsets as the interval [r_start:r_end]
    [r_start, r_end] = set_offset_search_interval(r_fixed_L, N, s, length_of_likelihood, closestRecLeft, closestRecRight , closestRecLeft_index, closestRecRight_index);
    phi_L = zeros(1,length_of_likelihood); % The initializatin of log likelihood
    
    % Move each sample with index n to the next resolution level
    for n = 1:length(R)                
    
        % Proposed alignment values for intermediate distributions between resolutions
        [r_proposed, len] = propose_offsets_from_higher_resolution(R,n,r_start,r_end, P, H);
    
        % All the likelihood values phi_L(r_proposed) are computed for the proposed alignment values r_proposed.
        % This is more efficient because at each forward move of a sample, same likelihood values are computed
        % as shown in Manuscript: Section 3.2 Forward Markov Kernel Design
        phi_L = compute_phi_L_of_proposed_offsets(phi_L, length_of_likelihood, r_proposed, temp_S1, temp_S2, N, s, S, step, F, L, w);            

        % Sample with index n is moved to higher resolution
        [R(n), weight(n)] = forward_kernel(L, r_proposed, phi_L, Pi_1, P, H, weight, R, n);
    end
    
    if L == 1 % In the highest resolution compute phi_L(r=1) to compare all the results with a non-overlapping
        phi_L(1) = compute_non_overlapping_probability(temp_S1, temp_S2, N, s, S, step, L, F, w);
        R = find(phi_L~=0); % Samples at the original resolution
    else
        weight = check_effective_sample_size(weight, NumParticles, R);    
    end
    
    [mx,ind] = max(phi_L(R));
    r_max = R(ind);
    disp(['Optimum alignment for L=' num2str(L) ' is r=', num2str(R(ind)) , ' or ' num2str((R(ind)-N(s)-1)*L*0.02) , ' in secs'])
end

Num_overlapping_frames = find_amount_of_overlap(temp_S1, r_max, N, s);


function [S1, S2, temp_S1, temp_S2, length_of_likelihood] = compute_intermediate_values_for_likelihood(T, F, Cluster, r_fixed_L, N, s, S, step)
    % Intermediate values for computing likelihood function
    S1 = zeros(1,T); % Holds the number of sequence coefficients aligned at each time instant in the interval [1:T]
    S2 = zeros(F,T); % Holds the sum of aligned coefficients aligned at each time instant in the interval [1:T]
    for k = 1:length(Cluster)
        S1(r_fixed_L(k):r_fixed_L(k)+N(Cluster(k))-1) = S1(r_fixed_L(k):r_fixed_L(k)+N(Cluster(k))-1) + 1;        
        S2(:,r_fixed_L(k):r_fixed_L(k)+N(Cluster(k))-1) = S2(:,r_fixed_L(k):r_fixed_L(k)+N(Cluster(k))-1) + squeeze(S{step,Cluster(k)});
    end
    temp_S1 = circshift(S1,[0,N(s)]); 
    temp_S2 = circshift(S2,[0,N(s)]);
    
    length_of_likelihood = N(s)+sum(S1~=0); % Length of search space
    


function Num_overlapping_frames = find_amount_of_overlap(temp_S1, r_max, N, s)
    temp1 = zeros(size(temp_S1));
    temp1(temp_S1>0)= 1;
    temp1(r_max:r_max+N(s)-1) = temp1(r_max:r_max+N(s)-1) + 1;
    Num_overlapping_frames =  sum(temp1>1); 
    
function weight = check_effective_sample_size(weight, NumParticles, R)
    W = weight/sum(weight);
    ESS = 1/sum(W.^2); % Effective sample size
    % Check the effective sample size
    if NumParticles/2 > ESS 
        % Resampling
        tempR = randgen(R,NumParticles,W);
        R = sort(tempR);
        weight = ones(1,NumParticles)*1/NumParticles;
        disp(['Resampling at stage L=' num2str(L)]);
    end
    
function [r_start, r_end] = set_offset_search_interval(r_fixed_L, N, s, length_of_likelihood, closestRecLeft, closestRecRight , closestRecLeft_index, closestRecRight_index)

    % The offset search boundaries are set with interval [r_start:r_end] in the current resolution level
    if isempty(closestRecLeft)
        r_start = 2;
    else
        r_start = r_fixed_L(closestRecLeft_index) + N(closestRecLeft) + 1;
    end
    if isempty(closestRecRight)
        r_end = length_of_likelihood;
    else
        r_end = r_fixed_L(closestRecRight_index) - N(s) ;        
    end
    
function [r_proposed, len] = propose_offsets_from_higher_resolution(R,n,r_start,r_end,P,H)
    % Proposed alignment values for intermediate distributions between
    % resolutions      
    len = P(1)+H(1); % The length of the region to be smoothed
    r_proposed = (1:len)-len/2+(2*R(n)-1);
    if r_proposed(1)< r_start
        r_proposed = r_proposed + r_start - r_proposed(1);
    elseif r_proposed(end) > r_end
        r_proposed = r_proposed-(r_proposed(end)-r_end);
    end

function phi_L = compute_phi_L_of_proposed_offsets(phi_L, length_of_likelihood, r_proposed, temp_S1, temp_S2,N , s, S, step, F, L, w)
    % All the likelihood values phi_L(r_proposed) are computed for the proposed alignment values.
    % This is more efficient because at each forward move of a sample, same likelihood values are computed
    % as shown in Manuscript: Section 3.2 Forward Markov Kernel Design
    for i=1:length(r_proposed)
        r = r_proposed(i);   
        if phi_L(r)==0
            temp1 = temp_S1;
            temp2 = temp_S2;
            temp1(r:r+N(s)-1) = temp1(r:r+N(s)-1) + 1;
            temp2(:,r:r+N(s)-1) = temp2(:,r:r+N(s)-1) + S{step, s};

            ind = find(temp1>0);
            F1 = temp2(:,ind); % Number of ones
            F2 = (temp1(ind)'*ones(1,F)*L)' - temp2(:,ind); % Number of zeros        

            phi_L(r) = sum(sum(log(0.5 * w.^(F1) .* (1-w).^(F2) + 0.5 * w.^(F2) .* (1-w).^(F1))));                                        
        end
    end

% Values of the log likelihood into a distribution
function Fp = log_2_distribution(r, f)
    Fp = f(r);
    Fp = Fp-max(Fp);
    Fp = exp(Fp);
    Fp = Fp/sum(Fp);
    
function [r_n, weight_n] = forward_kernel(L, r_proposed, phi_L, Pi_1, P, H, weight, R, n)
    
    % Apply the procedure for all the resolution levels except the highest (original) level
    if L ~= 1       
        Fp = log_2_distribution(r_proposed, phi_L);
        rs = r_proposed;

        number_of_smooth_distributions = length(P);
        % Move the sample between same resolution smooth distributions
        for i=1:number_of_smooth_distributions 
            win1 = Fp(1:P(i)); % The probabilities of the samples from the first window
            rs1 = rs(1:P(i)); % The samples from the first window
            win2 = Fp(H(i)+1:end); % The probabilities of the samples from the second window
            rs2 = rs(H(i)+1:end); % The samples from the second window

            % Move the sample to one of the windows with a probability proportional 
            % to the sum of probabilities in each window
            rt1 = sum(win1);
            rt2 = sum(win2);
            pr = rt1/(rt1+rt2);
            if rand(1) < pr                                
                % Weight update
                weight(n) = weight(n) * rt1/(Pi_1(n));
                % Update the initial distribution for the next sample move                
                Pi_1(n) = rt1; 
                % Assign the new set of offsets for the next sample move 
                rs = rs1; 
                % Compute the probabilities of the samples                
                Fp = log_2_distribution(rs, phi_L); 
            else                
                % Weight update
                weight(n) = weight(n) * rt2/(Pi_1(n));
                % Update the initial distribution for the next sample move                                    
                Pi_1(n) = rt2;
                % Assign the new set of offsets for the next sample move                 
                rs = rs2;
                % Compute the probabilities of the samples                
                Fp = log_2_distribution(rs, phi_L);                
            end
        end
        R(n)  = rs; % Accepted move of the proposed sample
    end
    r_n = R(n);
    weight_n =weight(n);

function pr_non_overlap = compute_non_overlapping_probability(temp_S1, temp_S2, N, s, S, step, L, F, w)
    r = 1;
    temp1 = temp_S1;
    temp2 = temp_S2;
    temp1(r:r+N(s)-1) = temp1(r:r+N(s)-1) + 1;
    temp2(:,r:r+N(s)-1) = temp2(:,r:r+N(s)-1) + S{step, s};

    ind = find(temp1>0);
    F1 = temp2(:,ind); % Number of ones
    F2 = (temp1(ind)'*ones(1,F)*L)' - temp2(:,ind); % Number of zeros        

    pr_non_overlap = sum(sum(log(0.5 * w.^(F1) .* (1-w).^(F2) + 0.5 * w.^(F2) .* (1-w).^(F1))));                                        
