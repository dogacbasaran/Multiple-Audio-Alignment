%% Feature Extraction Module:  %%
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
%
%  Feature extraction follows the method in Haitsma et. al. "A highly robust audio fingerprinting system"
    %   Preprocessing
    %       - Read file
    %       - Stereo to Mono
    %       - Resample to 16kHz
    %       - Normalize the signal
    %       - Remove the silence from start and end of the file
    %       - Save the file in ss{k} cell
    %   Extraction
    %       - Divide the signal into overlapping frames
    %       - Take square of the magnitude spectrum to obtain spectrogram
    %       energy
    %       - Find subband energies
    %       - First difference through time on spectrogram energies
    %       - First difference through frequency
    %       - Thresholding with 0 to obtain bit values that represents signs ->
    %       0 = "-" and   1 = "+", the result is saved in S{k}        
function dataset_features = feature_extraction_module(load_path)

    filenames = dir([load_path '*.wav']);
    
    % The number of audio sequences in the dataset
    K = length(filenames); 

    % Holds microphone number and recording number in that microphone
    % MicRec is of the form  [1 1 ; 1 2 ; 1 3; .. ; 2 1 ; 2 2 ; 2 3 ;...]
    % See Manuscript:
    MicRec = zeros(K,2); 

    % Set the initial parameters for Sequential Monte Carlo Sampler
    smc_parameters = initialize_SMC_parameters();
    
    % Set the parameters for feature extraction
    feature_parameters = initialize_feature_parameters(smc_parameters.numSteps, K);

    for k = 1:K
        fprintf('\nSequence %s\n',filenames(k).name);

        % Read the file 
        [s, fs] = audioread([load_path filenames(k).name]);
        
        % Preprocessing
        s = preprocess_audio(s, fs, feature_parameters.Fs);

        % Extract microphone and recording number
        MicRec(k,:) = extract_mic_rec_number(filenames(k).name);
        % Save the audio sequences in their raw form for displaying purposes
        feature_parameters.ss{k} = s;

        % Extract the binary features for audio source
        [feature_parameters.S, feature_parameters.N] = extract_features(s, k, smc_parameters, feature_parameters); 
    end

    % Sort the binary sequences according to their length N in decreasing order
    [feature_parameters.S, feature_parameters.ss, Nsteps, Nsorted, MicRec_sorted, sorting_index] = sort_sequences(feature_parameters.S, feature_parameters.ss, feature_parameters.N, MicRec);    
    
    % Print the ordered sequences on the screen
    filenames_sorted = print_ordered_files(filenames, sorting_index);
    
    % The structure holds necessary information for SMC procedure
    dataset_features = struct('S', {feature_parameters.S}, 'ss', {feature_parameters.ss}, 'Fs',feature_parameters.Fs,'F',feature_parameters.F,...
                              'Nsteps',Nsteps, 'wsteps1',smc_parameters.wsteps1,'wsteps2',smc_parameters.wsteps2, 'K', K,...
                              'MinNum_overlapping_frames', feature_parameters.Minimum_reliable_overlap_in_frames,...
                              'Lsteps',smc_parameters.Lsteps,'MicRec_sorted',MicRec_sorted, 'MicRec', MicRec,...
                              'minNumberOfSamples', smc_parameters.minNumberOfSamples, 'Ws', feature_parameters.Ws,...
                              'N', feature_parameters.N, 'Nsorted', Nsorted, 'numSteps', smc_parameters.numSteps);

function smc_parameters = initialize_SMC_parameters(numSteps,w_min,w_max,minNumberOfSamples)

    switch nargin
        case 0
            numSteps = 11;
            w_min = 0.51;
            w_max = 0.63;
            minNumberOfSamples = 100;
        case 1
            w_min = 0.51;
            w_max = 0.63;
            minNumberOfSamples = 100;
        case 2
            w_max = 0.63;
            minNumberOfSamples = 100;
        case 3
            minNumberOfSamples = 100;
        otherwise
    end

    Lsteps = 2.^((numSteps):-1:0); % The resolution levels  ...256 128 64 32 16 8 4 2 1
    
    % Hyperparameter w at each low resolution level
    % There two choices
    %       1- wsteps1 = w changes linearly from w_min to w_max
    %       2- wsteps2 = w changes non-linearly from w_min to w_max
    wsteps1 = linspace(w_min,w_max,numSteps+1);     
    coeff = 2; % Determines how sharp the non-linear change is
    t = 0:1/numSteps:1;        
    wsteps2 = (exp(coeff*t)-1)/(exp(coeff*t(end))-1)*(w_max-w_min) + w_min; 
    
    smc_parameters = struct('numSteps', numSteps,'Lsteps', Lsteps, 'wsteps1', wsteps1, 'wsteps2',wsteps2,'minNumberOfSamples',minNumberOfSamples);

function feature_parameters = initialize_feature_parameters(numSteps, K, Fs, Ws, F, minimum_reliable_overlap_in_secs)

    if nargin < 3
        Fs = 16000; % Sampling rate
        Ws = 0.064; % Frame length in seconds
        F = 32; % Number of subbands
        minimum_reliable_overlap_in_secs = 10; % In seconds
    elseif nargin < 4
        Ws = 0.064; % Frame length in seconds
        F = 32; % Number of subbands
        minimum_reliable_overlap_in_secs = 10; % In seconds
    elseif nargin <5
        F = 32; % Number of subbands
        minimum_reliable_overlap_in_secs = 10; % In seconds
    elseif nargin <6
        minimum_reliable_overlap_in_secs = 10; % In seconds
    end
    
    frameLen = round(Ws*Fs); % Frame length in samples
    frameInc = round(frameLen/2); % Hopsize is half of the window size
    audioWindow = hamming(frameLen); % Hamming window is used in framing 
    
    subbands = round(logspace(log10(100),log10(Fs/2),F+2)); % Subbands with logarithmic spacing   
    
    % Minimum reliable overlap - See Manuscript: Sec. 5.1 Evaluation criterion
    Minimum_reliable_overlap_in_frames = minimum_reliable_overlap_in_secs*Fs/frameInc; % In frames

    N = zeros(numSteps+1,K); % The length of sequences in frames at each resolution
    % The form of N is;  ( Assuming Nk is the length of k'th sequence)       
    
    %         .           .            .                   . 
    %         .           .            .                   .
    %       N1/4    N2/4     N3/4    ...     NK/4
    %       N1/2    N2/2     N3/2    ...     NK/2
    %       N1      N2       N3      ...     NK
       
    S = cell(numSteps+1,K); % Feature sequences all in cell 'S' at each resolution
    % The form of S is; (* Assuming xk is the k'th sequence and xk_L is the low resolution of xk at resolution L)  
    
    %         .            .             .                       . 
    %         .            .             .                       .
    %       x1_4     x2_4      x3_4     ...        xK_4
    %       x1_2     x2_2      x3_2     ...        xK_2
    %       x1       x2        x3       ...        xK
            
    ss = cell(1,K); % Original signals in cell 'ss'
    
    feature_parameters = struct('Fs', Fs, 'Ws', Ws, 'F', F, 'frameLen', frameLen, 'frameInc', frameInc, 'audioWindow', audioWindow,...
                                'Minimum_reliable_overlap_in_frames',Minimum_reliable_overlap_in_frames, ...
                                'subbands', subbands, 'N', N, 'S', {S}, 'ss', {ss});

function s = preprocess_audio(s, fs, Fs)    
    % Stereo - Mono conversion
    if size(s,2) == 2 % if signal is stereo, turn into mono
        s = (s(:,1) + s(:,2))/2;
    end

    % Set the sampling frequency to Fs 
    if (fs ~= Fs)
        s = resample(s,Fs,fs); % Resample the signal to have chosen sampling rate i.e., 16KHz
    end
    
    % Normalization
    s = s/std(s);

    % Remove the silence parts from the beginning and end of sequences if exists
    ind1 = find(s~=0);
    ind2 = find(s(end:-1:1)~=0);

    start_sample = ind1(1);
    end_sample = length(s)-ind2(1)+1;

    s = s(start_sample:end_sample);

function [S, N] = extract_features(s, k, smc_parameters, feature_parameters)
    audioWindow = feature_parameters.audioWindow;
    frameInc = feature_parameters.frameInc;
    N = feature_parameters.N;
    Fs = feature_parameters.Fs;
    F = feature_parameters.F;
    S = feature_parameters.S;
    subbands = feature_parameters.subbands;
    
    numSteps = smc_parameters.numSteps;
    
    % Framing
    frames = enframe(s,audioWindow,frameInc)';
    
    % Length of k'th sequence in frames
    N(end,k) = size(frames,2)-1; 

    % Energy of each coefficient in spectrogram
    try % To avoid the memory problems with long files
        power_spectrum = abs(fft(frames,Fs/2)).^2;
    catch me 
        disp(me.message)
        power_spectrum = abs(fft(frames(:,1:floor(size(frames,2)/2)),Fs/2)).^2;
        power_spectrum = [power_spectrum abs(fft(frames(:,floor(size(frames,2)/2)+1:end),Fs/2)).^2];
    end    

    % Remove all zeros from spectrogram
    power_spectrum(power_spectrum==0) = 1;

    % Finding the energy at each subband
    energy_in_subbands = zeros(F+1,N(end,k)+1);
    for f=1:F+1
        energy_in_subbands(f,:) = sum(power_spectrum(subbands(f):subbands(f+1),1:N(end,k)+1));            
    end
    
    % First difference through time on the spectrogram
    first_difference_in_time = energy_in_subbands(2:end,:) - energy_in_subbands(1:end-1,:);
    % First difference through frequency on the spectrogram
    S{end,k} = first_difference_in_time(:,2:end) - first_difference_in_time(:,1:end-1);
    % Thresholding with 0 to obtain bit values that represents signs ->
    % 0 = "-" and   1 = "+"
    S{end,k} = double(S{end,k}>0);

    % Compute low resolution sequences for each level. See Manuscript: sec. 2 - Model
    % Each column of S has different resolution versions of the same
    % sequence
    for step = numSteps:-1:1 % Steps from high to low
        L = 2^(numSteps - step + 1);
        for l = 1:L:N(end,k)-L+1
            S{step,k}(:,(l-1)/L + 1) = sum(S{end,k}(:,l:l+L-1),2);                
        end
        N(step,k) = size(S{step,k},2);
    end          
    
function micrec = extract_mic_rec_number(filename)
    % Extract the microphone number and record number from the filenames
    MicInd = strfind(filename,'mic'); 
    RecInd = strfind(filename,'rec');
    DotInd = strfind(filename,'.');
    
    micNum = str2num(filename(MicInd+3:RecInd-2)); % Microphone number
    recNum = str2num(filename(RecInd+3:DotInd-1)); % Recording number from the same microphone
    
    % MicRec holds the microphone and record number of the file. 
    micrec = [micNum recNum];

function [S, ss, Nsteps, Nsorted, MicRec_sorted, ind] = sort_sequences(S, ss, N, MicRec)
    % Sort the songs according to their lengths in decreasing order
    [Nsorted,ind] = sort(N(end,:),'descend');
    Nsteps = N(:,ind);
    ss = ss(ind);
    S = S(:,ind);
    MicRec_sorted = MicRec(ind,:);

function filenames_sorted = print_ordered_files(filenames, ind)
    % Print the sorted filenamess on the screen
    filenames_sorted = filenames(ind);
    fprintf('\nThe sorted filenames: \n');
    for i=1:length(filenames)
        fprintf('     %d -  %s\n',i, filenames_sorted(i).name);
    end                  
