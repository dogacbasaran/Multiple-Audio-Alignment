tic
current_folder = pwd;
if isempty(strfind(current_folder,'/'))==1 % OS is Windows
    parent_folder = current_folder(1:strfind(current_folder,'\SMC_based_alignment'));
    load_path = [parent_folder 'audio_data\'];    
else % OS is Linux
    parent_folder = current_folder(1:strfind(current_folder,'/SMC_based_alignment'));
    load_path = [parent_folder 'audio_data/'];    
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
    clk = clock;
    date_of_experiment = [num2str(clk(3)) '_' num2str(clk(2)) '_' num2str(clk(1))];
    fileID = fopen(['offset_estimation_fingerprinting_thr_' num2str(thr) '_result_' date_of_experiment '.txt'],'w');

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
