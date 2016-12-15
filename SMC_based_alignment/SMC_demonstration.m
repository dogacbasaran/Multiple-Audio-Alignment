% SMC Demonstration
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

current_folder = pwd;
if isempty(strfind(current_folder,'/'))==1 % OS is Windows
    parent_folder = current_folder(1:strfind(current_folder,'\SMC_based_alignment'));
    load_path = [parent_folder 'audio_data\'];
    
else % OS is Linux
    parent_folder = current_folder(1:strfind(current_folder,'/SMC_based_alignment'));
    load_path = [parent_folder 'audio_data/'];    
end

% Extract features from audio dataset
dataset_features = feature_extract_module(load_path);

% Align the unsyncronized audio signals (might take a while..)
[Clusters, r_clusters, time_elapsed] = SMC_main_module(dataset_features);

% The results are written in a text file under "SMC_offset_estimation_results" folder
% to be further used in the evaluation codes.
offset_estimates_for_evaluation(Clusters, r_clusters, dataset_features);

