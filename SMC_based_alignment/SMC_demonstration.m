% SMC_DEMONSTRATION - Demonstration of multiresolution alignment based on SMC samplers
%   Demonstrates the SMC based multiresolution, multiple audio alignment
%   method

% Copyright (C) 2016  Author: Dogac Basaran

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

% Choose the directory where the audio dataset is present
load_path = uigetdir(current_folder,'Select the path to the audio dataset');

% Extract features from audio dataset
dataset_features = feature_extract_module(load_path);

% Align the unsyncronized audio signals (might take a while..), the results
% are written in a text file.
[Clusters, r_clusters, time_elapsed] = SMC_main_module(dataset_features);

