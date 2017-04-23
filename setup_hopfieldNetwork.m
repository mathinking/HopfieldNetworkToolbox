% Run this file to bring the folders to the path. This is not needed if you
% install the Toolbox directly.

myPath = fileparts(which('setup_hopfieldNetwork'));
addpath(fullfile(myPath));
addpath(genpath(fullfile(myPath,'chn')));
addpath(genpath(fullfile(myPath,'data','TSPFiles')));
addpath(fullfile(myPath,'help'));
disp('Continuous Hopfield Network Toolbox Installation completed.')
disp('Read the documentation for further instructions.')

savepath;
