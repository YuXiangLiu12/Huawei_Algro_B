% Run this script from any working directory. Inputs stay beside this file.
projectRoot=fileparts(mfilename('fullpath'));
addpath(fullfile(projectRoot,'q1'));
out=q1_main(projectRoot,'revised');
