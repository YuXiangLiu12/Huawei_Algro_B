function out=q1_main(root,mode)
% Reproducible complete workflow. Old v0.4 results are never overwritten.
% q1_main(projectRoot) reruns both calibration stages and verification.
% q1_main(projectRoot,'resume') verifies saved v0.5 stage fits.
if nargin<1,root=fileparts(fileparts(mfilename('fullpath')));end
if nargin<2,mode='fit';end
assert(isfile(fullfile(root,'Q1','q1_config.m')),'Pass the project root, not its Q1 subfolder');
if strcmp(mode,'fit')
    q1_stage1(root);q1_stage2(root);
    q1_fit_transport(root,'direct');q1_fit_transport(root,'local');
elseif ~strcmp(mode,'resume')
    error('mode must be fit or resume');
end
out=q1_finalize(root);
end
