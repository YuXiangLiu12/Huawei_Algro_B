function out=q1_main(root,mode)
% Reproducible complete workflow. Old v0.4 results are never overwritten.
% q1_main(projectRoot) reruns both calibration stages and verification.
% q1_main(projectRoot,'resume') verifies saved v0.5 stage fits.
% q1_main(projectRoot,'revised') calibrates and verifies the local phase model.
if nargin<1,root=[];end
root=q1_project_root(root);
if nargin<2,mode='fit';end
assert(isfile(fullfile(root,'q1','q1_config.m')),'Pass the project root, not its q1 subfolder');
if strcmp(mode,'revised')
    out=q1_continue(root);return;
elseif strcmp(mode,'revised_resume')
    out=q1_continue(root,0);return;
elseif strcmp(mode,'fit')
    q1_stage1(root);q1_stage2(root);
    q1_fit_transport(root,'direct');q1_fit_transport(root,'local');
elseif ~strcmp(mode,'resume')
    error('mode must be fit, resume, revised, or revised_resume');
end
out=q1_finalize(root);
end
