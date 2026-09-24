function root=q1_project_root(root)
% Resolve project-relative input/output paths, independent of MATLAB pwd.
if nargin<1 || isempty(root)
    root=fileparts(fileparts(mfilename('fullpath')));
end
root=char(root);
if isfile(fullfile(root,'q1_config.m'))
    root=fileparts(root); % Also accept the q1 source directory.
end
assert(isfolder(root),'Project directory does not exist: %s',root);
[ok,info]=fileattrib(root);assert(ok,'Cannot resolve project directory');
root=info.Name;
assert(isfile(fullfile(root,'q1','q1_config.m')), ...
    'Cannot find q1 source directory under: %s',root);
assert(isfile(fullfile(root,'附件2.xlsx')), ...
    'Cannot find 附件2.xlsx under: %s',root);
end
