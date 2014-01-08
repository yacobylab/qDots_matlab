function gitinfo = get_git_info(wd)
%function gitinfo = get_git_info(wd)
% get the git information from the working directory specified in wd
%wd defaults to the working directory
if ~exist(which('git'),'file')
    error('cannot find git.m on you path');
end
if nargin==0
    wd = pwd;
end
gitinfo = [];
cd(wd);
gitinfo.SHA1=sscanf(git('show-branch --sha1-name working'),'[%d],%s');
gitstr=sscanf(git('show-branch working'),'%s]%s');%this will be [branch_name]
if strfind(gitstr,'fatal')
    gitinfo.branch = sprintf('%s is not a git repo',wd);
else
    gitinfo.branch = gitstr(2:end-1); %get rid of brackets
end
end