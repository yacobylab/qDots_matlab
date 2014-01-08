function blank_function_template(file,config)
%function jcal = jcal_fine(file, config)
% this is a template for writing nice functions
% Nice comments go here
%
% config is struct with lots of good stuff. it has:
%
%  opts: options, can be... defaults to...


if ~exist('file','var') || isempty(file)
  file=uigetfile('sm*.mat');
end

if ~exist('config','var')
    config=struct();
elseif ischar(config)
    config=struct('opts',config);
elseif iscell(config)
    config=struct(config{:});
end
config = def(config,'opts','defaults');
config = def(config,'blah',[1 Inf]);


%lots of nice code goes here
 
end %function ends here

% Apply a default.
function s=def(s,f,v)
  if(~isfield(s,f))
      s=setfield(s,f,v);
  end
return;
end

function b=isopt(config,name)
  b=~isempty(strfind(config.opts,name));
return;
end