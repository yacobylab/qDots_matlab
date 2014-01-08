function pd = pdload(name, opts)
%function pd = pdload(name, opts)
% Load most recent entry in a pulse dictionary.  Load the entire 
% dictionary if opts='all'
% if opts is a number, load the most recent dictionary before that.
%   where opts is a time as returned by 'now' or 'getscantime'

% (c) 2010 Hendrik Bluhm.  Please see LICENSE and COPYRIGHT information in plssetup.m.

global plsdata;

if isstruct(name) % This allows pdload(pload('x')) to be equiv to pdload('x')
    pd=name;
    return;
end

if ~exist('opts','var')
    opts = '';
end


if isempty(strfind(opts,'all')) && (isempty(opts) || ~isnumeric(opts))
  load([plsdata.grpdir, 'pd_', name,'_last']);
else
  load([plsdata.grpdir, 'pd_', name]);
  if isnumeric(opts) && ~isempty(opts)
     times=cellfun(@(x) getfield(x,'time'),pd);
     i=find(times < opts,1,'last');
     if isempty(i)
        error('cannot find dictionary from specified time'); 
     else
        pd=pd{i};
     end
  end  
end
    
return
