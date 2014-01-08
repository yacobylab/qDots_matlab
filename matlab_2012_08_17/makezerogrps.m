function makezerogrps(npls, pulselength, dict, chan)
%function makezerogroups(npls, pulselength)
% makes a pulse group full of zeros.
% npls is number of pulses in group
% pulselength in the length of the pulse (in microseconds).
% dict and chan could be 'right' and [3 4], for example
% if dict and chan are left blank then groups for both sides will be made
% if only one of dict or chan is given, then you're stupid and it throws an
% error

global plsdata;

if ~exist('dict','var') || isempty(dict)
   dict = {'left','right'}; 
end
if ischar(dict)
   dict = {dict}; 
end

if ~exist('chan','var') || isempty('chan')
   chan = {[2 1], [3 4]}; 
end

if isnumeric(chan)
   chan = {chan}; 
end

if length(chan) ~= length(dict)
   error('length of chan and length of dict do not match \n'); 
end

namepat = 'zeros_%02d_%d_%s';

pg.pulses=25;
pg.ctrl='loop pack'; 

for j=1:length(dict)
   pg.name = sprintf(namepat,pulselength, npls,upper(dict{j}(1))); 
   pg.dict = dict{j};
   pg.chan = chan{j};
   pg.params = pulselength;
   pg.varpar = pulselength*ones(npls,1);
   try 
       plsupdate(pg)
       fprintf('group %s already exists',pg.name);
   catch
      plsdefgrp(pg);
      fprintf('made group named %s \n',pg.name);
   end
   
end




end