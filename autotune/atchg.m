function atchg(ind, dv)

% (c) 2010 Hendrik Bluhm.  Please see LICENSE and COPYRIGHT information in plssetup.m.

if(abs(dv) > 25e-3)
   fprintf('Crazy big change.  I don''t believe you; ignoring.\n');
   return;
end
global tunedata;

if ischar(ind)
    indchr = ind;
    ind = strmatch(ind, tunedata.basenames);
    if isempty(ind)
       error('cannot find basis direction named %s',indchr); 
    end
end

v = smget(tunedata.gatechan);
smset(tunedata.gatechan, [v{:}] + tunedata.basis(:, ind)' * dv);
return;
