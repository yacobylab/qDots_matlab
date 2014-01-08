function delta = atdelta(gates, r)
% atdelta(gates, r)

% (c) 2010 Hendrik Bluhm.  Please see LICENSE and COPYRIGHT information in plssetup.m.


global tunedata;

if nargin < 2
    r = ones(1, length(gates)-2);
else
    r(end+1:length(gates)-2) = 1;
end
gates = smchanlookup(gates);

[p, gind] = ismember(gates, tunedata.gatechan);
if ~all(p)
    error('No gradient info for some gate(s).\n');
end

grad = tunedata.runs(end).grad(:, 1:2);
delta = [r, - r*grad(gind(1:end-2), :)/grad(gind(end-1:end), :)];
