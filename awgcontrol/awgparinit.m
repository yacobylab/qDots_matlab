function awgparinit(file, pulses, plsdef, plsinf, pardef, varpar, fixdef, trafofn, ind)
% awgparinit(file, pulses, plsdef, plsinf, pardef, varpar, fixdef, trafofn)
% file: name of file to store data
% pulses: pulse indices (or first index)
% plsdef: pulse definition to be modified and passed to awgmakeinf
% plsinf: pinf struct to be passed to awgmakeinf
% pardef: n x 2 matrix with pulsedef and data entry indices, 
%         -ve second indices refer to time, +ve ones val.
% varpar: length(pulses) x m matrix with parameters varying between pulses,  
%         nan entries are ignored, used for last m parameters.
% fixdef: vector with default parameters, length determines number of parameters.
% trafofn (optinal): transforms given parameters into pulse parameters. Can be a function or matrix.
% ind: index to pulse group. Each group uses the same parameters, but different
%      pulse definitions and transformations.

% (c) 2010 Hendrik Bluhm.  Please see LICENSE and COPYRIGHT information in plssetup.m.


if length(pulses) == 1 && ~isempty(varpar)
    pulses = pulses + (0:size(varpar, 2)-1);
end

if nargin < 8
    trafofn = [];
end

file = ['pp_', file];

if exist(file, 'file') || exist([file, '.mat'], 'file')
    load(file);
    if length(fixdef) ~= nparams
       error('Number of parameters not consistent.');
    end
    if nargin < 9
        ind = length(plsdata)+1;
    end
else
    nparams = length(fixdef);
    logdata = struct([]);
    plddata = struct([]);
    ind = 1;
end

if isempty(varpar) % mostly useful for single pulse.
    varpar = zeros(length(pulses), 0);
end

plsdata(ind).pulses = pulses;
plsdata(ind).plsdef = plsdef;
plsdata(ind).plsinf = plsinf;
plsdata(ind).pardef = pardef;
plsdata(ind).varpar = varpar; 
plsdata(ind).fixdef = fixdef;
plsdata(ind).trafofn = trafofn;

save(file,  'plsdata', 'nparams', 'logdata');
