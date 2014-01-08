function awgupdate(file, params, logind)
% awgupdate(file, params, logind)
% Update parametrized pulse group.
% Missing or nan entries or params are taken from default values.
% non-nan variable defaults always supersede other parameters.

% (c) 2010 Hendrik Bluhm.  Please see LICENSE and COPYRIGHT information in plssetup.m.


load(['pp_', file]);


for k = 1:length(plsdata)
    % constant defaults
    ind = [find(isnan(params)), length(params)+1:nparams];
    params(ind) = plsdata(k).fixdef(ind);


    for i = 1:length(plsdata(k).pulses)
        ind = isnan(plsdata(k).varpar(i, :));

        % defaults varying between pulses
        params(end-size(plsdata(k).varpar, 2)+find(~ind)) = plsdata(k).varpar(i, ~ind);
        params(end-size(plsdata(k).varpar, 2)+find(ind)) = plsdata(k).fixdef(end-size(plsdata(k).varpar, 2)+find(ind));

        if ~isempty(plsdata(k).trafofn)
            if isreal(plsdata(k).trafofn)
                params2 = plsdata(k).trafofn * params';
            else
                params2 = plsdata(k).trafofn(params);
            end
        else
            params2 = params;
        end
        plsdef = plsdata(k).plsdef;
        
        for j = 1:size(plsdata(k).pardef, 1)
            if plsdata(k).pardef(j, 2) < 0
                plsdef(plsdata(k).pardef(j, 1)).time(-plsdata(k).pardef(j, 2)) = params2(j);
            else
                plsdef(plsdata(k).pardef(j, 1)).val(plsdata(k).pardef(j, 2)) = params2(j);
            end
        end
        pinf = awgmakeinf(plsdef, plsdata(k).plsinf);
        awgmakepulse('', pinf, plsdata(k).pulses(i));
        %pause
    end
end

if nargin < 3
    logind = length(logdata) + 1;
elseif logind <= 0
    logind = length(logdata) + logind;
end
    

logdata(logind).time = now;
logdata(logind).params = params;

save(['pp_', file],  'plsdata', 'nparams', 'logdata');
logentry('Updated pulses in %s.', file);
