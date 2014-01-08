function atprint(ctrl, ind)
% atprint(ctrl, ind)
% ctrl: run index | 'grad' | 'basis'| 'delta' | 'deltavirt' | 'gboth'
% ind run index (if not given in first argument). <=0 means relative to
% end, default = 0

% (c) 2010 Hendrik Bluhm.  Please see LICENSE and COPYRIGHT information in plssetup.m.


global tunedata;
global smdata;

if nargin < 1
    ctrl = 0;
end

if isnumeric(ctrl)
    ind = ctrl;
    ctrl = '';
end
if(~exist('ind','var'))
    ind=0;
end
if ind <= 0 
    ind = ind + length(tunedata.runs);
end

ng = length(tunedata.gatechan);

switch ctrl
    case 'gboth'
        atprint('',ind);
        atprint('grad',ind);
        
    case ''

    fprintf('%10s', '');
    fprintf('%20s', 'LL triple (V)', 'leads (mus)', 'T1 (mus)', 'Tload (ns)', 't_c (mV)')
    fprintf('\n------------------------------------------------------------------------------------------------------------------------\n')
    for i = ind
        fprintf('%10s  (%.4f, %.4f)     (%6.3f %6.3f)%20.2f%20.0f%20.3f\n', '', tunedata.runs(i).vals .* [ones(1, 5), 1e3 1e3] )
    end
    fprintf('------------------------------------------------------------------------------------------------------------------------\n')


    case 'grad'
        if ~isempty(tunedata.runs(ind).grad)
        for i = 1:ng
            
            if isempty(smdata)
                fprintf('%-10s:', '')
            else
                fprintf('%-10s:', smdata.channels(tunedata.gatechan(i)).name);
            end
            fprintf(' (%.4f, %.4f)     (%6.3f %6.3f)%20.2f%20.0f%20.3f\n', tunedata.runs(ind).grad(i, :) .* [ones(1, 5), 1e3 1e3] )
        end
        else
            warning('tunedata.runs(%i).grad is empty\n',ind);
        end
    
    case 'basis'   
        if isempty(smdata)
            fprintf(['%-10s', repmat('%8d', 1, ng), '\n'], '', 1:ng);
            %fprintf('%-10s:', '')

        else
            fprintf(['%-10s',repmat('%-8s', 1, ng+1), '\n'], '',smdata.channels(tunedata.gatechan).name);
            fprintf('\n');
        end
        fprintf('------------------------------------------------------------------------------------------------------------------------\n')
        for i = 1:size(tunedata.basis,2)
            fprintf(['%-9s:', repmat('%8.3g', 1, ng), '\n'], tunedata.basenames{i}, tunedata.basis(:, i));
        end

    case 'delta'
        d = tunedata.runs(ind(1)).gates - tunedata.runs(ind(2)).gates;
        if isempty(smdata)
            fprintf([repmat('%9d', 1, ng), '\n'], 1:ng);
            %fprintf('%-10s:', '')

        else
            fprintf([repmat('%-9s', 1, ng), '\n'], smdata.channels(tunedata.gatechan).name);
            fprintf('\n');
        end
        fprintf('------------------------------------------------------------------------------------------------------------------------\n')

        fprintf([repmat('%9.3g', 1, ng), '\n'], d);
        
    case 'deltavirt'
        d = inv(tunedata.basis) *(tunedata.runs(ind(1)).gates - tunedata.runs(ind(2)).gates)';
        
        fprintf([repmat('%11s', 1, ng), '\n'], tunedata.basenames{:});
        fprintf('\n');
        fprintf('------------------------------------------------------------------------------------------------------------------------------------\n')

        fprintf([repmat('%11.3g', 1, ng), '\n'], d);

end
