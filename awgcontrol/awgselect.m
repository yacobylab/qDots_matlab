function awgselect(type, ind, targch, srcch)
% awgselect(type, ind, targch, srcch)
% program output channels.
% 
% type: 'pulse' (from database. Selects continuous run mode.)
%       'file' (from AWG hard drive, sets continuous run mode)
%       'seq'  (same as file, but  selects enhanced run mode).
% ind:  pulse/sequence index or file name(s)
% targch: output channel(s). Default = 1 or 1:2
% srcch: file channel (for type = 'pulse'). Default = 1 or 1:2.

% (c) 2010 Hendrik Bluhm.  Please see LICENSE and COPYRIGHT information in plssetup.m.


global awgdata;

switch type
    case 'pulse'
        load(awgdata.datafile);

        name = pulsedata(ind).name;
        if ~isempty(name)
            name = ['_', name];
        end

        if nargin < 4
            srcch = 1:size(pulsedata(ind).pulsetab, 1) - 1;
        end

        if nargin < 3
            targch = srcch;
        end

        for i = 1:length(srcch)
            %fprintf(awgdata.awg, '%s\n', sprintf('SOUR%d:FUNC:USER "/wfm/p%06d%s_%d.wfm"', ...
            fprintf(awgdata.awg, '%s\n', sprintf('SOUR%d:FUNC:USER "p%06d%s_%d.wfm"', ...
                targch(i), ind, name, srcch(i)));
            pause(.1);
        end
        fprintf(awgdata.awg, 'AWGC:RMOD CONT');
        logentry('AWG: Loaded pulse %d, %s.', ind, name(2:end));
        
    case {'file', 'seq'}
        if ~ischar(ind) 
            fprintf('Second argument must be (a) filename(s). Aborting.\n');
            return;
        end
        
        if nargin < 3
            targch = 1:size(ind, 1);
        end
        
        if strcmp(type, 'seq')
            fprintf(awgdata.awg, 'AWGC:RMOD ENH');
            %fprintf(awgdata.awg, 'AWGC:RMOD SEQ', 'ASYNC'); % for 7000 series
        else
            fprintf(awgdata.awg, 'AWGC:RMOD CONT');
        end

        
        for i = 1:size(ind, 1);
            fprintf(awgdata.awg, '%s\n', sprintf('SOUR%d:FUNC:USER "%s"', targch(i), ind(i, :)));       
        end
        ind = cellstr(ind);
        logentry('AWG: Loaded file(s) %s.', sprintf('%s ', ind{:}));

       
    otherwise
        fprintf('Invalid first argument. Use "pulse", "file" or "seq"\n'); 
end
