function awgmakeseq(cntrl, pulselist, filename, chans,filenames,zerofmt)
% awgmakeseq(cntrl, pulselist, filename, chans,filenames,zerofmt)
% cntrl: file, logjump
% pulselist: pulse number(1) and optional further parameters for each line:
% (one row for each line)
% repeat count(2) (default = 0, i.e. infinity)
% wait trigger(3) (default = 0)
% goto-1 (default 0)(4)
% goto line(5) (default -1 = next line)
% per-group map of what channels are zero, as int.  0011 means ch3, ch4 are
%   zeros.  Default 11,
%
% filenames, if set, overrides default filenames
% zerobase, if set, is format string for zero pulse names.  Defaults to zero_%08d_1.wfm 
% If pulselist is a row vector, it is assumed to contain pulse numbers and 
% transposed first. Thus, it is not possible to generate single line
% sequence files with parameters other than the default.

% (c) 2010 Hendrik Bluhm.  Please see LICENSE and COPYRIGHT information in plssetup.m.



global awgdata;

if ~exist('chans') || isempty(chans)
    chans = awgdata.chans;
end

if(~exist('zerofmt') || isempty(zerofmt))
   zerofmt='zero_%08d_1.wfm'; 
end
%load(awgdata.datafile);

if size(pulselist, 1) == 1
    pulselist = pulselist';
end

if ~isfield(awgdata, 'nchan') || awgdata.nchan == 0 
    nchan = max(chans); %last HW channel
else
    nchan = awgdata.nchan;
end

ndef = 6 - size(pulselist, 2);
nlines = size(pulselist, 1);

%default(end, 3) = 1; % go to first line. % does not combine with inf loop
default = [zeros(nlines, 3), -ones(nlines, 1), zeros(nlines,1)];
pulselist = [pulselist, default(:, end-ndef+1:end)];

if isfield(awgdata, 'plsgen');
    plsgen = awgdata.plsgen;
else
    plsgen = '';
end

    
buf =  sprintf('MAGIC 300%1d\r\nLINES %04d\r\n', nchan, nlines);

for i = 1:nlines
    
    name = awgdata.pulsedata(pulselist(i, 1)).name;
    if ~isempty(name)
        name = ['_', name];
    end

    if pulselist(i, 1)==awgdata.trigp
        zch=11;
    else
        if(pulselist(i,6) ~= 0) % Cheap hack to pack channel list into 1 column
            zch = pulselist(i,6);
        else                
            zch=sum(10.^(chans-1));
        end
    end
    for j = 1:nchan % loop over HW channels
        if mod(floor(zch/10^(j-1)),10)==1 
            buf = [buf, sprintf('"%s_%1d.wfm",', sprintf(fixfn(filenames{i})), j)];
        else
            %buf = [buf, ','];
            pd = awgdata.pulsedata(pulselist(i, 1));
            buf = [buf, sprintf(['"' fixfn(zerofmt) '",'], round(pd.pulsetab(1, end) * pd.tbase))];
        end
    end
    
    buf = [buf, sprintf('%05d,%1d,%1d,%04d\r\n', pulselist(i, 2:5))];
end


if strfind(cntrl, 'logjump')
    buf = [buf, sprintf('JUMP_MODE LOGIC\r\nJUMP_TIMING SYNC\r\nSTROBE 0\r\n')];
else
    buf = [buf, sprintf('JUMP_MODE SOFTWARE\r\nJUMP_TIMING SYNC\r\nSTROBE 0\r\n')];
end

out = fopen([awgdata.datadir, filename], 'w');
fwrite(out, buf);
fclose(out);
return

function f=fixfn(f)
f=strrep(f,'/','\\');
return
    
