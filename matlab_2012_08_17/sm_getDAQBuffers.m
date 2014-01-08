function out = sm_getDAQBuffers(chan, samprate, npts)
%function out = sm_getDAQBuffer(chan, samprate, npts,opts)
% reads DAQ channel chan, at a sampling rate of samprate for npts
% chan is like chan for smget. it can be vector of chan numbers, cell array
% of channel names, or single channel name

if isnumeric(chan)
    for j = 1:length(chan)
        chan2{j} = ['DAQ', num2str(chan(j))];
    end
    chan = chan2;
end

if ischar(chan)
    chan = {chan};
end

if ~iscell(chan)
    error('chan must be numerical array of channel numbers or cell array of channel names \n');
end

smset('samprate',samprate); %not sure if this is needed?
inst = sminstlookup('ATS660'); %hardcoded, bad :(

smcATS660v4([inst, 1, 5],npts,samprate);%configure DAQ card
smcATS660v4([inst, 1, 4]); %arm
smdata.inst(inst).data.mask = [];

%warning, this does a software and hardware trigger
%might cause problem in teh future
smcATS660v4([inst, 1, 3]); %force trigger
smatrigAWG(16);
out=smget(chan);

