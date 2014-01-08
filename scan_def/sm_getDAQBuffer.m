function out = sm_getDAQBuffer(chan, samprate, npts)
%function out = sm_getDAQBuffer(chan, samprate, npts)
% reads DAQ channel chan, at a sampling rate of samprate for npts


if isnumeric(chan)
    chan = ['DAQ', num2str(chan)];
end
if ~ischar(chan)
   error('chan must be number or channel name'); 
end
smset('samprate',samprate)
inst = sminstlookup('ATS660');
channum = str2num(chan(end));
smcATS660v4([inst, channum, 5],npts,samprate,1);
smcATS660v4([inst, channum, 4]); %arm
smatrigAWG(16); %trigger
out=smget(chan);

