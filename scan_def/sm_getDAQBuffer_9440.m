function out = sm_getDAQBuffer_9440(chan, samprate, npts, opts)
%function out = sm_getDAQBuffer(chan, samprate, npts)
% reads DAQ channel chan, at a sampling rate of samprate for npts
% opts can be 'softwaretrig' to use a software trigger of DAQ card.
% opts can be 'timer' to throw tic toc around the step

if ~exist('opts','var')
   opts = ''; 
end

opts = lower(opts);

force = ~isempty(strfind(opts,'softwaretrig'));

global smdata;
if isnumeric(chan)
    chan = ['DAQ', num2str(chan)];
end
if ~ischar(chan)
   error('chan must be number or channel name'); 
end
smset('samprate',samprate)
inst = sminstlookup('ATS9440');  
func = smdata.inst(inst).cntrlfn;
channum = str2num(chan(end));
func([inst, channum, 5],npts,samprate,1);

if force
  smdata.inst(inst).data.nrec = [0 0];
end
func([inst, channum, 4]); %arm

if force
  func([inst, channum, 3]); %software trigger
else
  smatrigAWG(16); %trigger
end
out=smget(chan);

