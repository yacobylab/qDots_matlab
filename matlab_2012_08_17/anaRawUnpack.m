function rawdata=anaRawUnpack(scan, data)
%data=anaRawUnpack(scan, data)
% Unpack raw data in a scan into a more usable format.
% Output data is in data{channel}(group, pulse, rep)
if nargin == 0
   s = load(uigetfile('sm*.mat'));
   scan = s.scan;
   data = s.data;
end

nchan = length(scan.loops(1).getchan);
if iscell(scan.data.conf.datachan)
  dchan = length(scan.data.conf.datachan);
else
  dchan=1;
end
rawchans = nchan+(1:dchan);
if ndims(data{rawchans(1)}) == 2
    for l=1:dchan
        data{rawchans(l)} = permute(data{rawchans(l)},[1 3 2]);
    end
end
  rawdata=cell(dchan,1);
npulse=scan.data.pulsegroups(1).npulse(1);
for l=1:dchan
    rawdata{l}=permute(data{rawchans(l)},[3,1,2]);
    ppg=size(rawdata{l},1)*size(rawdata{l},2); % points-per-group
    rawdata{l}=reshape(rawdata{l},npulse,ppg/npulse,size(rawdata{l},3));
    rawdata{l}=permute(rawdata{l},[3 1 2]);
end
    