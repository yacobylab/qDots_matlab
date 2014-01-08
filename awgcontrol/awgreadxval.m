function awgreadxval
% read xvals from file to awgdata.xval 

% (c) 2010 Hendrik Bluhm.  Please see LICENSE and COPYRIGHT information in plssetup.m.


global awgdata;

load(awgdata.datafile, 'pulsedata');
%awgdata.xval = [pulsedata.xval];
awgdata.xval = zeros(1, length(pulsedata));

for i = 1:length(pulsedata)
    awgdata.xval(i) = pulsedata(i).xval(1);
end


