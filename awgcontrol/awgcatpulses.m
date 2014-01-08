function pinf = awgcatpulses(pi)

% (c) 2010 Hendrik Bluhm.  Please see LICENSE and COPYRIGHT information in plssetup.m.


global awgdata;
if ~isstruct(pi)
    pi = awgdata.pulsedata(pi);
end
pinf = pi(1);
pinf.pulsetab = [pi.pulsetab];
pinf.marktab = [pi.marktab];

k = size(pi(1).pulsetab, 2);
l = size(pi(1).marktab, 2);

for i = 2:length(pi)

    n = size(pi(i).pulsetab, 2);
    pinf.pulsetab(1, k+(1:n)) = pinf.pulsetab(1, k+(1:n)) + pinf.pulsetab(1, k);

    m = size(pi(i).marktab, 2);
    pinf.marktab(1, l+(1:m)) = pinf.marktab(1, l+(1:m)) + pinf.pulsetab(1, k);
    
    k = k + n;
    l = l + m;
end
