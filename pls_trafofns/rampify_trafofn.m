function wfo=rampify_trafofn(wf, channel, pts)
% function wf=rampify_trafofn(wf, channel, pts)
% This function takes a waveform and replaces each point with a running
% average of the previous pts points,

% (c) 2010 Hendrik Bluhm.  Please see LICENSE and COPYRIGHT information in plssetup.m.
if pts == 0
    wfo=wf;
    return;
end
wfo=wf;
for i=1:(pts-1)
    wfo=wfo+wf(min((1:end)+i,end));
end
if pts>0
  wfo=wfo/pts;
end



