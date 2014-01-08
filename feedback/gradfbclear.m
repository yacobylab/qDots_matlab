function gradfbinit(t)

% (c) 2010 Hendrik Bluhm.  Please see LICENSE and COPYRIGHT information in plssetup.m.


global fbdata;

if nargin < 1
    t = .001;
end
%mask = (fbdata.time(end) - fbdata.time)/60 >  abs(t);
mask = (now*24*3600 - fbdata.time)/60 >  abs(t);
mask(end) = 0;
fbdata.time(mask) = nan;
fbdata.fbval(mask) = nan;
fbdata.intval(mask) = nan;
fbdata.ctrlval(mask)  = nan;
