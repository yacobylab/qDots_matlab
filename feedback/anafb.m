function [t, fbval, ctrlval, corr, dt] = anafb(file, trng, step)

% (c) 2010 Hendrik Bluhm.  Please see LICENSE and COPYRIGHT information in plssetup.m.


global awgdata;

load(file, 'scan', 'data')

t = (data{2}(:, 1)-data{2}(1, 1))*24 * 3600;


if nargin < 2
    trng = [0, Inf];
end


if nargin < 3
    if 0 && scan.data.pulsegroups(1).npulse > size(data{1}, 2) * 3
        step = 3;
    else
        step = 1;
    end
end

mask = t >= trng(1) & t <= trng(2) & isfinite(data{3}(:, 1));

fbval = data{3}(mask, 1);
ctrlval = data{3}(mask, 3);

if std(ctrlval) < 1e-10 && size(data{3}, 2) >=7 ;
    ctrlval = data{3}(mask, 7);
end

mask2 = mask & data{3}(:, 6)==0;
mask = mask & data{3}(:, 6) == 1;

  
f0 = 1000/(diff(awgdata.xval(scan.data.pulsegroups(1).pulses(step*(2:3)+1))));
% convert to mus^-1

figure(5);

subplot(311);
plot(t(mask), data{3}(mask, 1)*f0, 'b.-', t(mask), data{3}(mask, 2)*f0, 'r', ...
    t(mask2), data{3}(mask2, 1)*f0, 'b.--', t(mask2), data{3}(mask2, 2)*f0, 'r--')
sd = std(fbval)*f0;
title(sprintf('\\sigma = %.2f MHz, T_2^* = %.1f ns', sd, 1000/(sqrt(2) * pi * sd)));
ylabel('f (MHz)');

subplot(312)
plot(t(mask), data{3}(mask, 3), 'b', t(mask2), data{3}(mask2, 3), 'b--')
ylabel('cntrlval');

subplot(313)
if size(data{3}, 2) >= 7
    plot(t(mask), data{3}(mask, 4:5), t(mask|mask2), data{3}(mask|mask2, 7:end))
else
    plot(t(mask), data{3}(mask, 4:5))
end
ylabel('pulseind')
xlabel('t (s)');
t = t(mask|mask2);

if ~exist('xcorr')
    return;
end

ml = 500;
dt = polyfit(1:length(t), t', 1);
dt = dt(1);
corr = xcorr(fbval- mean(fbval), ml, 'unbiased')*f0^2;

if nargout == 0
    figure(6);
    plot((-ml:ml)*dt, (corr(ml+1)-corr), '.-');
    xlabel('\deltat (s)')
    ylabel('<(f(t+\Deltat)-f(t))^2> (MHz^2)')
end


