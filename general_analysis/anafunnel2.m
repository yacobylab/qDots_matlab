function [fp, b, fit] = anafunnel2(col, ds, rng)

% (c) 2010 Hendrik Bluhm.  Please see LICENSE and COPYRIGHT information in plssetup.m.


global dvdata;
global awgdata;


d = dvplot(col, ds);


b = d.y;
load(dvdata.collections{col}.datasets(ds).file{1}, 'scan');
x = awgdata.xval(scan.data.pulsegroups(end).pulses);

if nargin >= 3
    d.data = d.data(:, rng(1):rng(2));   
end


n = size(d.data, 1);

[m, mi] = max(d.data, [], 2);
fp = fitwrap('plinit plfit ', x, d.data, [mean(d.data, 2), repmat(1e-3, n, 1), x(mi)', repmat([.05, 0], n, 1)], @fitfn);

fit = polyfit(b, fp(:, 3), 3);

figure(10)
plot(b, fp(:, 3), b, polyval(fit, b));


function y = fitfn(beta, x) 
y = beta(1) + beta(2) * exp(-((x-beta(3))./ beta(4)).^2 ) + beta(5)*x;
