function [y, b, c, fitpar] = anafunnel(col, ds, lim, cntrl)
%[y, b] = anafunnel(col, ds, lim)

% (c) 2010 Hendrik Bluhm.  Please see LICENSE and COPYRIGHT information in plssetup.m.


global imagedata;

if nargin < 3;
    lim = [-Inf, Inf];
end

if nargin < 4
    cntrl = '';
end
    
data = implot(col, ds);
rng = data.range;
data = data.data;
npix = fliplr(size(data));


x = linspace(rng(1, 1), rng(1, 2), npix(1))*1e3;
bplot = linspace(rng(2, 1), rng(2, 2), npix(2));
scan = imagedata.cache{imagedata.collections{col}.images(ds).index}.params.scan;

if strfind(cntrl, 'comb')
    loop = 3;
else
    loop = 2;
end
if isempty(scan.loops(loop).npoints)
    b = scan.loops(loop).rng;   
elseif isempty(scan.loops(loop).rng)
    b = 1:scan.loops(loop).npoints;
else
    b = linspace(scan.loops(loop).rng(1), scan.loops(loop).rng(2), scan.loops(loop).npoints);
end

mask = x > lim(1) & x < lim(2);
x = x(mask);
data = data(:, mask);


figure(40);
clf
imagesc(x, rng(2, :), data);
set(gca, 'ydir', 'normal');

[m, mind] = max(data, [], 2);
fitpar(:, 1) = median(data, 2);
fitpar(:, 2) = m - fitpar(:, 1);
fitpar(:, 3) = x(mind);
fitpar(:, 4) = .1;
if strfind(cntrl, 'slp')
    fitpar(:, 5) = 0;
end
if strfind(cntrl, 'skew')
    fitpar(:, 6) = 0;
end

figure(41)

for i = sum(~any(isnan(data), 2)):-1:1
    clf
    %plot(x, data(i, :), x, fitfn(fitpar(i, :), x));
    hold on
    fitpar(i, :) = nlinfit(x, data(i, :), @fitfn, fitpar(i, :));
    plot(x, data(i, :), x, fitfn(fitpar(i, :), x));
    if strfind(cntrl, 'pause')
        pause;
    end
end

figure(40);
hold on;
plot(fitpar(:, 3), bplot, '.');

y = fitpar(~any(isnan(data), 2), 3);
c = fitpar(~any(isnan(data), 2), 2);
b = b(~any(isnan(data), 2));

function y = fitfn(beta, x)
if length(beta) < 5
    beta(5) = 0;
end
if length(beta) < 6
    beta(6) = 0;
end

y = beta(1) + (beta(2) + + beta(6)*x).* exp(-((x-beta(3))./ beta(4)).^2 ) + beta(5)*x;
