function fitpar = anajunction(col, ds, lim, trc)
%[y, b] = anafunnel(col, ds, lim)

% (c) 2010 Hendrik Bluhm.  Please see LICENSE and COPYRIGHT information in plssetup.m.


global imagedata;

if nargin < 3;
    lim = [-Inf, Inf];
end

data = implot(col, ds);
rng = data.range;
data = data.data*1e6;
npix = fliplr(size(data));

if nargin < 4
    trc = 1:npix(2);
end

x = linspace(rng(1, 1), rng(1, 2), npix(1))*1e3;
scan = imagedata.cache{imagedata.collections{col}.images(ds).index}.params.scan;


mask = x > lim(1) & x < lim(2);
x = x(mask);
data = mean(data(trc, mask), 1);




[m, mind] = max(abs(gradient(data)));
med(1) = median(data(1:mind));
med(2) = median(data(mind:end));

fitpar(1) = mean(med);
fitpar(2) = diff(med)/2;
fitpar(3) = x(mind);
fitpar(4) = .3;

figure(40)
clf
plot(x, data, x, fitfn(fitpar, x));
hold on

fitpar = nlinfit(x, data, @fitfn, fitpar);
plot(x, fitfn(fitpar, x), 'r');


if abs(fitpar(3)) < 100
    fprintf('t_c = %.1f mueV\n', fitpar(4) * .2/2.8 * 1e3);
else    
    fprintf('t_c = %.1f mueV\n', fitpar(4) * .2/3.4 * 1e3);
end


function y = fitfn(beta, x)

x = x-beta(3);
y = beta(1) + beta(2) * x./sqrt(x.^2 + beta(4).^2);
if length(beta) >= 5 
    y = y + beta(5) * x;
end
    
