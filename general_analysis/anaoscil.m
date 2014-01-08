function fitpar = anaoscil(col, ds, tmax, startpar)
% fitpar = anaoscil(col, ds, tmax, startpar)
% parameters: offset, amplitude, freq, phase, T2*^2

% (c) 2010 Hendrik Bluhm.  Please see LICENSE and COPYRIGHT information in plssetup.m.


if nargin < 3
    tmax = inf;
end
    
data = implot(col, ds);

fitpar = zeros(length(ds), 5);


fitpar(:, 5) = 400;  
fitpar(:, 4) = 0;

for i = 1:length(ds)
    x = imsetgraphprop(col, ds(i), 'data', 'XData');
    mask = x >= 0 & x <= tmax;
    x = x(mask);
    y = data(i).data(mask)';

    
    %[m, mi] = max(y);
    %fitpar(i, 3) = pi./x(mi);
    fitpar(i, 3) = pi./x(find(diff(y(3:end))< 0, 1)+3);
    fitpar(i, 2) = -std(y) * sqrt(2);
    %fitpar(i, 1) = y(1) - fitpar(i, 2);
    fitpar(i, 1) = mean(y);
    if nargin >= 4
        fitpar(i, isfinite(startpar)) = startpar(isfinite(startpar));
    end
    
    
    x0 = linspace(0, max(x));
    figure(40);
    fp0 = fitpar(i, :);
    fitpar(i, :) = nlinfit(x, y, @fitfn, fitpar(i, :));
    
    plot(x, y, '.', x0, fitfn(fp0, x0),':', x0, fitfn(fitpar(i, :), x0));
    
    
%pause;
    %scan = imagedata.cache{imagedata.collections{col}.images(ds).index}.params.scan;
end
mask = fitpar(:, 2) < 0;
fitpar(mask, 4) =  mod(fitpar(mask, 4) + pi, 2*pi);
fitpar(mask, 2) = -fitpar(mask, 2);


function y = fitfn(beta, x)

y = beta(1) + beta(2) * cos(beta(3) * x + beta(4)) .* exp(-x.^2/beta(5));
