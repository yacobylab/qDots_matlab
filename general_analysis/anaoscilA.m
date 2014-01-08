function fitpar = anaoscilA(col, ds, tmax, startpar)
% fitpar = anaoscilA(col, ds, tmax, startpar)
% fit damped oscillations to (avearged) curves, using ref pulse as zero
% level. (anaoscil is similar, but fits ref level)
% parameters: amplitude, freq, phase, T2*^2

% (c) 2010 Hendrik Bluhm.  Please see LICENSE and COPYRIGHT information in plssetup.m.


if nargin < 3
    tmax = inf;
end
    
data = implot(col, ds);

fitpar = zeros(length(ds), 5);


fitpar(:, 4) = 400;  
fitpar(:, 3) = 0;

for i = 1:length(ds)
    x = imsetgraphprop(col, ds(i), 'data', 'XData');
    mask = x >= 0 & x <= tmax;
    y = data(i).data(mask)' - mean(data(i).data(x < 0));
    x = x(mask);

    
    %[m, mi] = max(y);
    %fitpar(i, 3) = pi./x(mi);
    fitpar(i, 2) = pi./x(find(diff(y(3:end))< 0, 1)+3);
    fitpar(i, 1) = mean(y(x>0));%std(y) * sqrt(2);
    %fitpar(i, 1) = y(1) - fitpar(i, 2);
    if nargin >= 4
        fitpar(i, isfinite(startpar)) = startpar(isfinite(startpar));
    end
    
    
    x0 = linspace(0, max(x), 1000);
    figure(40);
    fp0 = fitpar(i, :);
    if max(x) > 50        
        fitpar(i, 5) = 0;
        fitpar(i, :) = nlinfit(x, y, @fitfn2, fitpar(i, :));
        plot(x, y, '.', x0, fitfn(fp0, x0),':', x0, fitfn2(fitpar(i, :), x0));
    else
        fitpar(i, :) = nlinfit(x, y, @fitfn, fitpar(i, :));
        plot(x, y, '.', x0, fitfn(fp0, x0),':', x0, fitfn(fitpar(i, :), x0));
    end
    
%pause;
    %scan = imagedata.cache{imagedata.collections{col}.images(ds).index}.params.scan;
end
mask = fitpar(:, 1) < 0;
fitpar(mask, 3) =  mod(fitpar(mask, 3) + pi, 2*pi);
fitpar(mask, 1) = -fitpar(mask, 1);
fitpar = fitpar(:, [1 1:end]);

function y = fitfn(beta, x)

y = beta(1) * (1 - cos(beta(2) * x + beta(3)) .* exp(-x.^2/beta(4)));

function y = fitfn2(beta, x)

y = beta(1) * (1 - cos(beta(2) * x + beta(3)) .* exp(-x.^2/beta(4))) + beta(5)*x;
