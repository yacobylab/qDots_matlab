function fitpar = anaoscil2(col, ds, step, startpar, ctrl)
% fitpar = anaoscil(col, ds, step, startpar, ctrl)
% parameters: offset, in-phase, out-phase, freq
% ctrl: usefirst, plot, useref

% (c) 2010 Hendrik Bluhm.  Please see LICENSE and COPYRIGHT information in plssetup.m.


if nargin < 5
    ctrl = '';
end

if strfind(ctrl, 'usefirst')
    firstp = 1;
else
    firstp = 2;
end

data = implot(col, ds);

if strfind(ctrl, 'useref')
    data.data = data.data(:, 2:end) - repmat(data.data(:, 1), 1, size(data.data, 2)-1);
    firstp = 1;
end

if length(step) == 1
    nfit = floor(size(data.data, 1)/step);
else
    nfit = length(step);
end
fitpar = zeros(nfit, 4);
x = 0:size(data.data, 2)-firstp;

if strfind(ctrl, 'plot')
    figure(40);
    pause
    clf;
end

for i = 1:nfit
    if length(step) == 1
        y = mean(data.data((i-1)*step+(1:step), firstp:end), 1);
    else
        y = data.data(step(i), firstp:end);
    end
    %[m, mi] = max(y);
    %fitpar(i, 3) = pi./x(mi);
    if i == 1 || 1 
       
        ft = ifft(y);
        ft = ft(1:round(end/2));
        [m, mi] = max(abs(ft(2:end)));
        
        
        fitpar(i, 4) = 2* pi* mi/length(x);
        fitpar(i, 2) = 2*real(ft(mi+1));
        fitpar(i, 3) = 2*imag(ft(mi+1));
        %fitpar(i, 1) = y(1) - fitpar(i, 2);
        fitpar(i, 1) = ft(1);

        if nargin >= 4 && ~isempty(startpar)
            fitpar(i, isfinite(startpar)) = startpar(isfinite(startpar));
        end
    else
        fitpar(i, :) = fitpar(i-1, :);
    end
    
    x0 = linspace(0, max(x));
    
    fp0 = fitpar(i, :);
    if strfind(ctrl, 'useref')
        fp0(1) = sqrt(sum(fitpar(i, 2:3).^2));
        fitpar(i, 2:4) = nlinfit(x, y, @fitfn2, fitpar(i, 2:4));
        fitpar(i, 1) = sqrt(sum(fitpar(i, 2:3).^2));
    else        
        fitpar(i, :) = nlinfit(x, y, @fitfn, fitpar(i, :));
    end
    
    if strfind(ctrl, 'plot')
        plot(x, y, '-.', x0, fitfn(fp0, x0),':', x0, fitfn(fitpar(i, :), x0));
        title(sprintf('%d: freq = %f', i, fitpar(i, 4)));
        drawnow;
        if strfind(ctrl, 'pause');            
            pause%(.05);
        end
    end
end

%fitpar(:, 4) = mod(fitpar(4, pi));


function y = fitfn(beta, x)

y = beta(1) + beta(2) * cos(beta(4) * x) + beta(3) * sin(beta(4) * x) ;


function y = fitfn2(beta, x)
y = sqrt(sum(beta(1:2).^2)) + beta(1) * cos(beta(3) * x) + beta(2) * sin(beta(3) * x) ;
