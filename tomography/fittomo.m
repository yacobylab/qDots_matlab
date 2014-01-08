function [beta1, f, y2] = fittomo(data, xscale, beta0, mask)

% (c) 2010 Hendrik Bluhm.  Please see LICENSE and COPYRIGHT information in plssetup.m.


if nargin < 4
    mask = true(size(beta0));
end

tmin = 7;
tmax = size(data, 2)-1;
data(1, :) = (data(1, :)-1)*xscale+1;
x = [tmin:tmax, -(tmin:tmax)];
y = [data(2, tmin+1:end), data(1, tmin+1:end)];

x0 = linspace(tmin, tmax);
y0 = fitfn(beta0(mask), [x0 -x0]);

beta1 = beta0;
beta1(mask) = nlinfit(x, y, @fitfn, beta0(mask));

y1 = fitfn(beta1(mask), [x0 -x0]);

% use Z parameters for X.
beta2 = beta1;
%beta2(3:2:7) = beta2(2:2:6); % for parametrization with two angles
beta2([3, 7]) = beta2([2, 6]); % for parametrization with mean and difference in theta
beta2(5) = -beta2(5);
y2 = fitfn(beta2(mask), [x0 -x0]);

figure(7)
clf
subplot(311);

w1 = sqrt(beta1(2)^2 + beta1(1)^2);
w2 = sqrt(beta1(3)^2 + beta1(1)^2);

plot(0:tmax, data(1, :), '.',  x0, y1(end/2+1:end), 'r')%, x0, y0(end/2+1:end), ':')%, ...
%    (x0 * w1 + beta1(6) - beta1(7))/w2, y2(end/2+1:end), '--')
title('X')

subplot(312);
plot(0:tmax, data(2, :), '.', x0, y1(1:end/2), 'r')%, x0, y0(1:end/2), ':')
    %x0, y2(1:end/2), '--', x0, y2(end/2+1:end), '-.')
%, (x0 * w1 + beta1(7) - beta1(6))/w2, y2(end/2+1:end), '--')
title('Z')

if size(data, 1) >= 3
    a1 = sign(diff(data(3, 1:2))) * sin(beta1(4) + .5 * beta1(5));
    a2 = sign(diff(data(3, 1:2))) * sin(beta1(4) - .5 * beta1(5));
    yy1 = a1 * sin(x0*w1 + beta1(6)) .* exp(-x0.^2./beta1(8)^2);
    %yy2 = a2 * sin(x0*sqrt(beta1(1)^2 + beta1(3)^2) + beta1(7)) .* exp(-x0.^2./beta1(8)^2);

    f = nlinfit(tmin:tmax, data(3, tmin+1:end), @fitfn2, [w1 beta1(6), a1]);
    f(3) = a1;
    yy2 = fitfn2(f, x0);
    subplot(313)
    plot(0:tmax, data(3, :), '.', x0, yy2, 'r')%, x0, yy1);
    title('Y')
  
end
if 0
    figure(8)
    clf
    hold on;
    
    x1 = linspace(0, 30);
    x2 = 0:size(data, 2)-1;
    h = plot3(interp1((x2*w2 + beta1(7)- beta1(6))./w1, data(1, :), x1, 'spline', nan), ...
        interp1((x2*f(1) - beta1(6)+ f(2))./w1, data(3, :), x1, 'spline', nan), ...
        interp1(x2, data(2, :), x1, 'spline', nan), 'g');
  
    set(h(1), 'linewidth', 1.5);
        
    %plot(x1, interp1((x2*w2 + beta1(7)- beta1(6))./w1, data(1, :), x1, 'spline', nan), ...
    %x2, data(2, :), 'r.', x1, interp1(x2, data(2, :), x1, 'spline', nan), 'r', ...
    %    x1, interp1((x2*f(1) - beta1(6)+ f(2))./w1, data(3, :), x1, 'spline', nan), ...
    %    x1, interp1(x2, data(3, :), x1, 'spline', nan));
    
    % this is Z freq and amplitudes
    h = plot3(y2(end/2+1:end), a1 * sin(x0*w1 + beta1(6)) .* exp(-x0.^2./beta1(8)^2), y2(1:end/2), 'b', ...
        [0 beta1(1)/w1], [0 0], [0, beta1(2)/w1], 'r.-');
    
    set(h(1), 'linewidth', 2);
    %data(1, :), d3(2, :), d3(ypls, :),..
    %    [0 mean(d3(1, :))], [0 mean(d3(2, :))], [0 mean(d3(ypls, :))],
    %    'r');

    [phi, theta] = meshgrid(linspace(0, 2*pi), linspace(0, pi/2, 4));
    plot3(cos(phi)'.* cos(theta)', sin(phi)'.* cos(theta)', sin(theta)', 'color', .8 * ones(1, 3));
    [phi, theta] = meshgrid(linspace(0, 2*pi, 13), linspace(0, pi/2, 25));
    plot3(cos(phi).* cos(theta), sin(phi).* cos(theta), sin(theta), 'color', .8 * ones(1, 3));
    set(gcf, 'color', .99*ones(3, 1))
    
    xlabel('X');
    ylabel('Y');
    zlabel('Z');
    axis off;
    text(0, 0, 1.2, '|S\rangle');
    text(1, -.2, -.15, '|\uparrow\downarrow\rangle');
    text(-1.1, 0, -.15, '|\downarrow\uparrow\rangle');
    h = plot3([-1 0 1], [0 0 0], [0 1 0], 'k.');
    set(h, 'MarkerSize', 20);
    set(gca, 'DataAspectRatio', [1 1 1], 'CameraPosition', [14.4075   -7.3699    7.6063], ...
        'CameraViewAngle', 5.5);
    
end

    function y = fitfn(beta, x)
        % x>0: Z-data, x<0: X-data
        % beta = [dBz, Jz, Jx, thetaz, thetax, phiz, phix, t2*];

        beta([find(mask), find(~mask)]) = [beta, beta0(~mask)];
        beta(4:5) = beta(4:5) * [1 1; .5, -.5];

        wz = sqrt(beta(2)^2 + beta(1)^2);
        wx = sqrt(beta(3)^2 + beta(1)^2);

        a = beta(1)/wz * sin(beta(4));% * exp(-x(x>0).^2./beta(8)^2));
        b = beta(3)/wx * sin(beta(5));% * exp(-x(x<0).^2./beta(8)^2));

        %y(x>0) = beta(2)/wz * cos(beta(4) * exp(-x(x>0).^2./beta(8)^2)) - a .* cos(wz*x(x>0)+beta(6));% .* exp(-x(x>0).^2./beta(8)^2);
        %y(x<0) = beta(1)/wx * cos(beta(5) * exp(-x(x<0).^2./beta(8)^2)) + b .* cos(wx*x(x<0)+beta(7));% .* exp(-x(x<0).^2./beta(8)^2);
        y(x>0) = beta(2)/wz * cos(beta(4)) - a .* cos(wz*x(x>0)+beta(6)) .* exp(-x(x>0).^2./beta(8)^2);
        y(x<0) = beta(1)/wx * cos(beta(5)) + b .* cos(-wx*x(x<0)+beta(7)) .* exp(-x(x<0).^2./beta(8)^2);

    end

    function y = fitfn2(beta, x)
        y = beta(3) * sin(beta(1) * x + beta(2)) .* exp(-x.^2./beta1(8)^2);
    end
end
