function blochplot(U, psi_0, varargin)

% (c) 2010 Hendrik Bluhm.  Please see LICENSE and COPYRIGHT information in plssetup.m.


if length(psi_0) == 2;
    psi_0(3:4) = 0;
end

psi = vmmult2(psi_0, U);

x = 2 * real(psi(:, 1) .* conj(psi(:, 2)));
y = 2 * imag(psi(:, 1) .* conj(psi(:, 2)));
z = abs(psi(:, 1)).^2 - abs(psi(:, 2)).^2;

figure(10);
%clf;
hold on;
h = plot3(x, y, z, 'r');
set(h(1), 'linewidth', 1, varargin{:});

[phi, theta] = meshgrid(linspace(0, 2*pi), linspace(-pi/2, pi/2, 7));
plot3(cos(phi)'.* cos(theta)', sin(phi)'.* cos(theta)', sin(theta)', 'color', .8 * ones(1, 3));
[phi, theta] = meshgrid(linspace(0, 2*pi, 13), linspace(-pi/2, pi/2, 50));
plot3(cos(phi).* cos(theta), sin(phi).* cos(theta), sin(theta), 'color', .8 * ones(1, 3));
set(gcf, 'color', 255/256*ones(3, 1))

campos = [14.4075   -7.3699    7.6063; 0, -20, 0; 0 0 -20];
set(gca, 'DataAspectRatio', [1 1 1], 'CameraPosition', campos(1, :), ...
    'CameraViewAngle', 6.5, 'cameraTarget', [0 0 0]);

xlabel('X');
ylabel('Y');
zlabel('Z');
axis off;
