function plotpulse(pulsetab, scale, offset);
% plotpulse(pulsetab, scale, offset);
% plot pulse on top of a dataset.
% pulsetab can be a pulse table, pulse info struct, or pulse index

% (c) 2010 Hendrik Bluhm.  Please see LICENSE and COPYRIGHT information in plssetup.m.


if isstruct(pulsetab)
    pulsetab = pulsetab.pulsetab;
elseif length(pulsetab) == 1
    pulsetab = awggetpulseinf(pulsetab);
    pulsetab = pulsetab.pulsetab;
end

if length(scale) == 1;
    scale = [scale, scale];
end

hold on;
plot(pulsetab(2, :) * scale(1) + offset(1), pulsetab(3, :) * scale(2) + offset(2), '.-k');
hold off;
