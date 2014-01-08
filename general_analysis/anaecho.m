function fp = anaecho(col, ds, lim, fitmask, freq)
% fp = anaecho(col, ds, lim, fitmask)

% (c) 2010 Hendrik Bluhm.  Please see LICENSE and COPYRIGHT information in plssetup.m.


global imagedata;

if nargin < 3 || isempty(lim)
    lim = [-inf, inf];
end

if nargin < 4 || isempty(fitmask);
    fitmask = logical([0 1 1 1]);
end

if nargin < 5 
    freq = .7;
end

if length(lim) == 1
    lim = [-lim, lim];
end

x = imagedata.collections{col}.images(ds(1)).graphobjs{1}.values{1};
mask = x > lim(1) & x < lim(2);
d = implot(col, ds);
data = [d.data]'*1e6;


fitfn = @(beta, x)beta(2)*cos(beta(1) * x) + beta(3)*sin(beta(1) * x)+beta(4);
fp = fitwrap('plfit', x(mask), data(:, mask), [freq, -1 0 0], fitfn, fitmask);

fp(:, 2:3) = [sqrt(sum(fp(:, 2:3).^2, 2)), atan2(fp(:, 2), fp(:, 3))];

if nargout == 0
    figure(501);
    subplot(211)
    plot(1:length(ds), fp(:, 2));
    subplot(212)
    plot(1:length(ds), fp(:, 3));
end
