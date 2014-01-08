function [trip, slp] = chargeana(d, rng)

% (c) 2010 Hendrik Bluhm.  Please see LICENSE and COPYRIGHT information in plssetup.m.


d = diff(d, [], 2);
d = d- median(d(:));

thresh = 1.5 * std(d(:));

ix = 1:size(d, 2);
jy = 1:size(d, 1);
x = linspace(rng(1, 1), rng(1, 2), size(d, 2)+1);
x = 0.5 * (x(2:end)+x(1:end-1));
y = linspace(rng(2, 1), rng(2, 2), size(d, 1));

% find junction as (negative) minima of derivative
[m, imin] = min(filter2(.5*ones(2), d), [], 2);
mm = min(m);
%mind2 = find(m < .3 * mm);
jmin = find(m < -thresh); % consider only sufficiently steep regions
cntr = round([mean(imin), mean(jmin)]); % rough estimate of center of junction

[ix2, jy2] = meshgrid(ix, jy);

% edges in UR sector
[m, imax1] = max(d .* (ix2 - cntr(1) > -jy2 + cntr(2)) , [], 2);
jmax1a = find(m' >  thresh & jy > jmin(end) & imax1' < imin(jmin(end))); % left and above UR end of junction
jmax1b = find(m' >  thresh & jy < jmin(end) & imax1' > imin(jmin(end))); % right and below UR end of junction
%jmax1a = find(m' >  thresh & jy > jmin(end));
%jmax1b = find(m' >  thresh & imax1' > imin(jmin(end)));

% LL sector
[m, imax2] = max(d .* (ix2 - cntr(1) < -jy2 + cntr(2)) , [], 2);
jmax2a = find(m' >  thresh & jy > jmin(1) & imax2' < imin(jmin(1))); % left and above LL end of junction
jmax2b = find(m' >  thresh & jy < jmin(1) & imax2' > imin(jmin(1))); % right and below LL end of junction

figure(1);
clf;
imagesc(rng(1, :), rng(2, :), d);
set(gca, 'ydir', 'normal')
axis image;
box on;
hold on;

% fit transitions
if all([length(jmax1a) length(jmax1b) length(jmax2a) length(jmax2b)] > 3)
    fit1a = polyfit(y(jmax1a), x(imax1(jmax1a)), 1);
    fit1b = polyfit(y(jmax1b), x(imax1(jmax1b)), 1);
    fit2a = polyfit(y(jmax2a), x(imax2(jmax2a)), 1);
    fit2b = polyfit(y(jmax2b), x(imax2(jmax2b)), 1);


    % find intersections;
    df1 = fit1a-fit1b;
    df2 = fit2a-fit2b;

    trip = [[fit1a(1) , 1] * -df1(2)/df1(1) + [fit1a(2), 0], ...
        [fit2a(1) , 1] * -df2(2)/df2(1) + [fit2a(2), 0]];

    slp = [fit1a(1), fit2b(1), fit1b(1), fit2a(1)];
    

    h = plot(x(imin(jmin)), y(jmin), '.r', x(cntr(1)), y(cntr(2)), 'x' , x(imax1(jmax1a)), y(jmax1a), '^k', ...
        x(imax1(jmax1b)), y(jmax1b), 'vk', x(imax2(jmax2a)), y(jmax2a), '>k', x(imax2(jmax2b)), y(jmax2b), '<k', ...
        polyval(fit1a, y([jy(end), jmin(end)])), y([jy(end), jmin(end)]), 'k', ...
        polyval(fit1b, y([jmin(end) jy(1)])), y([jmin(end), jy(1)]), 'k',...
        polyval(fit2a, y([jmin(1) jy(end)])), y([jmin(1), jy(end)]), 'k', ...
        polyval(fit2b, y([jy(1), jmin(1)])), y([jy(1), jmin(1)]), 'k', trip([1 3]), trip([2 4]), '.k');


    str = input('Accept/manual? (y/n/m)', 's');

    while 1
        switch str
            case {'y', 'Y'}
                return
            case {'n' 'N'}
                trip = [];
                slp = [];
                return;

            case {'m' 'M'}
                delete(h);
                break
        end
    end
end

fprintf('Click on UR and LL triple point.\n');
trip2 = ginput(2);
cntr = mean(trip2);

[x2, y2] = meshgrid(x, y);

%UR sector
[m, imax1] = max(d .* (x2 - cntr(1) > -y2 + cntr(2)) , [], 2);
jmax1a = find(m' >  thresh & y > trip2(1, 2) & x(imax1) < trip2(1, 1)); % left and above UR end of junction
jmax1b = find(m' >  thresh & y < trip2(1, 2) & x(imax1) > trip2(1, 1) & imax1' < ix(end-4)); % right and below UR end of junction

% LL sector
[m, imax2] = max(d .* (x2 - cntr(1) < -y2 + cntr(2)) , [], 2);
jmax2a = find(m' >  thresh & y > trip2(2, 2) & x(imax2) < trip2(2, 1)); % left and above LL end of junction
jmax2b = find(m' >  thresh & y < trip2(2, 2) & x(imax2) > trip2(2, 1)); % right and below LL end of junction


fit1a = polyfit(y(jmax1a), x(imax1(jmax1a)), 1);
fit1b = polyfit(y(jmax1b), x(imax1(jmax1b)), 1);
fit2a = polyfit(y(jmax2a), x(imax2(jmax2a)), 1);
fit2b = polyfit(y(jmax2b), x(imax2(jmax2b)), 1);


% find intersections;
df1 = fit1a-fit1b;
df2 = fit2a-fit2b;

trip = [[fit1a(1) , 1] * -df1(2)/df1(1) + [fit1a(2), 0], ...
    [fit2a(1) , 1] * -df2(2)/df2(1) + [fit2a(2), 0]];

slp = [fit1a(1), fit2b(1), fit1b(1), fit2a(1)];


plot(x(imin(jmin)), y(jmin), '.r', cntr(1), cntr(2), 'x' , x(imax1(jmax1a)), y(jmax1a), '^k', ...
    x(imax1(jmax1b)), y(jmax1b), 'vk', x(imax2(jmax2a)), y(jmax2a), '>k', x(imax2(jmax2b)), y(jmax2b), '<k', ...
    polyval(fit1a, [y(jy(end)), trip2(1, 2)]), [y(jy(end)), trip2(1, 2)], 'k', ...
    polyval(fit1b, [trip2(1, 2), y(jy(1))]), [trip2(1, 2), y(jy(1))], 'k',...
    polyval(fit2a, [trip2(2, 2), y(jy(end))]), [trip2(2, 2), y(jy(end))], 'k', ...
    polyval(fit2b, [y(jy(1)), trip2(2, 2)]), [y(jy(1)), trip2(2, 2)], 'k', trip([1 3]), trip([2 4]), '.k');
