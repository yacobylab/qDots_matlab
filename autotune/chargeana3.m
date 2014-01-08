function [trip, slp] = chargeana2(d, rng)

% (c) 2010 Hendrik Bluhm.  Please see LICENSE and COPYRIGHT information in plssetup.m.


d = diff(d, [], 2);

% Added to reduce outlier density by OD.
sp=1;  % Smoothing amount
if 1
 for l=1:size(d,1)
     d(l,:)=smooth(d(l,:),3);
 end
end
d = d- median(d(:));
d = d .* sign(mean(d(:)));

thresh = 2.5 * std(d(:)); %changed from 3 sigma 09/20/10

ix = 1:size(d, 2);
jy = 1:size(d, 1);
x = linspace(rng(1, 1), rng(1, 2), size(d, 2)+1);
x = 0.5 * (x(2:end)+x(1:end-1));
y = linspace(rng(2, 1), rng(2, 2), size(d, 1));

% Updated 10/24/2010 by OD to allow fine-tuning fit; skew distribution
% test seems to frequently guess wrong on the data sign.
trip2 = [];
while(size(trip2,1) < 2)
  figure(1); 
  clf;
  imagesc(rng(1, :), rng(2, tun:), d);
  set(gca, 'ydir', 'normal')
  axis image;
  box on;
  hold on;

  fprintf('Click on R and L triple point.  n to negate data.\n');
  while(size(trip2,1) < 2)
    [xp, yp, kp] = ginput(1);
    switch(kp)
        case 1
            trip2 = [trip2 ; [xp, yp]];
        case 'n'
            d=-d;
            break;
    end
  end
end

if(trip2(1,1) < trip2(2,1)) % User is a doofus and clicked l-r rather than r-l
  trip2 = circshift(trip2,1);
end

cntr = mean(trip2);

slp = -.585;
[x2, y2] = meshgrid(x, y);

%UR sector
imax1 = find(d .* ((slp * (x2 - cntr(1)) < y2 - cntr(2))) > thresh);
jmax1a = find(y2(imax1) > trip2(1, 2) & x2(imax1) < trip2(1, 1)); % left and above R end of junction
jmax1b = find(y2(imax1) < trip2(1, 2) & x2(imax1) > trip2(1, 1)); % right and below R end of junction

% LL sector
imax2 = find(d .* ((slp * (x2 - cntr(1)) < y2 - cntr(2))) > thresh);
jmax2a = find(y2(imax2) > trip2(2, 2) & x2(imax2) < trip2(2, 1)); % left and above LL end of junction
jmax2b = find(y2(imax2) < trip2(2, 2) & x2(imax2) > trip2(2, 1)); % right and below LL end of junction


fit1a = polyfit(y2(jmax1a), x2(jmax1a), 1);
fit1b = polyfit(y2(jmax1b), x2(jmax1b), 1);
fit2a = polyfit(y2(jmax2a), x2(jmax2a), 1);
fit2b = polyfit(y2(jmax2b), x2(jmax2b), 1);

% find intersections;
df1 = fit1a-fit1b;
df2 = fit2a-fit2b;

%trip = [[fit1a(1) , 1] * -df1(2)/df1(1) + [fit1a(2), 0], ...
%    [fit2a(1) , 1] * -df2(2)/df2(1) + [fit2a(2), 0]];
trip = trip2([1 3 2 4]);

slp = [fit1a(1), fit2b(1), fit1b(1), fit2a(1)];


plot(cntr(1), cntr(2), 'x' , x2(jmax1a), y2(jmax1a), '^k', ...
    x2(jmax1b), y2(jmax1b), 'vk', x2(jmax2a), y2(jmax2a), '<k', x2(jmax2b), y2(jmax2b), '>k', ...
    polyval(fit1a, [y(jy(end)), trip2(1, 2)]), [y(jy(end)), trip2(1, 2)], 'k', ...
    polyval(fit1b, [trip2(1, 2), y(jy(1))]), [trip2(1, 2), y(jy(1))], 'k',...
    polyval(fit2a, [trip2(2, 2), y(jy(end))]), [trip2(2, 2), y(jy(end))], 'k', ...
    polyval(fit2b, [y(jy(1)), trip2(2, 2)]), [y(jy(1)), trip2(2, 2)], 'k', trip([1 3]), trip([2 4]), '.k');

return;
