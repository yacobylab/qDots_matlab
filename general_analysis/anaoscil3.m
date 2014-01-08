function [fp, labt] = anaoscil3(col, ds, nt, bgslp)
% anaoscil3(col, ds, nt, bgslp)
% fit oscialltions with decay, from dvdata.

% (c) 2010 Hendrik Bluhm.  Please see LICENSE and COPYRIGHT information in plssetup.m.


global dvdata;

load(dvdata.collections{col}.datasets(ds).file{1}, 'data');
labt = data{2}(1, :, 1);
labt = (labt-labt(1)) * 24 * 3600;


t = 0:nt-1;

d = dvplot(col, ds);

y = d.data;
n3 = size(y, 1);

y = y(:, 2:nt+1) - repmat(y(:, 1), 1, nt) - repmat(t * bgslp, n3, 1);


fp(:, 2) = pi./t(sum(cumsum(diff(y(:, 3:nt), [], 2) < 0, 2)==0, 2)+3);
fp(:, 1) = mean(y, 2);
fp(:, 4) = 10;  
fp(:, [3 5 6]) = 0;


fp = fitwrap('plinit plfit woff', t, y, fp, @fitfn);


function y = fitfn(beta, x) 
y = beta(1) * (1 - cos(beta(2) * x + beta(3)) .* exp(-x.^2/beta(4)^2) + beta(5) + beta(6)*x);
