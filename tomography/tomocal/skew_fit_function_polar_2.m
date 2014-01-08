function [yf ff] = skew_fit_function_polar_2(ps, opts)
% function [yf ff] = skew_fit_function_se(ps,b0)
% return model.fn and model.yfn for fitting axes skews in decays
% includes sensor errors in a reasonable way.
% yfn takes both y and variances, and transforms both.
% paramerers 1:9 are theta0, phi0, r0, theta1, phi1, r1, theta2, phi2, r2
% parameters 10,11 are correct sensor values for singlet and triplet (nominally -1 and 1).
% Parameters ps(1), ps(2), ps(3) are the amplitude and tau and offsetfor decay part of signal.
% 12,13,14 might be smart to use.  b0 is the initial basis.
if exist('opts','var') && ~isempty(strfind(opts,'basis'))
    yf=trafop(ps);
    ff=0;
    return;
end
[yf ff] = skew_fit_function_se(ps, opts);
yf=@(p,y,ys) yf(trafop(p),y,ys);
ff=@(p,x) ff(trafop(p),x);
end

function p=trafop(p)
p(1:3) = bvec(p(1:3));
p(4:6) = bvec(p(4:6));
p(7:9) = bvec(p(7:9));
p(1:9) = inv(reshape(p(1:9),3,3));
end

function b=bvec(p)
%b = p(3)*([cos(p(1)); sin(p(1))*cos(p(2)) ; sin(p(1)) * sin(p(2))]);
b = p(3)*([sin(p(1))*cos(p(2)); sin(p(1))*sin(p(2)) ; cos(p(1))]);
end