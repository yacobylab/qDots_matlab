function [yf ff] = skew_fit_function_se(ps, opts)
% function [yf ff] = skew_fit_function_se(ps,b0)
% return model.fn and model.yfn for fitting axes skews in decays
% includes sensor errors in a reasonable way.
% yfn takes both y and variances, and transforms both.
% parameters 1:9 are the basis, with basis vectors in rows (reshaped)
% parameters 10,11 are correct sensor values for singlet and triplet (nominally -1 and 1).
% Parameters ps(1), ps(2), ps(3) are the amplitude and tau and offsetfor decay part of signal.
% 12,13,14 might be smart to use.  b0 is the initial basis.
if ~exist('opts','var')
   opts='';
end

if ~isempty(strfind(opts,'gauss'))
  ff=@(p,x) fff_gauss(ps,p,x);
elseif ~isempty(strfind(opts,'quad'))
  ff=@(p,x) fff_quad(ps,p,x);
else    
  ff=@(p,x) fff(ps,p,x);
end
 yf = @(p,y,ys) yff(p,y,ys);
end

function [y,ys]=yff(p,y,ys)
  % get input vector into a reasonable form
  y=reshape(y,length(y)/3,3);
  ys=reshape(ys,length(ys)/3,3);
  pm=reshape(p(1:9),3,3);
  
  % First apply sensor correction
  s=p(10);
  t=p(11);
  a=-s/2 + t/2;
  b=s/2 + t/2;
  yp = y * a + b;
  
  % Propagate errors for sensor correction
  ys = ys * a * a;
  
  % Basis change and magnitude
  yp=yp*pm;
  y=sqrt(sum(yp.^2,2))'; 
  
  % propagate errors for basis change
  if 0 % neglect covariances
    pms = pm.^2;
    ys=ys*pms;  % now std. deviations for variables, neglecting covariances.
    ys=sum((2*yp).^2 .* ys,2)'; % variance of r^2 
    ys=(ys./(4*(y.*y))); % variance of sqrt(r^2);
  else  % include covariances
    ys = (sum( ((pm*((yp)')).^2)' .* ys ,2))'./(y.^2);
  end
  end

function y=fff(ps,p,x)
  y=sqrt((p(ps(1))*exp(-x./(100*p(ps(2))))).^2 + p(ps(3)).^2 );
end

function y=fff_gauss(ps,p,x)
  y=sqrt((p(ps(1))*exp(-(x./(100*p(ps(2)))).^2)).^2 + p(ps(3)).^2 );
end

function y=fff_quad(ps,p,x)
  y=p(ps(1))-p(ps(2))*(x/100).^2;
end
