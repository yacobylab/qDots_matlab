function [yf ff] = skew_fit_function_m(ps, opts)
% function [yf ff] = skew_fit_function_m(ps,b0)
% return model.fn and model.yfn for fitting axes skews in decays
% yfn takes both y and variances, and transforms both.
% parameters 1:9 are the basis, with basis vectors in rows (reshaped)
% parameters 10-12 are the offsets to the basis sensors.
% Parameters ps(1), ps(2), ps(3) are the amplitude and tau and offsetfor decay part of signal.
% 13,14,15 might be smart to use.  b0 is the initial basis.
if ~exist('opts','var')
   opts='';
end
ps
if ~isempty(strfind(opts,'gauss'))
  ff=@(p,x) fff2(ps,p,x);
else
  ff=@(p,x) fff(ps,p,x);
end
 yf = @(p,y,ys) yff(p,y,ys);
end

function [y,ys]=yff(p,y,ys)
  y=reshape(y,length(y)/3,3);
  ys=reshape(ys,length(ys)/3,3);
  pm=reshape(p(1:9),3,3);
  yp=y*pm;
  % propagate erros
  pms = pm.^2;
  ys=ys*pms;  % now std. deviations for variables, neglecting covariances.
  yp=yp + repmat(p([10 11 12]),size(y,1),1)*pm;
  ys=sum((2*yp).^2 .* ys,2)'; % variance of r^2 
  y=sqrt(sum(yp.^2,2))'; 
  ys=(ys./(4*(y.*y)));
end

function y=fff(ps,p,x)
  y=p(ps(1))*exp(-x./(100*p(ps(2)))) + p(ps(3));
end

function y=fff2(ps,p,x)
  y=p(ps(1))*exp(-(x./(100*p(ps(2)))).^2) + p(ps(3));
end
