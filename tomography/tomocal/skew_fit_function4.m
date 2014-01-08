function [yf ff] = skew_fit_function4(ps,b0, opts)
% function [yf ff] = skew_fit_function2(ps,b0)
% return model.fn and model.yfn for fitting axes skews in decays
% parameters 1:6 are the x,y,z coordinates of the y and z axes
% [1 0 0 1] would be perfect axes.  Shorter vectors correspond to decays.
% Parameters ps(1), ps(2), ps(3) are the amplitude and tau and offsetfor decay part of signal.
% 10,11,12 might be smart to use.  b0 is the initial basis.
if ~exist('opts','var')
   opts='';
end
ps
if ~isempty(strfind(opts,'gauss'))
  ff=@(p,x) fff2(ps,p,x);
else
  ff=@(p,x) fff(ps,p,x);
end
yf = @(p,x) yff(p,x,b0);
end

function y=yff(p,y,b0)
  y=reshape(y,length(y)/3,3);
  pm=b0;
  pm(1:3,2:3)=reshape(p(1:6),3,2);  
  yp=(y+repmat(p([7 8 9]),size(y,1),1))*pm;
  y=sqrt(sum(yp.^2,2))'; 
end
function y=fff(ps,p,x)
  y=p(ps(1))*exp(-x./(100*p(ps(2)))) + p(ps(3));
end

function y=fff2(ps,p,x)
  y=p(ps(1))*exp(-(x./(100*p(ps(2)))).^2) + p(ps(3));
end
