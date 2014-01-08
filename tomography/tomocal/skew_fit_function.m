function [yf ff] = skew_fit_function(ps)
% function [yf ff] = skew_fit_function(ps)
% return model.fn and model.yfn for fitting axes skews in decays
% parameters 1:4 are the y,z coordinates of the y and z axes
% [1 0 0 1] would be perfect axes.  Shorter vectors correspond to decays.
% Parameters ps(1), ps(2) are the amplitude and tau for decay part of signal.
% 5,6 might be smart to use.
ff=@(p,x) fff(ps,p,x);
yf = @yff;
end

function y=yff(p,y)
  y=reshape(y,length(y)/3,3);
  pm(2:3,2:3)=reshape(p(1:4),2,2);
  pm(1,1)=1;
  yp=y*pm;
  y=sqrt(sum(yp.^2,2))'; 
end
function y=fff(ps,p,x)
  y=p(ps(1))*exp(-x./p(ps(2)));
end
