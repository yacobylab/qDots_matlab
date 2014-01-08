function f=mlogspace(low,high,npts, mult)
% function f=mlogspace(low,high,npts, mult)
%  Return npts log spaced points betwixt low and high.
%  if mult is specified, round the points to be a multiple of mult.
% ie
% mlogspace(1,10,4) =   1.0000    2.1544    4.6416   10.0000
% mlogspace(1,10,4,1) = 1 2 5 10
if exist('mult','var')
  f=mult*round(logspace(log10(low),log10(high),npts)/mult);
else
  f=logspace(log10(low),log10(high),npts);
end

