function freq = guessfreq(x, y)
%function freq = guessfreq(x, y)
% guess the frequency of oscillations of y as a function of x

% ctrl = 1: offset, cos, sin, freq, shift = 0, decay prefac = 1/mean(abs(x));

% (c) 2010 Hendrik Bluhm.  Please see LICENSE and COPYRIGHT information in plssetup.m.

if nargin ==1
   y = x;
   x = 1:length(x);
end

dx = mean(diff(x));
%window=0.5*(1+cos(2*pi*(0:length(x)-1)/(length(x)-1)));
window=1;
ft = fft((y-mean(y)).*window) .* exp(1i * x(1) * (0:length(x)-1) * 2*pi /(dx*length(x)))/length(x);
% hack to get phase shift, factors determined experimentally
ft = ft(1:round(end/2));
[m, mi] = max(abs(ft(2:end)));
xi=max(1,mi-2):min(length(ft)-1,mi+2);
mia=sum(abs(ft(xi+1)).*xi)/sum(abs(ft(xi+1)));
freq = 2* pi* (mia)/((length(x)+1) * dx);
% f=ft(mi+1);
% fp(3) = -3*real(f);
% fp(2) = -3*imag(f);
% fp(1) = mean(y);
% fp(5) = 0;
% fp(6) = .5/mean(abs(x));

