function wf=rc_trafofn(wf, channel, tau)
% This function takes a waveform sampled at 1ns/point, and precompensates
% it for an RC rise time.  The argument is the rise time in ns.

% (c) 2010 Hendrik Bluhm.  Please see LICENSE and COPYRIGHT information in plssetup.m.

if tau == 0 || abs(tau) < 1e-3
    return
end
downsamp=4;
% Double the sampling rate of the data
y=wf;
for i=1:downsamp
  y=([ y(1) y(floor(1+.5:0.5:end+.5)) ] + [ y(floor(1:0.5:end)) y(end) ])/2;
end

nyquist=(2^downsamp)*pi*1e9;  %1 Ghz nyquist frequency on the resampled data
k= ([ 0:(length(y)/2-1) length(y)/2 -((length(y)/2-1):-1:1) ])*nyquist/(length(y)/2);

f=fft(y);
if(tau > 0)
  wf=real(ifft(f./rc(k,tau),'symmetric'));
else
  wf=real(ifft(f.*rc(k,tau),'symmetric'));
end
for i=1:downsamp
  wf = (wf(1:2:end)+wf(2:2:end))/2;
end
return

function j=rc(k,tau)
% Skin effect transfer function.  k=omega
   %j=exp(-sqrt((1i*k+alpha).*(1i*k))).*exp(1i*k);  % Second term mostly
   %cancels phase   
   %kmax=pi*1e9;
   tau=abs(tau)*1e-9;
   j=exp(1i*tau*k/1.464)./(1i*k*tau+1);
%   j=1./(1i*k*tau+1);
 return

