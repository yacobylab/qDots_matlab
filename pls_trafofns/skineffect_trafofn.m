function wf=skineffect_trafofn(wf, channel, atten, clock)
% This function takes a waveform sampled at 1ns/point, and precompensates
% it for skin effect damping.  Atten specifies the amount of damping in
% dB at 1 Ghz.
% This function

if ~exist('clock','var') || isempty(clock)
    clock=1e9;
end

% (c) 2010 Hendrik Bluhm.  Please see LICENSE and COPYRIGHT information in plssetup.m.

alpha=atoa(atten,2*pi*clock);
downsamp=1;
% Double the sampling rate of the data
y=wf;
for i=1:downsamp
  y=([ y(1) y(floor(1+.5:0.5:end+.5)) ] + [ y(floor(1:0.5:end)) y(end) ])/2;
end

nyquist=(2^downsamp)*pi*clock;  %1 Ghz nyquist frequency on the resampled data
k= ([ 0:(length(y)/2-1) length(y)/2 -((length(y)/2-1):-1:1) ])*nyquist/(length(y)/2);


f=fft(y);

if(atten > 0)
  wf=real(ifft(f./skin(k,alpha,clock),'symmetric'));
else
  wf=real(ifft(f.*skin(k,alpha,clock),'symmetric'));
end
for i=1:downsamp
  wf = (wf(1:2:end)+wf(2:2:end))/2;
end
return

function j=skin(k,alpha,clock)
% Skin effect transfer function.  k=omega
   %j=exp(-sqrt((1i*k+alpha).*(1i*k))).*exp(1i*k);  % Second term mostly cancels phase
  kmax=max(abs(k));
  kmax=2*pi*clock;
  j=exp(-(1+1i)*sqrt(k*alpha)) .* exp(1i*(k/kmax)*sqrt(kmax*alpha));
return

function alpha=atoa(atten,w)
  alpha = ((log(10)*atten/(20))^2)/w;    
return
