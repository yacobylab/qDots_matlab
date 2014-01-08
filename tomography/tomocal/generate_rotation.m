function simdata=generate_rotation(initial, n, theta, q, steps)
% function simdata=generate_rotation(initial, n, theta, steps)
%   generate perfect simulated bloch rotation data around axis n.

gt=zeros(steps, 3);
t=linspace(0,theta,steps);
for i=1:length(t)
  gt(i,:) = rmat(n,t(i))* initial * exp(-t(i)/(2*pi*q));
end
simdata.t=t;   % rotation angles.
simdata.gt=gt; % ground truth
return;
