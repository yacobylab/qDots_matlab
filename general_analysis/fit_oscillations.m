function [fn fp]=fit_oscillations(x,y,config)
%function [fitfn fp]=fit_oscillations(file,config)
% fits 1-D oscillations
% returns the fit function and fit parameters
% config has fields:
%    opts: can be gauss, exp, nodecay, nocenter, noplot; defaults to gauss
% fp = [offset, amplitude, freq, phase, decay_const, decay_center]

% if ~exist('file','var') || isempty(file)
%   file=uigetfile('sm*.mat');
% end

if ~exist('config','var')
    config=struct();
elseif ischar(config)
    config=struct('opts',config);
elseif iscell(config)
    config=struct(config{:});
end
config = def(config,'opts','guass ');
%config = def(config,'blah',[1 Inf]);

% prepare data to be right shape for fitwrap
x = (x(:))';
y = y(:)';

   
   
if isopt(config,'exp')
   fitfn = @(p,x) p(1)+p(2)*cos(p(3)*x+p(4)).*exp(-(x-p(6))./p(5));  
end

if isopt(config,'gauss')
    fitfn = @(p,x) p(1)+p(2)*cos(p(3)*x+p(4)).*exp(-((x-p(6))./p(5)).^2);  
end

if isopt(config,'nodecay')
   fitfn = @(p,x) p(1)+p(2)*cos(p(3)*x+p(4))*p(5)*p(6);
end

p=fioscill(x,y,1);  
init=[mean(y), range(y)/2, p(4),pi,range(x)/10, x(end/2)];
mask = [1 1 1 1 1 1];

if isopt(config,'nocenter') % for ramsey, to force decay to start at t=0
   init(6) =0;
   mask(6) = 0;
end

if isopt(config,'nodecay')
   init(5:6) = [1 1]; 
   mask(5:6) = [0 0];
end

init=fitwrap('',x,y,init,fitfn,[1 0 0 1 0 0]);
p=fitwrap('fine',x,y,init,fitfn,mask);  

if ~isopt(config,'noplot')
   figure(502); clf; hold on;
   plot(x,y, 'rx');
   plot(x,fitfn(p,x),'b');
end

fn = fitfn;
fp = p;


%lots of nice code goes here
 
end %function ends here

% Apply a default.
function s=def(s,f,v)
  if(~isfield(s,f))
      s=setfield(s,f,v);
  end
return;
end

function b=isopt(config,name)
  b=~isempty(strfind(config.opts,name));
return;
end