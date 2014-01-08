function [rawdata, scalefuncs]=anaRawScale(rawdata,t1s,grps)
%data=anaRawUnpack(scan, data)
% Rescale raw data.  t1s is a vector of t1 time estimtes.
% ASSUMES NO CROSSTALK
% Output data is in data{channel}(group, pulse, rep)

if ~exist('grps','var') || isempty(grps)
    grps=1:size(rawdata{1},1);
end

for i=1:length(t1s)
   [fitfn, initfn] = getfn(t1s(i));
   fd=rawdata{i}(grps,:,:);
   [n,v] = hist(fd(:),200);   % HIST USES CENTERS!
   figure(500+i);
   fp=fitwrap('plinit plfit samefig',v,n,initfn,fitfn,[1 1 1 1 0 0 1]);
   % Mean voltages for singlet, triplet
   mv = ((1-exp(-abs(fp(5:6))))./abs(fp(5:6))-.5).*[-1 1]./fp(2) + fp(1);
   scalefuncs{i}=makescalefunc(1/diff(mv),-mv(1)/diff(mv));
   rawdata{i} = scalefuncs{i}(rawdata{i});   
end

end
    
function f=makescalefunc(scale,off)
  f=@(x) x*scale+off;
end

% t1 is the ratio of t1 to the relevant measurement time 
function [fitfn, initfn] = getfn(t1)

distfn = @(a, x) exp(-a(1)) * exp(-(x-1).^2./(2* a(2)^2))./(sqrt(2 * pi) * a(2)) + ...
     a(1)/2 * exp(a(1)/2 * (a(1) * a(2)^2 - 2 * x)) .* ...
     (erf((1 + a(1) * a(2)^2 - x)./(sqrt(2) * a(2))) + erf((-a(1) * a(2)^2 + x)./(sqrt(2) * a(2))));
% parameters: [t_meas/T1, rms amp noise/peak spacing]
 
fitfn = @(a, x) a(3) * distfn(abs(a([5 7])), .5-(x-a(1)).*a(2)) + a(4) * distfn(abs(a([6 7])), .5+(x-a(1)).*a(2));
%parameters: [ center between peaks, 1/spacing, coeff left peak, coeff right peak, t_m/T1,left, t_m/T1,right, rms amp noise/peak spacing]
% sum(coefficients) = 1 corresponds to a PDF for unity peak spacing.
% If fitting raw histograms, # samples = sum(fp(:, 3:4), 2) ./(fp(:, 2) * diff(d.x(1:2)));

initfn.fn = @(x, y)[sum(x.*y)/sum(y), 3/range(x), max(y), max(y), 1e-4, t1, .2];
initfn.args = {};
%fifn.vals = [nan(1, 4), -10, 0];

end