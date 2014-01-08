function [data, scalefuncs, mv,fp]=anaHistScale(scan, data,t1s,grps)
%[data,scalefuncs, meanvalues,fp]=anaHistScale(scan, data,t1s,grps)
% Rescale histogrammed data.  t1s is a vector of t1 time estimtes.
% ASSUMES NO CROSSTALK
% currently only works w/ only one channel of data.
 
if length(size(data{end})) == 2 % only 1 group
    data{end}=permute(data{end},[1 3 2]);
end

if ~exist('grps','var') || isempty(grps)
    grps=1:size(data{end},2); %second dim of histogram is # groups. 
end
nds=floor(length(data)/2); %nds = number of data sets. 
for i=1:length(t1s)
   [fitfn, initfn] = getfn(t1s(i)); %fit function for histograms, in terms of peak spacing, noise, t1. 
   ii = [];
   for j = 1:length(scan.loops(1).procfn) %find the procfn that does histogramming
       if length(scan.loops(1).procfn(j).fn)==1 && strcmpi(func2str(scan.loops(1).procfn(j).fn.fn),'histc')
           ii  = [ii, j];
       end
   end
   if length(ii) ~= 1
       error('cannot determine correct procfn');
   end
   %v=scan.loops(1).procfn(end-(nds-i)).fn.args{1};
   v=scan.loops(1).procfn(ii).fn.args{1}; %scan.loops(1).procfn(3).fn=histc, %scan.loops(1).procfn(3).args=set of histogram vals, from fbdata and fConfSeq2, scan.loops(1).procfn(3).dim=180
   v=v+(v(2)-v(1))/2;    % HistC gives edges, not centers.
   v(end)=[];
   data{nds+i+1}(isnan(data{end})) = 0; %any nans in histogrammed set to 0, nds+i+1 is histogram associated w/ data set i. 
   if all(data{nds+1+1}==0) %it was all populated with nans, means histogramming didin't work
      error('histogram data was all 0 or NaN. anaHistScale wont work'); 
   end
   n=squeeze(sum(sum(data{nds+i+1}(:,grps,:),2),1)); %averages over all reps and groups to find the number of elements at each voltage. 
   n(end)=[];
   n=n/mean(n); % Make fitwrap happy.
   figure(400+i);
   clf; hold on;
   fp=fitwrap('plfit plinit samefig fine',v,n',initfn,fitfn,[1 1 1 1 0 0 1]);    
   %fprintf('Peak distance %g\n',1/fp(2));
   % uncomment below to fit t1.
   %fp=fitwrap('plfit samefig resid',v,n',fp,fitfn,[1 1 1 1 1 1 1]);
   % Mean voltages for singlet, triplet
   mv = ((1-exp(-abs(fp(5:6))))./abs(fp(5:6))-.5).*[-1 1]./fp(2) + fp(1); 
   scalefuncs{i}=makescalefunc(1/diff(mv),-mv(1)/diff(mv)); %takes inverse peak spacing and offset w.r.t inverse peak spacing to rescale from 0 to 1.  
   data{i} = scalefuncs{i}(data{i});   
end

end
    
function f=makescalefunc(scale,off)
  f=@(x) x*scale+off;
end

% t1 is the ratio of t1 to the relevant measurement time 
function [fitfn, initfn] = getfn(t1)
if isnan(t1)
    t1=1e-5;
end
distfn = @(a, x) exp(-a(1)) * exp(-(x-1).^2./(2* a(2)^2))./(sqrt(2 * pi) * a(2)) + ...
     a(1)/2 * exp(a(1)/2 * (a(1) * a(2)^2 - 2 * x)) .* ...
     (erf((1 + a(1) * a(2)^2 - x)./(sqrt(2) * a(2))) + erf((-a(1) * a(2)^2 + x)./(sqrt(2) * a(2))));
% parameters: [t_meas/T1, rms amp noise/peak spacing]
 
fitfn = @(a, x) a(3) * distfn(abs(a([5 7])), .5-(x-a(1)).*a(2)) + a(4) * distfn(abs(a([6 7])), .5+(x-a(1)).*a(2));
%parameters: [ center between peaks, 1/spacing, coeff left peak, coeff right peak, t_m/T1,left, t_m/T1,right, rms amp noise/peak spacing]
% sum(coefficients) = 1 corresponds to a PDF for unity peak spacing.
% If fitting raw histograms, # samples = sum(fp(:, 3:4), 2) ./(fp(:, 2) * diff(d.x(1:2)));

initfn.fn = @(x, y)[sum(x.*y)/sum(y),  1/sqrt(sum(y.*((x-sum(x.*y)/sum(y)).^2)/sum(y))), max(y), max(y), 1e-5, t1, .2];
initfn.args = {};
%fifn.vals = [nan(1, 4), -10, 0];

end