function [data, scalefuncs, mv,fp]=anaHistScaleFGrpwise(scan, data,t1s,grpsw)
%[data,scalefuncs, meanvalues,fp]=anaHistScaleGrpwise(scan, data,t1s,grpsw)
% Rescale histogrammed data.  t1s is a vector of t1 time estimtes.
% scale data frame-by-frame and group-by-group.
% ASSUMES NO CROSSTALK
% currently only works w/ only one channel of data.

if length(size(data{end})) == 2 % only 1 group
    data{end}=permute(data{end},[1 3 2]);
end

if ~exist('grpsw','var') || isempty(grpsw)
    grpsw=1:size(data{end},2);
end
nds=floor(length(data)/2);
for i=1:length(t1s)
    for j=1:size(data,1)
        for grps=grpsw
            [fitfn, initfn] = getfn(t1s(i));
            v=scan.loops(1).procfn(end-(nds-i)).fn.args{1};
            % HistC gives edges, not centers.
            v=v+(v(2)-v(1))/2;
            v(end)=[];
            data{nds+i+1}(isnan(data{end})) = 0;
            n=squeeze(sum(sum(data{nds+i+1}(j,grps,:),2),1));
            n(end)=[];
            figure(400+i);
            clf; hold on;
            fp=fitwrap('plfit plinit samefig fine',v,n',initfn,fitfn,[1 1 1 1 0 0 1]);
            fprintf('Peak distance %g\n',1/fp(2));
            % uncomment below to fit t1.
            %fp=fitwrap('plfit samefig resid',v,n',fp,fitfn,[1 1 1 1 1 1 1]);
            % Mean voltages for singlet, triplet
            mv = ((1-exp(-abs(fp(5:6))))./abs(fp(5:6))-.5).*[-1 1]./fp(2) + fp(1);
            scalefuncs{i}=makescalefunc(1/diff(mv),-mv(1)/diff(mv));
            data{i}(j,grps,:) = scalefuncs{i}(data{i}(j,grps,:));
        end
    end
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

initfn.fn = @(x, y)[sum(x.*y)/sum(y),  1/sqrt(sum(y.*((x-sum(x.*y)/sum(y)).^2)/sum(y))), max(y), max(y), 1e-5, t1, .2];
initfn.args = {};
%fifn.vals = [nan(1, 4), -10, 0];

end