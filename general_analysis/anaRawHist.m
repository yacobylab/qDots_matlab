function [spcnt tpcnt] = anaRawHist(col, ds, t1)

% (c) 2010 Hendrik Bluhm.  Please see LICENSE and COPYRIGHT information in plssetup.m.

for i=1:length(ds)  
  pd=dvplot(col,ds(i));
  rol = pd.y(3)
  
  dvdisplay(col, ds(i), 'curve');
  [fitfn, initfn] = getfn(rol/t1);
  dvfit(col, ds(i),  'plinit plfit', initfn, fitfn, [1 1 1 1 0 0 1]);
end
  fp = dvfit(col, ds, 'getp');  
  spcnt=fp(:,3)./(fp(:,3)+fp(:,4));
  tpcnt=fp(:,4)./(fp(:,3)+fp(:,4));  
for i=1:length(ds)
    dvdisplay(col,ds(i),'title','string',sprintf('%.1f %% S, %.1f %% T', 100*spcnt(i), 100*tpcnt(i)));    
end
  dvplot(col,ds);
end

function [fitfn, initfn] = getfn(t1)

distfn = @(a, x) exp(-a(1)) * exp(-(x-1).^2./(2* a(2)^2))./(sqrt(2 * pi) * a(2)) + ...
     a(1)/2 * exp(a(1)/2 * (a(1) * a(2)^2 - 2 * x)) .* ...
     (erf((1 + a(1) * a(2)^2 - x)./(sqrt(2) * a(2))) + erf((-a(1) * a(2)^2 + x)./(sqrt(2) * a(2))));
% parameters: [t_meas/T1, rms amp noise/peak spacing]
 
fitfn = @(a, x) a(3) * distfn(abs(a([5 7])), .5-(x-a(1)).*a(2)) + a(4) * distfn(abs(a([6 7])), .5+(x-a(1)).*a(2));
%parameters: [ center between peaks, 1/spacing, coeff left peak, coeff right peak, t_m/T1,left, t_m/T1,right, rms amp noise/peak spacing]
% sum(coefficients) = 1 corresponds to a PDF for unity peak spacing.
% If fitting raw histograms, # samples = sum(fp(:, 3:4), 2) ./(fp(:, 2) * diff(d.x(1:2)));

initfn.fn = @(x, y)[sum(x.*y)/sum(y), 125, max(y), max(y), 1e-3, t1, .2];
initfn.args = {};
%fifn.vals = [nan(1, 4), -10, 0];

end
