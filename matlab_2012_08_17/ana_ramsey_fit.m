function [ output_args ] = ana_ramsey_fit( out )
%Fits ramsey fringes
%   Input struct is output of ana_rescale_ramsey
offset=round(.05*size(out.newtime(1,:),2));
figure(1); clf;
clear tau;
clear period;
clear const;
clear amp;
for k=1:length(out.detunings)
    %setup the mfit
    data = [];
    model = [];
        
    fitfn = @(p,x) p(1) +(p(2)*sin(2*pi*x/p(3) +p(6))).*exp(-((x)/p(4)).^p(5));
    pars = [.45, .3, 200, 500, 0, 1.6,pi/2];
    lb = [.1 .1 75 50 1.4  0];
    ub=  [.5 .5 500 2000 1.6 2*pi];
        
    mask2 = offset:(size(out.newdata,2)-offset);
    data(end+1).y = out.newdata(k,mask2);
    data(end).x = out.newtime(k,mask2);
    options = optimset('Display','off','MaxIter',1000);
    pars=lsqcurvefit(fitfn,pars,data(1).x,data(1).y,lb,ub,options);
    fit=fitfn(pars,data(1).x);
    tau(k)=abs(pars(4));
    const(k)=pars(5);
    period(k)=pars(3);
    amp(k)=pars(2);
    figure(1); clf; hold on;
    subplot(4,1,1);
    plot(data(1).x,data(1).y, '.-', data(1).x, fit, 'k')
    subplot(4,1,2); plot(tau);
    xlabel('Detuning index');
    ylabel('T2* (ns)');
    subplot(4,1,3); plot(amp);
    xlabel('Detuning index');
    ylabel('Amplitude');
    subplot(4,1,4); plot(const);
    xlabel('Detuning index');
    ylabel('Decay exponent');


end

