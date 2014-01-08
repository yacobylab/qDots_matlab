file=get_files;
d=load(file{1}); 
allData=double(d.scan.data.FPGA.data); 
allData=allData(:);

bins=100;
cen = mean(allData);
rng = 6*std(allData);
Vt=linspace(cen-rng,cen+rng,bins); 
hist_data = histc(allData, Vt);

scantime=getscantime(d.scan,d.data); 
%t1p=att1('right',scantime,'after'); 
t1p=1e-5; 
measT=1.2e-6; 
t1=measT/t1p; 
st1=0; 
figure(11); 
plot(Vt,hist_data); 
distfn = @(a, x) exp(-a(1)) * exp(-(x-1).^2./(2* a(2)^2))./(sqrt(2 * pi) * a(2)) + ...
    a(1)/2 * exp(a(1)/2 * (a(1) * a(2)^2 - 2 * x)) .* ...
    (erf((1 + a(1) * a(2)^2 - x)./(sqrt(2) * a(2))) + erf((-a(1) * a(2)^2 + x)./(sqrt(2) * a(2))));
% parameters: [a(1): t_meas/T1, a(2): rms amp noise/peak spacing x: voltage]

fitfn = @(a, x) a(3) * distfn(abs(a([5 7])), .5-(x-a(1)).*a(2)) + a(4) * distfn(abs(a([6 7])), .5+(x-a(1)).*a(2));
%parameters: [a(1): center between peaks, a(2): 1/spacing, a(3): coeff left peak, a(4): coeff right peak, a(5): t_m/T1,left
% a(6): t_m/T1,right, a(7): rms amp noise/peak spacing]
% sum(coefficients) = 1 corresponds to a PDF for unity peak spacing.
% If fitting raw histograms, # samples = sum(fp(:, 3:4), 2) ./(fp(:, 2) * diff(d.x(1:2)));

initfn.fn = @(x, y)[sum(x.*y)/sum(y), 5/range(x), max(y), max(y), st1, t1, .2];
initfn.args = {};
%fifn.vals = [nan(1, 4), -10, 0];


fitpar2=fitwrap('plinit plfit',Vt,hist_data',initfn,fitfn,[1 1 1 1 0 0 1]);
samp_num=length(allData);


    Sfit=fitpar2; Sfit(3)=1; Sfit(4)=0; % here, we set the triplet peak coeff. to 0, singlet coeff. to 1
    fitSing=(fitpar2(3)+fitpar2(4))*fitfn(Sfit,Vt)/samp_num; %Then recalculate the curves
    Tfit=fitpar2; Tfit(3)=0; Tfit(4)=1;
    fitTrip=(fitpar2(3)+fitpar2(4))*fitfn(Tfit,Vt)/samp_num;
    
   SfidFit=cumsum(fitSing); %Sum number in bins under each voltage.
    TfidFit2=cumsum(fitTrip(end:-1:1));
    TfidFit=TfidFit2(end:-1:1);

    
    fid=(SfidFit*fitpar2(4)+TfidFit*fitpar2(3))/(fitpar2(3)+fitpar2(4)); 
    [mfid indfid]=max(fid);
    


    fprintf('Max fidelity is %f at Vt equal to %f \n',mfid,Vt(indfid));