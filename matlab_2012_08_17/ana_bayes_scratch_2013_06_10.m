data=-1+2*indata;
data = 2*(data>0)-1;
domega = 1/abs(diff(xv(1:2))); %in GHz
omega = .08+.5*linspace(-domega,domega,2000);
best_omega = zeros(1,size(data,1));
figure(1); clf
tic
for j = 1:10%size(data,1)
    p_omega = ones(size(omega));
    for k = 1:size(indata,2)
        p_omega = p_omega.*(1+data(j,k)*cos(pi*xv(k)*omega/domega));
        plot(omega,p_omega); pause
    end
    p_omega(1:end/2) = 0;
    p_omega(end-10:end) = 0;
    %plot(omega,p_omega); pause     
    [m mi] = max(p_omega);
    best_omega(j) = omega(mi);
end

% for j = 100
%     %p_omega = ones(size(omega));
%     for k = 1:length(omega)
%         p_omega(k) = sum(log(1+data(j,:).*cos(pi*xv(1:1000)*omega(k)/domega)));
%         %plot(omega,p_omega); pause
%     end
%     %p_omega(1:end/2) = 0;
%     %p_omega(end-10:end) = 0;
% end

%plot(omega,p_omega);
toc

%% load file and set up two scales
s=load('sm_dBz_phase_bayes_R_4219.mat');
data=s.data;
scan=s.scan;

scantime=getscantime(scan,data);
[t1t t1] = att1('right',scantime,'after'); % FIX ME
data_all=anaHistScale(scan,data,t1); % FIX ME
data=squeeze(data_all{1});
mmask = find(any(isnan(data)),2);
data(mmask,:)=[];
xv=plsinfo('xval',s.scan.data.pulsegroups.name,[],scantime);
data1 = data(:,1:1000); data2 = data(:,1001:2000);
%data1 = data(:,1:2:2000); data2 = data(:,2:2:2000); 
evos = xv(1+2*size(data1,2):end);
Tevo1 = xv(1);
Tevo2 = xv(501);

nevo = length(evos); rrng = 1:nevo; evo_inds = (1+size(data1,2)+size(data2,2)):length(xv);
%evos=  xv(evo_inds);
Nrep = size(data,1);
scale1 = mean(data1,2);
scale2 = mean(data2,2);


%% fit raw data
fitfn = @(p,x) p(1) +(p(2)*sin(2*pi*x/p(3) +p(4))).*exp(-((x+p(6))/p(5)).^2);
c=[];
%first without scaling, to get dbz_avg and deltaDBZ (1/t2*)
beta0 = [.3, .3, 10, .001, 80 0];
pars = fitwrap('plinit plfit',evos,nanmean(data(:,evo_inds),1),beta0,fitfn, [1 1 1 1 1 0]);
dbz_avg = 2*pi/pars(3);
d_dbz = 1/abs(pars(5));
out.dbzmean=pars; 

%% fit first scale
t=linspace(1,nevo,nevo);
ti=linspace(1,nevo+50,nevo);
rep=1:length(scale1);
alpha=linspace(-0.25,.25,200);
clear dBx dataI tsV
for i=1:length(alpha);
    for j=1:Nrep;
        ts=(1+alpha(i)*scale1(j))*t;
        tsV(j,:)=ts;
        dataI(j,:)=qinterp1(ts,data(j,evo_inds),ti);
    end
    dBx(:,i)=nanmean(dataI);
end
[tV,repV]=meshgrid(t,rep);
[tiV,repV]=meshgrid(ti,rep);
[alphaV,talpha]=meshgrid(alpha,ti);

%now fit all scales
for j = 1:length(alpha)
    beta0 = pars;%[.3, .3, 14, .001, 80 0 1.1];
    mask = ~isnan(dBx(rrng,j));
    %pars = fitwrap('plinit plfit',(talpha(rrng,j)./(1+alphaV(rrng,j)))',dBx(rrng,j)',beta0,fitfn, [1 1 1 1 1 0]);
    pars = fitwrap('plinit plfit',(talpha(mask,j))',dBx(mask,j)',beta0,fitfn, [1 1 1 1 1 0]);
    c(j) = abs(pars(5));
    %pause
end
c(c>5000)=0;
[T2,I]=max(c);
out.T2 = T2;
alphabest=alpha(I);
figure(1);  clf;
plot(alpha,c,'x')
xlabel('scale factor');
ylabel('T_2^* (ns)');


%% second scale factor
beta=linspace(-.001,.001,200);
omegas = 1./(alphabest*scale1);
clear dBx dataI tsV
for i=1:length(beta);
    for j=1:Nrep
        ts=(1+alphabest*scale1(j)+(scale1(j)^-1)*beta(i)*scale2(j)^1)*t;
        %if abs(1/scale2(j))<10
        %ts=(1+alphabest*scale1(j)+beta(i)*scale2(j)^-1)*t;
        %else
         %   ts=(1+alphabest*scale1(j))*t;
        %end
        tsV(j,:)=ts;
        dataI(j,:)=qinterp1(ts,data(j,evo_inds),ti);
    end
    dBx(:,i)=nanmean(dataI);
end
[tV,repV]=meshgrid(t,rep);
[tiV,repV]=meshgrid(ti,rep);
[betaV,tbeta]=meshgrid(beta,ti);



clear c
%now fit all scales
for j = 1:length(beta)
    beta0 = [.3, .3, 10, .001, 80 0];
    %pars = fitwrap('plinit plfit',(talpha(rrng,j)./(1+alphaV(rrng,j)))',dBx(rrng,j)',beta0,fitfn, [1 1 1 1 1 0]);
    pars = fitwrap('plinit plfit',(tbeta(rrng,j))',dBx(rrng,j)',beta0,fitfn, [1 1 1 1 1 0]);
    c(j) = abs(pars(5));
    %pause
end
c(c>5000)=0;
[T2,I]=max(c);
out.T2 = T2;
betabest=beta(I);
figure(1);  clf;
plot(beta,c,'x')
xlabel('scale factor');
ylabel('T_2^* (ns)');

for j = I
    beta0 = [.3, .3, 10, .001, 80 0];
    %pars = fitwrap('plinit plfit',(talpha(rrng,j)./(1+alphaV(rrng,j)))',dBx(rrng,j)',beta0,fitfn, [1 1 1 1 1 0]);
    pars = fitwrap('plinit plfit',(tbeta(rrng,j))',dBx(rrng,j)',beta0,fitfn, [1 1 1 1 1 0]);
    bestT2star = abs(pars(5));
    %pause
end

%% second scale factor
gamma=linspace(-.1,.1,200);
clear dBx dataI tsV
for i=1:length(gamma);
    for j=1:Nrep;
        ts=(1+alphabest*scale1(j)+betabest*scale2(j)^3+gamma(i)*scale1(j)^2)*t;
        tsV(j,:)=ts;
        dataI(j,:)=qinterp1(ts,data(j,evo_inds),ti);
    end
    dBx(:,i)=mean(dataI);
end

[tV,repV]=meshgrid(t,rep);
[tiV,repV]=meshgrid(ti,rep);
[gammaV,tgamma]=meshgrid(gamma,ti);

clear c
%now fit all scales
for j = 1:length(beta)
    beta0 = [.3, .3, 10, .001, 80 0];
    %pars = fitwrap('plinit plfit',(talpha(rrng,j)./(1+alphaV(rrng,j)))',dBx(rrng,j)',beta0,fitfn, [1 1 1 1 1 0]);
    pars = fitwrap('plinit plfit',(tgamma(rrng,j))',dBx(rrng,j)',beta0,fitfn, [1 1 1 1 1 0]);
    c(j) = abs(pars(5));
    %pause
end
c(c>5000)=0;
[T2,I]=max(c);
out.T2 = T2;
bestgamma=gamma(I);
figure(1);  clf;
plot(gamma,c,'x')
xlabel('scale factor');
