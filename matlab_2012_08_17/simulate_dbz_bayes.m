function out = simulate_dbz_bayes(config)

if ~exist('config','var')|| isempty(config)
   config = struct(); 
end
if iscell(config)
   config = struct(config{:}); 
end

config = def(config,'opts','profile'); 
config = def(config,'omega_dbz',.082*2*pi); omega_dbz = config.omega_dbz;
config = def(config,'omega_max',.2*2*pi); omega_max = config.omega_max;
config = def(config,'omegas',linspace(0,omega_max,1e5)); omegas = config.omegas;
tau = pi/omega_max;
config = def(config,'N_reps',100);N_reps = config.N_reps; % number of times to repeat the simulation
config = def(config,'big_N',500); big_N = config.big_N;
%big N in Sergeevich et al. Number of total measurements we want to make to estimate dbz

config = def(config,'evotime',(1:100)*tau); evotime = config.evotime;
n_meas = big_N*ones(size(evotime))/length(evotime);
big_M = big_N/length(evotime);
evolist = [];
for j = 1:length(evotime)
   evolist = [evolist, ones(1,n_meas(j))*evotime(j)]; 
end
simple = 0;%length(evotime)==1;
config = def(config,'figind',200); figind = config.figind;

%set up p(omega)
config = def(config,'d_dbz',.01*2*pi); d_dbz = config.d_dbz;
if d_dbz==0 %delta fn
    p_omega = 0*omegas; [m mi] = min(abs(omegas-omega_dbz)); p_omega(mi) = 1;
    sigma_prior = .005*2*pi;
    prior_omega = exp(-(omegas-omega_dbz).^2/(2*sigma_prior^2)); prior_omega = prior_omega/sum(prior_omega);
    %prior_omega = zeros*omegas;
    %prior_omega(omegas>omega_dbz-.005*2*pi & omegas < omega_dbz+.005*2*pi) = 1; prior_omega/sum(prior_omega);
else
    p_omega = exp(-(omegas-omega_dbz).^2/(2*d_dbz^2)); p_omega = p_omega/sum(p_omega);
    prior_omega = p_omega;
end

if isopt(config,'flat')
   prior_omega = ones(size(omegas)); prior_omega = prior_omega/sum(prior_omega);
end

out.evotime = evotime;
out.n_meas = n_meas;
out.p_omega = p_omega;
out.omegas = omegas;

%sample an omegan for fun
real_omega = sample_pdf(omegas,p_omega);
fprintf('real omega = %d \n omega_dbz = %d \n',real_omega,omega_dbz);
%fprintf('quantization error = %d\n', real_omega-omega_dbz);
figure(figind); clf; hold on; figind = figind+1;
plot(omegas,p_omega);
xlabel('\omega');ylabel('P(\omega)');
YL = get(gca,'YLim');
plot(real_omega*[1 1],YL,'r');

%plot one rep of experiment (n_meas measurements)
meas = binornd(1,cos(real_omega*evotime(1))^2,1,n_meas(1));
figure(figind); clf; hold on; figind=figind+1;
plot(meas,'.'); 
title('Singles Shot Values');
xlabel('Triplet Probability');
plot([0, length(meas)],mean(meas)*[1 1],'r');
plot([0, length(meas)],cos(real_omega*evotime(1))^2*[1 1],'g');
plot(cumsum(meas)./(1:length(meas)),'c');
legend({'Single Shot Vals','Mean of shots', 'Actual \Omega_{dBz}','Running Average'});

%sample the data
results = struct();

tic
for j = 1:N_reps
   real_omega = sample_pdf(omegas,p_omega);
   results(j).real_omega = real_omega;
   results(j).meas = [];
   for k = 1:length(evotime)
       results(j).meas = [results(j).meas, binornd(1,cos(real_omega*evotime(k))^2,1,n_meas(k))];
   end   
end
fprintf('sampling the data took %f seconds \n',toc);

if isopt(config,'compress') %set up compressive sensing matrix
    L = length(evotime);
    T_samp = abs(diff(evotime(1:2))); F_samp = 1/T_samp;
    NFFT = 2^nextpow2(L);
    ffreqs = omega_dbz+2*pi*F_samp/2*linspace(-.5,.5,NFFT/2+1);
    fine_freqs = linspace(ffreqs(1),ffreqs(end),10^4);
    t  = linspace(0,2*pi,1000);
    x  = linspace(0,2*pi,length(t));
    for j = 1:length(ffreqs)
    for k = 1:length(fine_freqs)
        y = cos(ffreqs(j)*t).*cos(fine_freqs(k)*t);
        A(j,k) = trapz(y,x);
    end  
    end
end

if ~simple
    %plot distribution of sampled omegas
    if N_reps>1
        [h b]=hist([results.real_omega],50);
        figure(figind); clf; hold on; figind=figind+1;
        plotyy(omegas,p_omega,b,h);
        xlabel('\omega (Grad/s)');
        ylabel('P(\omega)');
        legend('original distribution', 'sampled');
        %plot(omegas,p_omega);
        %plot(b,h/sum(h),'r')
    end
    tic
    for j = 1:N_reps
        P_MLE = prior_omega; %prior distribution
        %figure(1); plot(omegas,P_MLE);pause
        if isopt(config,'meantimes')
            f=mean(reshape(results(j).meas,big_M,length(evotime)));
            eee = evotime;
        else
            f=results(j).meas;
            eee = evolist;
        end
        results(j).P_k = zeros(size(f));
        if isopt(config,'bayes')
        for k = 1:length(f)
            P_MLE = P_MLE.*(1+(-2*f(k)+1)*cos(omegas*eee(k)))/sum(P_MLE);
            [m mi]=max(P_MLE);
            if isopt(config,'video')
                figure(1); clf; plot(omegas,P_MLE); hold on;
                YL = get(gca,'YLim'); plot(results(j).real_omega*[1 1],YL,'r');
                title(sprintf('MSE = %d',(omegas(mi)-results(j).real_omega)^2/(results(j).real_omega)^2));
                pause;
            end   
            results(j).var(k) = std(P_MLE)^2;
            results(j).P_k(k) = ((omegas(mi)-results(j).real_omega)/results(j).real_omega)^2;
            %title(sprintf('MSE = %d',(omegas(mi)-results(j).real_omega)^2/(results(j).real_omega)^2));
            %pause%keyboard
            [m mi]=max(P_MLE);
            results(j).P_MLE = P_MLE;
            results(j).MLE = omegas(mi);
            results(j).MSE_bayes = (omegas(mi)-results(j).real_omega)^2/(results(j).real_omega)^2;
        end
        end
        if isopt(config,'fourier')
            L = length(eee);
            T_samp = abs(diff(eee(1:2))); F_samp = 1/T_samp;
            NFFT = 2^nextpow2(L);
            freqs = omega_dbz+2*pi*F_samp/2*linspace(-.5,.5,NFFT/2+1);
            ft = fft(f-mean(f),NFFT)/L;
            P_MLE = 2*abs(ft(1:NFFT/2+1));%.*(prior_omega);
            if isopt(config,'compress')
                P_MLE = compressive_sample(P_MLE,A);
            end
            
            [m mi]=max(P_MLE);
            results(j).P_MLE = P_MLE;
            if 1
                f2 = freqs(mi)+F_samp*linspace(-1,1,length(freqs));
                prior = ones(size(f2)); prior= prior/sum(prior);
                P_MLE = prior;
                for k = 1:length(f)
                    P_MLE = P_MLE.*(1+(-2*f(k)+1)*cos(f2*eee(k)))/sum(P_MLE);
                    figure(1); clf; plot(f2,P_MLE); pause
                end
                [m mi]=max(P_MLE);
                results(j).var(k) = std(P_MLE)^2;
                
                [m mi]=max(P_MLE);
                results(j).P_MLE = P_MLE;
                results(j).MLE = omegas(mi);
                results(j).MSE_bayes = (omegas(mi)-results(j).real_omega)^2/(results(j).real_omega)^2;
                
                
            end
            results(j).MLE = freqs(mi);
            results(j).MSE_bayes = (freqs(mi)-results(j).real_omega)^2/(results(j).real_omega)^2;
            %keyboard
        end
        %keyboard
    end
    if isopt(config,'profile'); fprintf('analysis took %f seconds \n',toc); end;
    %out.mean_MSE =(mean(vertcat(results.P_k)-repmat([results.real_omega],big_N,1)')).^2;
    out.mean_MSE =(mean(vertcat(results.P_k)));
    figure(figind); clf; hold on; figind = figind+1;
    plot(out.mean_MSE);
    xlabel('number of measurements');
    ylabel('mean MSE');
    
    try
        fitfn = @(p,x) p(1)*x.^(p(2));
        pars = fitwrap('plinit plfit',1:length(out.mean_MSE),out.mean_MSE,[1 -2],fitfn);
        figure(figind); clf; figind = figind+1;
        loglog(out.mean_MSE,'bx'); hold on;
        loglog(1:length(out.mean_MSE),fitfn(pars,1:length(out.mean_MSE)),'r');
        xlabel('number of measurements');
        ylabel('mean MSE');
    catch
        fprintf('power law fitting of MSE failed \n');
    end
    
    [h b]=hist([results.MLE],50);
    figure(figind); clf; hold on; figind=figind+1;
    plotyy(omegas,p_omega,b,h);
    xlabel('\omega (Grad/s)');
    ylabel('P(\omega)');
    legend('original distribution', 'MLE disribution');
    
    figure(figind); clf; figind = figind+1; hold on;
    plot([results.real_omega],[results.MLE],'.');
    xlabel('real \Omega'); 
    ylabel('Extracted \Omega');
    
end % end not simple

if simple
    tic;
    if length(results)>1
        %now plot results of many simulations
        for j = 1:length(results)
            MSE_tot(j,:) =  abs((cumsum(results(j).meas)./(1:length(results(j).meas)))-cos(results(j).real_omega*evotime)^2)/cos(results(j).real_omega*evotime)^2;
        end
        MSE_all = mean(MSE_tot);
        %fit a power law
        fitfn = @(p,x) p(1)*x.^(p(2));
        pars = fitwrap('plinit plfit',1:length(results(1).meas),MSE_all,[1 -2],fitfn);        
        figure(figind); clf; figind=figind+1;
        loglog(MSE_all,'bx'); hold on;
        loglog(1:length(results(1).meas),fitfn(pars,1:length(results(1).meas)),'r');
        xlabel('Number of Measurements');
        ylabel('Mean Square Error');
        title(sprintf('Average of %i reps',N_reps));
        legend({'simulation',sprintf('MSE ~ N^{%.1f}',pars(2))});
    end
    if isopt(config,'profile'); fprintf('analysis took %f seconds \n',toc); end;
end

if 0
%look at a bayesian estimator
for j = 1:N_reps
    P_MLE = p_omega; %prior distribution
    figure(1); plot(omegas,P_MLE);
%     for k = 1:length(omegas)
%        P_MLE(k) = sum(log(1+(2*results(j).meas-1)*cos(omegas(k)*evotime)));
%     end

%     for k = 1:length(omegas)
%        P_MLE(k) = prod(1+(2*results(j).meas-1)*cos(omegas(k)*evotime));
%     end
%     keyboard
    for k = 1:length(results(1).meas)
        P_MLE = P_MLE.*(1+(2*results(j).meas(k)-1)*cos(omegas*evotime))/sum(P_MLE);
        figure(1); plot(omegas,P_MLE); keyboard
    end
    [m mi]=max(P_MLE);
    results(j).MLE = omegas(mi);
    results(j).MSE_bayes = (omegas(mi)-results(j).real_omega)^2/(results(j).real_omega)^2;
    keyboard
end

fprintf('mean MSE from bayes = %d',mean([results.MSE_bayes]));
end


out.results = results;
out.config = config;
end

function out = compressive_sample(ft_old, A)
opts = spgSetParms('verbosity',0);
out = spg_bp(A, ft_old', opts);
end

function val = sample_pdf(xv,prob)
     [~,ii]= find(cumsum(prob)>rand,1,'first');
    val = xv(ii);
end

% Apply a default.
function s=def(s,f,v)
  if(~isfield(s,f))
      %s=setfield(s,f,v);
      s.(f) = v;
  end
end

function b=isopt(config,name)
  b=~isempty(strfind(config.opts,name));
end 

%omega_dbz= .082*2*pi; %(80MHz);
%omega_max = .25*2*pi; %400MHz
%omegas = linspace(0,omega_max,100000);
%N_reps = 1000; 
%big_N = 500; 
%set up evo times
%evotime =59; simple = 1;% one evo time
%evotime(1:10)*tau; %many evo times
%figind = 200; 

%generate p(omega)
%delta function
%p_omega = 0*omegas; [m mi] = min(abs(omegas-omega_dbz)); p_omega(mi) = 1;

%normally distributed
%d_dbz = .005*2*pi; 
%p2_omega = exp(-(omegas-omega_dbz).^2/(2*d_dbz^2)); p2_omega = p2_omega/sum(p2_omega);
%p_omega = normpdf(omegas,omega_dbz,d_dbz); p_omega = p_omega/sum(p_omega);


% %lots of samples at one evo time
% results = struct();
% for j = 1:N_reps
%    real_omega = sample_pdf(omegas,p_omega);
%    results(j).meas = binornd(1,cos(real_omega*evotime)^2,1,n_meas);
%    results(j).real_omega = real_omega;
% end