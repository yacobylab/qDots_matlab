file=uigetfile('');
load(file)

%1D scan: data{1} has form reps x eps x time
%the two pulses are given end to end in time part. 
data=data{1};
[reps epssteps timesteps]=size(data);
   
    
data=squeeze(nanmean(data));   %Take mean across reps
data_shaped = reshape(data,epssteps,timesteps/2,2);    %Split the two pulses into different rows
data_diff = -squeeze(diff(data_shaped,[],3)); %Take the difference between two pulses 

samprate = scan.consts.val(1);
x = (1:size(data_diff,2))./samprate * 1e6;
 mask = round(3.5*samprate*1e-6):round(size(data_diff, 2)-.5*samprate*1e-6);
%% 1d fitting and plotting. 
figure(100)
x0=x(1);
fitfn=@(p, x)p(1)*exp(-(x-x0)./p(2)); %decaying exp.
for eps_ind=1:epssteps
   beta0 = [max(data_diff(eps_ind,mask)),1]; 
   fits(eps_ind,1:2)=fitwrap('plinit plfit fine samefig', x(mask),data_diff(eps_ind, mask), beta0,fitfn);        
   fprintf('Readout T1: %g usec\n',fits(eps_ind,2));
end

%Plot T1 vs. eps
eps=1e3*linspace(scan.loops(1).rng(1),scan.loops(1).rng(2),epssteps);
plot(eps,fits(:,2)) 
xlabel('Epsilon (mV)')
ylabel('T1 (us)')

%% 2d data processing 
file=uigetfile('');
load(file)
data=data{1};
[epsy epsx timesteps]=size(data);
       
data_shaped = reshape(data,epsy,epsx,timesteps/2,2);    %Split the two pulses into different rows
data_diff = -squeeze(diff(data_shaped,[],4)); %Take the difference between two pulses 

samprate = scan.consts.val(1);
x = (1:size(data_diff,3))./samprate * 1e6;
 mask = round(3.5*samprate*1e-6):round(size(data_diff,3)-.5*samprate*1e-6);
 %% 2d fitting
 fits=[];
 figure(100)
 x0=x(1);
 fitfn=@(p, x)p(1)*exp(-(x-x0)./p(2));
for eps_indx=1:epsx
    for eps_indy=1:epsy
    beta0 = [max(data_diff(eps_indy,eps_indx,mask))*exp(min(x(mask))),1]; 
   fits(eps_indx,eps_indy,1:2)=fitwrap('plinit plfit fine samefig', x(mask),squeeze(data_diff(eps_indy,eps_indx, mask))', beta0,fitfn);        
   fprintf('Readout T1: %g usec\n',fits(eps_indx,eps_indy,2));
    end
end
%% 2d plotting
epsyval=linspace(scan.loops(1).rng(1),scan.loops(1).rng(2),epsy)
epsxval=linspace(scan.loops(2).rng(1),scan.loops(2).rng(2),epsx)
imagesc(epsxval, epsyval,fits(:,:,2))
ylabel('PlsRamp2')
xlabel('PlsRamp1')


