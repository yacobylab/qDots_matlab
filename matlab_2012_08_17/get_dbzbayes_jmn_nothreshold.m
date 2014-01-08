function [meandbz] = get_dbzbayes_jmn_nothreshold(fdata,times)
%function dbzs = get_dbzbayes(fdata,times,mdbz)
    % estimates the frequency of dbz based on fourier data. inputs are
    % fdata: the fourier data: N_reps x number of different evo times used 
    %           N_reps different frequencies will be estimated
    % times: the different evolution times used to estimate it
    % B: bandwidth of estimation (used to figure out how many times
    %           aliased)
    
    %version uses thresholded data
    %version implements Bayes inversion
    
    %assume that data is measuring each time once
    
meandbz=zeros(size(fdata,1),1);
%fdata=2*(fdata<0.4)-1; %threshold comes from histogram analysis in histogramqubitdetection.m, output is now ones and -1's, 1 corresponds to a Singlet outcome which is consistent with convention in our notes

alpha=.25;
beta=.67;

Nfreqs=256;
indices=linspace(1,Nfreqs,Nfreqs);

Nmeas=length(times);
Nexps=size(fdata,1);

B=1/(2*abs(times(1)-times(2)));

clear avg;
clear fftavg;
clear fftfreq;

for k=1:Nexps
prior=ones(1,Nfreqs);
post=ones(1,Nfreqs);
for j=1:Nmeas
    rk=fdata(k,j);
    p=.5*(1+rk.*(alpha+beta.*cos(pi*j.*[1:Nfreqs]/Nfreqs)));
    post=post.*p;
    post=post./max(post);
    avg(k,j)=sum(indices.*post./(sum(post)));  
end
    meandbz(k)=(Nfreqs-avg(k,Nmeas))/(Nfreqs)*B+B; %Hard coded alias correction
end



