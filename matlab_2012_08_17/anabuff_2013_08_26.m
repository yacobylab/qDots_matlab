%this is for looking at the output of the daq. 
%we are interesting in several different types of data, so I will leave this as cells to copy
%and paste between. 

%filename=uigetfile('','MultiSelect','on');
%if ~iscell(filename)
%    filename={filename}; 
%end


%%
    for i=1:8
        figure(i)
    end
for i=1:16 
    filename=sprintf('BuffDAQ_%02d',i);    
     load(filename); 
     cv = CDSvals(:,1:99); dv = DAQvals(:,2:100);
     %figure(i)
     
     scind=mod(i,4); 
     if scind==0
         scind=4; 
     end
     pdind=ceil(i/4); 
     figure(scind); 
     subplot(2,2,pdind)        
     plot(dv(:),cv(:),'.')   
          title(sprintf('Sample Count %d, PostDelay %d',scind*2,pdind+1));     

     fitfn=@(p,x) p(1).*x+p(2); 
     beta0=[0.5 0.001];
     call=cv(:); 
     dall=dv(:);
     [pars res] = fitwrap('plinit plfit',dall',call',beta0,fitfn, [1 1]);       
     err(i)=sqrt(sum(res.^2))/length(res); 
     sl(i)=pars(1); 
     off(i)=pars(2);
     %axis([5e-3 27e-3 -12e-3 22e-3]); 
     
      figure(pdind+4); 
     subplot(2,2,scind)        
     plot(dv(:),cv(:),'.')   
     title(sprintf('Sample Count %d, PostDelay %d',scind*2,pdind+1));     
     
     
end


%Samplecount 2 4 6 8 
%Post Delay 2 3 4 5

%% For looking at the DMM Output and the average DAQ values that cause it.  
load('DMMDAQ_01')
%cv are the DMM values, dv the DAQ values. 
%this used Ns=4, Npulse=999;

%it might be interesting to try to fit a line, find the variance, etc. Get a sense of the size
%of the cloud. 

%First plot it. 
%using 4 samplecount. 999 samples. 
Ns=4; Sc=999;
fitfn=@(p,x) p(1).*x+p(2); 
beta0=[250 4];
figure(9); clf; 

for i=1:4  %iterate over post from 2 to 5    
    for k=1:4 %iterate over pre from 2 to 5. 
        figure(9);
        subplot(4,4,sub2ind([4 4],k,i)); %%note!! this is weird! sub2ind and subplot index differently
        plot(dv{i,k}(:),cv{i,k}(:),'.'); 
        title(sprintf('Post %d, Pre %d',i+1,k+1))      
        call=cv{i,k}(:); 
        dall=dv{i,k}(:);
        [pars res] = fitwrap('plinit plfit',dall',call',beta0,fitfn, [1 1]);       
        %pause
        err(i,k)=sqrt(sum(res.^2))/length(res); 
        sl(i,k)=pars(1); 
        off(i,k)=pars(2);        
    end
end

%remember if we run multiple scans, because DAC Gain may not be constant, we should think
%largely of  the relative size across a single set of scans. not 


%% For looking at the DMM Output and the average DAQ values that cause it.  
load('DMMDAQ_04')
%cv are the DMM values, dv the DAQ values. 
Npulse=999;

%it might be interesting to try to fit a line, find the variance, etc. Get a sense of the size
%of the cloud. 

Sc=999;
fitfn=@(p,x) p(1).*x+p(2); 
beta0=[250 4];
figure(9); clf; 
err=[]; sl=[]; off=[];
n=length(dv); 
for i=1:n  %iterate over post from 2 to 5    
        subplot(1,n,i)
        plot(dv{i}(:),cv{i}(:),'.');     
        title(sprintf('Post Delay %d Sample Count %d',i,n+1-i))      
        call=cv{i}(:); 
        dall=dv{i}(:);
        %[pars res] = fitwrap('plinit plfit',dall',call',beta0,fitfn, [1 1]);       
        %pause
        %err(i)=sqrt(sum(res.^2))/length(res); 
        %sl(i)=pars(1); 
        %off(i)=pars(2);        
end

figure(10); clf; hold on; 
load('DMMDAQ_04_0')
plot(dv{4}(:),cv{4}(:),'b.');
load('DMMDAQ_06')
plot(dv{4}(:),cv{4}(:)+0.4818,'r.');



%remember if we run multiple scans, because DAC Gain may not be constant, we should think
%largely of  the relative size across a single set of scans. not 
%% Then, we are interested in the average input versus the output. 
Ns=4; Sc=999;

vin=mean(vins); 
gtrans=@(x) x.*2^(-34).*500.*Ns.*Sc;
vout=gtrans(g).*vin +off; %this is linear, so we can just take the average. 
%6.4404 versus 6.4321 


