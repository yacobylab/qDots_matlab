function ana_cds(ctrl,csdata)
%Controls: gain, off decide if we are plotting a gain or offset sweep
%dcy lets you fit the decay. FIX ME
%This function will process offset and gain sweeps from the CDS box. 
%It will also check that the average DMM voltage and average DAQ voltage is constant 
%Requires that you log dBz and parameters of the CS data. If you don't, can add manually: 
%csdata(1)=NSamps csdata(2)=SampleCount csdata(3)=offset csdata(4)=dbz_ave
%to plot the DMM vals, use dmmshot

if ~isempty(strfind(ctrl,'gain'))
    sweepstr='Box Gain';     
    if ~isempty(strfind(ctrl,'act'))
        sweepstr='Actual Gain';       
    end
elseif ~isempty(strfind(ctrl,'off'))
    sweepstr='Box Offset'; 
elseif ~isempty(strfind(ctrl,'nsamps'))
    sweepstr='NSamps'; 
end
 

 rabi_amp=@(p,x) p(1).^2./(p(1).^2+(x-p(2).^2)/4);             
 rabi_freq=@(p,x) (p(1).^2+(x-p(2).^2)/4).^0.5; 
%we can measure amplitude and frequency and find w21 from that, yes?
%so we would fit the amplitude to this, and or we owuld     
filename=uigetfile('','MultiSelect','on');
if ~iscell(filename)
    filename={filename}; 
end
decay=2; 
fitfn = @(p,x) p(1) +(p(2)*sin(2*pi*x/p(3) +p(4))).*exp(-((x+p(6))/p(5)).^p(7));
beta0 = [.35, .35, 140, -pi/2, 600 0 2];
figure(7); clf; figure(8); clf; 
for i=length(filename):-1:1    
    load(filename{i});        
    if isfield(scan.data,'post_dbz') && ~isempty(scan.data.post_dbz)
        dbz_ave(i)=(scan.data.pre_dbz+scan.data.post_dbz)./2;  %gives it in MHz. 
    elseif isfield(scan.data,'pre_dbz') && ~isempty(scan.data.pre_dbz)
        dbzave(i)=scan.data.pre_dbz;
    else
        dbzave(i)=csdata(4); 
    end

    scantime=getscantime(scan,data);
    [t1t t1] = att1('right',scantime,'after'); % FIX ME
    [data_all,scalefuncs,mv,fp]=anaHistScale(scan,data,t1); % FIX ME
    data=squeeze(data_all{1});
    
    xv=plsinfo('xval',scan.data.pulsegroups.name,[],scantime); 
    t=xv(2,:); %vector of times for each pulse. 
    shottime=t(1); %the time at which we sample repeatedly. 
    ind=find(diff(t)~=0,1,'first'); %to sort the dbZ times from the shot times. 
    evo_inds=(ind+1):length(data); 
    evos=xv(2,evo_inds); 
    shots=2:ind; 
    beta0(1)=mean(nanmean(data(:,evo_inds),1)); beta0(2)=range(nanmean(data(:,evo_inds),1));
%     cosfn5 = '@(y, x)y(1)+y(2)*cos(y(4)*x+y(3)).*exp(-(x-y(5)).^2 * y(6).^2)';
%     cosfn2 = '@(y, x)y(1)+(y(2)*cos(y(4)*x) + y(3) * sin(y(4)*x)).*exp(-(x-y(5)).^2 * y(6).^2)';
%     fifn.fn = @fioscill;
%     fifn.args = {2};
    try        
%         fp=fitwrap('plinit plfit',evos,nanmean(data(:,evo_inds),1),fifn, cosfn5, [1 1 1 1 0 1]);
%         fp = [fp(1), fp(2)*cos(fp(3)), -fp(2)*sin(fp(3)), fp(4:6)];
%         pars = fitwrap('plinit plfit',evos,nanmean(data(:,evo_inds),1),fp,cosfn2,[1 1 1 1 0 1]); 
        pars=[];
        pars = fitwrap('plfit plinit',evos,nanmean(data(:,evo_inds),1),beta0,fitfn, [1 1 1 1 1 0 0]);
        pars = fitwrap('plfit plinit',evos,nanmean(data(:,evo_inds),1),pars,fitfn, [1 1 1 1 1 0 1]);%fit oscillations 
    catch
        if isempty(pars); 
            pars=nan(1,6);
        end
    end
    %w/ cosfn 2, 
    % y(1)=offset sqrt(y(2)^2+y(3)^2)=amp t2*=1/y(6) freq=2*pi*y(4) (in Ghz)
    % phase=arctan(-y(3)/y(4); 
%     off_all(i)=pars(1); % in V
%     amp_all(i)=(pars(2).^2+pars(3).^2).^0.5.*diff(mv); %in V
%     t2(i)=abs(1./pars(6)); %in ns
%     rabi_avg(i)=1e3/(2*pi)*pars(4); %in MHz
%     phase(i)=atan(-pars(3)./pars(2)); %in radians
    

     rabi_avg(i) = 1000./pars(3); %given in MHz     
     t2(i) = abs(pars(5)); %Given in ns      
     off_all(i)=pars(1); 
     amp_all(i)=pars(2)*diff(mv);     
     phase(i)=pars(4)+pi/2; %we expect oscillations to start at -pi/2 = singlet
    
    v_ave=mean(nanmean(data(:,shots))); %the average input voltage, gives a sense of ddBz
    fprintf('Actual val %0.3f, Optimal val %0.3f \n',v_ave,off_all(i)); %optimal val is at '0' pt of sine function. 
    
    %beta0=pars; %update guess.  
    d3(i,:)=nanmean(data(:,evo_inds,1)); %This keeps track of data for the color plot. 
    
    if (isfield(scan.data,'CSNSamps') && ~isempty(scan.data.CSNSamps)) Ns=scan.data.CSNSamps; 
    else  Ns=csdata(1);  end
    
    if isfield(scan.data,'CSSampleCount') && ~isempty(scan.data.CSSampleCount) Sc=scan.data.CSSampleCount; 
    else  Sc=csdata(2);  end
    
    if isfield(scan.data,'CSOffset') && ~isempty(scan.data.CSOffset) off(i)=scan.data.CSOffset;
    else off(i)=csdata(3);  end    
    gtrans=@(x) x.*2^(-34).*500.*Ns.*Sc;
    gain(i)=scan.data.CSGain;
    if ~isempty(strfind(ctrl,'gain'))
        sweepvar(i)=gain(i);
        if ~isempty(strfind(ctrl,'act'))            
            sweepvar(i)=gtrans(sweepvar(i));
        end
    elseif ~isempty(strfind(ctrl,'off'))
        sweepvar(i)=off(i);
    elseif ~isempty(strfind(ctrl,'nsamps'))
        sweepvar(i)=Ns;
    end
    figure(7);
    imagesc(evos,sweepvar,d3);
    ylabel(sweepstr);
    set(gca,'YDir','normal')
    xlabel('Time (ns)');
    
      
    daqave(i)=v_ave*diff(mv)+mv(1); 
    
    if length(scan.loops(1).getchan>2)
        dmmave(i)=nanmean(data{3});
        dmmmax(i,:)=max(data{3});
        dmmmin(i,:)=min(data{3});
        dmmstd(i)=nanstd(data{3});        
        figure(8)
        if ~isempty(strfind(ctrl,'dmmshot'))            
            nfls=length(filename);
            r=ceil(nfls/4);
            dmmshots=data{3}(:);
            if isfield(scan.data,'setpt') && ~isempty(scan.data.setpt)
                dmmshots=(dmmshots-scan.data.setpt)./gtrans(gain(i));
            end
            subplot(4,r,i)
            plot(nanmean(data(:,shots),2).*diff(mv),dmmshots,'.')
        end
    end
     v_ave_all(i)=v_ave; 

end

    daqvar=daqave-mean(daqave);
    if isfield(scan.data,'ampl') && ~isempty(scan.data.ampl)
        ampl=scan.data.ampl; 
    else
        ampl=nanmean(amp_all); 
    end
    df=daqvar./(nanmean(amp_all).*shottime*1e-9*2*pi)/1e6;
    figure(10); clf;
    subplot(1,2,2);    
    plot(sweepvar,df,'.-')
    xlabel(sweepstr)        
    ylabel('Error in Average dBz (MHz)')

    Nreps=size(data,1);
    
    descr1=sprintf('Actual val %0.3f, Optimal val %0.3f \n',mean(v_ave_all),mean(off_all));
    descr2=sprintf('NSamps: %d. SampleCount: %d PreDelay: %d. PostDelay: %d NReps: %d \n',Ns,Sc,scan.data.CSPreDelay,scan.data.CSPostDelay,Nreps);
    descr3=sprintf('Average Rabi (Mhz): %2.3f, Shot Time: %d \n Phase: %1.2f',mean(rabi_avg),shottime,mean(phase)/pi);
    
if length(scan.loops.getchan)>2

    %Now, print the variation in daq frequencies.
    figure(10); hold on;
    subplot(1,2,1);
    
    %save=8./1e6;
    v=[4 5 6 7 8 9]; %These are the frequencies of the VCO at a set of voltages.
    f=[117.4 124.2 131.4 139.1 147.3 155.9].*1e6;
    %s=[6.75 7.09 7.59 8.1 8.53 8.82].*1e6;
    freqave=interp1(v,f,dmmave)./1e6;
    
    freqmax=interp1(v,f,dmmmax)./1e6;
    freqmin=interp1(v,f,dmmmin)./1e6;
    freqmaxave=interp1(v,f,dmmave+dmmstd)./1e6;
    freqminave=interp1(v,f,dmmave-dmmstd)./1e6;
%     
%     freqactmin=save.*(dmmmin'-dmmave)+freqave;
%     freqactmaxave=save.*(dmmstd)+freqave;
%     freqactminave=-save.*(dmmstd)+freqave;
%     freqactmax=save.*(dmmmax'-dmmave)+freqave;
    plot(sweepvar,freqave,'.-',sweepvar,freqmax,'.-',sweepvar,freqmin,'.-',sweepvar,freqmaxave,'.-',sweepvar,freqminave,'.-')%,gain,freqactmax,'.-',gain,freqactmin,'.-',gain,freqactmaxave,'.-',gain,freqactminave,'.-');
    %legend('Average','Max','Min','1 std','1 std')
    xlabel(sweepstr); ylabel('Frequency Output by Specs')
end    
figure(9); clf;

[ax h1 h2]=plotyy(sweepvar,rabi_avg,sweepvar,t2);
xlabel(sweepstr)
set(get(ax(1),'Ylabel'),'String','Rabi (MHz)'); 
set(get(ax(2),'Ylabel'),'String','T2* (ns)'); 
set(h1,'LineStyle','.'); 
set(h2,'LineStyle','.');
set(h1,'LineStyle','-'); 
set(h2,'LineStyle','-');
figure(11); clf; 
[ax2 h21 h22]=plotyy(sweepvar,rabi_avg,sweepvar,1./amp_all);
set(get(ax2(1),'Ylabel'),'String','Rabi (MHz)'); 
set(get(ax2(2),'Ylabel'),'String','Inverse Amplitude (1/V)'); 
set(h21,'LineStyle','.'); 
set(h21,'LineStyle','-'); 
set(h22,'LineStyle','.');
set(h22,'LineStyle','-');
xlabel(sweepstr); 

figure(12); clf; 
plot(sweepvar,phase); 
xlabel(sweepstr); 
ylabel('Phase (pi)')
ppt=guidata(pptplot);
prettyfile1=regexprep(filename{end},'(sm_)|(\.mat)','');
prettyfile2=regexprep(filename{1},'(sm_)|(\.mat)','');
prettyfile=sprintf('%s to %s',prettyfile1,prettyfile2);
set(ppt.e_title,'String',prettyfile);
set(ppt.e_figures,'String',['[',sprintf('%d ',7:10),']']);
set(ppt.exported,'Value',0);
set(ppt.e_body,'String',[descr1 descr2 descr3]);
set(ppt.e_file,'String',filename{end});

end


