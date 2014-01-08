function [sf config alldata alltau]=ana_fast_pt1(files,config) 
if ~exist('files','var') || isempty(files)
    files=uigetfile('','MultiSelect','on');
end
if ~exist('config','var')
    config=struct();
end
tic
config=def(config,'fit',1);
config=def(config,'window',.7);
config=def(config,'step',8);
if ischar(files)
    files={files};
end
alldata=[];
alltau=[];
for f=1:length(files)     
    sf{f}.file = load(files{f}); %each separate file is a different struct? in sf.  
    s=sf{f}.file;    
    figure(200+f);
    frames=1:size(s.data{1},1); %number of repetitions   
    imagesc(permute(nanmean(s.data{1}(frames,:,:),1),[2 3 1]));
    scantime=getscantime(s.scan,s.data);
    %d=anaHistScale(s.scan,s.data,att1('right',scantime,'before'));
    d=anaHistScaleGrpwise(s.scan,s.data,att1('right',scantime,'before')); %fixes contrast? d same size as data. 
    d=d{1}; 
    d=permute(d,[2 3 1]); % so d is  double of size [time epsilon reps]
    if 1
        lpcnt=.9; % if this percentage is below the bottom tpcnt of data, cull
        hpcnt=.9; % or hpcnt is in the top tpcnt of data.
        tpcnt=.7; % threshold percentage
        npts=size(d,2) * lpcnt;
        npts2=size(d,2) * hpcnt;
        for i=1:size(d,3)
            pts=d(:,:,i);
            pts=sort(pts(~isnan(pts)));
            if isempty(pts)
                continue;
            end
            thresh=pts(floor(end*tpcnt)); %this yields the maximum value
            badrows=(find(sum(d(:,:,i) > thresh,2) > npts))'; % gets rid of all above maximum. 
            thresh2=pts(floor(end*(1-tpcnt))); % minimum. 
            badrows=[badrows (find(sum(d(:,:,i) < thresh2,2) > npts2))'];
            if ~isempty(badrows)
                figure(f*10+i);
                clf
                subplot(1,2,1);
                imagesc(d(:,:,i));
                subplot(1,2,2);
                d(badrows,:,i)=1;
                imagesc(d(:,:,i));
                d(badrows,:,i)=nan; %and then we set all of the bad rows to be nan, which will get ignroed when we do nanmean. 
            end
        end
    end
    nodbz=1:size(d,1);
    eps = plsinfo('xval',s.scan.data.pulsegroups(2).name,3,scantime);
    tau=[];
    for k=nodbz
        if k == 1 || ~strcmp(s.scan.data.pulsegroups(k-1).name,s.scan.data.pulsegroups(k).name)
            tmp= plsinfo('xval',s.scan.data.pulsegroups(k).name,1,scantime); %all these come out as 1. hmm. 
        end
        tau(k) =tmp(1);
    end
    sf{f}.eps = eps;
    sf{f}.tau = tau;
    sf{f}.d = d;
    sf{f}.nodbz = nodbz;
    tau=tau(1);
    sf{f}.t = 1e9 * sf{f}.tau ./ s.scan.data.freqs(nodbz); % time in ns. what are freqs?
    figure(2);
    imagesc(eps,sf{f}.t,nanmean(d(nodbz,:,frames),3));
    dt=sf{f}.t;
    [overlap, oa, ob] = intersect(alltau, dt);
    dt(ob)=[];
    nodbz(ob)=[];
    size(dt)
    size(d)
    alldata((end+1):(end+length(nodbz)),:)=nanmean(d(nodbz,:,frames),3);
    alltau=[alltau dt];
end
[alltau ati] = sort(alltau);
alldata=alldata(ati,:);
[alltau ati] = unique(alltau);
alldata=alldata(ati,:);
figure(5);
clf;
imagesc(eps,alltau,alldata);

figure(6);
clf;
[x,y]=meshgrid(eps,alltau);
surf(x,y,alldata);
shading interp
view(2)
axis tight;
toc
end


function b=isopt(config,name)
b=~isempty(strfind(config.opts,name));
return;
end 

% Apply a default.
function s=def(s,f,v)
if(~isfield(s,f))
    s=setfield(s,f,v);
end
return;
end 