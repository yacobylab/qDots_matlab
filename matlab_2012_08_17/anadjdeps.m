%initialize
j={};
eps={};

%% collect many j 

%[a b c]=ana_echo('',struct('opts','ramsey','rng',[4 inf],'frames',[1:23 34:50])); %Ramsey_R_3235
%[a b c]=ana_echo('',struct('opts','ramsey','rng',[4 inf],'frames',[1:50])); %Ramsey_R_3234
%[a b c]=ana_echo('',struct('opts','ramsey','rng',[4 inf],'frames',[1:41])); %Ramsey_R_3233
%[a b c]=ana_echo('',struct('opts','ramsey','rng',[4 inf],'frames',[1:30])); %Ramsey_R_3222
[a b c]=ana_echo('',struct('opts','ramsey','rng',[4 inf],'frames',[1:31])); %Ramsey_R_3214




%jsing = 1e3*b.params(:,4)./(2*pi);
%omega_dbz=b.omega_dbz; 
jsing=1e3*sqrt( abs(b.params(:,4)./(2*pi)).^2-((b.omega_dbz/((2*pi)))^2)) .* sign(b.params(:,4) - b.omega_dbz)
epsing = b.freq(1,:);

%%
j{end+1}=jsing'
eps{end+1}=epsing



%% plot many j
figure(15); clf; 
colors='rgbcmk';
color=@(x) colors(mod(x,end)+1);
for i=1:size(j,2)
    plot(eps{i}(:),j{i}(:),[color(i) '.-']); 
    hold on 
end

figure(13); clf; 
for i=1:size(j,2)
     semilogy(eps{i}(:),j{i}(:),[color(i) '.-']); 
     hold on 
end

%%
%get dj/deps 
dj={};
j2={};
mask={};
for k=2:5
    mask{k}=eps{k}>1.46;
end
mask{2}=eps{2}>1.5
mask{1}=eps{1}>1.85;
colors='rgbcmk';
color=@(x) colors(mod(x,end)+1);
figure(12); clf; 
figure(14); clf; 
for i=1:size(j,2)
    dj{i}=diff(j{i}(mask{i}));
    dj{i}=dj{i}./diff(eps{i}(mask{i}));
    j2{i}=j{i}(mask{i}); 
    figure(12) 
    plot(-dj{i}(:),j2{i}(2:end),[color(i) '.-']); 
    hold on
    figure(14)
    plot(j2{i}(2:end),-j2{i}(2:end)'./dj{i}(:),[color(i) '.-']); 
    hold on 
end


%% FUCK. IDIOT. 
figure(10); clf; params=[];
for i=1:size(jdata.j,2)    
    if i==1
        mask=eps{i}>1.85;
    elseif i==2
        mask=eps{i}>1.5;
    else
        mask=eps{i}>1.46;
    end
    jdum=j{i}(mask); 
    epsdum=eps{i}(mask);
    %fitfcn = @(p,x)p(3)+p(1).*exp(-x./p(2)+p(4)./abs(x-p(5)).^p(6));
    fitfcn=@(p,x)p(3)+p(1).*exp(-x./p(2));
    init=[1 1 min(jdum)];
    slp=-(log(jdum(8))-log(jdum(2)))./(epsdum(8)-epsdum(2));
    init(2)=1/slp; 
    slp2=exp(log(jdum(5))+slp*epsdum(5));
    init(1)=slp2; 
    params1 = fitwrap('plinit plfit',epsdum,jdum,init,fitfcn,[1 1 1]);
    %params2 = fitwrap('plinit plfit',epsdum,jdum,params1,fitfcn,[1 1 1 1 0 0]);
    params(i,:)=params1; 
    %syms x
%for i=1:size(jdata.j,2)   
    %difffunc=diff(fitfcn(params(i,:),x),x);
    djdeps=@(p,x) -p(1)./p(2).*exp(-x./p(2)); 
    dx=.01; m1=min(epsdum); m2=max(epsdum); 
    epsplot=m1:dx:m2;
    %djdeps=subs(difffunc,x,epsplot);
        figure(10)
    plot(-djdeps(params1,epsplot), fitfcn(params1,epsplot),color(i))
    hold on
    plot(-smooth(diff(jdum(:))./diff(epsdum(:)),1),(jdum(2:end)+jdum(1:end-1))./2,[color(i) '.']); 
end
