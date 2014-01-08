s=load('sm_PLine_L_0422');
x=plsinfo('xval',s.scan.data.pulsegroups(1).name);
d=s.data{1};
at=s.data{2};
t=[];
fp=[];
% slope, mean, tc, temp, ss, center
ff=@(p,x) p(1)*x+p(2) + p(5)*(x-p(6))./(sqrt((x-p(6)).^2+p(3)^2)).*tanh(sqrt((x-p(6)).^2 + p(3)^2)/(2*p(4)))
for i=1:size(d,1)
    y=d(i,:);
    if any(isnan(y))
        continue;
    end
    t(end+1)=at(i);
    figure(1);
    clf;
    plot(x,y);
    fp(end+1,:)=fitwrap('plfit fine',x,y,[-1e-3 mean(y) .012 .1 range(y) mean(x)],ff,[1 1 0 1 1 1]);
end

%% Load the temp log and assign it to the data

f=fopen('z:/qDots/fridge/matlab_fridge_log/MCLog_2011_12_19.txt');
tlog=[];
ind = 0;
while 1 
    tl = fgets(f);
  if ~ischar(tl);
      break;
  end
  if ind >1
  d=regexp(tl,'\s*','split');
  if ~isempty(d)
    tmp=sscanf(d{3},'%f');
    if ~isempty(tmp)
      tm=datenum(sprintf('%s %s',d{1},d{2}));
      tlog(end+1,:)=[tm ,1e3*tmp];
    end
  end
  end
  ind = ind+1;
end
fclose(f);
%% Assign it to data.  Column 1 of fp is now temp
for i=1:length(t)
    [mtd,ind] = min(abs(tlog(:,1)-t(i)));    
    fp(i,1)=tlog(ind,2);        
end
%%   
x=fp(:,1);
y=1e3*fp(:,4)*1e-3/8.617e-5;  % mK
figure(1);
clf;
ind = find(x < 200 & x > 1);
plot(x,y)
tf=fitwrap('plinit plfit',[0; x(ind)]',[0; y(ind)]',[1 0],@(p,x) p(1)*x+p(2),[1 0])
tf