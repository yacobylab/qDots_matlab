%function rtd=ana_scrub(file, config)
% function ana_scrub(rtfile, config)
% Try to estimate mean sensor values based on experimental data
clear file;
clear config;

if ~exist('file','var') || isempty(file)
  [file dir]=uigetfile('sm*.mat','Scrub Set');
  file = [dir file];    
end
allsides=1;
rt=load(file);
scantime=getscantime(rt.scan,rt.data);
t1=[];
for sides=allsides
  figbase=10*sides;
  figures=[];
  % Make some color plots for switch checks
  figure(1+figbase);
  figures = unique([figures 1+figbase]);
  clf;
  d=permute(rt.data{sides},[1 3 2]);
  imagesc(reshape(d,size(d,1),size(d,2)*size(d,3)));
  figure(2+figbase);
  clf;
  d=permute(rt.data{sides},[1 3 2]);
  plot(mean(reshape(d,size(d,1),size(d,2)*size(d,3)),1));
  switch rt.scan.loops(1).getchan{sides}
      case 'DAQ1'
          side='left';
      case 'DAQ2'
          side='right';
      otherwise
          error('ugh');
  end
  
  [t1lt t1l] = att1(side,scantime,'before' );
  t1(sides)=t1l;
end
grp=find(cellfun(@(p) ~isempty(p),regexp({rt.scan.data.pulsegroups.name},'^ref_')));
[rtd rtff mv fp]=anaHistScale(rt.scan,rt.data,t1,size(rt.data{1},2));
%rtd=rt.data;
for sides=allsides
    xv(sides,:)=plsinfo('xval',rt.scan.data.pulsegroups(grp).name,length(allsides)-(sides-1));
end
% make a nice data matrix.
[ds dsi dsa]=unique(xv','rows');
avgs={};

for sides=allsides
  avgs{sides}=nan(max(dsa),max(dsa));
  data=squeeze(mean(rtd{sides}(:,grp,:),1));          
  for i=1:max(dsa)     
      i
      pts=find(dsa == i)
      d=mean(data(pts));
      fprintf('%d: %g\n',i,d);
      if size(ds,2) > 1{1}
        avgs{sides}(ds(i,1),ds(i,2))=d;
      else
        avgs{sides}(ds(i,1),ds(i,1))=d;
      end
  end  
end

if length(sides) == 2
% Overall plot
figure(1);
clf;
hold on;
for i=1:length(avgs{1}(:))  
   if i < 12
      plot(avgs{1}(i),avgs{2}(i),'rx');
   else
      plot(avgs{1}(i),avgs{2}(i),'gx');
   end
end
%
xlabel(rt.scan.loops(1).getchan{1}); ylabel(rt.scan.loops(1).getchan{2});

% Deduced level plots
figure(2);
clf;
lx(2,1:2)=avgs{1}(2,1:2);
lx(1,1:2)=2*avgs{1}(1,1:2)-avgs{1}(2,1:2);
ly(1:2,2)=avgs{2}(1:2,2);
ly(1:2,1)=2*avgs{2}(1:2,1)-avgs{2}(1:2,2);
subplot(223);
plot(lx(2,2),ly(2,2),'rx',0,0,'kx');
xlabel(rt.scan.loops(1).getchan{1}); ylabel(rt.scan.loops(1).getchan{2});
subplot(224)
plot(lx(1,2),ly(1,2),'rx',1,0,'kx');
xlabel(rt.scan.loops(1).getchan{1}); ylabel(rt.scan.loops(1).getchan{2});
subplot(221);
plot(lx(2,1),ly(2,1),'rx',0,1,'kx');
xlabel(rt.scan.loops(1).getchan{1}); ylabel(rt.scan.loops(1).getchan{2});
subplot(222);
plot(lx(1,1),ly(1,1),'rx',1,1,'kx');
xlabel(rt.scan.loops(1).getchan{1}); ylabel(rt.scan.loops(1).getchan{2});

% Useful output
fprintf('Triplet=1,Singlet=-1 Sensor %s levels are [%g,%g]\n',rt.scan.loops(1).getchan{1},(mean(lx,2)*2-1));
fprintf('Triplet=1,Singlet=-1 Sensor %s levels are [%g,%g]\n',rt.scan.loops(1).getchan{2},(mean(ly,1)*2-1));
%
s=mean(lx,2)*2-1;
fprintf('%s: Map -1,1 to %g,%g\n', rt.scan.loops(1).getchan{1},-(2+s(1) + s(2)) / (s(1)-s(2)),-(-2+s(1) + s(2)) / (s(1)-s(2)));
s=mean(ly,1)*2-1;
fprintf('%s: Map -1,1 to %g,%g\n', rt.scan.loops(1).getchan{2},-(2+s(1) + s(2)) / (s(1)-s(2)),-(-2+s(1) + s(2)) / (s(1)-s(2)));
else    
  lx=[];
  lx(1)=avgs{1}(1,1)*2-avgs{1}(2,2);
  lx(2)=avgs{1}(2,2);
  % Useful output
  fprintf('Triplet=-1,Singlet=1 Sensor %s levels are [%g,%g]\n',rt.scan.loops(1).getchan{1},-(lx*2-1));
  tl=-(lx(1)*2-1);
  sl=-(lx(2)*2-1);
  %
  s=lx*2-1;
  fprintf('%s: Map 1,-1 to %g,%g\n', rt.scan.loops(1).getchan{1},(2+s(1) + s(2)) / (s(1)-s(2)),(-2+s(1) + s(2)) / (s(1)-s(2)));    
  ref_fun=@(d) -1 + 2*(d-tl)/(sl-tl);
end
    
