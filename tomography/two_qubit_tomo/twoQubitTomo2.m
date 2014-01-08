%% load data and scan
filename=uigetfile('sm*Tomo*.mat')
s=load(filename);
rdr=anaRawUnpack(s.scan,s.data);
scantime=getscantime(s.scan,s.data);

%% histogram it
[t1lt t1l] = att1('left',scantime,'before');
[t1rt t1r] = att1('right',scantime,'before');

rd=anaRawScale(rdr,[t1l t1r],2:length(s.scan.data.pulsegroups));
for i=1:2
   rd{i}=-(rd{i}*2-1);
end

%rd  = 10 x 51 x 5000
%       10 x 51 x 5000

%% Apply sensor correction
%load file
ssL = load(uigetfile('sm*_mat*'));
ssR = load(uigetfile('sm*_mat*'));
% apply sensor correction
sensorcorrection = @(x,po) po(1)+.5*(po(2)-po(1))*(x+1);

rd{1} = sensorcorrection(rd{1},ssL.po);
rd{2} = sensorcorrection(rd{2},ssR.po);
%% unpack it nicely
%tomo order {'ST', 'UD', 'Y'}; R inner, L outer.  ST=X, Y=Y, UD = Z 
%  ST-ST, ST-UD, ST-Y
%  UD-ST, UD-UD, UD-Y
%  Y- ST, Y-UD, Y-Y
%dbz1 = mean(squeeze(rd{1}(1,:,:)),2); dbz= mean(squeeze(rd{2}(1,:,:)),2); 

%rd{1} = 2*(rd{1}-mean(dbz1))/(range(dbz1));
%rd{2} = 2*(rd{2}-mean(dbz2))/(range(dbz2));

% y-readouts are flipped by our dBz rotation.
Ymult = 1; % multiplier for ydata. use 1 if tomocal will be used, -1 if not

bs = {'XI','YI','ZI','IX','IY','IZ', 'XY', 'XZ','YX', 'YZ', 'ZX', 'ZY', 'XX', 'YY', 'ZZ'};



IX=  squeeze(nanmean(nanmean(rd{2}([2 5 8],:,:),3),1))';
IY=  Ymult*squeeze(nanmean(nanmean(rd{2}([4 7 10],:,:),3),1))';
IZ=  squeeze(nanmean(nanmean(rd{2}([3 6 9],:,:),3),1))';
XI=  squeeze(nanmean(nanmean(rd{1}([2 3 4],:,:),3),1))';
YI=  Ymult*squeeze(nanmean(nanmean(rd{1}([8 9 10],:,:),3),1))';
ZI=  squeeze(nanmean(nanmean(rd{1}([5 6 7],:,:),3),1))';
XX = nanmean(squeeze(rd{1}(2,:,:).*rd{2}(2,:,:)),2);
XY = Ymult*nanmean(squeeze(rd{1}(4,:,:).*rd{2}(4,:,:)),2);
XZ = nanmean(squeeze(rd{1}(3,:,:).*rd{2}(3,:,:)),2);
YX = Ymult*nanmean(squeeze(rd{1}(8,:,:).*rd{2}(8,:,:)),2);
YY = nanmean(squeeze(rd{1}(10,:,:).*rd{2}(10,:,:)),2);
YZ = Ymult*nanmean(squeeze(rd{1}(9,:,:).*rd{2}(9,:,:)),2);
ZX = nanmean(squeeze(rd{1}(5,:,:).*rd{2}(5,:,:)),2);
ZY = Ymult*nanmean(squeeze(rd{1}(7,:,:).*rd{2}(7,:,:)),2);
ZZ = nanmean(squeeze(rd{1}(6,:,:).*rd{2}(6,:,:)),2);

dstring=['[' sprintf('%s(dtau,:) ',bs{:}) '];'];
groups={1:3,4:6,7:length(bs)};
colors={'r','b','m'};

% plot it
dtau =25; % which index of dtau.
  
figure(11); clf; hold on;
for i=1:length(groups)
    data=eval(dstring);
    H=bar(groups{i},data(groups{i}),colors{i});
end
set(gca,'XTick',1:length(bs));
set(gca,'XTickLabel',bs);

%% rotate things with tomocal data
%load file
%ssL = load(uigetfile('sm*_mat*'));
%ssR = load(uigetfile('sm*_mat*'));
% apply sensor correction

ML = ssL.pm;
MR = ssR.pm;

DDD = [XI YI ZI IX IY IZ XX, XY,XZ, YX, YY, YZ, ZX, ZY, ZZ];
newD = [XI YI ZI IX IY IZ XX, XY,XZ, YX, YY, YZ, ZX, ZY, ZZ];

L = [ML' zeros(3,12)];
R = [zeros(3,3) MR' zeros(3,9)];
MM = kron(ML',MR');
M = [L; R; [zeros(9,6) MM]];

for i = 1:length(XX)
   newD(i,:)  = M*DDD(i,:)';
end

%% cell to export data nicely to Bell state finder
estring=['[' sprintf('%s(j) ',bs{:}) '];'];
dataexp = zeros(15,size(XX,1));
for j = 1:size(XX,1)
   dataexp(:,j) = eval(estring); 
    
end
%% Use corrected data
  dataexp=newD';

%% Use Oliver's bell state finder. 
clear bestangles;
clear beststates;
clear bestfids;
for j=1:size(dataexp,2)
    [ba,bs,bf,statefunc] = BellStateFinder(dataexp(:,j));
    bestangles(j,:) = ba;
    beststates(j,:) = bs;
    bestfids(j) = bf;
    norms(j)=(sum(dataexp(:,j).^2)+1)/4;
end
%%
figure(55)
clf;
x=1:length(norms);
plot(x,bestfids,'rx-',x,norms,'bx-');
%% rotate it using tomocal matricies, ssL and SSR are tomocal files
nrms = zeros(1,size(dataexp,2));
for j = 1:size(dataexp,2)
   dataexp(:,j) = tomoCalRot2(dataexp(:,j),ssL.pm, ssR.pm); 
   nrms(j) = (1+sum(dataexp(:,j).^2))/4;
end
figure(1); clf;
plot(nrms);
%% hack at a movie
clear out;
for dtau = 1:size(XX,1)
figure(11); clf; hold on;
data=[eval(dstring)];
for i=1:length(groups)
    data=eval(dstring);
    H=bar(groups{i},data(groups{i}),colors{i});
end
set(gca,'XTick',1:length(bs));  
set(gca,'XTickLabel',bs);
set(gca,'YLim',[-1.1 1.1]);  
fid=(sum(abs(data(end-3:end)))+1)/4;
norm=(1+sum(abs(data.^2)))/4;
enorm=sum(data(7:end).^2);
title(sprintf('dt = %d, norm=%g, enorm=%g, fid=%g',dtau,(1+sum(abs(data.^2)))/4,sum(data(7:end).^2),fid));
out(dtau,1:3)=[fid norm enorm];
pause;
end

%% hack at a horribly complicated figure
figure(666); clf;
for dtau = 1:size(XX,1)
figure(666); hold on;
subplot(8,8,dtau);  hold on;
data=[eval(dstring)];for i=1:length(groups)
    data=eval(dstring);
    H=bar(groups{i},data(groups{i}),colors{i});
   set(gca,'XTick',1:length(bs));  
   set(gca,'XTickLabel',bs);
   set(gca,'YLim',[-1.1 1.1]);
   set(gca,'XLim',[1 length(bs)]);
   %title(sprintf('dt = %d, norm=%g, enorm=%g',dtau,(1+sum(abs(data.^2)))/4,sum(data(7:end).^2)));   
end


end