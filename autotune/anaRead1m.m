function [mvis mtime vthresh t1] = anaRead1m(scan, data)
%function [mvis mtime vthresh t1] = anaRead1(scan, data)
%  Given a mp. varying scan dBzT1m, analyze T1 at each measurement point.
%  This is done by reformatting the data to fool anaRead1

% (c) 2010 Hendrik Bluhm.  Please see LICENSE and COPYRIGHT information in plssetup.m.


gd=plsinfo('gd',scan.data.pulsegroups(1).name);
dt = 1/scan.configfn.args{3}(2); %integration time bin

pulselength = plsinfo('zl', scan.data.pulsegroups(1).name);
pulselength = abs(pulselength(1));

data1 = data{1};

data2 = reshape(data1',round(1e-9*size(gd.varpar,1)*pulselength/dt),...
                       round(size(data1,1)*size(data1,2)/(1e-9*size(gd.varpar,1)*pulselength/dt)));
% gd now has all pulses in each row, reps in columns.  
nmp=size(gd.varpar,1)/2;
plen=size(data2,1)/nmp;
for i=1:nmp
  data=data2((1+plen*(i-1)):(plen*i),:);
  data={data'};
  [mvis(i) mtime(i) vthresh(i) t1(i)] = anaRead1(scan,data,79+i);  
end
figure(79+nmp+1);
clf;
mx=min(gd.varpar(:,1));
Mx=max(gd.varpar(:,1));
rx=Mx-mx;
my=min(gd.varpar(:,2));
My=max(gd.varpar(:,2));
ry=My-my;
axis([mx-rx/2 Mx+rx/2 my-ry/2 My+ry/2]);
for i=1:nmp
    ind=1+2*(i-1);
    text(gd.varpar(ind,1),gd.varpar(ind,2),sprintf('%.2f\n%.2f\n%.2f',mvis(i),t1(i),mtime(i)*1e6));
end
return;
