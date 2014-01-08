%% load data and scan

% (c) 2010 Hendrik Bluhm.  Please see LICENSE and COPYRIGHT information in plssetup.m.

load sm_RamseyETomo_R_4254 scan;
load sm_RamseyETomo_R_4254 data;

%% 
data1 = squeeze(mean(data{1},1));
c = 'rgbcmyk';
figure(12); clf; hold on;
for i = 2:4%size(data1,1)
    plot(data1(i,:), c(i-1));
end
xlabel('d\tau (ns)');
ylabel('V_{RF}');
legend('ST', 'UD', 'Y');
%% pirate olivers code to make bloch sphere
rad = .5*range(data1(1,:));
cntrx1 = mean(data1(2,:));
cntry1 = mean(data1(3,:));
cntrz1 = mean(data1(4,:));
cntrx2 = mean(data1(5,:));
cntry2 = mean(data1(6,:));
cntrz2 = mean(data1(7,:));
[x y z] = sphere(100);

c=[255/255 139/255 29/255];
cmds={};
cmds=[cmds struct('type','disc','val',rad*[1 0 0],'color',[c .2],'color2',c)];
cmds=[cmds struct('type','disc','val',rad*[0 1 0],'color',[c .2],'color2',c)];
cmds=[cmds struct('type','disc','val',rad*[0 0 1],'color',[c .2],'color2',c)];
%cmds=[cmds struct('type','sphere','color',[c .05])];
cmds=[cmds struct('type','label','val',[1.15 0 0],'label','|S>')];
cmds=[cmds struct('type','label','val',[-1.15 0 0],'label','|T>')];
cmds=[cmds struct('type','label','val',[0 0 1.15],'label','|\uparrow\downarrow>')];
cmds=[cmds struct('type','label','val',[0 0 -1.15],'label','|\downarrow\uparrow>')];

scale = 5;
dbz=[0 0 1/32];
j1 =[1/9 0 0];
atan2(j1(3),j1(1))
j2 =[dbz(1) 0 0]+dbz;
cmds=[cmds struct('type','vector','val',dbz*scale,'size',1,'color',[0.5 0.5 0.1],'label','\delta B_z')];
cmds=[cmds struct('type','vector','val',j1*scale,'size',1,'color',[1 0 0],'label',' J')];
cmds=[cmds struct('type','vector','val',j2*scale,'size',1,'color',[0 0 1],'label',' J''')];

opts=struct('view',[.2 1 .1],'opts',' ','dt',.1);

figure(16); clf; hold on;
opengl software;
dbloch(cmds, opts);
rotate3d on;
%% plot in 3D

data2 = data1-mean(data1(1,:));
%surf(rad*x+cntrx,rad*y+cntry,rad*z+cntrz);
plot3(data2(3,:)/rad, data2(2,:)/rad, data2(4,:)/rad, 'b');
plot3(data2(6,:)/rad, data2(5,:)/rad, data2(7,:)/rad, 'g');
xlabel('ST'); ylabel('Yhat'); zlabel('Up Down');
