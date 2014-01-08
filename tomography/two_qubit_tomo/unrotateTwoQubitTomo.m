function val = unrotateTwoQubitTomo(file, config)
%function val = unrotateTwoQubitTomo(file, config)
%unrotates data from the results file "file"
%config is a struct with the following fields:
%fitopts: gets passed to the twoqubittomo code 
%mfit: gets passed to the twoqubittomo code
%rotopts: gets passed to paulijunk2, decides how to rotate things, i.e.
%'STonly'
%plotopts: decides what gets plotted: 
% can have angles, fullplot, pauli and maybe more later...
%val returns a struct with the rotated data, raw data, fidelity, and angles

if ~exist('file','var') || isempty(file)
   file = uigetfile('sm*results*.mat');
end

if ~exist('config','var')
  config=struct();
end

config = def(config, 'fitopts',' ana noppt norm fid dqwd fidfit fullplot concurrence'); %default fit options for ana_twoQubitTomo_2
config = def(config, 'mfit','fid norm');% default mfit for ana_twoQubitTomo_2
config = def(config, 'opts','');%
config = def(config, 'rotopts','STonly noplot');% opts to pass to paulijunk2 (that undoes rotations)
config = def(config, 'plotopts','angles fullplot pauli'); %what should get plotted

%fit the data
fitdata2 = ana_twoQubitTomo_2(file, struct('opts', config.fitopts,'mfit', config.mfit));

%set up dephasing:
pfunc=@(p,x) [cos((x(:)+p(2))*p(3)) * p(1), ...  % <Y,I> - like coefficent
   cos((x(:)+p(2))*p(3)) * p(1), ...  % <I,Y> - like coefficent
   repmat(p(1)^2,length(x),1),   ...  % <Y,Y> - like coefficent
   sin((x(:)+p(2))*p(3)) * p(1), ... % <Y,ST> - like coefficent
   sin((x(:)+p(2))*p(3)) * p(1)]; % <ST,Y> - like coefficent

 n = 3; 
 alpha=1;
 pfunc1=@(p,x) pfunc(p,x) .* [exp(-x(:)/p(n+1)), exp(-x(:)/p(n+2)),...
     exp(-((x(:)*(1/p(n+1))).^alpha+(x(:)/p(n+2)).^alpha)),exp(-(x(:)/p(n+1)).^alpha),exp(-(x(:)/p(n+2)).^alpha) ];
 %pfunc1=@(p,x) pfunc(p,x) .* [exp(-x(:)/p(n+1)), exp(-x(:)/p(n+1)), exp(-2*x(:)/p(n+1)),exp(-x(:)/p(n+1)),exp(-x(:)/p(n+1)) ];
 %n=4;
 %pfunc2=@(p,x) pfunc1(p,x) .* [exp(-(x(:)/p(n+1)).^2),exp(-(x(:)/p(n+1)).^2), ones(length(x),1), ...
     %exp(-(x(:)/p(n+1)).^2),exp(-(x(:)/p(n+1)).^2) ];
v = [2 5 14 8 11]; 
pvecbank = [];
pvbanknodec = [];
pvbanknodec =[];

Tlen = length(fitdata2.dataexp);

%make the unrotated, dephase+entangled pauli vectors
for m = 1:Tlen
    paulivec = [0 1 0 0 1 0 0 0 0 0 0 0 0 1 0]; %initial state = Y direction in each qubit
    rr =pfunc1(fitdata2.p, fitdata2.mfdata(1).x(m));
    paulivec(v) = rr;
    pvecbank(:,end+1) = paulivec;
    paulivec = [0 1 0 0 1 0 0 0 0 0 0 0 0 1 0]; %initial state = Y direction in each qubit
    rr =pfunc1([fitdata2.p(1:end-2) Inf Inf], fitdata2.mfdata(1).x(m));
    paulivec(v) = rr;
    pvbanknodec(:,end+1) = paulivec;
end

aa = [pi/6 0.1 0.1 pi/6 0.1 0.1]'; %initial angle guess
rdata = []; rfid = [];
for m=1:Tlen
       aa(:,end+1) = paulijunk2(fitdata2.dataexp(:,m),pvecbank(:,m)', aa(:,end), config.rotopts);
       rdata(:,m) =PauliRotator(fitdata2.dataexp(:,m), (aa(:,end)'.*[-1 1 1 -1 1 1])',''); %rotated data
       rfid(m) = .25*(1+sum(abs(rdata([8 11 14],end)))); % (rotated) fidelity
end

close(999);

%decide what is X and what is Z axis
if ~isempty(strfind(config.plotopts,'ST=z'))|| ~isempty(strfind(config.plotopts,'ST=Z'))
pscramble = [3 2 1 6 5 4 12 11 10 9 8 7 15 14 13];
else
   pscramble = 1:15; 
end

rdata = rdata(pscramble,:);
fitdata2.dataexp = fitdata2.dataexp(pscramble,:);

val.rotdata = rdata;
val.rfid = rfid;
val.angles = aa(:,2:end);
val.rawdata = fitdata2.dataexp;
val.params = fitdata2.p;
val.targetstate = pvecbank;

%now lets plot stuff

fignum = 777;
c='rgbcmyk';
yoff = 1;
if strfind(config.plotopts, 'fullplot')
      figure(fignum); clf; hold on;
      fignum = fignum+1;
      %      1    2     3    4    5    6     7    8    9    10     11    12    13
      bs = {'XI','YI','ZI','IX','IY','IZ', 'XY', 'XZ','YX', 'YZ', 'ZX', 'ZY', 'XX', 'YY', 'ZZ'};
      colors=[1 0 0; 0 0 1 ; .7 0 .7];
      cind=[1 1 1 2 2 2 3 3 3 3 3 3 3 3 3 ];
      l=size(fitdata2.dataexp,1);
        for i=1:l
            area(fitdata2.mfdata(1).x, rdata(i,:)+yoff*(l-i), yoff*(l-i),'FaceColor',colors(cind(i),:));
        end
      set(gca,'YTick',yoff*(0:1:(l-1)));
      set(gca,'YTickLabel',bs(end:-1:1));
      set(gca,'YLim',[-1 (l-1)*yoff+1]);
      xlabel('T_{entangle} (\mus)');
end

if strfind(config.plotopts, 'fullplotunrot')
      figure(fignum); clf; hold on;
      fignum = fignum+1;
      %      1    2     3    4    5    6     7    8    9    10     11    12    13
      bs = {'XI','YI','ZI','IX','IY','IZ', 'XY', 'XZ','YX', 'YZ', 'ZX', 'ZY', 'XX', 'YY', 'ZZ'};
      colors=[1 0 0; 0 0 1 ; .7 0 .7];
      cind=[1 1 1 2 2 2 3 3 3 3 3 3 3 3 3 ];
      l=size(fitdata2.dataexp,1);
        for i=1:l
            area(fitdata2.mfdata(1).x, fitdata2.dataexp(i,:)+yoff*(l-i), yoff*(l-i),'FaceColor',colors(cind(i),:));
        end
      set(gca,'YTick',yoff*(0:1:(l-1)));
      set(gca,'YTickLabel',bs(end:-1:1));
      set(gca,'YLim',[-1 (l-1)*yoff+1]);
      xlabel('T_{entangle} (\mus)');
end

if strfind(config.plotopts, 'pauli') || strfind(config.plotopts, 'Pauli')
   [maxfid Tmaxfid] = max(val.rfid);
   T = [2 Tmaxfid Tmaxfid+1 2*Tmaxfid];
   for j=T
      figure(fignum); clf; hold on;
      fignum = fignum+1;
      subplot(2,1,1); hold on;
      PauliPlot(rdata(:,j),'samefig');
      title('Rotated data');
      if j==Tmaxfid
         %text(1,-.5,sprintf('Bell State Fidelity = %.2f',maxfid)); 
      end
      subplot(2,1,2); hold on;
      PauliPlot(pvecbank(:,j), 'samefig');
      title('Expected Quantum State');
       
   end
   
end

if strfind(config.plotopts, 'fid')
   figure(fignum); clf; hold on;
   fignum = fignum+1;
   p=plot(fitdata2.mfdata(1).x,rfid, 'rx');
   set(p,'MarkerSize',8);
   set(p,'LineWidth',2);
   xlabel('\tau(\mus)');
   ylabel('Bell State Fidelity');
   g=get(gca);
   ll = line(g.XLim,[.5 .5]);
   set(ll,'Color', 'g');
   %legend('Bell State Fidelity', 'Product States Forbidden');
    
end

if strfind(config.plotopts, 'angles')
    mask = 3:Tlen;
    figure(fignum); clf; hold on;
    fignum = fignum+1;

    for j = 1:3
        plot(fitdata2.mfdata(1).x(mask),(aa(j,mask)), [c(j) '.']);
    end
    legend('angle', 'theta', 'phi');
    xlabel('T_{entangle}');
    figure(fignum); clf; hold on;
    fignum = fignum+1;
    for j = 1:3
        plot(fitdata2.mfdata(1).x(mask),aa(3+j,mask), [c(j) '.']);
    end
    xlabel('T_{entangle}');
    legend('angle', 'theta', 'phi'); 
end
 
if strfind(config.plotopts, 'norms')
   wfid1 = []; wfid2 = []; stnrm = [];

    for m = 1:Tlen
       stnrm(end+1) = .25*(1+sum(pvbanknodec(:,m).*rdata(:,m)));
        wfid1(end+1) = .25*(1+sum(rdata([8 11 14],m)));
        wfid2(end+1) = .25*(1+sum(rdata([8 11 14],m).*[-1 -1 1]')); 
    end
mask = 2:Tlen;

if strcmp(config.plotopts, 'fit_w_ent')|| strcmp(config.plotopts, 'fit_W_ent')
 
end

figure(fignum); clf; hold on; fignum = fignum + 1;
p=plot(fitdata2.mfdata(1).x(mask),stnrm(mask),'rx');
%set(p,'MarkerSize',8); 
set(p,'LineWidth',2.5);
p=plot(fitdata2.mfdata(1).x(mask),wfid1(mask),'b+');
set(p,'MarkerSize',8); set(p,'LineWidth',2);
p=plot(fitdata2.mfdata(1).x(mask),wfid2(mask),'gs');
set(p,'MarkerSize',8); set(p,'LineWidth',2);
legend('State Norm', 'Bell Fidelity 1', 'Bell Fidelity 2');
xlabel('\tau (\mus)');

dummy =pfunc1(fitdata2.p,fitdata2.mfdata(1).x(mask));
%pfunc1=@(p,x) pfunc(p,x) .* [exp(-x(:)/p(n+1)), exp(-x(:)/p(n+2)),...
     %exp(-x(:)*(1/p(n+1)+1/p(n+2))),exp(-x(:)/p(n+1)),exp(-x(:)/p(n+2)) ];
normfunc = @(p,x) .25*(1+sum(pfunc1(p,x).*pfunc([p(1:end-2)],x),2));
p=plot(fitdata2.mfdata(1).x(mask), normfunc(fitdata2.p,fitdata2.mfdata(1).x(mask)),'r');
LW = 2; %line width
set(p,'LineWidth',LW);
p=plot(fitdata2.mfdata(1).x(mask),(1+sum(dummy(:,3:5),2))*.25,'b');
set(p,'LineWidth',LW);
p=plot(fitdata2.mfdata(1).x(mask),(1+(dummy(:,3:5)*[1 -1 -1]'))*.25,'g');
set(p,'LineWidth',LW);
end
    
end





% Apply a default.
function s=def(s,f,v)
if(~isfield(s,f))
  s=setfield(s,f,v);
end
return;
end
function b=isopt(config,name)
b=~isempty(strfind(config.opts,name));
return;
end

function val = fit_W_cphase(data, pguess)

%paulifunc = @(p,x)

pfunc=@(p,x) [cos((x(:)+p(2))*p(3)) * p(1), ...  % <Y,I> - like coefficent
   cos((x(:)+p(2))*p(3)) * p(1), ...  % <I,Y> - like coefficent
   repmat(p(1)^2,length(x),1),   ...  % <Y,Y> - like coefficent
   sin((x(:)+p(2))*p(3)) * p(1), ... % <Y,ST> - like coefficent
   sin((x(:)+p(2))*p(3)) * p(1)]; % <ST,Y> - like coefficent

 n = 3;    
 pfunc1=@(p,x) pfunc(p,x) .* [exp(-x(:)/p(n+1)), exp(-x(:)/p(n+2)),...
     exp(-x(:)*(1/p(n+1)+1/p(n+2))),exp(-x(:)/p(n+1)),exp(-x(:)/p(n+2)) ];
 
 fidfunc1 = @(p,x) .25*(1+pfunc1(p,x)*[0  0 1 1 1]');
 fidfunc2 = @(p,x) .25*(1+pfunc1(p,x)*[0 0 1 -1 -1]');
 normfunc1 = @(p,x) .25*(1+sum(pfunc1(p,x).*pfunc([p(1:end-2)],x),2));
 mmodel = {fidfunc1,fidfunc2,normfunc1};
 [p, chisq] = mfitwrap(data,mmodel,pguess,'plfit plinit',[0 0 1 0 0]);
 val.p = p;
 val.chisq = chisq;

end