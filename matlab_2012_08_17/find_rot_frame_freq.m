function goodfreq = find_rot_frame_freq(file)

if ~exist('file','var')|| isempty(file)
    file = uigetfile('sm*.mat');
end
d= ana_avg(file,'noppt');
data = squeeze(nanmean(d.data{1}));
fitfn = @(p,x) p(1) +(p(2)*sin(2*pi*x/p(3) +p(4))).*exp(-((x+p(6))/p(5)).^1);
freq = []; ampl=[];
pp=plsinfo('xval',d.scan.data.pulsegroups(1).name,[],d.scantime);
xv = pp(2,:);
fvals = linspace(d.scan.loops(1).rng(1),d.scan.loops(1).rng(2),size(data,1));

cosfn2 = '@(y, x)y(1)+(y(2)*cos(y(4)*x) + y(3) * sin(y(4)*x)).*exp(-(x-y(5)).^2 * y(6).^2)';
cosfn5 = '@(y, x)y(1)+y(2)*cos(y(4)*x+y(3)).*exp(-(x-y(5)).^2 * y(6).^2)';
fifn.fn = @fioscill;
fifn.args = {2};

for j = 1:size(data,1)
     fifn.args={2};
     fp=fitwrap('fine woff',xv,data(j,:),fifn, cosfn5, [1 0 1 1 0 0]);
     fp = [fp(1), fp(2)*cos(fp(3)), -fp(2)*sin(fp(3)), fp(4:6)];
     pars=fitwrap('plinit plfit fine woff',xv,data(j,:),fp, cosfn2, [1 1 1 1 0 1]);
     freq(j) = pars(4); %abs(1/pars(3));
     ampl(j) = sqrt(pars(2)^2+pars(3)^2);%abs(pars(2));
end
mask = ampl>0.9 | freq <.005 | freq > 1;
ampl(mask)=0;
FOM = smooth(ampl./freq);
[m mi] = max(FOM);
goodfreq = fvals(mi);
%freq(mask) = NaN;
figure(1); clf;
subplot(1,2,1); hold on;
plot(fvals,freq);
YL = get(gca,'YLim');
plot(goodfreq*[1 1],YL,'r');
xlabel('generator frequency'); ylabel('rabi frequency');
subplot(1,2,2); hold on;
YL = get(gca,'YLim');
plot(goodfreq*[1 1],YL,'r');
plot(fvals,ampl);
xlabel('generator frequency'); ylabel('rabi amplitude');

figure(2); clf; plot(fvals,FOM); hold on;
YL = get(gca,'YLim');
plot(goodfreq*[1 1],YL,'r');
figure(3); imagesc(xv,fvals,data);

end