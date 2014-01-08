function [xp yp] = at_correlate(runnumber,tsize)
% function [xp yp] = at_correlate(runnumber,tsize)
% function [xp yp] = at_correlate(scan,data)
%   guess where the triple points are.
global tunedata;
if isstruct(runnumber) % we were passed scan,data
    d.scan = runnumber;
    d.data=tsize;
else
    if ~exist('runnumber','var')
        runnumber=length(tunedata.runs);
    end
    if ~exist('tsize','var')
        tsize=11;
    end
    runname = tunedata.name; %name of current set
    if ~isempty(runname)
        runname = ['_', runname];
    end
    file = sprintf('%s/sm_chrg%s_%03i', tunedata.dir, runname, runnumber);
    d=load(file);
end
[x y] = meshgrid(linspace(d.scan.loops(1).rng(1),d.scan.loops(1).rng(2),size(d.data{1},2)), ...
                 linspace(d.scan.loops(2).rng(1),d.scan.loops(2).rng(2),size(d.data{1},1)));
             
figure(55);
clf;
c=normxcorr2(tunedata.chrg.imgr,d.data{1});
[max_c, imax] = max(c(:));
[yrpeak, xrpeak] = ind2sub(size(c),imax(1));
subplot(311);
imagesc(x(1,:),y(:,1),c);
axis xy;

c=normxcorr2(tunedata.chrg.imgl,d.data{1});
[max_c, imax] = max(c(:));
[ylpeak, xlpeak] = ind2sub(size(c),imax(1));
xlpeak=xlpeak-round(size(tunedata.chrg.imgl,2)/2);
ylpeak=ylpeak-round(size(tunedata.chrg.imgl,1)/2);
xrpeak=xrpeak-round(size(tunedata.chrg.imgr,2)/2);
yrpeak=yrpeak-round(size(tunedata.chrg.imgr,1)/2);

subplot(312);
imagesc(x(1,:),y(:,1),c);
axis xy;
subplot(313);
imagesc(d.scan.loops(1).rng,d.scan.loops(2).rng,d.data{1});
axis xy;
hold on;
plot(x(yrpeak,xrpeak),y(yrpeak,xrpeak),'kx');
plot(x(ylpeak,xlpeak),y(ylpeak,xlpeak),'k+');
xp=[x(yrpeak,xrpeak), x(ylpeak,xlpeak)];
yp=[y(yrpeak,xrpeak), y(ylpeak,xlpeak)];
end

