function at_template(runnumber,tsize)
%function at_template(runnumber,tsize)
%Generate templates for automated triple point detection.
%
global tunedata;
if ~exist('runnumber','var')
    runnumber=length(tunedata.runs);
end
if ~exist('tsize','var')
    tsize=10;
end
tunedata.chrg.imgl=[];
tunedata.chrg.imgr=[];
runname = tunedata.name; %name of current set
if ~isempty(runname)
    runname = ['_', runname];
end

file = sprintf('%s/sm_chrg%s_%03i', tunedata.dir, runname, runnumber);
d=load(file);
%coeff = fit_plane(d.data{1});
%[mx,my]=meshgrid(1:size(d.data{1},2),1:size(d.data{1},1));
%d.data{1}=d.data{1}-mx*coeff(1)-my*coeff(2)-coeff(3);
figure(1);
clf;
imagesc(d.scan.loops(1).rng,d.scan.loops(2).rng,d.data{1});
axis xy;
fprintf('Please click on right, then left triple points.\n');
[xp yp] = ginput(2);
[x y] = meshgrid(linspace(d.scan.loops(1).rng(1),d.scan.loops(1).rng(2),size(d.data{1},2)), ...
                 linspace(d.scan.loops(2).rng(1),d.scan.loops(2).rng(2),size(d.data{1},1)));
  spx=find((abs(x(1,:)-xp(1)) < tsize*abs(x(1,2)-x(1,1))));
  spy=find((abs(y(:,1)-yp(1)) < tsize*abs(y(2,1)-y(1,1))))';  
  tunedata.chrg.imgr=d.data{1}(spy,spx);
  spx=find((abs(x(1,:)-xp(2)) < tsize*abs(x(1,2)-x(1,1))));
  spy=find((abs(y(:,1)-yp(2)) < tsize*abs(y(2,1)-y(1,1))))';  
  tunedata.chrg.imgl=d.data{1}(spy,spx);
  
end


function coeff=fit_plane(data)
data=data(~any(isnan(data),2),:);
[gx,gy] = gradient(data);
sm=2;
for l=1:size(gx,1)
    gx(l,:)=smooth(gx(l,:),3);
end
for l=1:size(gy,2)
    gy(:,l)=smooth(gy(:,l),3);
end
coeff(1)=median(gx(cull(gx)));
coeff(2)=median(gy(cull(gy)));
coeff(3)=mean(mean(data));
end



function se=cull(data)
m = median(data(:));
s = median(abs(data(:)-m));
se = find(abs(data(:)-m) < 2*s);
m = mean(data(se));
se = find(abs(data(:)-m) < 2*s);
end