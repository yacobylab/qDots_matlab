function val = ana_switch(files, config)
%function val = ana_switch(files, config);
% config is a config struct describing what to do.
% leaving files blank will promt uigetfile, which can accept multiple files
%files can also be a cell of file names
% defaults to plot gradient and flattened data with figsize = 2 x 2
%defaults to not opening pptplot. passing opts 'pptplot' will prepare the
%figure for pptplot- for now, pptplot will plot all of the files under the
%same file/scan name, which, for now, defaults to the first file.
% config.opts; any of noppt, nograd, norng, nocolorbar, noflat.
% config.smoothing; how much to smooth gradient.  Default 1.

if ~exist('files','var') || isempty(files)
    files = uigetfile('sm*.mat', 'MultiSelect', 'on');
end

if ~iscell(files)
    files = {files};
end

if ~exist('config','var')
    config=struct();
end
%set up defaults for config
config = def(config,'opts','');   % Random boolean options
config = def(config,'smoothing',1); % default smoothing for gradient plot;
config = def(config,'plotsize',[2 2]);
count = 1; %counter for figure placement purposes
fignum = 769;
np=prod(config.plotsize);
for j = 1:length(files)
    load(files{j}, 'scan', 'data');
    d = data{1};
    val.d{j}=d;
    %calculate and plot 2d gradient data
    grad_data=gradnorm(filter(d,config.smoothing));
    val.gd{j}=grad_data;
    xv = linspace(scan.loops(1).rng(1), scan.loops(1).rng(2), scan.loops(1).npoints);
    yv = linspace(scan.loops(2).rng(1), scan.loops(2).rng(2), scan.loops(2).npoints);
    figure(ceil(fignum + count/np));
    if mod(count,np)==1
        clf; hold on
    end
    if ~isopt(config, 'nograd')
        sp = (mod(count,np)==0)*np+mod(count,np)*(mod(count,np)~=0);
        subplot(config.plotsize(1),config.plotsize(2),sp); hold on;
        imagesc(xv,yv,grad_data); axis image;
        t = title(files{j});
        set(t, 'Interpreter', 'none');
        count = count+1;
    if ~isopt(config,'norng')
       cd=grad_data(:);
       cd=sort(cd(~isnan(cd)));
       set(gca,'clim',cd([ceil(end*0.05) floor(end*0.95)]))
    end
    if isopt(config,'colorbar')
        colorbar;
    end
    end
    
    %calculate and plot flattened data
    coeff=fit_plane(d);
    [mx,my]=meshgrid(1:size(d,2),1:size(d,1));
    cdata=d-mx*coeff(1)-my*coeff(2)-coeff(3);
    if ~isopt(config, 'noflat')
        sp = (mod(count,np)==0)*np+mod(count,np)*(mod(count,np)~=0);
        subplot(config.plotsize(1),config.plotsize(2), sp); hold on;
        imagesc(xv,yv,cdata); axis image;
        t = title(files{j});
        set(t, 'Interpreter', 'none');
        count = count+1;
    if ~isopt(config,'norng')
       cd=cdata(:);
       cd=sort(cd(~isnan(cd)));
       set(gca,'clim',cd([ceil(end*0.05) floor(end*0.95)]))
    end
      if isopt(config,'colorbar')
        colorbar;
    end
    end
end

%set up pptplot
%the default is to only grab the first file. not ideal
if ~isopt(config,'noppt')
     figs = fignum + (1:floor(count/np));
     ppt=guidata(pptplot);     
     set(ppt.e_file,'String',files{1});     
     set(ppt.e_figures,'String',['[',sprintf('%d ',figs),']']);
     set(ppt.e_title,'String','ana_switch');
     set(ppt.e_body,'String','finding those damned switches');     
end
  
val.files = files;
end


%lots of internal functions below

function out=filter(data, sigma)
if (~exist('sigma','var'))
    sigma=3;
end
%fprintf('Sigma is %g\n',sigma);
wid=sigma;
ks=min(5,floor(sigma*3/2)*2+1);
kw=floor(ks/2);
kc=ceil(ks/2);
[x y]=meshgrid(-kw:kw,-kw:kw);
kernel=exp(-(x.*x+y.*y)/(2*wid*wid));
kernel=kernel / sum(kernel(:));
out=filter2(kernel,data);
remove=3;
out([(1:remove) (end-remove:end)],:)=nan;
out(:,[(1:remove) (end-remove:end)])=nan;
end
function out=gradnorm(data)
[gx,gy]=gradient(data);
out=sqrt(gx.^2 + gy.^2);
end

function b=isopt(config,name)
b=~isempty(strfind(config.opts,name));
return;
end
% Apply a default.
function s=def(s,f,v)
if(~isfield(s,f))
    s=setfield(s,f,v);
end
return;
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