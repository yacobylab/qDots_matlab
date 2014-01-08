function [figs pars sdata ]=ana_echo(file,config)
%function [figs pars sdata] = ana_echo(file,config) 
% config is a config struct describing what to do.  File can be blank
% ts is the total evolution times as a vector, ie. 0.05:0.05:2.  If ts
%   is nan(default), try to guess from group xvals.
% opts can be ramsey, noscale, fitdecay, phase, period, frq, amp, guessxval, even,
%  afitdecay -- fit decay only if amplitude is bigger than 5%
% odd, ramq, pdiff, nocenter, nocolor, ramt2, fitoffset (for frq decay)
% t1 for right is ~=0.05
% t1 for left is ~0
% dbz is a list of which channels have dbz data. defaults to any pulsegroup w/ dbz in the name.
% rng is points to fit. [1 10] means first 10 xvals.  [5 inf] means 5 onwared.
% grng is groups to fit.  [ .5 inf] means .5 onward.  uses group xvals
% channel is a list of channels to fit.
% nodbz -- skip dbz fit
figs=[];
pars = [];
plotnum=1;
fitdescr='';
persistent lastname;
persistent sdata_cache;
persistent fileinfo;
if ~exist('file','var') || isempty(file)
  if ~isstr(lastname)
      lastname='';
  end
  file=uigetfile('sm*.mat','ana_echo',lastname);
end

if isempty(file)
    return;
end
cache=0;
if strcmp(lastname,file) && ~isempty(sdata_cache) && ~isempty(fileinfo)
    st=dir(file);
    if st.bytes == fileinfo.bytes && st.datenum == fileinfo.datenum
        cache=1;
    end
    fileinfo=st;
end
lastname=file;

if ~exist('config','var')
    config=struct();
end

% Load the data file
if cache
  if ~isopt(config,'quiet')
    fprintf('Using cached copy of file %s\n',file);
  end
  sdata=sdata_cache;
else
  sdata=load(file);  % Load the scan, data.. This lets us auto-generate some options.
  fileinfo=dir(file);
  sdata_cache=sdata;
end
  data=sdata.data;
  scan=sdata.scan;

scantime=getscantime(scan,data);
pars.scantime = scantime;
% Parse options
config = def(config,'opts','');   % Random boolean options
config = def(config,'rng',[]);    % Range of points to fit, ie [10 inf];
config = def(config,'grng',[]);   % Group range points to fit, ie [1.1 1.5];
config = def(config,'channel',1); % Range of channels to fit. fixme to be smarter
config = def(config,'t1',nan);   % Default ratio of t1 to tmeas.  fixme to be smarter.
config = def(config,'frames',[]); % Default frames to fit.
config = def(config,'xlabel','T (\mus)'); % Default xlabel for inter-group series
config = def(config,'dxlabel','T (ns)');  % Default xlabel for intra-group series
config = def(config,'fitopts','interp');  % Interpolate xbals for fits\
config = def(config,'fb',100);    % Base number of figures to output.
config = def(config,'spsize',[2 2]); % Number of subplots for 'minor' figures.
% cellfun(@(p) ~isempty(p),regexp({scan.data.pulsegroups.name},'[dD][bB][zZ]'))
config = def(config,'dbz', find(cellfun(@(p) ~isempty(p),regexp({scan.data.pulsegroups.name},'[dD][bB][zZ]')))); % Is there a dBz group?
config = def(config,'side',[]); % Side of dot examined.
config = def(config, 'ts', nan); % Group xvals.
config = def(config,'acut',nan);  % Cutoff amplitude

if isopt(config,'ramsey') || isopt(config,'Ramsey')
    if isopt(config, 'singlegroup')
        config.opts = [config.opts 'fitdecay nocenter'];
    else
        config.opts=[config.opts 'amp per frq fitdecay nocenter ramq ramt2'];
    end
end
       
if isopt(config,'echo')
   config.opts=[config.opts 'freq amp per'];
end

if isempty(config.frames)
    config.frames=1:size(data{config.channel},1);
end
if isempty(config.side)
    switch scan.loops(1).getchan{1}
        case 'DAQ2'   
          config.side='right';
        case 'DAQ1'
          config.side='left';
        otherwise
          error('Unable to determine side');
    end
end
   
if isnan(config.t1)
    [t1t config.t1] = att1(config.side,scantime,'before');
end

notdbz=setdiff(1:length(scan.data.pulsegroups),config.dbz);

%======================
% Find the dt's for individual scan lines.  
for i=1:length(scan.data.pulsegroups)
    xvt=plsinfo('xval', scan.data.pulsegroups(i).name, [],scantime);
    dxvt=diff(xvt,[],2) ~= 0;
    [mc ind]=max(sum(dxvt,2));    
    dt(i,:)=xvt(ind,:);
    if ismember(i,notdbz)
      xv(:,i)=xvt(:);
    end
end
if isnan(config.ts)
% Guess the group xval from the params
if isopt(config,'guessxval') 
    pulseparams=[];    
    for i = notdbz
        p=plsinfo('params', scan.data.pulsegroups(i).name, [],scantime);
        if isempty(p)
            fprintf('Warining; no parameters on group %s\n',scan.data.pulsegroups(i).name);
            fprintf('Probable logging error.  Throwing out data.\n');
            p=repmat(nan,size(pulseparams,2),1);
        end
        pulseparams(:,i)=p;
    end
    % Find the varpar that changes
    [r c]=find(diff(pulseparams(:,notdbz),[],2) ~= 0);
    switch mode(c)
        case 2
            xlab = 'T $(\mu{}s)$';
        case 3
            xlab = '$\epsilon$ (mV)';
        otherwise
            xlab = '?';
    end
    ts=pulseparams(mode(c),:)';
else
    ts=[];
    [r c] = find(diff(xv,[],2) ~= 0);
    if ~isempty(r)        
      ts = xv(mode(r),:)';
    else
      fprintf('No xval changes from group to group.  Try guessxval\n');           
      ts = xv(1,:);
    end
end
else
    ts=config.ts;
    if (length(ts) ~= length(scan.data.pulsegroups))
        error('The length of TS (%d) must be the same as the number of groups (%d)\n',length(ts),length(scan.data.pulsegroups));
    end
end

ts=ts(:);
if isopt(config,'even')
  notdbz=notdbz(2:2:end);
end
if isopt(config,'odd')
  notdbz=notdbz(1:2:end);
end
% Scale the data
if ~isopt(config,'noscale')
  data_all=anaHistScale(scan,data,config.t1);
  yl='P(T)';
  offset=.3;
else
  data_all=data;
  yl='V_{rf} (mV)';
  offset=4e-4;
end

%pars.pulseparams = pulseparams;
omega_dbz=2*pi/32; % Assume this if there is no dbz reference.

for i=1:length(config.channel)
  data=data_all{config.channel(i)};
  fb=config.fb+100*(i-1);
 
  % Concatinate giant data struct for mfitwrap
for j = 1:length(notdbz) 
   mdata(j).x = dt(j,:);
   mdata(j).y = squeeze(nanmean(data(config.frames,notdbz(j),:),1))';    
end
  
  %=================
  % Plot the main data
  o=0;
  figure(fb); figs=[figs gcf];
  clf;  
  colors='rgbcmyk';
  color=@(x) colors(mod(x,end)+1);
  params=[];
  for j=1:length(notdbz) % Fit all the rest of the data.
     ind=notdbz(j);
     rdata=squeeze(nanmean(data(config.frames,ind,:),1))'+o;
     o=o+offset;
     if ~isopt(config, 'noplot')
     plot(dt(ind,:),rdata,[color(j) '.']);
     end
     [fp,ff]=fitosc(dt(ind,:),rdata,['plot' config.opts],config.rng,[color(j) '-']);
     ylabel(yl);
     xlabel(config.dxlabel);
     params(j,:)=fp;     
     hold on;
  end
  pars.params=params;
  pars.ts=ts;
  pars.dt=dt;
  if length(params) > 3
      if size(params,1) > 10
        jbar=mean(params(3:end-3,4));
      else
        jbar=mean(params(:,4));
      end
  else
      jbar=nan;
  end
  title(sprintf('|J|=%g Mhz',1e3*jbar/(2*pi)));
  fitdescr = [ fitdescr sprintf('|J|=%g Mhz\n',1e3*jbar/(2*pi)) ];
  
  if ~isopt(config,'nodbz')
      for i=1:length(config.dbz)
          dbzdata=squeeze(nanmean(data(config.frames,config.dbz(i),:),1))';
          [plotnum figs] = nextfig(config,plotnum,fb,figs);
          plot(dt(1,:),dbzdata,'b.');
          xlabel('T (ns)');
          ylabel(yl);
          
          if ~isopt(config,'nofitdbz')
              [fp,ff]=fitosc(dt(1,:),dbzdata,['fitdecay nocenter plot ' config.fitopts],[]);
              hold on;
              str=sprintf('T_2^*=%.3g ns, V=%.3f, T=%.3f, phi=%f',1./fp(6),2*sqrt(fp(2)^2+fp(3)^2),2*pi/fp(4),atan2(fp(3),fp(2))-pi);
              title(str);
              pars.dbzt2=1./fp(6);
              fitdescr = [ fitdescr sprintf('dBz_%d: ',config.dbz(i)) str sprintf('\n') ];
              omega_dbz = fp(4);
              pars.omega_dbz=omega_dbz;
          else
              title('dBz reference');
          end
      end
  end
  
  
  ts=ts(notdbz);
  % Plot various handy things.
  % amplitude vs. xval
  if isopt(config,'amp')
    [plotnum figs] = nextfig(config,plotnum,fb,figs);
    %ampfunc=@(x) 2*((abs(x(:,4)) > 0.01) .* (abs(x(:, 6)) < .2).*(sqrt(x(:, 2).^2 + x(:, 3).^2)));  
    ampfunc=@(x) 2*(sqrt(x(:, 2).^2 + x(:, 3).^2));  
    a=ampfunc(params);
    plot(ts,a,'b.');   
    if ~isnan(config.acut)       
       cut=config.acut*median(a(1:3));
       mask=(a > cut);
    else
       mask=~isnan(a);
    end
    if isopt(config, 'rmoutlier')
        mask = mask & (a <1.2);
    end
    plot(ts(mask),ampfunc(params(mask,:)),'b.');   
    if isopt(config,'linfit')
      ind = find(ts' > config.grng(1) & ts' < config.grng(2) & mask);
      fpp = fitwrap('plfit',ts(ind)', ampfunc(params(ind,:))', [-1 0], @(p,x) p(1)*x + p(2));
      str = [str sprintf('Amp %g t + %g',fpp(1),fpp(2))];
      fp=[];
      fp(1) = fpp(2); fp(2) = (fpp(2)+fpp(1)*mean(ts(ind)))/fpp(1);
    elseif isopt(config,'logfit')
      ind = find(ts' > config.grng(1) & ts' < config.grng(2) & mask(:)');
      fpp = fitwrap('plfit',ts(ind)', log(ampfunc(params(ind,:)))', [-1 0], @(p,x) p(1)*x + p(2));
      str = [str sprintf('Amp exp^(%g t + %g)',fpp(1),fpp(2))];
      fp=[];
      fp(1) = fpp(2); fp(2) = 1/fpp(1);
    else
      [fp,junk,fstr]=fitdecay(ts(mask)',ampfunc(params(mask,:))',['plot' config.opts],config.grng);
      str=['Amp' fstr];
      str= [str sprintf('Amp=%.3g T_2^{echo}=%.3g Q=%.1f T=%.3g',fp(1),fp(2),fp(2)/(1e-3*2*pi/jbar),2*pi/mean(jbar))];
    end
    title(str);
    fitdescr = [ fitdescr 'Amp: ' str sprintf('\n') ];
    ylabel(yl);
    xlabel(config.xlabel);
    pars.amp = fp(1);
    if ~isempty(config.grng)
      ind = find(ts > config.grng(1) & ts < config.grng(2));
    else
      ind=1:length(ts);
    end 
    pars.maxamp=max(ampfunc(params(ind,:)));
    pars.T2 = fp(2);
    pars.Q = fp(2)/(1e-3*2*pi/jbar);
    pars.Jbar = jbar;
    pars.T= 2*pi/mean(jbar);
    pars.afp = fp;
  end
  
  % Plot various handy things.
  % amplitude vs. xval
  if isopt(config,'per')
    [plotnum figs] = nextfig(config,plotnum,fb,figs);
    perfunc=@(params) 2*pi./params(:,4);
    plot(ts,perfunc(params),'b.');
    title('Period');
    ylabel('T (ns)');
    xlabel(config.xlabel);
  end
  
  if isopt(config,'frq')
    [plotnum figs] = nextfig(config,plotnum,fb,figs);
    freqfunc=@(params) real(sqrt((params(:,4)./(2*pi)).^2-((omega_dbz/((2*pi)))^2)));   
    plot(ts,1e3*freqfunc(params),'b.');
    %fp=fitdecay(ts',ampfunc(params)',['plot' opts]);
    %title(sprintf('Amp=%.3g T_2^{echo}=%.3g Q=%.1f T=%.3g',fp(1),fp(2),fp(2)/(1e-3*2*pi/jbar),2*pi/mean(jbar)));
    ylabel('J (MHz)');
    xlabel(config.xlabel);
    if isopt(config,'fitoffset')
      [decp decf] = fitdecay(ts',1e3*freqfunc(params)','plot fitoffset',config.grng);
    else
      [decp decf] = fitdecay(ts',1e3*freqfunc(params)','plot',config.grng);
    end
    str=sprintf('decay const = %.2d', decp(2));    
    title(str);
    fitdescr = [ fitdescr 'Freq: ' str sprintf('\n') ];
    
    % 2nd form is more useful for quickly guestimating J's.
    %fitdescr = [ fitdescr sprintf('J(eps) = %.3g * exp(-eps/%.3g) + %.3g Mhz',decp(1:3)) ];
    fitdescr = [ fitdescr sprintf('J(eps) = 100 * exp(-(eps-%.4g)/%.3g) + %.3g Mhz',log(decp(1)/100)*decp(2),decp(2),decp(3))] ;
    pars.freq = [ts';1e3*freqfunc(params)'];
  end
    if isopt(config,'ramt2')
       [plotnum figs] = nextfig(config,plotnum,fb,figs);
       plot(ts,1./params(:,6));
       xlabel(config.xlabel);
       ylabel('T_2^*, ns');
   end
  if ~isopt(config,'nocolor')
     [plotnum figs] = nextfig(config,plotnum,fb,figs);
     rdata=reshape(permute(data,[1 3 2]),size(data,1),size(data,2)*size(data,3));     
     imagesc(rdata(config.frames,:));
  end
  if isopt(config,'mean')
     [plotnum figs] = nextfig(config,plotnum,fb,figs);
     rdata=nanmean(data(config.frames,:,:),1);     
     imagesc(squeeze(rdata));
  end
  
end
%pars = params;
figs = unique(figs);

  if ~isopt(config,'noppt')
     prettyfile=regexprep(file,'(sm_)|(\.mat)','');
     indentdescr= regexprep(fitdescr,'^(.)','\t$1','lineanchors');
     ppt=guidata(pptplot);     
     set(ppt.e_file,'String',file);     
     set(ppt.e_figures,'String',['[',sprintf('%d ',figs),']']);     
     set(ppt.e_title,'String',prettyfile);
     set(ppt.e_body,'String',fitdescr);     
     clipboard('copy',sprintf('%s\n%s\n\n',['===' prettyfile],indentdescr));
     set(ppt.exported,'Value',0);
  end
return;


function [fp,ff]=fitosc(x,y,opts,rng,style)
   % initialization function
   fig=gcf;
   fifn.fn = @fioscill;
   fifn.args = {1};
   
   
   
   cosfn = '@(y, x)y(1)+y(2)*cos(y(4)*x) + y(3) * sin(y(4)*x)';
   cosfn2 = '@(y, x)y(1)+(y(2)*cos(y(4)*x) + y(3) * sin(y(4)*x)).*exp(-(x-y(5)).^2 * y(6).^2)';
   cosfn3 = '@(y,x)y(1)+(y(2)*cos(y(4)*x) + y(3) * sin(y(4)*x)).*(y(5)*x)'; %linear decay
   cosfn4 = '@(y, x)y(1)+y(2)*cos(y(4)*x) + y(3) * sin(y(4)*x)';
       
   if exist('rng','var') && ~isempty(rng)
       pts=x>rng(1) & x <rng(2);
       x=x(pts);
       y=y(pts);
   end
   
   fdata.x = x;
   fdata.y=y;
   
   %fprintf('%g\n',std(y)); %why is this here?
   if isempty(strfind(opts,'fitdecay')) || (~isempty(strfind(opts,'afitdecay')) && std(y) < 1e-2)      
     fp=mfitwrap(fdata, cosfn, fifn,'fine', [1 1 1 1 0 0]);
     ff=str2func(cosfn);    
   elseif ~isempty(strfind(opts,'nocenter'))          
     fp=fitwrap('fine',x,y,fifn, cosfn2, [1 1 1 1 0 0]);
     fp=mfitwrap(fdata, cosfn2, fp,'fine', [1 1 1 1 1 1]);
     ff=str2func(cosfn2);
   else     
     fp=mfitwrap(fdata,cosfn2,fifn,'fine', [1 1 1 1 0 1]);
     ff=str2func(cosfn2);
   end
   if ~isempty(strfind(opts,'plot'))
       figure(fig);
       hold on;
       if ~isempty(strfind(opts,'interp'))
           x=linspace(x(1,1),x(1,end),512);
       end
       if exist('style','var') && ~isempty('style')
         plot(x,ff(fp,x),style);
       else
         plot(x,ff(fp,x),'r-');
       end
   end
return;

function [fp,ff,fitstring]=fitdecay(x,y,opts,rng,style)
   % initialization function
   fig=gcf;

   if isempty(strfind(opts,'fitoffset'))
      mask=[1 1 0];
   else
      mask=[1 1 1];
   end
      
   if exist('rng','var') && ~isempty(rng)
       pts=x>rng(1) & x <rng(2);
       x=x(pts);
       y=y(pts);
   end

   if ~isempty(strfind(opts,'gauss'))
     ff=@(p,x) p(1)*exp(-(x/p(2)).^2)+p(3);
     init=[1 max(x)/3 min(y)];
     init(1)=range(y)/(ff(init,min(x)-init(3)));     
     fmt='Decay: %.3f exp(-(t/%3f)^2)+%g\n';
     fperm=[1 2 3];
   elseif ~isempty(strfind(opts,'both'))
     ff=@(p,x) p(1)*exp(-(x/p(2)).^2-x/p(4))+p(3);
     init=[1 max(x)/3 min(y) max(x)/10];
     init(1)=range(y)/(ff(init,min(x))-init(3));
     mask(4)=1;
     fperm=[1 2 3];
     fmt='Decay: %.3f exp(-t/%3f -(t/%3f)^2)+%g\n';
   elseif ~isempty(strfind(opts,'power'))
     ff=@(p,x) p(1)*exp(-(x/p(2)).^p(4))+p(3);       
     init=[1 max(x)/3 min(y) 1];
     mask(4)=1;
     init(1)=range(y)/(ff(init,min(x))-init(3));   
     fperm=[1 2 4 3];
     fmt='Decay: %.3f exp(-(t/%3f)^{%f})+%g\n';
   else       
     ff=@(p,x) p(1)*exp(-x/p(2))+p(3);       
     init=[1 max(x)/3 min(y)];
     init(1)=range(y)/(ff(init,min(x))-init(3));
     fmt='Decay: %.3f exp(-t/%3f)+%g\n';
     fperm=[1 2 3];
   end

   if isempty(strfind(opts,'fitoffset'))      
      init(3)=0;
   end

   fp = fitwrap('plinit plfit',x,y,init,ff,mask);
   fitstring=sprintf(fmt,fp(fperm));
   if ~isempty(strfind(opts,'plot'))
       figure(fig);
       hold on;
       if ~isempty(strfind(opts,'interp'))
           x=linspace(x(1,1),x(1,end),512);
       end
       if exist('style','var') && ~isempty('style')
         plot(x,ff(fp,x),style);
       else
         plot(x,ff(fp,x),'r-');
       end
   end
return;

% Apply a default.
function s=def(s,f,v)
  if(~isfield(s,f))
      s=setfield(s,f,v);
  end
return;

function b=isopt(config,name)
  b=~isempty(strfind(config.opts,name));
return;

function [plotnum figs] = nextfig(config, plotnum, fb, figs)   
   figure(1+fb+floor((plotnum-1)/prod(config.spsize))); 
   if(mod(plotnum-1,prod(config.spsize)) == 0)
       clf;
   end   
   subplot(config.spsize(1),config.spsize(2), mod(plotnum-1,prod(config.spsize))+1); 
   figs=unique([figs gcf]); 
   plotnum=plotnum+1;
return