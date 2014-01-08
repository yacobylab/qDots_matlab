function [grad out]=sm_getgradient(grpname,opts)
% function [grad out]=sm_getgradient(grpname,opts)
%  grpname; name/number of group to use.
%     empty/does not exist: use first group whose name starts with dBz_
%  opts.datachan: default to DAQ2.
%  opts.opts: nopol, nodisp, reget
%  opts.ampthresh: 0.1; If visiblity is less than this, assume fit is bad and
%     make grad_dev huge.
%  opts.chisqthresh: 10; If chi^2 is more than this, assume fit is bad and
%     make grad, grad_dev hughe.

%
%     reget assumes the most recent previous scan to run was sm_getgradient,
%     and skips re-initializing the DAQ card.
% Return the current magnetic field gradient. 
% Sign of the gradient is only trustworthy if the gradient is locked.
% Negative means triplet-side.
% out contains 
%    norm_sig: normalized signal; 1 indicated 100% visibility
%        scan: scan used to take the data
%    grad_dev: estimate of error bar on gradient measurement
%       xd,yd: raw data
%       ff,fp: fit funtion, fit parameters.



% fill in default options, make the scan
  global fbdata;
  global awgdata;
  global smdata;
  if ~exist('opts','var')
      opts=struct();
  end  
  opts=def(opts,'datachan',{'DAQ2'}); 
  opts=def(opts,'chisqthresh',5);
  opts=def(opts,'ampthresh',.2); 
  opts=def(opts,'fitwrapopts','');
  opts=def(opts,'opts','');
  opts=def(opts,'nloop',200);
  if ~exist('grpname','var') || isempty(grpname)
     dbzgrps=find(~cellfun('isempty',regexp({awgdata(1).pulsegroups.name},'^dBz_')));
     grpname=dbzgrps(1);  % Default to first DBZ group.
  end
  
  if exist('fConfSeq2_v4')==2 && ~smdata.inst(sminstlookup('ATS660')).data.extclk
      scanfunc = @fConfSeq2_v4; 
  else
     scanfunc = @fConfSeq2; 
  end
  if ~iscell(opts.datachan)
      opts.datachan={opts.datachan};
  end
  
  if isstruct(grpname)
      scan=grpname;
      grpname=scan.data.pulsegroups(1).name;
  else
      if isopt(opts,'nopol')
        scan=scanfunc(grpname,struct('nloop',opts.nloop,'nrep',1,'opts','raw ampok','datachan',opts.datachan));
      else
        scan=scanfunc(grpname,struct('nloop',opts.nloop,'nrep',1,'opts','pol raw ampok','datachan',opts.datachan));
      end
      scan.datachan=opts.datachan; % Protect cell arrays.
  end
  if isopt(opts,'nodisp') 
      scan.disp=[];
  else
      for i=1:length(scan.datachan)
        scan.disp=struct('loop',1,'channel',i,'dim',1);
      end
  end
      
  scan.configch=[];
  
  % Take the data  
  if isopt(opts,'reget')
    smatrigfn(1,smchaninst(scan.loops(1).getchan{1}),4);
    for i=2:length(scan.loops(1).prefn)
      if ischar(scan.loops(1).prefn(i).fn)
        f=str2func(scan.loops(1).prefn(i).fn);
        f(1,scan.loops(1).prefn(i).args{:});          
      else
        scan.loops(1).prefn(i).fn(1,scan.loops(1).prefn(i).args{:});
      end
    end
    d=smget(scan.loops(1).getchan{1});
    smset('PulseLine',1);
    d{3}=d{1};
    d{4}=d{1};
    d{1}=mean(reshape(d{1},scan.loops(1).procfn(1).fn(3).args{1}),2)';
    d{4}=histc(d{4},scan.loops(1).procfn(4).fn.args{1})';
  else
    d=smrun(scan);
  end
  if any(isnan(d{1}))
      grad=nan;
      return;
  end 

  % Work out x,y, sigma_y
  %data_std=std(reshape(d{3},length(d{1}),opts.nloop),0,2)/sqrt(opts.nloop);
  %changed to this line 2013/04/08
  data_std=std(reshape(d{3},length(d{1}),opts.nloop),0,2)/sqrt(opts.nloop);
  
  % Cache the xvals for faster repeated scans.
  if isfield(scan.data,'xv')
      xv=scan.data.xv;
  else      
      xv=plsinfo('xval',grpname);
      scan.data.xv = xv;
  end
  out.scan=scan;
  out.xd=xv;
  out.yd=squeeze(d{1});
  % Hack alert!  This will need to get fixed in a big way for multiple dots.
  fbdata.refval(str2num(opts.datachan{1}(end))) = out.yd(1);
  
  % Initial guess for fit
  fp=fioscill(xv,squeeze(d{1}),1);
  fp(5:6)=0;
  cosfn2 = '@(y, x)y(1)+(y(2)*cos(y(4)*x) + y(3) * sin(y(4)*x)).*exp(-(x-y(5)).^2 * y(6).^2)';
  
  % Actual fit
  [fp,chisq,cov]=mfitwrap(struct('x',xv,'y',d{1},'vy',(data_std.^2)'),struct('fn',str2func(cosfn2)),fp,'',[1 1 1 1 0 0]);

  % Work out the visibility of the oscillation
  vc=scan.loops(1).procfn(end).fn.args{1};
  vc=vc+(vc(2)-vc(1))/2;
  n=sum(d{end});
  xbar = sum(d{end} .* vc)/n;
  xxbar = sum(d{end} .* vc .* vc)/n;
  sigstd=sqrt(xxbar-xbar*xbar);
  out.norm_sig=sqrt(fp(2)^2+fp(3)^2)/sigstd;

  % Use derivative of fit function to guess singlet or triplet side 
  dcosfn = @(y,x) -y(2)*sin(y(4)*x) + y(3) * cos(y(4)*x);
  fbtime=fbdata.taus(2);  
  grad = 1e3*(fp(4)/(2*pi))*sign(dcosfn(fp,fbtime));  
  out.grad_dev = 1e3*sqrt(cov(4,4))/(2*pi);
  out.fp=fp;
  out.ff=str2func(cosfn2);
  
  if ~isnan(opts.ampthresh) && out.norm_sig < opts.ampthresh
     out.grad_dev=1000; 
     grad=0;
  end
  if ~isnan(opts.chisqthresh) && chisq > opts.chisqthresh
      grad=0;
      out.grad_dev=1000;
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