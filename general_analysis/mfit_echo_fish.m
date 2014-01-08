function [out figures]=mfit_echo(fname,opts)
% function [out figures]=mfit_echo(fname,opts)
% opts.grps ([2 inf] default)
% opts.opts ('residuals colorplot plotdbz' default )
% opts.taurange ([-inf,inf] default)
% opts.mfitopts ('plinit plfit plotiterr err' default)
% opts.frames ; usual
% Simultaneous constrained fit on echo data.
% model: echo signal: y0_i+a*cos(w_i * x + phi_i) exp(-((x-x0_i)/t2*)^2) exp(-(tau/t2)^alpha)
%         + fish_a * cos(w_i * x + fish_phi_i) exp(-((tau/2 + x)/t2*)^2)
% p(1)=t2
% p(2)=alpha
% p(3)=a
% p(4)=fish_a
% p(5)=t2*
% p(6+5i)=w_i
% p(7+5i)=phi_i
% p(8+5i)=fish_phi_i
% p(9+5i)=y0_i
% p(10+5i)=x0_i
if ~exist('opts','var')
    opts=struct();
end
opts=def(opts,'grps',[2 inf]);
opts=def(opts,'opts','residuals colorplot plotdbz'); % nocenter also allowed
opts=def(opts,'taurange',[-inf inf]);
opts=def(opts,'mfitopts','plinit plfit plotiter err');
opts=def(opts,'frames',[]);

fignum=1024;
figures=[];
%  parameters with _i are allowed to vary from curve to curve.
if ~exist('grps','var')
    grps=[2 inf];
end

%% Load the data.
s=ana_avg(fname,struct('grps',opts.grps,'opts','noplot noppt','frames',opts.frames));
good=find(s.tv > opts.taurange(1) & s.tv < opts.taurange(2));
s.tv=s.tv(good);
s.grps=s.grps(good);

opts=def(opts,'dbz', find(cellfun(@(p) ~isempty(p),regexp({s.scan.data.pulsegroups.name},'[dD][bB][zZ]')))); % Is there a dBz group?

x=squeeze(nanmean(s.data{1}(:,s.grps,:),1));
if any(any(isnan(s.data{1})))
    fprintf('NAN data: variances will be wrong\n');
end
% denominator here gets variance of mean, only if no NANs in data
xs=squeeze(nanstd(s.data{1}(:,s.grps,:),1)).^2/(size(s.data{1},1));

%% Make a 2D color plot of the data and variances
if isopt(opts,'colorplot')
    figure(fignum);
    figures=[figures fignum];
    fignum=fignum+1;
    
    clf;
    [mx,my]=meshgrid(s.xv{opts.grps(1)},s.tv);
    subplot(121);
    pcolor(mx,my,x);
    shading flat;
    colorbar;
    subplot(122);
    pcolor(mx,my,xs);
    shading flat;
    colorbar;
end

%% Plot and fit the dbz reference.
if ~isempty(opts.dbz)
    for i=1:length(opts.dbz)
        dbzdata =squeeze(nanmean(s.data{1}(:,opts.dbz(i),:),1));
        dbzdatas=squeeze(nanstd(s.data{1}(:,opts.dbz(i),:),1)).^2/(size(s.data{1},1));

        % initialization function
        fifn.fn = @fioscill;
        fifn.args = {1};
        cosfn2 = @(y, x)y(1)+(y(2)*cos(y(4)*x) + y(3) * sin(y(4)*x)).*exp(-(x-y(5)).^2 * y(6).^2);        
        initial=fioscill(s.xv{opts.dbz(i)}, dbzdata', 1);
        [fp chisq cov]=mfitwrap(struct('x',s.xv{opts.dbz(i)},'y',dbzdata','vy',dbzdatas'),struct('fn',cosfn2),initial,[1 1 1 1 0 1]);
        out.dbz(i)=1e3*fp(4)/(2*pi);
        out.dbzt2(i)=1./fp(6);
        out.dbzs(i)=sqrt(cov(4,4))*1e3/(2*pi)
        out.dbzt2s(i)=sqrt(cov(6,6)/(fp(6)^4));
        
        if isopt(opts,'plotdbz')
            figure(fignum);
            figures=[figures fignum];
            fignum=fignum+1;
            clf;
            hold on;
            plot(s.xv{opts.dbz(i)},dbzdata,'kx');
            dbzx=linspace(min(s.xv{opts.dbz(i)}),max(s.xv{opts.dbz(i)}),512);
            plot(dbzx,cosfn2(fp,dbzx),'r-');
            str=sprintf('T_2^*=%.3g ns, V=%.3f, T=%.3f, phi=%f',1./fp(6),2*sqrt(fp(2)^2+fp(3)^2),2*pi/fp(4),atan2(fp(3),fp(2))-pi);
            title(str);            
        end
    end
end

%% Stuff the data structure and initial guess
if isopt(opts,'fishless')
    initial=[max(s.tv) 1 1 0 30];
else
  initial=[max(s.tv) 1 1 0.5 30];
end
for j=1:length(s.grps)
    mdata(j).y=x(j,:);
    mdata(j).vy=xs(j,:);
    mdata(j).x=s.xv{s.grps(j)};
    
    mmod(j).fn=@echo_fit;
    off=5*(j-1);
    mmod(j).pt=@(p) [ p(1:5) s.tv(j)  p((6:10)+off)];
    p=fioscill(mdata(j).x,mdata(j).y,1);
    
    % propagate forward frequency guess when signal is very small
    if (nanstd(mdata(j).y) < nanmean(sqrt(mdata(j).vy))) && j > 1
        p(4)=initial(end-4);
        fprintf('Propagating omega forward on group %d\n',j);
    end
    initial=[initial, p(4), atan2(p(2),p(3)), 0, mean(mdata(j).y), 0];
end

% Median filter the initial frequency guess.
for i=4:5:length(initial)
    initial(i)=median([initial(max(4,i-5)), initial(i), initial(min(end,i+5))]);
end

out.s=s;
out.mdata=mdata;
out.mmod=mmod;
out.initial=initial;
out.p=initial;
for j=1:length(s.grps)
    mask=initial*0;
    mask(7+5*(j-1))=1;
    if ~isopt(opts,'fishless')
      mask(8+5*(j-1))=1;
    end
    if isopt(opts,'nocenter')
      mask(10:5:end)=0;
    end
    out.p=mfitwrap(mdata(j),mmod(j),out.p,opts.mfitopts,mask);
%    if j==2
%        mask=mask*0;
%        mask(5)=1;
%        out.p=mfitwrap(mdata(1:2),mmod(1:2),out.p,opts.mfitopts,mask);
%    end
end
mask=initial*0+1;
mask(2)=0;
if isopt(opts,'fishless')
    mask(8:5:end)=0;
    mask(4)=0;
end
if isopt(opts,'nocenter')
    mask(10:5:end)=0;
end

[out.p out.chisq, out.cov]=mfitwrap(mdata,mmod,out.p,opts.mfitopts,mask);
mask=initial*0+1;
if isopt(opts,'fishless')
    mask(8:5:end)=0;
    mask(4)=0;
end
if isopt(opts,'nocenter')
    mask(10:5:end)=0;
end

[out.p out.chisq, out.cov]=mfitwrap(mdata,mmod,out.p,opts.mfitopts,mask);


%as=sqrt(diag(out.chisq));
%mj=1e3*sum(out.p(:,6:5:end)./((as:,6:5:end).^2),2)./sum(1./as(:,6:5:end).^2,2)/(2*pi);
%%mj=squeeze(1e3*aa(:,4+10)/(2*pi));
%mjs=(1e3/(2*pi))^2*size(aa(:,6:5:end),2)./sum(1./as(:,6:5:end).^2,2);
%mjo=mj;
%mj=sqrt(mj.^2-([a.nu_dbz].^2)')

if isopt(opts,'residuals')
    % Make a pretty plot of fit qualities
    offset=0.5;
    ac=autocolor;
    figure(fignum);
    figures=[figures fignum];
    fignum=fignum+1;
    clf;
    hold on;
    for j=1:length(out.mdata)
        subplot(121);
        hold on;
        c=ac();
        if isopt(opts,'noerrorbars')
          h=plot(out.mdata(j).x,out.mdata(j).y+offset*(j-1),'rx','Color',c);
        else            
          h=errorbar(out.mdata(j).x,out.mdata(j).y+offset*(j-1), sqrt(out.mdata(j).vy),'LineStyle','None','Color',c);
        end
        x=linspace(min(out.mdata(j).x),max(out.mdata(j).x),512);
        y=out.mmod(j).fn(out.mmod(j).pt(out.p),x);
        plot(x,y+offset*(j-1),[c '-']);
        x=out.mdata(j).x;
        y=out.mmod(j).fn(out.mmod(j).pt(out.p),x);
        %plot(x,y+offset*(j-1),[c '.']);
        set(gca,'XLim',[min(x),max(x)])
        subplot(122);
        hold on;
        em=.2;
        plot(x,offset*em*(j-1)+x*0,'k:');
        errorbar(x,out.mdata(j).y-y+offset*em*(j-1), sqrt(out.mdata(j).vy),'Color',c);
        title(sprintf('%s: \\chi^2=%3f',fname,out.chisq));
        set(gca,'XLim',[min(x),max(x)])
    end
    subplot(121);
    x0=out.p(8:5:end);
    plot(x0,((1:length(x0))-1)*offset,'k-');
end
end
% model: echo signal: y0_i+a*cos(w_i * x + phi_i) exp(-((x-x0_i)/t2*)^2) exp(-(tau/t2)^alpha)
%         + fish_a * cos(w_i * x + fish_phi_i) exp(-((tau/2 + x)/t2*)^2)
% p(1)=t2
% p(2)=alpha
% p(3)=a
% p(4)=fish_a
% p(5)=t2*
% p(6)=tau (fake)

% p(7)=w_i
% p(8)=phi_i
% p(9)=fish_phi_i
% p(10)=y0_i
% p(11)=x0_i

function y=echo_fit(p,x)
y=p(10)+p(3)*cos(p(7)*x+p(8)) .* exp(-abs((x-p(11))/p(5)).^2 - abs(p(6)/p(1))^p(2)) + ...
       p(4) * cos(p(7)*x+p(9)) .* exp(-((1e3*p(6)/2 + x)/p(5)).^2);
end


% Apply a default.
function s=def(s,f,v)
if(~isfield(s,f))  
      s=setfield(s,f,v);    
end
end
function b=isopt(config,name)
b=~isempty(strfind(config.opts,name));
end