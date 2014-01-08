function [good]=sm_setgradient(tgt,opts)
% function [good]=sm_setgradient(tgt,opts)
%  Lock the gradient to tgt.  tgt defaults to -35.
%  options:
%     datachan   default DAQ2
%      grpname   Group name or number for dBz measurement
%                   default: first group that starts with dBz_
% capturerange   Range in Mhz to consider locked. Default 5
%     pumptime   Time to pump when trying to lock.  Default 0.002
%    finaltime   Time to pump when testing lock.  Default 0.3
%     attempts   How many times to try.  Default 10
% pumpattempts   How many times to try pumpng.  Default 1000
%    locktests   How many times the gradient must stay fixed to call it 
%                    good.  Default 5
%       updown   fbdata.buttonpls indices for up, down pump (singlet, triplet).  Default [4 3]
%       figure   Figure to use for status display.  Default 1035

% The heart of this code is a Kalman filter.  See wikipedia page for similar notation.
  global fbdata;
  global awgdata;
  flipcount=0;
  signcount=0;
  quiet = awgdata.quiet;
  if ~exist('opts','var')
      opts=struct();
  end
  good = 0;
  awgcntrl('on start wait err');
  
  opts=def(opts,'datachan','DAQ2');  
  opts=def(opts,'grpname','');
  opts=def(opts,'pumptime',[]);
  opts=def(opts,'capturerange',5);
  opts=def(opts,'finaltime',0.3);    
  opts=def(opts,'attempts',10);
  opts=def(opts,'pumpattempts',1000);
  opts=def(opts,'locktests',5);
  opts=def(opts,'opts','');
  opts=def(opts,'gopts','nopol nodisp');
  opts=def(opts,'updown',[4 3]);  % triplet then singlet pump
  opts=def(opts,'figure',1035);
  
  if ~exist('tgt','var') || isempty(tgt)
      tgt=fbdata.gradtgt(str2num(opts.datachan(end)));
  end
  if isempty(opts.pumptime)
      opts.pumptime=fbdata.pumptime(str2num(opts.datachan(end)));
  end
  
  scan=opts.grpname;
  gopts=opts;
  gopts.opts=gopts.gopts;
  fb_ind=awgseqind(-fbdata.plsgrp);
  
  x_hist=[];
  obs_hist=[];
  pump_hist=[];
  if ~awgcntrl('ison')
      error('awg is off. fool!')
  end
  for i=1:opts.attempts
    % Lock the gradient to closest fixed point
    smset('PulseLine',fb_ind+fbdata.buttonpls(5));          % turn on feedback
    pause(opts.finaltime);                           % wait for lock.
    
    % Measure gradient.  Sign might be trustworthy here.
    [grad out] = sm_getgradient(scan,gopts); % measure gradient. We can trust this sign-wise.
    if ~isopt(gopts,'reget')
      gopts.opts=[gopts.opts ' reget'];
    end
    scan=out.scan;
        
    % 5,5 is a pump rate of 5Mhz/pump period.  Not a bad starting guess.
    x=[grad 5 5]';                       % estimated state
    % We have no idea on the pump rate to start with, assume sigma=100
    P=diag([out.grad_dev^2, 1,1]);         % estimated covariance matrix.  We start with no knowledge of pump rates.
    % The gradient fluctuates a lot, but the pump rate should be pretty stable.
    Q=diag([100 1 1]);                 % Process noise, half-assed guess.
    while 1
        if abs(x(1) - tgt) < opts.capturerange  % Are we close enough to lock?
            break;
        end
        % gradient is not in capture range.  Pump a bit
        if grad > tgt
            smset('PulseLine',fb_ind+fbdata.buttonpls(opts.updown(2)));
            F=[ 1 -1 0; 0 1 0 ; 0 0 1];
            pump_hist=[pump_hist -1];
            prate=x(2);
        else
            smset('PulseLine',fb_ind+fbdata.buttonpls(opts.updown(1)));
            F=[ 1 0 1 ; 0 1 0 ; 0 0 1];
            pump_hist=[pump_hist 1];
            prate=x(3);
        end
            
        if 1 % adaptive pump time
            ideal_time=.75*abs(grad-tgt)/prate;
            max_accel=10;
            ideal_time=max(opts.pumptime,min(max_accel*opts.pumptime,ideal_time));
        else
            ideal_time=opts.pumptime;
        end
        tic;
        mpause(ideal_time);
        smset('PulseLine',1);
        t2=toc;
        b=t2/opts.pumptime;
        F(1,2:3)=F(1,2:3)*b; % Correct process matrix for actual pump time.
        
        % Predict the new state
        x_priori = F * x ;
        P_priori = F * P * (F') + Q;
        
        % Measure the new state.
        [grad out] = sm_getgradient(scan,gopts);                
        smset('PulseLine',1);
        grad = abs(grad)*sign(x_priori(1));  % Guess on sign of gradient.
        obs_hist=[obs_hist , [grad ; out.grad_dev]];
        
        % Perform a correction step.  
        H = [1 0 0];
        y = grad - x_priori(1);   % Innovation
        S = P_priori(1,1) + out.grad_dev^2; % Innovation covariance
        K = P_priori * (H') / S;
        x = x_priori + K * y; % a-postori value
        P = (eye(size(P)) - K * H) * P_priori;
        x_hist=[x_hist x];
        
        % Check for sign error; if the apparent pump rate is negative, we probably have the wrong
        %   sign on the gradient.
        if flipcount >= 3 && ((x(2) < -0.5*sqrt(P(2,2))) || (x(3) < -0.5*sqrt(P(3,3))))  % We appear to have
                                                             % misidentified gradient
            x(1)=-x(1);
            x(2:3)=abs(x(2:3));
            if ~quiet
               fprintf('FLIP!\n');
            end
            flipcount = 0;
        else
            flipcount=flipcount+1;
        end
        
        if x(2) > 20 || x(3) > 20
           x(2:3)=[5 5];           
           P(2,2)=1e4;
           P(3,3)=1e4;
        end
        % Update the plots
        if opts.figure
           opts=plot_grad(opts,out);

           xvals=1:size(x_hist,2);
           pumpon=find(pump_hist > 0);

           if ~isfield(opts,'filterhistory')               
               subplot(223);
               cla;
               opts.filterhistory=plot(x_hist(1,:),'k-'); 
               hold on;
               opts.pumppoints=plot(xvals(pumpon),x_hist(1,pumpon),'rx');
               set(gca,'YLim',[-100,100]);               
               opts.minband=plot([xvals(1) xvals(end)],[tgt+opts.capturerange tgt+opts.capturerange],'k-');
               opts.maxband=plot([xvals(1) xvals(end)],[tgt-opts.capturerange tgt-opts.capturerange],'k-');
               opts.obs=errorbar(1:size(obs_hist,2),obs_hist(1,:),obs_hist(2,:),'b.');
               xlabel('Iteration');
               ylabel('Gradient');
               subplot(224);
               cla;
               opts.uprate=plot(x_hist(2,:),'gx-');
               hold on;
               opts.downrate=plot(x_hist(3,:),'bx-');
               xlabel('Iteration');
               ylabel('S(g)/T(b) pump rate\nMHz/pump cycle');
               set(gca,'YLim',[-2 10]);
           else
               set(opts.filterhistory,'XData',xvals,'YData',x_hist(1,:));
               set(opts.pumppoints,'XData',xvals(pumpon),'YData',x_hist(1,pumpon));
               set(opts.minband,'XData',[xvals(1) xvals(end)]);
               set(opts.maxband,'XData',[xvals(1) xvals(end)]);
               set(opts.obs,'XData',(1:size(obs_hist,2))','YData',obs_hist(1,:)','UData',obs_hist(2,:)','LData',obs_hist(2,:)');
               set(opts.uprate,'XData',xvals,'YData',x_hist(2,:));
               set(opts.downrate,'XData',xvals,'YData',x_hist(3,:));
           end
           if get(opts.figure, 'CurrentCharacter') == char(27)
               set(opts.figure, 'CurrentCharacter', char(0));
               sleep;
               return
           end
        end
    end
    
    smset('PulseLine',fbdata.buttonpls(5)+fb_ind);
    tests=0;
    while 1
       pause(opts.finaltime);
       [grad out] = sm_getgradient(scan,gopts);        
       smset('PulseLine',fbdata.buttonpls(5)+fb_ind);
       opts=plot_grad(opts,out);
       if abs(grad - tgt) > 2*opts.capturerange 
          fprintf('Test %d failed: %g\n',tests,grad)
          break;
       end
       tests = tests + 1;
       %fprintf('Test %d passed: %g\n',tests,grad)
       if tests > opts.locktests
           good = 1;
           return;
       end
    end    
  end
end

function opts=plot_grad(opts,out)
  if ~opts.figure
      return;
  end
  if ~isfield(opts,'graddata')
      figure(opts.figure);
      subplot(211);
      cla;
      opts.graddata=plot(out.xd,out.yd,'rx'); hold on;
      opts.gradfit=plot(out.xd,out.ff(out.fp,out.xd),'r-');
      xlabel('ns');
      ylabel('V_{rf}');
  else
      set(opts.graddata,'YData',out.yd);
      set(opts.gradfit,'YData',out.ff(out.fp,out.xd));
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