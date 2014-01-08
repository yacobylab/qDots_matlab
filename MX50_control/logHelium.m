% logHelium(mode)
% Handle logging of helium level in qDot fridge.  Prerequisites in
% getHelium.m must be met for this to do anything.  Only 3 weeks of data
% is kept.
% mode is one of:
%   start    - Start logging every few minutes (see logInterval)
%   stop     - Cease logging
%   delete   - Completely kill the timer
%   clear    - Clear the log
%   plot     - Plot the entire log
%   plotweek - Plot the last week of data.
%   fitweek  - Fit the last week of data.
%   fortune  - Try to predict the future
%
% returns:  The log in x-y format.

% (c) 2010 Hendrik Bluhm.  Please see LICENSE and COPYRIGHT information in plssetup.m.

function [log] = logHelium(mode)
% We store log data in the userdata for the timer;
% a 2 x n matrix with the first column "datenums" (1=1 day), the 2nd column
% percentages.
  logInterval = 15*60; %seconds
  timerTag = 'heliumLogTimer_qDot';  % must be unique, so make it weird.

  % Find the timer.  Make a new one if need be.
  tl = timerfind('Tag',timerTag);
  if(isempty(tl)) % Not timer was found
    disp('Helium timer not found..  Making a new one');
    t = timer('Tag',timerTag);
    set(t,'ExecutionMode','fixedRate');
    set(t,'Period',logInterval);
    set(t,'TimerFcn',@timerCallback);
    set(t,'StartDelay',1);
    % Run the callback once to put some data in the array.
    timerCallback(t,'');
  else
    assert(length(tl) == 1);
    t = tl(1);
  end
  
  ud = get(t,'UserData');
  switch(mode)
      case 'fortune'
          fittime=8.0/24.0; % Fit the last 8 hours
          pts=find(ud(:,1) > now()-fittime);
          xo=mean(ud(pts,1));
          fit = polyfit(ud(pts,1)-xo, ud(pts,2),1);
          disp(sprintf('The current boiloff is %g%% per day',-fit(1)));
          disp(sprintf('You have %g days to 40%%, %g days to 20%%',...
              (40-fit(2))/fit(1)-(now()-xo), (20-fit(2))/fit(1)-(now()-xo)));
          figure(9090);
          clf;
          plot(ud(pts,1),ud(pts,2),'rx',ud(pts,1),fit(1)*(ud(pts,1)-xo)+fit(2),'b-');
          datetick('x','mm/dd HH:MM');
          title('Helium Level');
          xlabel('Date');
          ylabel('Helium (%)');
    case 'fitweek'
          fittime=7; % Fit the last 8 hours
          pts=find(ud(:,1) > now()-fittime);
          xo=mean(ud(pts,1));
          fit = polyfit(ud(pts,1)-xo, ud(pts,2),1);
          disp(sprintf('The current boiloff is %g%% per day',-fit(1)));
          disp(sprintf('You have %g days to 40%%, %g days to 20%%',...
              (40-fit(2))/fit(1)-(now()-xo), (20-fit(2))/fit(1)-(now()-xo)));
          figure(9090);
          clf;
          plot(ud(pts,1),ud(pts,2),'rx',ud(pts,1),fit(1)*(ud(pts,1)-xo)+fit(2),'b-');
          datetick('x','mm/dd HH:MM');
          title('Helium Level');
          xlabel('Date');
          ylabel('Helium (%)');          
      case 'clear'
          ud = [];
          set(t,'UserData',ud);
          % Run the callback twice to put some data in the array.
          timerCallback(t,'');
          pause(1);
          timerCallback(t,'');
      case 'delete'
          delete(t);
      case 'start'
          set(t,'Period',logInterval);
          start(t);
      case 'stop'
          stop(t);
      case 'plot'
          plotLevel(ud);
      case 'plotweek'
          ind = find(ud(:,1)-now() > -7);
          plotLevel(ud(ind,:));
      otherwise
          disp(strcat('Unknown mode "', mode,'"'));
          help logHelium;          
  end  
  log = ud;
end

% Timer callback function.  Internal use only
function timerCallback(obj,data)
  ud=get(obj,'UserData');
  ud = [ud; [now(), getHelium()] ];
  dind = find(ud(:,1) < now - 21);
  if(~isempty(dind))
      ud(dind,:)=[];
  end
  set(obj,'UserData',ud);
end

% Plot levels.  Internal use.  Plot 80%/20% lines in red.
function plotLevel(ud)
  figure(9090); 
  clf;
  
  [r c] = size(ud);
  if(r < 2)
      return;
  end

  plot(ud(:,1),ud(:,2), 'b-', ...
      [ud(1,1) ud(end,1)] , [20 20], 'r-', ...
      [ud(1,1) ud(end,1)] , [80 80], 'r-');
  datetick('x','mm/dd HH');
  title('Helium Level');
  xlabel('Date');
  ylabel('Helium (%)');
end
