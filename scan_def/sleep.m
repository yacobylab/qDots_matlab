function varargout = sleep(varargin)
%set experiment in state where it can be left alone:
%needs to be updated (can't guess the right files); 

global smdata;
global tunedata;
global fbdata;
try
  smset('PulseLine',awgseqind('all_off_LR'));
catch
  fprintf('Error setting pulseline\n');
end
evalin('base','save stdscans_2012_08_22 scanPhase2 scanPhase1 scanSDR scanRFL scanRFR scanSDL scanCalibR scanCalibL scanQPC -append');
if exist('smdata','var') && ~isempty(smdata)
  save z:\qDots\sm_config\smdata_MX50_2012_08_22 smdata
else 
  error('IDIOT')
end

if exist('tunedata','var') && ~isempty(tunedata)
  save tunedata_2012_08_22 tunedata
else
 error('Finely tuned IDIOT');
end

if exist('fbdata','var') && ~isempty(fbdata)
  save z:/qDots/sm_config/fbdata fbdata;
else
 error('feedbacky IDIOT');
end

fprintf('zzzzzzz\n');

%fprintf('why: '); why

if(length(varargin)>0)
    varargout=varargin;
else
end
end