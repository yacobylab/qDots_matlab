function val=setMCSP(val,path)
%function val=setMCSP(val,[path])
% Set the mixing chamber set point on an oxford fridge to val.
% if path is specified, it is the path to the labview root for 
% oxford nonsense.  Defaults to c:\Our Labview\Oxford\

  e=actxserver('LabVIEW.Application');    
  if ~exist('path','var') 
      path='C:\Our Labview\Oxford\';
  end
  Kelvinox=[path 'KELVPNLS.LLB\KelvFrontPanel.vi'];
  KelvinoxVi= invoke(e,'GetVIReference',Kelvinox);
  KelvinoxVi.SetControlValue('Set T',1);
  pause(2);  % give the window time to show.
  SetT=[path 'Kelvpnls.llb\KelvPromptForT.vi'];
  SetTVi = invoke(e,'GetVIReference',SetT);
  SetTVi.SetControlValue('Set T', val);
  SetTVi.SetControlValue('OK',1);
 end