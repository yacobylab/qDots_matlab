% Reads He level on the ILM using activex.  ILMFrontPanel must be running
% in LabView

% (c) 2010 Hendrik Bluhm.  Please see LICENSE and COPYRIGHT information in plssetup.m.


function val=getHelium
    e=actxserver('LabVIEW.Application');    
    ILM='C:\Our Labview\Oxford\ILMFrontPanel.vi';
    ILMvi = invoke(e,'GetVIReference',ILM);
    val=str2num(ILMvi.GetControlValue('Display1'));
end
