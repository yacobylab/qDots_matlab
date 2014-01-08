function val=getIGH(ctrl)
%  function val=getIGH(ctrl)
%   sensible controls are:
%     M/C mixing chamber
%     G1  condenser pressure in mB
%     P1  condesner pressure in mB
   e=actxserver('LabVIEW.Application');    
   %IGH='C:\Oxford MX400\IGHSUBS.LLB\IGHFrontPanel.vi'; %from mx400 code
   IGH='C:\Our Labview\Oxford\IGHSUBS.LLB\IGHFrontPanel.vi';
   IGHvi = invoke(e,'GetVIReference',IGH);
   val=(IGHvi.GetControlValue(ctrl));
   switch ctrl
       case {'G1','P1'}
           val=val/100;   
   end
end