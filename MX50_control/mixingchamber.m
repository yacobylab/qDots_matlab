function val=mixingchamber
   e=actxserver('LabVIEW.Application');    
   %IGH='C:\Oxford MX400\IGHSUBS.LLB\IGHFrontPanel.vi'; %from mx400 code
   IGH='C:\Our Labview\Oxford\IGHSUBS.LLB\IGHFrontPanel.vi';
   IGHvi = invoke(e,'GetVIReference',IGH);
   val=(IGHvi.GetControlValue('G1'));
end