function scan=csdata(scan,setpt)

scan.data.CSOffset = CS.DACOffset.get();
scan.data.CSGain = CS.DACGain.get();
scan.data.CSPreDelay=CS.PreDelay.get();

       scan.data.CSPostDelay=CS.PostDelay.get();
       scan.data.CSNSamps=CS.NSamps.get();
       scan.data.CSSampleCount=CS.SampleCount.get();
if exist('setpt','var')
    scan.data.setpt=setpt; 
end
   