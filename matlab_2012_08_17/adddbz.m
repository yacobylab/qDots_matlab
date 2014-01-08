function adddbz(side)
  if ~exist('side','var')
      side='LR';
  end
  global awgdata;
  zl=plsinfo('zl',awgdata(1).pulsegroups(end).name);
  gen=1;
  dbzname=sprintf('dBz_%02d_%d_%d_%s',abs(zl(1)*1e-3),awgdata(1).pulsegroups(end).npulse(1),gen,side);
  
  try      
      awgadd(dbzname);      
  catch e
      if ~awgdata.quiet
         fprintf('making groups \n'); 
      end
      config = struct('gen',1,'side',side, 'len',abs(zl(1)*1e-3),...
          'npulse',awgdata(1).pulsegroups(end).npulse(1));
      make_dbz_groups(config);
      awgadd(dbzname);
      
  end
      if ~awgdata.quiet
          fprintf('Adding %s\n',dbzname);
      end
end