% function atcenter(opts)
% Magically step to the measurement point.
% opts - 'offset' to center w/ offset.
% (c) 2010 Hendrik Bluhm.  Please see LICENSE and COPYRIGHT information in plssetup.m.

function atcenter(opts)
global tunedata;
global awgdata;

if ~exist('opts','var')
    opts = '';
end

quiet = ~isempty(strfind(opts,'quiet'));
  if(isfield(tunedata,'xychan'))
      xychan=tunedata.xychan;
  else
      xychan={'X','Y'};
  end    
  cx = -(tunedata.measp(1,1)); %runs(end).chrg(3) + tunedata.offset(1,1));
  cy = -(tunedata.measp(2,1)); %runs(end).chrg(4) + tunedata.offset(2,1));
  if ~isempty(strfind(opts,'offset'))
      c1=str2num(tunedata.chrg.scan.loops(1).setchan{1}(end));
      c2=str2num(tunedata.chrg.scan.loops(2).setchan{1}(end));      
      newoffset=awgdata(1).offset;
      fprintf('Old awg offsets are [%s]\n',sprintf('%g ',newoffset));
      newoffset([c1 c2]) = newoffset([c1 c2]) + 1e3*tunedata.measp([1 2]);      
      fprintf('New awg offsets will be [%s]\n',sprintf('%g ',newoffset));
      if strcmp(input('Accept (yes/[no])? ','s'), 'yes') == 0
          fprintf('not centering \n');
          return;
      end

       for j = 1:length(awgdata)
               awgdata(j).offset(awgdata(j).chans)=newoffset(awgdata(j).chans);
       end
  else      
      if all([cx,cy] == 0)
          fprintf('No change\n');
          return;
      end
      if ~quiet
         fprintf('Step of %g mV on %s, %g mV on %s\n', cx *1e3, xychan{1},cy * 1e3,xychan{2});
      end
      if isempty(strfind(opts,'noconfirm')) && (strcmp(input('Accept (yes/[no])? ','s'), 'yes') == 0)
          fprintf('not centering \n');
          return;
      end
      if ~quiet
         fprintf('atchg(''%s'',%g)\n',xychan{1},cx);
         fprintf('atchg(''%s'',%g)\n',xychan{2},cy);
      end
      atchg(xychan{1},cx);
      atchg(xychan{2},cy);
      tunedata.measp(1:2, 1) =  [0 ; 0];
  end
end
