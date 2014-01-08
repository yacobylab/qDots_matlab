function [s1] = atfinechrg(opts,res,tmult)
%function [s1 s2 s3] = atfinechrg(opts,res,tmult)
% Take fine charge scans at 0,45,90,135 degrees.
% options: 
%           x,y,d,D; skip individual scans.
%           r; reverse next scan.
%           g; grow following scans by 25%.
% res: resolution of scan.
% tmult: multiply default ramp time by this.
  awgcntrl('amp');
  if ~exist('opts','var')
      opts='x';
  end
  
  global tunedata;
  side = upper(tunedata.name(1));
  scan=tunedata.chrg.scan;
  if exist('res','var') 
      if(length(res) == 1)
        scan.loops(1).npoints=res;
        scan.loops(2).npoints=res;
      else
        scan.loops(1).npoints=res(1);
        scan.loops(2).npoints=res(2);
      end        
  end
  
  if exist('tmult','var')
     scan.loops(1).ramptime=scan.loops(1).ramptime*tmult;
  end

  sc1=scan.loops(1).setchan;
  sc2=scan.loops(2).setchan;
  rx = 0;
  ry = 0;
  for i=1:length(opts)
      scan.figure=1000+i;
      switch opts(i)
          case 'g'
              scan.loops(1).rng=scan.loops(1).rng*1.25;
              scan.loops(2).rng=scan.loops(1).rng*1.25;
              continue;
          case 'r'
              rx = ~rx;
              continue;
          case 'R'
              ry = ~ry;
              continue;
          case 'x'
              ss=scan;
          case 'y'
              ss=scan;
              ss.loops(2).setchan=sc1;
              ss.loops(1).setchan=sc2;
          case 'd'
              ss=scan;
              ss.loops(1).setchan={sc1,sc2};
              ss.loops(2).setchan='count';
              ss.loops(1).trafofn{1}='@(x,c) (x(1)+x(2))/sqrt(2)';
              ss.loops(1).trafofn{2}='@(x,c) (x(1)-x(2))/sqrt(2)';
              ss.loops(1).rng=scan.loops(1).rng/sqrt(2);
              ss.loops(2).rng=scan.loops(2).rng/sqrt(2);
          case 'D'
              ss=scan;
              ss.loops(1).setchan={sc1,sc2};
              ss.loops(2).setchan='count';
              
              ss.loops(1).trafofn{1}='@(x,c) (x(1)+x(2))/sqrt(2)';
              ss.loops(1).trafofn{2}='@(x,c) (-x(1)+x(2))/sqrt(2)';
              ss.loops(1).rng=scan.loops(1).rng/sqrt(2);
              ss.loops(2).rng=scan.loops(2).rng/sqrt(2);
      end
    if ry    
        ss.loops(2).rng=ss.loops(2).rng([2 1]); 
    end
    if rx    
        ss.loops(1).rng=ss.loops(1).rng([2 1]); 
    end
    revstr='';
    if rx
        revstr=[revstr '_rx'];
    end
    if ry
        revstr=[revstr '_ry'];
    end

    s1=smrun(ss,smnext(sprintf('chrg_fn_%s%s%s',opts(i),revstr,side)));  
    if any(isnan(s1{1})); sleep; return; end
  end
  sleep;
end