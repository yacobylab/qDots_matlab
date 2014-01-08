function  make_dbz_groups(config)
%function  make_dbz_groups(config)
%   makes the appropriate dbz groups specified by config:
% config is a struct with the following fields:
% gen- generation, defaults to 1
% len- length of pulse, defaults to 8 (this is silly);
% scale- scale for pulse- usefull for awg7k, defaults to 1
% pulse- pulse in pulsedata to use. defaults to 21722, probably won't work
% with other pulses
% opts: defaults to nostag (don't make stag). can also have nocombine
% (won't make combined LR groups)
% sides: defaults to 'left right'
%npulse: number of pulses in the group: defaults to 64

global awgdata;
global plsdata;

if ~exist('config','var')
   config = struct(); 
end

config=def(config,'gen',1); 
config=def(config,'len',1);
config=def(config,'mtime',nan);
config=def(config,'scale',1);
config=def(config,'npulse',64);
config=def(config,'pulse',12);
config=def(config,'opts','nostag'); 
config=def(config,'sides','left right');


doleft = ~isempty(strfind(config.sides,'left'));
doright = ~isempty(strfind(config.sides,'right'));

len=config.len; pulses=config.npulse; mtime=config.mtime; 
generation=config.gen; scale=config.scale;

func=@plsdefgrp;

pg.pulses=config.pulse;

if isopt(config,'nostag')
   stagrng = 0;
else
    stagrng = (0:1);
end

for stag=stagrng

% make a 2-dot dbz group.
if stag
    stagstr='_stag';
else
    stagstr='';
end
namebase=sprintf('dBz_%02d_%d_%d%s',len,pulses,generation,stagstr);

if doleft
    pg.ctrl='loop pack';
    namepat=[namebase '_L']; pg.chan = [2 1]; pg.dict='left';
    if stag
        pg.dict = {'stagl', pg.dict};
    end
    pg.params=[len, pulses*scale+1, mtime, 0]; %pulselength, max dt, meas time, septime
    pg.varpar=(0:(pulses-1))'*scale;
    pg.name=namepat;
    func(pg);
end

if doright
    namepat=[namebase '_R']; pg.chan = [3 4]; pg.dict='right';
    if stag
        pg.dict = {'stagr', pg.dict};
    end
    pg.params=[len pulses mtime 0]; %pulselength, max dt, meas time, septime
    pg.varpar=(0:(pulses-1))';
    pg.name=namepat;
    func(pg);
end

if ~isopt(config,'nocombine')
    clear pg;
    load compmatrix;
    pg.pulses.groups={[namebase '_R'],[namebase '_L']};
    pg.matrix=eye(4);
    pg.pulseind(2,:) = [1:pulses];
    pg.pulseind(1,:) = [1:pulses];
    pg.ctrl='grp loop pack';
    pg.name=[namebase '_nc_LR'];
    func(pg);
    
    
    clear pg;
    load compmatrix;
    pg.pulses.groups={[namebase '_R'],[namebase '_L']};
    pg.matrix=compmatrix;
    pg.pulseind(2,:) = [1:pulses];
    pg.pulseind(1,:) = [1:pulses];
    pg.ctrl='grp loop pack';
    pg.name=[namebase '_LR'];
    func(pg);
end

end

end

% Apply a default.
function s=def(s,f,v)
  if(~isfield(s,f))
      s=setfield(s,f,v);
  end
return;
end

function b=isopt(config,name)
  b=~isempty(strfind(config.opts,name));
return;
end
