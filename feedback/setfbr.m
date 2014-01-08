function setfbr(nr,opts)
% setfbr(nrep)
% Set nreps and jumps for standard two-qubit feedback group
% nrep can be [LT+, LFB; RT+, RFB], or a 1 x 8 vector with the rep count
% for each of the following pulse combinations:
% T T O T S 0 S S
% T O T S T S 0 S
% n(1) = 0 or nr(8)=0 may not work.

% (c) 2010 Hendrik Bluhm.  Please see LICENSE and COPYRIGHT information in plssetup.m.

if ~exist('opts','var') 
   opts = '' ;
end

seq = ~isempty(strfind(opts, 'seq'));

global awgdata;
global fbdata;

fbdata.fbr = nr;

if size(nr, 1) == 2
    if seq
        nr = [0, nr(1), nr(2), 0, 0, nr(4), nr(3), 0];
    else
       ntot = max(sum(nr, 2));
       nst = [sum(nr([1 4])), sum(nr([2 3]))]-ntot;
    % reps counts for each of the eight combinations
       nr = max(0, [min(nr(:, 1)), nr(1:2)-nr([2 1])-max(nst, 0), nst, nr([4 3])-nr(3:4)-max(nst, 0), min(nr(:, 2))]);
    end
end
% see software.txt, 12/30/10
if seq
    jmp = [3, 7; 6, 2];
else
   jmp = [find(nr(1:end-1)>= 1 & nr(2:end)==0), 8;  % lines from where to jump
      find(nr(2:end)>= 1 & nr(1:end-1)==0)+1, 1]; % jump targets
end

if isnan(fbdata.plsgrp)
   error('fbdata.plsgrp = NaN. The pulsegroup is likely not in AWG memory \n'); 
end
pg.name = awgdata(1).pulsegroups(fbdata.plsgrp).name;
nps = floor(awgdata(1).pulsegroups(fbdata.plsgrp).npulse(1)/8); %# fb pulse combinations.

pg.jump = repmat(jmp, 1, nps) + reshape(repmat(0:nps-1, size(jmp, 2)*2, 1)*8, 2, size(jmp, 2) * nps);

fbdata.nrep = nr;
nr(nr==0) = 1; %avoid infs appaearing, cosmetics
pg.nrep = [repmat(nr, 1, nps), zeros(1, 4)];

plsupdate(pg);
awgadd(pg.name);
