function mpause(x)
% More accurate for short pauses than Matlab

% (c) 2010 Hendrik Bluhm.  Please see LICENSE and COPYRIGHT information in plssetup.m.

if(x > 0.01)
    pause(x);
else
    tic;while(toc < x) ; end; 
end
