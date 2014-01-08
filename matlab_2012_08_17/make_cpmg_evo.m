function evol=make_cpmg_evo(evo,pipls, npls, evotime)
%function make_cpmg_evo(evo,pipls, npls, evotime)

if ~exist('evotime','var') || isempty(evotime)
    evotime = evo.time(1);
end

evo.time(1) = evotime/npls;
evol = evo;
evol.time(1) = evo.time(1)/2;
for j = 2:npls 
   evol(end+1) = pipls;
   evol(end+1) = evo;
end
evol(end+1) = pipls;
evol(end+1) = evo;
evol(end).time(1) = evo.time(1)/2;



end