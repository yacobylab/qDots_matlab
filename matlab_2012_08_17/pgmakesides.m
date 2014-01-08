function varargout = pgmakesides(pg,opts)
% pgs= function pgmakesides(pg,opts)
% given a pulsegroup, make all the sides.
% todo; add ability to make a combined group here.
% opts; can be stag for staggered.
% 1 output arg; return left and right together.
% 2 output args; sep. them
sides={'right','left'};
chans={[3 4],[2 1]};
if ~exist('opts')
    opts='';
end

for i=1:length(sides)
    pgs{i}=pg;
    if isfield(pg,'dict') && ~isempty(strfind(opts,'stag'))
        pgs{i}.dict=[pgs{i}.dict, ['stag' sides{i}(1)]];
    else
        pgs{i}.dict=['stag' sides{i}(1)];
    end
    if isfield(pg,'dict')
        pgs{i}.dict=[pgs{i}.dict, sides{i}];
    else
        pgs{i}.dict=sides{i};
    end
    pgs{i}.chan = chans{i};
    pgs{i}.name = [pg.name '_' upper(sides{i}(1))]; 
end
if nargout <= 1
    varargout{1} = pgs;
else
    for i = 1:length(sides)
        varargout{i}=pgs{i};
    end
end
