function atrestore(ctrl, ind)


global tunedata;

if isempty(strfind(ctrl, 'ana'))
    ctrl = ['ana ', ctrl];
end

runname = tunedata.name; %name of current set
if ~isempty(runname)
    runname = ['_', runname];
end

for i = ind;
    
    if strfind(ctrl, 'gates')
        load(sprintf('%s/sm_chrg%s_%03i', tunedata.dir, runname, i), 'configvals', 'configch');
        [chfound, chind] = ismember(tunedata.gatechan, smchanlookup(configch));
        if ~all(chfound)
            error('Channel missing.');
        end
        tunedata.runs(i).gates = configvals(chind);
    end
   
    autotune(ctrl, i);
end
