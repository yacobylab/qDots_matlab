function atswap(setname,opts)
%  swap tunedata in/out
%    atswap()             -- list available
%    atswap('foo','save') -- save current tunedata as foo
%    atswap('bar')        -- save current tunedata, load bar

% (c) 2010 Hendrik Bluhm.  Please see LICENSE and COPYRIGHT information in plssetup.m.

global tunedata;
  nocopy={'sets','basis','basenames','dir'};
  if (~exist('setname'))
    fprintf('Available sets: (%s currently loaded)\n',tunedata.name);
    fprintf('\t%s\n',tunedata.sets.name);
  else
    if(~exist('opts','var'))
        opts='load';
    end    
    if isempty(setname)
        setname=tunedata.name;
    end
    i=find(strcmp({tunedata.sets.name},setname));        
    switch(opts)
        case 'rename'
            tunedata.name=setname;
        case 'save'           
            tunedata.name=setname;
            f=fields(rmfield(tunedata,nocopy));
            if(isempty(i))
                i=length(tunedata.sets)+1;
            end
            for l=1:length(f)
                tunedata.sets(i).(f{l})=tunedata.(f{l});
            end
        case 'load'
            atswap(tunedata.name,'save');
            if(isempty(i))
                error('Unable to load set "%s"\n',setname);
                atswap();
            end
            % Merge set data into tunedata so we don't nuke fields
            % we neither copy nor understand..
            ms=tunedata.sets(i);
            f=fieldnames(ms);
            for l=1:length(f)
               tunedata=setfield(tunedata,f{l}, getfield(ms,f{l}));
            end
            fprintf('Loaded tunedata set "%s"\n',tunedata.name);
            fprintf('Current run: %d\n',length(tunedata.runs));
    end
end
