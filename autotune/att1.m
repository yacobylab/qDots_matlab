function [ t1 ratio] =att1(side,scantime,opts )
%[t1 ratio] att1(side,[scantime],[opts]) Return t1 time appropriate to scantime.
%     side; left or right
% scantime; defaults to now
%     opts; 'before' t1 was measured before scan (default)
%           'after'  t1 was measured after scan.
%           'ask'    ask user.
if ~exist('opts','var')
    opts='before';
end

if ~isempty(strfind(opts,'ask'))
    [bt1 bt1r] = att1(side,scantime,'before');
    [at1 at1r] = att1(side,scantime,'after');
    if isnan(at1)
        t1=bt1;
        ratio=bt1r;
        return;
    end
    fprintf('T1 before: %g\n',bt1);
    fprintf('T1  after: %g\n',at1);
    a=input('[B]efore or [a]fter?', 's');
    if ~isempty(a) && (a=='a' || a=='A')
        t1=at1;
        ratio=at1r;
        fprintf('Using after\n');
    else
        t1=bt1;
        ratio=bt1r;
        fprintf('Using before\n');
    end
    return
end

before = isempty(strfind(opts,'after'));

if ~exist('scantime','var')
    scantime=now;
end

t1 = findt1(scantime,side,before);
    

if nargout > 1
    d=pdload(side,scantime);
    mt=d.meas(end); %end necessary for strange readout. 
    mt=mt.time(1)-(mt.time(4)+mt.time(5));
    ratio=1e-6*mt/t1;
end

end

function t1=findt1(scantime, side,before)
global tunedata;
t1=nan;
atswap(side);
    
for i=length(tunedata.runs):-1:1
    if isfield(tunedata.runs(i),'t1') && ~isempty(tunedata.runs(i).t1)         
        for j=length(tunedata.runs(i).t1):-1:1
            if isfield(tunedata.runs(i).t1(j),'t1') && isnan(tunedata.runs(i).t1(j).t1)
                continue;
            end            
            % find out the scan time.
            if isfield(tunedata.runs(i).t1(j),'time') % new style t1 measurement.
                time = tunedata.runs(i).t1(j).time;                       
            else
                scan = load(sprintf('%s/sm_t1_%s_%03i', tunedata.dir, side, i));
                time=getscantime(scan.scan,scan.data);                
            end           
            if before
              if time < scantime
                 t1 = tunedata.runs(i).t1(j).t1;                 
                 return;
              end
            else
                if time > scantime
                    t1 = tunedata.runs(i).t1(j).t1;
                else
                    return;
                end
            end
        end
    end
end
end