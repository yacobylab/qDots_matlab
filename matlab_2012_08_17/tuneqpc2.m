function out = tuneqpc2(dummy, rng, sweepgates, trafofn, qpcgate, tgt, err,k_p )
% function out = tuneqpc2(dummy, rng, sweepgates, trafofn, qpcgate, tgt, err,k_p )

% autocomensates by setting conductance to be correct at beginning of each line
% this is meant to be used as a prefn in the outter loop of a gate scan
% dummy is unused (it gets passed by smrun in the prefn)
% rng is the range of the inner loop (this function compensates at the
     % middle of the range
% sweepgates is really the setchans of the inner loop
% trafofn are the trafofn of the inner loop (if any)
% qpcgate is the gate to feedback on
% tgt is the target qpc conductance
% err is the allowable qpc error
% k_p is the proportional gain.
% this function assumes that the DAQ channel is T'ed into a dmm for easy
% readout of the qpc conductance. 
     
global smdata;
tic;
chatty = 0;
max_iters = 50;
momentum_ratio = .7; % step = ratio*error + (1-ratio)*old_error

if ~exist('k_p','var')|| isempty(k_p)
    k_p = .5; % proportional gain
end

big_step = .01; % biggest allowable qpc step

midvals = zeros(1,length(sweepgates));
if ~isempty(trafofn) % calculate what to set the gates to
    for j = 1:length(sweepgates)
        if isempty(trafofn(j))
            midvals(j) = mean(rng);
        else
            midvals(j) = trafofn(j).fn(mean(rng),trafofn(j).args);
        end
    end
else
    midvals = ones(1,length(sweepgates))*mean(rng);
end

smset(sweepgates, midvals);
pause(.01);

qpcval = getqpcval;
initval = qpcval;
if chatty
    fprintf('initial error = %.3f \n',qpcval-tgt);
end

count = 0;
while abs(qpcval-tgt) > err && count < max_iters
    if count < 2 % only taken one step
        if abs(k_p*(qpcval-tgt)) > big_step
            fprintf('crazy big change. taking a step of %.3f mV \n',big_step);
            sminc(qpcgate,big_step*sign(k_p*(qpcval-tgt)));
            keyboard
        elseif abs(qpcval-tgt)> abs(initval-tgt) %not making progress turn around
            sminc(qpcgate,-2.1*k_p*(qpcval-tgt));
            oldstep = -2.1*k_p*(qpcval-tgt);
            %keyboard
        else
            sminc(qpcgate,k_p*(qpcval-tgt));
            oldstep=k_p*(qpcval-tgt);
            if chatty
                fprintf('taking a step of %.3f mV \n',1e3*k_p*(qpcval-tgt));
            end
        end
        if chatty
            fprintf('iteration %.0f, error = %.4f \n',count,(tgt-qpcval))
        end
    else % we've already taken more than one step
        newstep = (1-momentum_ratio)*oldstep+momentum_ratio*k_p*(qpcval-tgt);
        if abs(newstep) > big_step
            fprintf('crazy big change. taking a step of %.3f mV \n',big_step);
            sminc(qpcgate,big_step*sign(newstep));
            keyboard % dangerous, so wait for user input
            oldstep = bigstep*sign(newstep);
        elseif abs(oldstep) < abs(newstep); % we're not making progress
            sminc(qpcgate,-1.3*oldstep); %turn around
            fprintf('turning around \n');
            oldstep = -1.3*oldstep;
        else
            sminc(qpcgate,newstep);
            if chatty
                fprintf('taking a step of %.3f mV \n',1e3*k_p*(qpcval-tgt));
            end
            oldstep = newstep;
        end
    end
    qpcval = getqpcval;
    count = count+1;
end

if count > 0 % write executive summary if we iterated
    fprintf('summary: time = %.2fsec, \n%.0f iterations\n total qpc change: %.2f mV \n \n',toc,count,1e3*(qpcval-initval));
end
if count < max_iters
    out =1;
else
    out = 0 ;
end

end

function val = getqpcval
% oldSR = smget('samprate');
% old_nrec =smdata.inst(sminstlookup('ATS660')).data.nrec;
% smdata.inst(sminstlookup('ATS660')).data.nrec = 0; % set nrec to 0
% val=smget(rdout);
% smdata.inst(sminstlookup('ATS660')).data.nrec = old_nrec;
val = cell2mat(smget('DMM'));
end