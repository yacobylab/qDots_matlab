function out = tuneqpc(dummy, rng, sweepgates, trafofn, qpcgate, tgt, err,k_p )

global smdata;
tic;
chatty = 0;
max_iters = 50;

if ~exist('k_p','var')|| isempty(k_p)
k_p = .1; % proportional gain    
end

big_step = .01; % biggest allowable qpc step

midvals = zeros(1,length(sweepgates));
if ~isempty(trafofn)
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
    if abs(k_p*(qpcval-tgt)) > big_step
        fprintf('crazy big change. taking a step of %.3f mV \n',big_step);
        sminc(qpcgate,big_step*sign(k_p*(qpcval-tgt)));
        keyboard
    elseif abs(qpcval-tgt)> abs(initval-tgt) %not making progress turn around
        sminc(qpcgate,-2.1*k_p*(qpcval-tgt));
        %keyboard
    else
      sminc(qpcgate,k_p*(qpcval-tgt));
      if chatty  
         fprintf('taking a step of %.3f mV \n',1e3*k_p*(qpcval-tgt));
      end
    end
    qpcval = getqpcval;
    count = count+1;
    if chatty
      fprintf('iteration %.0f, error = %.4f \n',count,(tgt-qpcval))
    end
end

if count > 0
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