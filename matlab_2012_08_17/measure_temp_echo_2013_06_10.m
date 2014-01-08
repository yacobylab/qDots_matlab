tList = [0 .065:.005:.13 .14:.01:.2 0];
%pause(30*60);
for j = tList
    fprintf('\n \n');
    rsetMCSP(j);
    fprintf('setting temperature to %.2f \n',j);
    if j>1e-3
        pause(900);
    end
    if sm_setgradient
        fprintf('dBz = %.2f MHz\n',(sm_getgradient))
       d=smrun(fConfSeq2([47 28:46],struct('nloop',75,'nrep',100, 'opts','pol','datachan','DAQ2')),smnext('ramseyE_R')); 
       if any(isnan(d{1}))
          rsetMCSP(0); 
          sleep;
           break; 
       end
       fprintf('dBz = %.2f MHz\n',(sm_getgradient))
    else
        fprintf('cannot lock gradient at temperature %.2f \n',j);
    end
end

rsetMCSP(0);
sleep;

%%
for j = 49
   if sm_setgradient
      sm_getgradient
      smrun(fConfSeq2([j],struct('nloop',1,'nrep',1200, 'opts','pol','datachan','DAQ2')),smnext('dBz_phase_bayes_R'));
      sm_getgradient
      fprintf('\n \n');
   end
    
end