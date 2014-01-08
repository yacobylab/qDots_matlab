function followMCSP(fname)

if ~exist('fname','var')
    fname='MCSP.txt';
end
ct=nan;
while 1
   try
       f=fopen('MCSP.txt');
       nt=fscanf(f,'%f');
       fclose(f);
       if isempty(nt) || isnan(nt) || nt < 0 || nt > .301
           fprintf('Unlikely temperature\n');
           fprintf(' --- %g\n', nt);
       elseif nt ~= ct
           setMCSP(nt);
           ct=nt;
           fprintf('Changing temperature to %f\n', nt);
       end      
   catch err
       fprintf('Error: %s.  Retry in 30 seconds\n',err.message);
   end
   pause(30);
   ds = datestr(now); 
   fid = fopen('z:/qDots/fridge/matlab_fridge_log/MCLog_2012_10_03.txt','a');
   TMC = getIGH('M/C');
   fprintf(fid,sprintf('%s %.3f %.3f \n',ds,TMC,ct));
   fclose(fid);     
end