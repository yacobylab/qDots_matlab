function rsetMCSP(temp, fname)

if ~exist('fname','var')
    fname='MCSP.txt';
end
       f=fopen(fname,'w');
       fprintf(f,'%f\n',temp);
       fclose(f);
end