function dbzs = get_dbz(fdata,times,mdbz)


L = length(times);
T_samp = abs(diff(times(1:2))); F_samp = 1/T_samp;
NFFT = 2^nextpow2(8*L); %padded with zeros
nyq = abs(F_samp/2);
n_alias = floor(mdbz/nyq);
ff = (n_alias+.5)*nyq+(-1)^n_alias*(F_samp/2)*linspace(-.5,.5,NFFT/2+1);
dbzs=zeros(size(fdata,1),1);
for j = 1:size(fdata,1)
   ft = fft(fdata(j,:)-mean(fdata(j,:)),NFFT)/L;
   ft = 2*abs(ft(1:NFFT/2+1));
   [~, mi] = max(ft(2:end));
   dbzs(j) = ff(mi+1); %+1 because looking at 2:end
end

end