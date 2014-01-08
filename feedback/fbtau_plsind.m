function plsind = fbtau_plsind(taus)
  global fbdata;
% function fbtau_plsind(taus)
% given a vector of taus in ns [ 21 , 21], return an appropriate pulseind.
% nb: fbdata.tau has a list of valid taus.  If an illegal tau is requested,
%   the first element in this array is substituted.    
    taus=min(taus,max(fbdata.tau));
    taus=max(taus,min(fbdata.tau));
    for i = 1:length(taus)
        [junk taus(i)] = find(taus(i)== fbdata.tau);
    end
    taus=taus-1;
    plsind = taus(1)*8 + taus(2)*8*length(fbdata.tau);
end