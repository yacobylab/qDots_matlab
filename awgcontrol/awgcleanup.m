function names= awgcleanup()
% Return a list of files in the pulse data directory w/out corresponding
% entries in the awgdata structure.  

% (c) 2010 Hendrik Bluhm.  Please see LICENSE and COPYRIGHT information in plssetup.m.


global awgdata
files=dir([awgdata.datadir '/ *.wfm']);

end
