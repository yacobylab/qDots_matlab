function cf = getconfig(files, chan)
% cf = getconfig(files, chan)

% (c) 2010 Hendrik Bluhm.  Please see LICENSE and COPYRIGHT information in plssetup.m.

if ischar(chan)
    chan = {chan};
end

if ischar(files)
    files = {files};
end

cv = zeros(length(files), length(chan));

for i = 1:length(files);
    load(files, 'configvals', 'configch');
    for j = 1:length(chan)
        cf(i, j) = configvals(strfind(chan{j}, configch));
    end
end
  
