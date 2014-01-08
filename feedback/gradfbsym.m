function gradfbsym(file, frames)

% (c) 2010 Hendrik Bluhm.  Please see LICENSE and COPYRIGHT information in plssetup.m.


global fbdata;

load(file, 'data');

data= data{fbdata.dataind}(frames, :, :);

data2{fbdata.dataind} = nan(size(data));
for i = 1:size(data, 1)

   data2{fbdata.dataind}(i, :, :) = data(i, :, :);
   gradfeedbfn([], data2);
   pause;
    
end
