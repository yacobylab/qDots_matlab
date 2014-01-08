function awgdatafileninit(filename)

% (c) 2010 Hendrik Bluhm.  Please see LICENSE and COPYRIGHT information in plssetup.m.


if exist(filename) == 2 ||exist([filename, '.mat']) == 2
    fprintf('File %s already exists. Overwrite? y/[n]', filename);
    in = '';
    while isempty(regexp('yYnN', in))
        in = input('', 's');
        if in == 'n' || in =='N'
            return
        end
    end
end
    

pulsedata.pulsetab = [];
pulsedata.name = [];
pulsedata.clk = [];
pulsedata.tbase = [];
pulsedata.taurc = [];

pulsedata(1) = [];

save(filename, 'pulsedata');
