%%
smset('1a',-.28)
smset('2a',-.45);
%%
sLcomp.loops(2).prefn(1).args{1} = sLcomp.loops(1).rng;
sLcomp.loops(2).prefn(1).args{2} = sLcomp.loops(1).setchan;
sLcomp.loops(2).prefn(1).args{5} = cell2mat(smget('DMM'));
fprintf('qpc setpoint = %.1f mV \n', 1e3*cell2mat(smget('DMM')));
sLcomp.loops(2).prefn(1).args{6} = .008; %error

%%
%%

smset('4a',-.2);
smset('3a',-.2);

%%

sRcomp.loops(2).prefn(1).args{1} = sRcomp.loops(1).rng;
sRcomp.loops(2).prefn(1).args{5} = cell2mat(smget('DMM'));
fprintf('qpc setpoint = %.1f mV \n', 1e3*cell2mat(smget('DMM')));
sRcomp.loops(2).prefn(1).args{6} = .007; %error

%%

sRcompwide.loops(2).prefn(1).args{1} = sRcompwide.loops(1).rng;
sRcompwide.loops(2).prefn(1).args{2} = sRcompwide.loops(2).getchan;
sRcompwide.loops(2).prefn(1).args{5} = cell2mat(smget('DMM'));
fprintf('qpc setpoint = %.1f mV \n', 1e3*cell2mat(smget('DMM')));
sRcompwide.loops(2).prefn(1).args{6} = .007; %error