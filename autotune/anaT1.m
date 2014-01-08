%%

% (c) 2010 Hendrik Bluhm.  Please see LICENSE and COPYRIGHT information in plssetup.m.

file = 'sm_T1DBz_R_4030';
load(file, 'scan'); 
load(file, 'data');

data1 = mean(data{1},1);
dt = 1/scan.configfn.args{3}(2); %time division


Smask = [];
Tmask = [1:120, 845:850];

singlet = data1(1:length(data1)/2);
triplet = data1(length(singlet)+1:end) - singlet;
triplet(Tmask) = [];
T = dt*(1:length(triplet));

beta0 = [min(triplet), range(triplet), 1/abs(range(triplet))*dt]; 
params = fitwrap('plinit plfit', T, triplet, beta0, @(p,x)p(1)+p(2).*exp(-1*x./p(3)));

fprintf(' Triplet T_{1} = %.2f usec \n', 1e6*params(3));
