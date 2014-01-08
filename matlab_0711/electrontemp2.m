%% Load, process, plot data for the warm and cold traces. 
%Warm, at 247 mK: 
load data0711\sm_jnct_1053 data;
load data0711\sm_jnct_1053 scan;

datah = mean(data{1});
epsilonval1  = (2)^.5*linspace(scan.loops(1).rng(1), scan.loops(1).rng(2), scan.loops(1).npoints);
figure; plot(epsilonval1, datah-mean(data{1}(:)), 'r'); hold; 
%figure; plot(epsilonval1, datah, 'r'); hold; 
datah2=datah-mean(data{1}(:));

%Cold, fridge at 44 mK
load data0711\sm_jnct_1055 data;
load data0711\sm_jnct_1055 scan;
datac = mean(data{1});
epsilonval2  = (2)^.5*linspace(scan.loops(1).rng(1), scan.loops(1).rng(2), scan.loops(1).npoints);
plot(epsilonval2, datac-mean(data{1}(:)), 'b'); 
datac2=datac-mean(data{1}(:));
%plot(epsilonval2, datac, 'b'); 

xlabel('detuning (V)'); ylabel('I_{QPC} (A)'); hold;
%% Find fit for warm electrons
s=1; % set to 0 for left and 1 for right
%scaled Kboltzmann = 8.617e-5 eV/kelvin
%T = 231mK = .231K, and T=Te, so beta0(3)=
%Don't know what the coupling is.
 
%Initial guesses: beta0(7), the lever arm, taken from previous measurement of electron temp. 
ts=(-1)^s*range(datah2)/1.5; %slope for tanh function, given range tanh about 1.6 over center
Th=1/(8.617e-5*2*.231);
mp=find(~logical(datah2-min(datah2))); %minplace 
sl=(datah2(end)-min(datah2))/(epsilonval(end)-epsilonval(mp));


tcv=1:50;
for i=1:50

beta0 = [0 ts Th sl 1e-5 tcv(i)*1e-6 .1];
fits = fitwrap('plinit plfit', epsilonval1, datah2, beta0, @qpcTraceFit2, logical([1 1 0 1 1 0 1])) ;
xlabel('detuning (V)'); ylabel('g'); 

levv(i)=fits(7);

end
plot(tcv, levv, 'b'); 
xlabel('Coupling (uV)'); ylabel('Lever Arm'); 
%% Cold electrons
s=1; % set to 0 for left and 1 for right

ts=(-1)^s*range(datac2)/1.6; %slope for tanh function, given range tanh about 1.6 over center
mp=find(~logical(datah2-min(datah2))) %minplace 
sl=(datac2(end)-min(datac2))/(epsilonval2(end)-epsilonval2(mp))

tcv=1:50;
for i=1:50

beta0 = [0 ts 8.9e4 sl 1e-5 tcv(i)*1e-6 levv(i)];
fits2 = fitwrap('plinit plfit', epsilonval2, datac2, beta0, @qpcTraceFit2, logical([1 1 1 1 1 0 0]));
xlabel('detuning (V)'); ylabel('V_{RF} (AU)'); hold;
Tnew(i) = 1/(2*8.617e-5*fits2(3));
Tc = fits2(6);

end

plot(tcv, Tnew, 'b'); 
xlabel('Coupling (uV)'); ylabel('Measured Temperature'); 

