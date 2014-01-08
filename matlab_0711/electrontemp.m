%% Load, process, plot data for the warm and cold traces. 
%Warm, at 247 mK: 
load data0711\tune_2011_07_21\sm_line_left_716 data;
load data0711\tune_2011_07_21\sm_line_left_716 scan;

datah = mean(data{1});
epsilonval1  = (2)^.5*linspace(scan.loops(1).rng(1), scan.loops(1).rng(2), scan.loops(1).npoints);
figure; plot(epsilonval1, datah-mean(data{1}(:)), 'r'); hold; 
%figure; plot(epsilonval1, datah, 'r'); hold; 
datah2=datah-mean(data{1}(:));
%Cold, at 44 mK
load data0711\tune_2011_07_21\sm_line_left_712 data;
 load data0711\tune_2011_07_21\sm_line_left_712 scan;

datac = mean(data{1});
epsilonval2  = (2)^.5*linspace(scan.loops(1).rng(1), scan.loops(1).rng(2), scan.loops(1).npoints);
plot(epsilonval2, datac-mean(data{1}(:)), 'b'); 
datac2=datac-mean(data{1}(:));
%plot(epsilonval2, datac, 'b'); 
xlabel('detuning (V)'); ylabel('I_{QPC} (A)'); hold;
%% Find fit for warm electrons

%scaled Kboltzmann = 8.617e-5 eV/kelvin
%T = 247mK = .2465K, and T=Te, so beta0(3)=
%Don't know what the coupling is. 
 
%Initial guesses: beta0(7), the lever arm, taken from previous measurement of electron temp. 
%(discuss this)
beta0 = [.008 6e-3 2.35e4 -.9 1e-5 5e-6 .159];
fits = fitwrap(' plfit', epsilonval1, datah2, beta0, @qpcTraceFit2, logical([1 1 0 1 1 1 1])) ;

%Solve for lever arm by assuming Tc=0
lev=fits(7)
xlabel('detuning (V)'); ylabel('G_{qpc} x 10^{10}'); %text(-7e-3,-9.76, 'aplha = .154       coupling = 31\muev');
%% new trace of junction: scan was 45 degrees (slope = -1)

 TcGuess = 5e-6;
 beta0 = [.008 6e-3 8.9e4 -.9 1e-5 TcGuess lev];

fits2 = fitwrap(' plfit', epsilonval2, datac2, beta0, @qpcTraceFit2, logical([1 1 1 1 1 1 0]));
Tnew = 1/(2*8.617e-5*fits(3))
Tc = fits(6)

xlabel('detuning (V)'); ylabel('V_{RF} (AU)'); hold;
