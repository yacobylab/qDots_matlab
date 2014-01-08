function out = loadPosAna(file, opts)

% (c) 2010 Hendrik Bluhm.  Please see LICENSE and COPYRIGHT information in plssetup.m.


load(file, 'data');
load(file, 'scan');

p = plsinfo('gd', scan.data.pulsegroups.name);
xv = p.varpar;
data1 = mean(data{1},1);

fitfn = @(p,x) p(7) + p(1)*(tanh((x-p(3))/p(2))+1)/2+ p(4)*(tanh((x-p(6))/p(5))+1)/2;

[m mi] = find(data1 ==min(data1));
mp = xv(mi);

%clever initial fit guess (well, not really clever, but it works)
ig = [range(data1), .1*range(xv),...
    xv((diff(data1) == max(diff(data1)))), range(data1),...
    -.1*range(xv),...
    xv((diff(data1) == min(diff(data1)))), min(data1)]; 

params = fitwrap('plinit plfit', xv', data1, ig, fitfn, [1 1 1 1 1 1 1]);

width = abs(params(3)-params(6));

fprintf('usable load position is %g uV wide \n', 1e3*width);

out = width;
