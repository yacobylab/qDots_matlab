function anaoneD(files, dind, scale)

% (c) 2010 Hendrik Bluhm.  Please see LICENSE and COPYRIGHT information in plssetup.m.


if ~iscell(files)
    files = {files}; 
end

if nargin < 2
    dind = 1;
end

if nargin < 3
    scale = ones(1, length(files));
end

if size(scale, 1) == 1
    scale(2, :) = 0;
end

figure(50);
clf;
hold on;

c = 'rgbcmyk';

for i = 1:length(files);
    load(files{i}, 'scan', 'data')
    
    rng = scan.loops(1).rng;  
    if isempty(rng)
        rng = 1:length(data{dind});
    else
        rng = linspace(rng(1), rng(2), scan.loops(1).npoints);
    end
    plot(rng+ scale(2, i), data{dind}*scale(1, i), c(mod(i-1,7)+1)); 
end
