function [o1,o2,o3,o4,o5,o6,o7,o8,o9,o10] = ldsmM(fn,i1,i2,i3,i4,i5,i6,i7,i8,i9,i10)
% Loads data files from MATLAB sm as for Labview SM
% format:  out = ldsmM(fn,i1,...,i10)
% fn - a string, the file name. 'ldsm' assumes that the extension is .dat.
% i1...i10 - optional inputs. these are the names of the channels you wish to load. 
%            useful if you don't want to read all of them.
% data - cell array of output arrays. the first of these correspond to the axis of the measurement.
% Written by Hendrik

% (c) 2010 Hendrik Bluhm.  Please see LICENSE and COPYRIGHT information in plssetup.m.



load(fn, 'data', 'scan');

% put fast ramp first
dim = length(scan.loops);

npoints = zeros(size(scan.loops));
for i = 1:length(scan.loops)
    if isempty(scan.loops(i).npoints)
        npoints(i) = length(scan.loops(i).rng);
    else
        npoints(i) = scan.loops(i).npoints;
    end
    if isempty(scan.loops(i).rng)
        scan.loops(i).rng = 1:npoints(i);
    end
end

switch dim
    case 1
        o1 = linspace(scan.loops(1).rng(1), scan.loops(1).rng(end), npoints(1));
        
    case 2
        x = linspace(scan.loops(1).rng(1), scan.loops(1).rng(end), npoints(1));
        y = linspace(scan.loops(2).rng(1), scan.loops(2).rng(end), npoints(2));
        [o1, o2] = ndgrid(x, y);
        
    case 3
        x = linspace(scan.loops(1).rng(1), scan.loops(1).rng(end), npoints(1));
        y = linspace(scan.loops(2).rng(1), scan.loops(2).rng(end), npoints(2));
        z = linspace(scan.loops(3).rng(1), scan.loops(3).rng(end), npoints(3));
        [o1, o2, o3] = ndgrid(x, y, z);

    otherwise 
        error('> 3D scans not supported yet.\n');
end

for i = 1:length(data)
    data{i} = permute(data{i}, max(2, dim):-1:1);
    eval(sprintf('o%d = data{i};', i+dim));    
end
