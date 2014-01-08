function [ simdata ] = skew_axes(simdata,axes)
% Skew the axes in simdata according to axes.
% Axes is a 3x3 matrix with the new axes in rows.
% ie, [1 0 0] 
%     [0 0 1]
%     [0 1 0]
% will swap y and z.

simdata.data = simdata.gt * axes;
% the help for this function is longer than the function.
end

