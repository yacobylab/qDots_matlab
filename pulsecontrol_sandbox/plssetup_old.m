function plssetup_old(plsfile)
%function plssetup_old(plsfile)
% plssetup for loading plsdata from a different file
% handy for analysis of old data

% (c) 2011 Hendrik Bluhm.  Please see LICENSE and COPYRIGHT information in plssetup.m.
% Copyright 2011 Hendrik Bluhm, Vivek Venkatachalam
% This file is part of Special Measure.
% 
%     Special Measure is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
% 
%     Special Measure is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
% 
%     You should have received a copy of the GNU General Public License
%     along with Special Measure.  If not, see <http://www.gnu.org/licenses/>.

global awgdata;
global plsdata;

if ~exist('plsfile','var') || isempty(plsfile)
   [plsfile dr] = uigetfile('*');
   plsfile = [dr plsfile];
end

tic;
% if strcmp(computer, 'GLNX86')
%     plsdata.datafile = '~ygroup/qDots/awg_pulses/plsdata_2012_08_22.mat';
% else
%     if exist('z:/qDots','file')
%       plsdata.datafile = 'z:/qDots/awg_pulses/plsdata_2012_08_22.mat';
%       %plsdata.datafile = 'z:/qDots/awg_pulses/plsdata_1110.mat';
%     else
%       plsdata.datafile = 'y:/qDots/awg_pulses/plsdata_2012_08_22.mat';  
%     end
% end
evalin('base','plsdata.datafile = plsfile;');
evalin('base','plssync(''load'');');
evalin('base','awgloaddata;');
return;

% using different interfaces
%awg = tcpip('140.247.189.142', 4000); % slow
%awg = visa('ni', 'TCPIP::140.247.189.142::INSTR');
%awg.OutputBufferSize = 2^18;
