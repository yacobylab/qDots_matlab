function [dataexp, covariances, scantime, s] = ana_TomoUnpack(file, config)
% function [dataexp, covariances, scantime, s] = ana_TomoUnpack(file, config)
% unpack two qubit tomo data
% file is a filename or cell array of filenames, left empty it prompts uigetfile
% multiselect is on. for now multiple files will just concatinate the data
% config is a struct with the usual format
% config.t1Time is a cell array that dictates whether T1 data is taken from
%  before or after. elements can be 'before' or 'after'. defaults to 
%  {'before', 'before');
% config.Ymult is a multiplier for the Y data- useful if not using tomocal
% config.opts cal (use cal files)
%   the function will prompt the user to load tomocal files (not the data
%   files but the results of tomocal (pm, po, etc) or to ignore it
% scantime and s are the last scantime and last (scan &data), respectively.
global tunedata;
%get file if needed
if ~exist('file','var') || isempty(file)
    file = uigetfile('sm*tomo*.mat', 'MultiSelect', 'on');
end

if ~iscell(file)
    file = {file};
end

if ~exist('config','var')
    config=struct();
end
% set up default config stuff
config = def(config,'opts','cal');   % Random boolean options
config = def(config,'t1Time',{'before', 'before'}); % defaults to try to get t1 from before {left, right}
config = def(config, 'Ymult', []); % default Y-multiplier = 1; can set to -1 if not calibrating tomo. empty=auto
config = def(config, 'calFile', '');

for j = 1:length(file)
fprintf('%s \n', file{j});
end

rd = {};

%load files, concatenate them together.
for j = 1:length(file)
    s=load(file{j});
    rdr=anaRawUnpack(s.scan,s.data);
    scantime=getscantime(s.scan,s.data);
    
    % histogram the data
    %first get t1
    [t1lt t1l] = att1('left',scantime,config.t1Time{1});
    [t1rt t1r] = att1('right',scantime,config.t1Time{2});
    
    %histogram it
    rdd=anaRawScale(rdr,[t1l t1r],2:length(s.scan.data.pulsegroups));
    for i=1:2
        rdd{i}=-(rdd{i}*2-1);
    end
    
    for k = 1:2;
        if j==1
            rd{k} = permute(rdd{k}, [3 2 1]);
        else
            rd{k} = [rd{k}; permute(rdd{k}, [3 2 1])];
        end
    end   
    
 end

 for k=1:2
       rd{k} = permute(rd{k}, [3 2 1]); 
 end
corrected =0;
figind = 0;
bs = {'XI','YI','ZI','IX','IY','IZ', 'XY', 'XZ','YX', 'YZ', 'ZX', 'ZY', 'XX', 'YY', 'ZZ'};

if isopt(config,'cal')
    if isempty(config.Ymult)
        config.Ymult=1;
    end
    
    if isempty(config.calFile)
        fprintf('find Left tomoCal file \n');
        ssL = uigetfile('sm*_L_*_mat*');
        if ~ischar(ssL) %protect against empty files: defaults to no calibration
           ssL.pm = eye(3);
           ssL.po = [-1 1];
           warning('no tomocal file. assuming no correction');
        else
            ssL = load(ssL);
        end
        fprintf('find Right tomoCal file \n');
        ssR = uigetfile('sm*_R_*_mat*');
        if ~ischar(ssR)
           ssR.pm = eye(3);
           ssR.po = [-1 1];
           warning('no tomocal file. assuming no correction');
        else
            ssR = load(ssR);
        end
    else
       ssL = load(config.calFile{1});
       ssR = load(config.calFile{2});
    end
    % apply sensor correction
    sensorcorrection = @(x,po) po(1)+.5*(po(2)-po(1))*(x+1);
    
    if ~isopt(config, 'nosensorcorrection')
        rd{1} = sensorcorrection(rd{1},ssL.po);
        rd{2} = sensorcorrection(rd{2},ssR.po);
    end
else
    if isempty(config.Ymult)
        config.Ymult=-1;
    end
end

% ST U/D Y in file.  Right sweeps faster.
% ST/ST(2) ST/UD(3) ST/Y(4)
% UD/ST(5) UD/UD(6) UD/Y(7)
%  Y/ST(8)  Y/UD(9)  Y/Y(10)
if isopt(config,'cal')
      %left and right rotation matricies
      ML = ssL.pm;
      MR = ssR.pm;
       
      L = [ML' zeros(3,12)];
      R = [zeros(3,3) MR' zeros(3,9)];
      MM = kron(ML',MR');
      M = [L; R; [zeros(9,6) MM]];
end

for j = 1:size(rd{2},2)
    IX =               nanmean(rd{2}([2 5 8 ],j,:),1);
    IY =  config.Ymult*nanmean(rd{2}([4 7 10],j,:),1);
    IZ =               nanmean(rd{2}([3 6 9 ],j,:),1);
    XI =               nanmean(rd{1}([2 3 4 ],j,:),1);
    YI =  config.Ymult*nanmean(rd{1}([8 9 10],j,:),1);
    ZI =               nanmean(rd{1}([5 6 7 ],j,:),1);
    XX =              rd{1}(2 ,j,:).*rd{2}(2 ,j,:);
    XY = config.Ymult*rd{1}(4 ,j,:).*rd{2}(4 ,j,:);
    XZ =              rd{1}(3 ,j,:).*rd{2}(3 ,j,:);
    YX = config.Ymult*rd{1}(8 ,j,:).*rd{2}(8 ,j,:);
    YY =              rd{1}(10,j,:).*rd{2}(10,j,:);
    YZ = config.Ymult*rd{1}(9 ,j,:).*rd{2}(9 ,j,:);
    ZX =              rd{1}(5 ,j,:).*rd{2}(5 ,j,:);
    ZY = config.Ymult*rd{1}(7 ,j,:).*rd{2}(7 ,j,:);
    ZZ =              rd{1}(6 ,j,:).*rd{2}(6 ,j,:);
    
    tomodata{j} = squeeze([ XI, YI, ZI, IX, IY, IZ, XX, XY, XZ, YX, YY, YZ, ZX, ZY, ZZ])';
    ndata=min(find(any(isnan(tomodata{j}),2)));   
    if isempty(ndata)
        ndata=size(XI,3);
    end
    covariances{j}=nancov(tomodata{j})/ndata;
    tomodata{j} = squeeze(nanmean(tomodata{j},1));
    if isopt(config,'cal')
        tomodata{j}=M*tomodata{j}'; 
        covariances{j}=M*covariances{j}*(M');
    end
    scramble=[1:6 8:10 12:14 7 11 15];
    dataexp(j,:)=tomodata{j}(scramble)'; % reshape back into order of bs
    covariances{j}=covariances{j}(scramble,scramble);
end
dataexp=dataexp';
end

% Apply a default.
function s=def(s,f,v)
if(~isfield(s,f))
    s=setfield(s,f,v);
end
return;
end

function b=isopt(config,name)
b=~isempty(strfind(config.opts,name));
return;
end