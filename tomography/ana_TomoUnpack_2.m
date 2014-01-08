function [dataexp, covariances, scantime, s,grps] = ana_TomoUnpack_2(file, config)
% function [dataexp, covariances, scantime, s] = ana_TomoUnpack_2(file, config)
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

rd = cell(1,2);

%load files, concatenate them together.
for j = 1:length(file)
    s=load(file{j});
    plsgrps = {s.scan.data.pulsegroups.name};
    dbzgrps = find(cellfun(@(p) ~isempty(p),regexp(plsgrps,'^([dD][bB][zZ])')));
    tomogrps = setdiff(1:length(plsgrps),dbzgrps);
    rdr=anaRawUnpack(s.scan,s.data);
    scantime=getscantime(s.scan,s.data);
    
    % histogram the data
    %first get t1
    [t1lt t1l] = att1('left',scantime,config.t1Time{1});
    [t1rt t1r] = att1('right',scantime,config.t1Time{2});
    
    %histogram it
    rdd=anaRawScale(rdr,[t1l t1r],tomogrps);
    for i=1:2 %just two sensors
        rdd{i}=-(rdd{i}*2-1);
    end
    
    for k = 1:2;    
        rd{k} = [rd{k}; rdd{k}];
    end   
    
end

%permute data to make it easy to work with. 
 for k=1:2
       rd{k} = permute(rd{k}, [3 2 1]); 
 end

tomos = {'ST','UD','Y'};
for j = tomogrps%1:length(tomogrps)
  gd{j} = plsinfo('gd',plsgrps{(j)},[],scantime);
  grpnames{j} = gd{j}.pulses.groups;
  if isempty(strfind(grpnames{j}{1},'_L'))
      grpnames{j}=fliplr(grpnames{j});
  end
  leftnames{j} = grpnames{j}{1};
  rightnames{j} = grpnames{j}{2};
end
STL = []; UDL = []; YL = [];
STR = []; UDR = []; YR = [];

STL = find(cellfun(@(p) ~isempty(p),regexp(leftnames,'ST')));
UDL = find(cellfun(@(p) ~isempty(p),regexp(leftnames,'UD')));
YL = find(cellfun(@(p) ~isempty(p),regexp(leftnames,'Y')));

STR = find(cellfun(@(p) ~isempty(p),regexp(rightnames,'ST')));
UDR = find(cellfun(@(p) ~isempty(p),regexp(rightnames,'UD')));
YR = find(cellfun(@(p) ~isempty(p),regexp(rightnames,'Y')));

xi = UDL;
yi = YL;
zi = STL;

if length([xi,yi,zi])~=9
   error('problems unpacking left groups'); 
end
 
ix = UDR;
iy = YR;
iz = STR;

if length([ix, iy, iz])~=9
    error('problem unpacking right groups');
end

xx = intersect(UDL,UDR); if length(xx)~=1; error('error unpacking groups'); end;
xy = intersect(UDL,YR); if length(xy)~=1; error('error unpacking groups'); end;
xz = intersect(UDL,STR); if length(xz)~=1; error('error unpacking groups'); end;
yx = intersect(YL,UDR); if length(yx)~=1; error('error unpacking groups'); end;
yz = intersect(YL,STR); if length(yz)~=1; error('error unpacking groups'); end;
zx = intersect(STL,UDR); if length(zx)~=1; error('error unpacking groups'); end;
zy = intersect(STL,YR);  if length(zy)~=1; error('error unpacking groups'); end;
yy = intersect(YL,YR);  if length(yy)~=1; error('error unpacking groups'); end;
zz = intersect(STL,STR); if length(zz)~=1; error('error unpacking groups'); end;


 for k=1:2
       rd{k} = permute(rd{k}, [3 2 1]); 
 end
corrected =0;
figind = 0;
bs = {'XI','YI','ZI','IX','IY','IZ', 'XY', 'XZ','YX', 'YZ', 'ZX', 'ZY', 'XX', 'YY', 'ZZ'};
grps.bs=bs;
grps.names = grpnames;
grps.inds = {xi,yi,zi,ix,iy,iz,xy,xz,yx,yz,zx,zy,xx,yy,zz};

if isopt(config,'cal')
    if isempty(config.Ymult)
        config.Ymult=1;
    end
    
    if isempty(config.calFile)
        fprintf('find Left tomoCal file \n');
        ssL = uigetfile('sm*_L_*_mat*');
        if ~ischar(ssL) %protect against empty files: defaults to no calibration
%            ssL.pm = eye(3);
%            ssL.po = [-1 1];
            ssL.data_corrector = @(x) x;
           warning('no tomocal file. assuming no correction');
        else
            ssL = load(ssL);
        end
        fprintf('find Right tomoCal file \n');
        ssR = uigetfile('sm*_R_*_mat*');
        if ~ischar(ssR)
%            ssR.pm = eye(3);
%            ssR.po = [-1 1];
           ssR.data_corrector = @(x) x;
           warning('no tomocal file. assuming no correction');
        else
            ssR = load(ssR);
        end
    else
       ssL = load(config.calFile{1});
       ssR = load(config.calFile{2});
    end
    % apply sensor correction
%    sensorcorrection = @(x,po) po(1)+.5*(po(2)-po(1))*(x+1);
    
%     if ~isopt(config, 'nosensorcorrection')
%         rd{1} = sensorcorrection(rd{1},ssL.po);
%         rd{2} = sensorcorrection(rd{2},ssR.po);
%     end
else
    if isempty(config.Ymult)
        config.Ymult=-1;
    end
end


if isopt(config,'cal')
      %left and right rotation matricies
      ML = ssL.data_corrector(eye(3));%ssL.pm;
      MR = ssR.data_corrector(eye(3));%ssR.pm;
       
      L = [ML' zeros(3,12)];
      R = [zeros(3,3) MR' zeros(3,9)];
      MM = kron(ML',MR');
      M = [L; R; [zeros(9,6) MM]];
end

for j = 1:size(rd{2},2)
    IX =               nanmean(rd{2}(ix,j,:),1);
    IY =  config.Ymult*nanmean(rd{2}(iy,j,:),1);
    IZ =               nanmean(rd{2}(iz,j,:),1);
    XI =               nanmean(rd{1}(xi,j,:),1);
    YI =  config.Ymult*nanmean(rd{1}(yi,j,:),1);
    ZI =               nanmean(rd{1}(zi,j,:),1);
    XX =              rd{1}(xx ,j,:).*rd{2}(xx ,j,:);
    XY = config.Ymult*rd{1}(xy ,j,:).*rd{2}(xy ,j,:);
    XZ =              rd{1}(xz ,j,:).*rd{2}(xz ,j,:);
    YX = config.Ymult*rd{1}(yx ,j,:).*rd{2}(yx ,j,:);
    YY =              rd{1}(yy,j,:).*rd{2}(yy,j,:);
    YZ = config.Ymult*rd{1}(yz ,j,:).*rd{2}(yz ,j,:);
    ZX =              rd{1}(zx ,j,:).*rd{2}(zx ,j,:);
    ZY = config.Ymult*rd{1}(zy ,j,:).*rd{2}(zy ,j,:);
    ZZ =              rd{1}(zz ,j,:).*rd{2}(zz ,j,:);
    
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
