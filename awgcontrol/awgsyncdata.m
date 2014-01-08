function awgsyncdata(cntrl, file)
% awgsyncdata(cntrl, file)
% transfer longer one of awgdata.pulesdata and puledata in 
% awgdata.datafile to shorter one.
% if cntrl=='clr' is given, awgdata.pulsedata is removed.
% if cntrl contains 'save', data from memory is always written to file
% % if cntrl contains 'laod', data from memory is always read from file
% if cntrl contains 'cpsrc', data is is written to srcdata;
% file overrides awgdata.datafile if given.

% (c) 2010 Hendrik Bluhm.  Please see LICENSE and COPYRIGHT information in plssetup.m.


global awgdata;

if nargin < 1 
    cntrl = '';
end

if nargin < 2
    file = awgdata.datafile;
end

if isempty(strfind(cntrl, 'save'))
    load(file, 'pulsedata');
end

if strfind(cntrl, 'cpsrc')
    awgdata.srcdata = pulsedata;
    return;
end

if isempty(strfind(cntrl, 'save')) && (~isfield(awgdata, 'pulsedata') || length(pulsedata) >= length(awgdata.pulsedata))   ...
        || ~isempty(strfind(cntrl, 'load'))
    awgdata.pulsedata = pulsedata;
    fprintf('Loading pulse data\n');
else
    pulsedata = awgdata.pulsedata;
    save(file, 'pulsedata');
    fprintf('Saving pulse data\n');
end

if strfind(cntrl, 'clr')
    awgdata = rmfield(awgdata, 'pulsedata');
end

% alternative version:
% awgadm(cntrl)
% cntrl: [load|save|sync] [mem|file];
%        load: datafile -> memory. Implies mem.
%        save: memory -> datafile
%        sync: update shorter one. Implies mem.
%         
%        mem:  use memory. 
%        file: use file. Deletes info in memory!

% global awgdata;
% 
% if strfind(cntrl, 'load')
%     load(awgddata.datafile, 'pulsedata');
%     awgdata.pulsedata = pulsedata;
% elseif strfind(cntrl, 'save') 
%     pulsedata = awgdata.pulsedata;
%     save(awgddata.datafile, 'pulsedata');
% elseif strfind(cntrl, 'sync')
%     load(awgddata.datafile, 'pulsedata');    
%     if ~ isfield(awgdata, 'pulsedata') || length(pulsedata) > awgdata.pulsedata
%         awgdata.pulsedata = pulsedata;
%     else
%         pulsedata = awgdata.pulsedata;
%         save(awgddata.datafile, 'pulsedata');
%     end
% end

