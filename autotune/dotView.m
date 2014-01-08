function dotView(cntl, file, gates)
% dotView(cntl, file, gates)
%
% cntl: init (takes argument file, gates), viewer asks user to click on
% position of the gates in the order that they are listed
% cntl: refresh: refresh after gate values are changed
% cntl: save: save viewdata
% cntl: reload: reload viewdata to avoid reinitializing 

% (c) 2010 Hendrik Bluhm.  Please see LICENSE and COPYRIGHT information in plssetup.m.

if ~exist('cntl','var')
    cntl='refresh';
end

global viewdata;
global smdata;
%global coords;
%filename = 'CNOpic';
%format = 'JPEG';


if strcmp(cntl, 'init')
    
    viewdata.gates = smchanlookup(gates);
    viewdata.file = file;
    set(figure(1010),'Name','Dot View');
    image(imread(viewdata.file))
    axis image;

    for i = 1:length(viewdata.gates)
        viewdata.txth(i) = text(20, 20, smdata.channels(viewdata.gates(i)).name);
        viewdata.coords(i, :) = ginput(1);
        set(viewdata.txth(i), 'position', viewdata.coords(i, :));        
    end
    return;
end


if strcmp(cntl, 'refresh')
    if isempty(viewdata)
       dotView('reload'); 
    end
    gatenum = length(viewdata.gates);
    if ~ishandle(1010)
        figure(1010);
        image(imread(viewdata.file));
        axis image
    else
        figure(1010); % Raise figure
    end
    
    for i = 1:gatenum
        gate = num2str(cell2mat(smget(viewdata.gates(i))), 3);
        if ~ishandle(viewdata.txth(i))
            viewdata.txth(i) = text(viewdata.coords(i,1), viewdata.coords(i,2), gate, 'color', 'r');
        else
            set(viewdata.txth(i), 'string', gate);        
        end
    end
elseif strcmp(cntl, 'save')
    save viewdata viewdata;
elseif strcmp(cntl, 'reload')
    load viewdata;    
else
    error('ivalid input arguments');
end
end

