function gradfbui

% (c) 2010 Hendrik Bluhm.  Please see LICENSE and COPYRIGHT information in plssetup.m.


global fbdata;
persistent hv;

if ishandle(1033);
    clf(1033);
else
    figure(1033);
    set(1033, 'position', [800, 300, 220, 480], 'MenuBar', 'none', ...
        'Name', 'Gradient feedback');
end

x1 = 10;
y1 = 480;
dx1 = 60;
dy1 = 25;

hv(1) = uicontrol('style', 'edit', 'position', [x1+dx1, y1+3, 100, 20], ...
    'HorizontalAlignment', 'right', 'string', sprintf('%g', fbdata.freq), 'userdata', fbdata.freq, 'callback', @editcheck);

hv(2) = uicontrol('style', 'edit', 'position', [x1+dx1, y1-dy1+3, 100, 20], ...
    'HorizontalAlignment', 'right', 'string', sprintf('%g', fbdata.setp), 'userdata', fbdata.setp, 'callback', @editcheck);
hv(3) = uicontrol('style', 'edit', 'position', [x1+dx1, y1-2*dy1+3, 100, 20], ...
    'HorizontalAlignment', 'right', 'string', sprintf('%g', fbdata.pgain), 'userdata', fbdata.pgain, 'callback', @editcheck);
hv(4) = uicontrol('style', 'edit', 'position', [x1+dx1, y1-3*dy1+3, 100, 20], ...
    'HorizontalAlignment', 'right', 'string', sprintf('%g', fbdata.igain), 'userdata', fbdata.igain, 'callback', @editcheck);
hv(5) = uicontrol('style', 'edit', 'position', [x1+dx1, y1-4*dy1+3, 100, 20], ...
    'HorizontalAlignment', 'right', 'string', sprintf('%g', fbdata.ctrlval(end)), 'userdata', fbdata.ctrlval(end), 'callback', @editcheck);

hv(6) = uicontrol('style', 'checkbox', 'position', [x1+dx1+20, y1-5*dy1+3, 20, 20], 'value', fbdata.fbon, 'callback', @fbenable);
hv(7) = uicontrol('style', 'checkbox', 'position', [x1+dx1+80, y1-5*dy1+3, 20, 20], 'value', fbdata.plot, 'callback', @plotenable);

hv(8) = uicontrol('style', 'pushbutton', 'position', [20, y1-160, 50, 30], 'string', 'Read', 'callback', @getfn);
hv(9) = uicontrol('style', 'pushbutton', 'position', [100, y1-160, 50, 30], 'string', 'Set', 'callback', @setfn);

t=uicontrol('style', 'edit', 'position', [x1, y1-270, 180, 90], 'HorizontalAlignment', 'left', 'Max', 10, 'KeyPressFcn', @runcmd);
uicontrol('style','pushbutton','position',[x1+185,y1-270,25,90],'string','^E','callback',@runcmd2,'userdata',t);
for l=1:8;
  t=uicontrol('style', 'edit', 'position', [x1, y1-295-25*(l-1), 180, 20], 'HorizontalAlignment', 'left', 'Max', 1, 'KeyPressFcn', @runcmd);
  uicontrol('style','pushbutton','position',[x1+185,y1-295-25*(l-1),25,20],'string','^R','callback',@runcmd2,'userdata',t);
end
%uicontrol('style', 'pushbutton', 'position', [140, y1-440, 50, 30], 'string', 'Run', 'callback', @runcmd);



uicontrol('style', 'text', 'position', [x1, y1, 50, 20], ...
    'HorizontalAlignment', 'left', 'string', 'Freq:', 'BackgroundColor', [.8 .8 .8]);
uicontrol('style', 'text', 'position', [x1, y1-dy1, 50, 20], ...
    'HorizontalAlignment', 'left', 'string', 'Setpoint:', 'BackgroundColor', [.8 .8 .8]);
uicontrol('style', 'text', 'position', [x1, y1-2*dy1, 20, 20], ...
    'HorizontalAlignment', 'left', 'string', 'P:', 'BackgroundColor', [.8 .8 .8]);
uicontrol('style', 'text', 'position', [x1, y1-3*dy1, 20, 20], ...
    'HorizontalAlignment', 'left', 'string', 'I:', 'BackgroundColor', [.8 .8 .8]);
uicontrol('style', 'text', 'position', [x1, y1-4*dy1, 60, 20], ...
    'HorizontalAlignment', 'left', 'string', 't_pol', 'BackgroundColor', [.8 .8 .8]);
uicontrol('style', 'text', 'position', [x1, y1-5*dy1, 80, 20], ...
    'HorizontalAlignment', 'left', 'string', 'Feedback on:', 'BackgroundColor', [.8 .8 .8]);
uicontrol('style', 'text', 'position', [x1+dx1+50, y1-5*dy1, 30, 20], ...
    'HorizontalAlignment', 'left', 'string', 'Plot:', 'BackgroundColor', [.8 .8 .8]);



function fbenable(h, eventdata, handles)
    fbdata.fbon = get(h, 'value');
end

function plotenable(h, eventdata, handles)
    fbdata.plot = get(h, 'value');
end


function getfn(h, eventdata, handles)    
    set(hv(1), 'string', sprintf('%g', fbdata.freq), 'userdata', fbdata.freq);
    set(hv(2), 'string', sprintf('%g', fbdata.setp), 'userdata', fbdata.setp);
    set(hv(3), 'string', sprintf('%g', fbdata.pgain), 'userdata', fbdata.pgain);
    set(hv(4), 'string', sprintf('%g', fbdata.igain), 'userdata', fbdata.igain);
    set(hv(5), 'string', sprintf('%g', fbdata.ctrlval(end)), 'userdata', fbdata.ctrlval(end));
    set(hv(6),  'value', fbdata.fbon);
    set(hv(7),  'value', fbdata.plot);
end    

function setfn(h, eventdata, handles)    
    fbdata.freq = get(hv(1), 'userdata');
    fbdata.setp = get(hv(2), 'userdata');
    fbdata.pgain = get(hv(3), 'userdata');
    fbdata.igain = get(hv(4), 'userdata');
    %if get(hv(5), 'userdata') > 0
    %    fbdata.ctrlval(end) = get(hv(5), 'userdata');
    %end
    fbdata.ctrlval(end) = get(hv(5), 'userdata');
    fbdata.intval(end) = fbdata.ctrlval(end) - fbdata.pgain * (fbdata.fbval(end) - fbdata.setp);
end    



function editcheck(h, eventdata, handles)    
    val = sscanf(get(h, 'string'), '%f');
    if ~isempty(val)
        set(h, 'userdata', val);
    end
    set(h, 'string', sprintf('%g', get(h, 'userdata')));
end
function runcmd2(h,eventdata,handles)
  cmd = get(get(h,'userdata'), 'string');
  cmd = cellstr([cmd, repmat(',', size(cmd, 1), 1)]);
  evalin('base', [cmd{:}]);
end
function runcmd(h, eventdata, handles)
    if length(eventdata.Modifier) == 1 && strcmp(eventdata.Modifier{1}, 'control') ...
            && strcmp(eventdata.Key, 'return')
        cmd = get(h, 'string');
        cmd = cellstr([cmd, repmat(',', size(cmd, 1), 1)]);
        evalin('base', [cmd{:}]);
    end
end
end
