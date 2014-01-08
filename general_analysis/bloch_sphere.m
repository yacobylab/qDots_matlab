function bloch_sphere(ecmds, opts)
% function bloch_sphere(cmds, opts)
% Draw a basic bloch sphere for sketching on

c=[255/255 139/255 29/255];
c=[139/255 255/255 29/255];
if ~exist('ecmds','var')
    ecmds=[];
end
if ~exist('opts','var')
    opts='raise';
end
cmds={};
cmds=[cmds struct('type','disc','val',[1 0 0],'color',[c .15],'color2',c)];
cmds=[cmds struct('type','disc','val',[0 1 0],'color',[c .15],'color2',c)];
cmds=[cmds struct('type','disc','val',[0 0 1],'color',[c .15],'color2',c)];
cmds=[cmds struct('type','sphere','color',[c .05])];
%cmds=[cmds struct('type','sphere','color',[1 0 0 .05],'radius',1/2)];
%cmds=[cmds struct('type','sphere','color',[1 0 0 .05],'radius',1/2)];
cmds=[cmds struct('type','label','val',[0 0 1.15 ],'label','|S>')];
cmds=[cmds struct('type','label','val',[0 0 -1.15],'label','|T>')];
cmds=[cmds struct('type','label','val',[1.15 0 0],'label','|\uparrow\downarrow>')];
cmds=[cmds struct('type','label','val',[-1.15 0 0],'label','|\downarrow\uparrow>')];
out.axes.x=[1 0 0];
out.axes.y=[0 1 0];
out.axes.z=[0 0 1];
if isfield(out.axes,'x')
  cmds=[cmds struct('type','vector','val',out.axes.x,'label','X','color',[0 0 0])];
end
if isfield(out.axes,'y')
  cmds=[cmds struct('type','vector','val',out.axes.y,'label','Y','color',[0 0 0])];
end
if isfield(out.axes,'z')
  cmds=[cmds struct('type','vector','val',out.axes.z,'label','Z','color',[0 0 0])];
end
dbloch([cmds ecmds],struct('opts', opts));

end

