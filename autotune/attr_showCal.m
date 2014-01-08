function data=attr_showCal(out,opts,cals)
% Nicely show calibration for an attr result.
if ~exist('cals','var') || isempty(cals)
    ccals=regexp(fieldnames(out),'^proccal_(.*)','tokens');    
    cals={};
    for i=1:length(ccals)
        if ~isempty(ccals{i})
         cals=[cals ccals{i}{1}];
        end
    end
end
if ~exist('opts','var')
    opts=struct();
end
if ischar(opts)
    opts=struct('opts',opts);
end
opts=def(opts,'opts','');

if ~iscell(cals)
   cals={cals}; 
end

% dbloch style
%
c=[255/255 139/255 29/255];
c=[139/255 255/255 29/255];
cmds={};
cmds=[cmds struct('type','disc','val',[1 0 0],'color',[c .15],'color2',c)];
cmds=[cmds struct('type','disc','val',[0 1 0],'color',[c .15],'color2',c)];
cmds=[cmds struct('type','disc','val',[0 0 1],'color',[c .15],'color2',c)];
cmds=[cmds struct('type','sphere','color',[c .05])];
cmds=[cmds struct('type','sphere','color',[1 0 0 .05],'radius',1/2)];
%cmds=[cmds struct('type','sphere','color',[1 0 0 .05],'radius',1/2)];
cmds=[cmds struct('type','label','val',[0 0 1.15 ],'label','|S>')];
cmds=[cmds struct('type','label','val',[0 0 -1.15],'label','|T>')];
cmds=[cmds struct('type','label','val',[1.15 0 0],'label','|\uparrow\downarrow>')];
cmds=[cmds struct('type','label','val',[-1.15 0 0],'label','|\downarrow\uparrow>')];
if isfield(out.axes,'x')
  cmds=[cmds struct('type','vector','val',out.axes.x,'label','X','color',[0 0 0])];
end
if isfield(out.axes,'y')
  cmds=[cmds struct('type','vector','val',out.axes.y,'label','Y','color',[0 0 0])];
end
if isfield(out.axes,'z')
  cmds=[cmds struct('type','vector','val',out.axes.z,'label','Z','color',[0 0 0])];
end

for i=1:length(cals)    
    f=out.(['proccal_' cals{i}]);
    if ~isfield(f,'data')
        fprintf('%-10s: Missing data\n',cals{i});
        continue;
    end
    af=f.data.approx_fid;
    data(i).name=cals{i};
    data(i).dict=f.dictentry;
    data(i).fid=f.data.fid;
    data(i).chisq=f.data.chisq;
    data(i).af = mean(af);
    
    fprintf('%-10s (%-8s): Fidelity %-6.4g Chisq %-4.2g, Apprx. Fid. %5.3g (%5.3g-%5.3g)\n', ...
        cals{i},f.dictentry,f.data.fid, f.data.chisq, mean(af), min(af), max(af));
    if isopt(opts,'detail')
      fprintf('   Axis: %-5.3g,%-5.3g,%-5.3g, angle %-5.3g pi radians\n', ...
        f.data.eigs.axis, f.data.eigs.angle/pi);    
      fprintf('  |Eigs|: %g,%g,%g\n',f.data.eigs.values);
      fprintf('  chi |Eigs|: %g,%g,%g,%g\n',sort(eig(MtoChi(f.data.pm))));
    end
    nm=regexprep(f.dictentry,'R(.)_','R_$1');
    nm=regexprep(nm,'pi','\\pi');
    nm=regexprep(nm,'pi2','pi/2');
    cmds=[cmds struct('type','vector','val',f.data.eigs.axis * f.data.eigs.angle/(pi),'label',nm,'labelsep',.3)];
    if ~isopt(opts,'noplot') && ~isopt(opts,'nopplot')
    figure(500+i);    
    clf;    
    ana_process_plot(f.data.pm,'cchibar samefig',f.target);        
    subplot(121);
    title(nm);
    subplot(122);
    title(sprintf('Fidelity %5.3g',f.data.fid));
    end
end 
if ~isopt(opts,'noplot') && ~isopt(opts,'nobplot')
figure(400);
clf;
dbloch(cmds,struct('opts',' resize raise','dt',.1));
end
end

function s=def(s,f,v)
if(~isfield(s,f))
    s.(f) = v;
end
end

function b=isopt(config,name)
b=~isempty(strfind(config.opts,name));
end