function at_grad_adjust2(bases,dv,side)

%This is meant to improve a specific set of basis vectors. By running charge 
%scans, we see how those bases move the triple points, and add the X/Y
%bases accordingly to stabilize. Iterate over increasingly large changes,
%as hopefully bases will move less and less. 
%Later write in a way to average. 
%opts can be 'short' (to use reduced number of iterations)
%bases can be 'all', 'leads', 'NTs', or any combination of gates. 

global tunedata;

[gates b]=getgates(bases); %puts selected bases into cell format. 
atswap('left'); atswap('right'); %to save tunedata
if ~exist('side','var')
    sides={'left','right'}; 
else 
    sides={side};
end


%For loop over the gates. For each gate, first run a charge scan to get original trip pts.
%Then increment by vstep, run scan for each side, and find how triple point
%moved. Update basis by adding in part of the X Y bases. 
for s=1:length(side) %over sides
    atswap(sides{s});
    chg0=mean(getnewtrip()); %the starting triple point 
    xy{1}=['X' upper(sides{s}(1))]; 
    xy{2}=['Y' upper(sides{s}(1))]; 
    basxy=basislookup(xy); 
    trippt=[];
    err=[];
    for i=1:2
       atchg(xy{i},dv)
       comp(1:2,i) =mean(getnewtrip(),1); % comp has format [stpX stpY; tlX tlY]
       atchg(xy{i},-dv)
    end        

        
    atchg(gates{i},-vstep(vn))

        dtrippt=trippt(:,:,vn)-trippt(:,:,1); %format is 1st row dxl dxr, 2nd dyl dyr
        %Now we know how much things have moved.
        grad=reshape(dtrippt./vstep,4,1);
        dbasis=tunedata.basis(1:4,9:12)*grad;
        err(vn)=mean(grad);
        fprintf('Add XL %.2f, YL %.2f, XR $.2f, YR %.2f to %5s? \n',grad(1), grad(2),grad(3),grad(4),gates{i});   
        info=input('(yes/no/[stop]:', 's');
        doit = strcmp(info, 'yes');
        stop=strcmp(info,'stop');
        if stop return; end 
        if doit tunedata.basis(1:4,b(i))=tunedata.basis(1:4,b(i))+tunedata.basis(1:4,9:12)*dbasis; end
        %the basis has format row 1:4 1a:4a, column 9:12 XL:YR
        %now, write code to change the basis. 
        end
    end
    fprintf('The error rate for %5s changed from %.2f to %.2f \n', gates{i},err(2),err(end)) %choose err(2) because the first is that for the charge scan. 
end

end

function [gates b]=getgates(opts)

bases=[];
if ~isempty(strfind(opts,'all'))
    bases={'Lead1', 'Lead2', 'Lead3', 'Lead4', 'T12', 'T34', 'N12', 'N34'};
end
if ~isempty(strfind(opts,'leads'))
    bases={'Lead1', 'Lead2', 'Lead3', 'Lead4'};
end
if ~isempty(strfind(opts,'NTs')) 
    bases= {'T12', 'T34', 'N12', 'N34'};
end

if ~isempty(strfind(opts,'Lead1')) bases{end+1}='Lead1'; end 
if ~isempty(strfind(opts,'Lead2')) bases{end+1}='Lead2';  end 
if ~isempty(strfind(opts,'Lead3')) bases{end+1}='Lead3';  end 
if ~isempty(strfind(opts,'Lead4')) bases{end+1}='Lead4';  end 
if ~isempty(strfind(opts,'T12')) bases{end+1}='T12';  end 
if ~isempty(strfind(opts,'T34')) bases{end+1}='T34';  end 
if ~isempty(strfind(opts,'N12')) bases{end+1}='N12';  end 
if ~isempty(strfind(opts,'N34')) bases{end+1}='N34';  end 

bases=unique(bases); 
b=basislookup(bases);
end        
        
function triple=getnewtrip()
global tunedata
        %Run the same program as in atxyfix2. 
        %Then, need to adjust basis mid scan: Find the amount the
        %triple points moved 
    tunedata.chrg.scan.consts(2).val=awgseqind(tunedata.chrg.pls); %PulseLine
    smset({tunedata.chrg.scan.consts.setchan}, [tunedata.chrg.scan.consts.val]); %un-obsolete
    scan = smscanpar(tunedata.chrg.scan, tunedata.cntr);
    data = smrun(scan); %set up and run scan 
    if any(isnan(data{1}(:))); smset(1:18, cell2mat(gates)); return; end
    figure(70); clf; hold on;
    imagesc(scan.loops(1).rng, scan.loops(2).rng, diff(data{1}));
    xlim([scan.loops(1).rng(1), scan.loops(1).rng(2)]);
    ylim([scan.loops(2).rng(1), scan.loops(2).rng(2)]);
    triple=ginput(2);
end

function b=basislookup(basis)
%give this either a single gate in string form. 
if ischar(basis)
    basis={basis}; 
end
if iscell(basis)
    for j = 1:length(basis); 
        b(j)=find(strncmp(basis{j},tunedata.basenames,length(basis{j}))); 
    end
else 
    fprintf('Please input a cell or a string')
    b=nan; 
end
end
