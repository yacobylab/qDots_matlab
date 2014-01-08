function at_grad_adjust2(opts,bases)

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
sides={'left','right'}; 
rfs={'PlsRamp1','PlsRamp2','PlsRamp3','PlsRamp4'}; 
if ~isempty(strfind(opts,'short'))
    dv=[1 4].*1e-4;
else
    dv=[1 2 4 8].*1e-4; 
end

%For loop over the gates. For each gate, first run a charge scan to get original trip pts.
%Then increment by vstep, run scan for each side, and find how triple point
%moved. Update basis by adding in part of the X Y bases. 
for s=1:length(sides); 
    atswap(sides{s});
    stp0=stpscan(rfs,[0 0 0 0],side{s},0); %No voltage applied on PlsRamps
    tl0=tlscan(rfs,[0 0 0 0],side{s},0);
    xy{1}=['X' upper(sides{s}(1))]; 
    xy{2}=['Y' upper(sides{s}(1))]; 
    basxy=basislookup(xy); 
    for i=1:2
       atchg(xy{i},dv(1))
       comp(1,i) = (stpscan(rfs,[0 0 0 0],side{s},i)-stp0)/dv(1); % comp has format [stpX stpY; tlX tlY]
       comp(2,i) = (tlscan(rfs,[0 0 0 0],side{s},i)-tl0)/dv(1);
       atchg(xy{i},-dv(1))
    end
    for p=1:length(gates); 
        cont=1; vn=1; vstep=dv; 
        while cont
            atchg(gates{p},vstep);
            stppt=stpscan(rfs,[0 0 0 0],side{s},p+2);
            tlpt=tlscan(rfs,[0 0 0 0],side{s},p+2);
            atchg(gates{p},-vstep)
            dvec=[stppt-stp0 ; tlpt-tl0]; 
            grad=-pinv(comp)*dvec./vstep(j);
            dbasis=tunedata.basis(1:4,basxy)*grad;
            err(vn)=mean(grad);
            fprintf('Add %s %.2f, %s %.2f to %5s? \n',xy{1},grad(1),xy{2},grad(2),gates{p});   
            info=input('(yes/no:', 's');
            doit = strcmp(info, 'yes');
            if doit tunedata.basis(1:4,b(p))=tunedata.basis(1:4,b(p))+dbasis; end
            %the basis has format row 1:4 1a:4a, column 9:12 XL:YR
            vn=vn+1; 
            info=input('Use a larger voltage? (yes/no): ','s');
            cont=strcmp(info,'yes'); 
            vstep=2*dv;           
        end
        fprintf('The error rate for %5s changed from %.2f to %.2f \n', gates{p},err(1),err(end)) 
    end
end

end

function [basis b]=getgates(opts)

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
