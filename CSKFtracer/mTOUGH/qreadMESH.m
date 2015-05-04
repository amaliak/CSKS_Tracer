function G = qreadMESH(folder,filename)
% function to read MESH file 
% faster version with sscanf

% example usage
% G = qreadMESH('./',[]);

cd(folder)

% change flag to 1 only if G needs to have values for material,volx and
% ahts
flag = 0;


if isempty(filename)
    txt = fileread('MESH'); 
else
    txt = fileread(filename); 
end

nhead = 6; % number of characters in header line
nchar = 80; % number of characters per line
nelem = (strfind(txt,'con00')-1-nhead)/(nchar+1);
if isempty(nelem)
    disp('Check MESH file - possibly wrong format');
end
lelem = 5; % length of element and material names
lspaces = 10;
lma = 5;

txt = txt(nhead+1:end); 
txt = txt(1:nelem*(nchar+1));
txt = reshape(txt', [nchar+1,nelem])';

% NOTE: ma, volx and ahtx are only needed because we use qreadMESH to
% create new files - to do: see how we can avoid reading ma, volx, ahtx, x,
% y and z every time since they don't change every time we need to modify
% MESH

elem = txt(:,1:lelem); % element name, 10 space characters, then material name
if flag; 
    ma = txt(:,lelem+lspaces+1:lelem+lspaces+lma+1);
    volx = sscanf(strjoin(cellstr(txt(:,21:30))'),'%f',[nelem 1]); 
    ahtx = sscanf(strjoin(cellstr(txt(:,31:40))'),'%f',[nelem 1]); 
end
% str2double is slower than sscanf
% pmx = str2double(cellstr(txt(:,41:50))); 
% x = str2double(cellstr(txt(:,51:60)));
% y = str2double(cellstr(txt(:,61:70)));
% z = str2double(cellstr(txt(:,71:80)));

pmx = sscanf(strjoin(cellstr(txt(:,41:50))'),'%f',[nelem 1]); 
x = sscanf(strjoin(cellstr(txt(:,51:60))'),'%f',[nelem 1]); 
y = sscanf(strjoin(cellstr(txt(:,61:70))'),'%f',[nelem 1]); 
z = sscanf(strjoin(cellstr(txt(:,71:80))'),'%f',[nelem 1]); 

% put data into structure G

G.elem = cellstr(elem);
if flag
    G.ma = cellstr(ma);
    G.volx = volx;
    G.ahtx = ahtx;
end
G.pmx = pmx;
G.x = x;
G.y = y;
G.z = z;


save('MESHmap.mat','G')

end