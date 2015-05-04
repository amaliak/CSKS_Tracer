function [string_before_pmx,string_after_pmx] = qsplitMESH(filename)
% function to split MESH file in three parts
% string up to PMX, PMX, and string after PMX

% create strings before and after PMX just once and use to create modified
% MESH file in qwriteMESH

% filename should include full path

if isempty(filename)
    txt = fileread('./MESH'); 
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

string_before_pmx = txt(:,1:40); % element name, 10 space characters, then material name
%string_pmx = txt(:,41:50);
%pmx = sscanf(strjoin(cellstr(txt(:,41:50))'),'%f',[nelem 1]); %gives pmx
%as number
string_after_pmx = txt(:,51:nchar);

end