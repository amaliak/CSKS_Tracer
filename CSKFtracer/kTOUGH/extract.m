function small=extract(big,fields,x,y,z,savearg)

% function that extracts the smaller domain with coordinates given by arg2
% from the bigger domain given by arg1

% Input = full solution from tough2 in vector form
%         x, y , z are the x,y,z coordinates in vector form, the number of
%         points will be x*y*z
% Output= small solution extracted from big solution in vector form for all
%         x,y,z pairs. Coordinates cycle first for z, then y and then x

% this function works for files created by TOUGH2 where the location of the
% block in the MESH file also corresponds to the location of the solution
% read from SAVE

% example usage: 
%x = [281:2:359];
%y = [281:2:359];
%z = [1019:-2:1001];
%fields = [1:5]; it is the fields of the structure that we extract from
%small = extract ( big , fields, x , y , z) 
% output will be a vector with 16000 elements (to be converted to a 20x20x10 volume)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% obtain line numbers from MESH file
% if x,y,z are param.x, param.y and param.z then savearg should be =1
if ~savearg
lineno=coord2line(x,y,z); 
else
load lineno.mat
end

% this lineno corresponds to the MESH file, which works when extracting small domain from big domain but does not work when extracting the measurements.

% repeat for all variables of structure big
names = fieldnames(big);

for n=1:fields
    eval(['value_big=big.',names{n},';'])
    value_small = value_big(lineno-1,:); % lineno in MESH file includes the header
    eval(['small.',names{n},'=value_small(:,:);']);
end
