function Gmod = qwriteMESH(infile,outfile,pmx_mod)

% function that reads an existing MESH file and modifies it to include
% the permeabilty modification given in vector K

% modified for performance (AK)

% infile: name of original MESH file --> string
% outfile: name of modified MESH file --> string
% pmx_mod: matrix with permeabilities --> matrix of reals

% example usage

% Gmod = qwriteMESH('./MESH','./MESHnew',9*ones(29805,1));

% NOTE: pmx_mod does not include con00!

s0 = dir(infile);

% either extract strings before and after pmx, or load them from memory
% extracting takes 0.154084 sec
% loading from memory takes 0.028 sec
if ~exist('MESHparts.mat')
    [string_before_pmx,string_after_pmx] = qsplitMESH(infile);
    string_before_pmx=cellstr(string_before_pmx);
    string_after_pmx=cellstr(string_after_pmx);
else
    load MESHparts
end

pmx = cellstr(num2str(pmx_mod,'%1.4E'));
pmx = pmx(1:end-1);
str = [string_before_pmx';pmx';string_after_pmx'];

fid  =  fopen('MESHtmp','wt') ;

fprintf(fid,'ELEME\n');
fprintf(fid,'%s%s%s\n',str{:});

%%
% read and copy the rest of the file
% disp('Copying CONNE...');

if exist('MESHCONN','file') %MESHCONN already includes con00
    disp('modifying MESH');
    % if MESHCONN exists it will be directly concatenated with MESH
else
    % can make this fast by using fileread
    disp('CONNE does not exist...');
    %fidconn  =  fopen('MESHCONN','wt') ;
    
    %fclose(fidconn);
end

fclose(fid);

% % combine MESH and MESHCONN
% 
system(['cat ','MESHtmp',' MESHCONN > ',outfile]);
delete('MESHtmp')
% % opentxt(infile);
% % opentxt(outfile);
% 
% % check size of new file to be sure
% 
s = dir(outfile);
% 
if abs(s.bytes - s0.bytes) > 0
    disp('Inconsistency in generated MESH file')
    keyboard;
end 

Gmod = qreadMESH('./',[]);

end


