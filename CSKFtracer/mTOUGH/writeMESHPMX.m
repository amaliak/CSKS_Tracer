function Gmod=writeMESHPMX(infile,outfile,Kmodv,G0)

% function that reads an existing MESH file and modifies it to include
% the permeabilty modification given in vector K

% infile: name of original MESH file --> string
% outfile: name of modified MESH file --> string
% K: matrix with permeabilities --> matrix of reals

%% figure out modifiers and put them in a vector

% Kbase=mean(mean(K));
% Kmod=K./Kbase;
% Kmodv=vec(Kmod);
% disp(['Mean value of K is',num2str(Kbase)])
% disp(['Use ',num2str(Kbase),' in  ROCKS '])

% disp('Now writing ELEME ...') % for speed issue
%% load initial MESH file and create Gmod
%G0       = readMESH(infile);
Gmod     = G0;
Gmod.pmx = Kmodv;
% if there are  negative Zs, the size of G0.Z has to be adjusted
% check if there are negative Zs



%% write to new MESH file

fid  =  fopen(outfile,'wt') ;

fprintf(fid,'%s\n','ELEME');
no_perm_block = 1; % blocks like con00 that do not have permeability values
for i=1:length(Kmodv)-no_perm_block
    % append G.elem (8c+7 spaces) and G.ma (5c) if needed.
    % function appstr (see below) will recognize how many 
    % characters exist, and will append as needed, no user intervention
    % required
    el0=appstr(G0.elem(i),8+7);
    ma0=appstr(G0.ma(i),5);
    %create line
    % if there are  negative Zs, the size of G0.Z has to be adjusted
    if min(Gmod.z)<0
        line = sprintf('%s%s%1.4E%1.4E%1.4E%1.4E%1.4E%1.3E',el0,ma0,G0.volx(i),...
        G0.ahtx(i),Kmodv(i),G0.x(i),G0.y(i),G0.z(i));
    else
        line = sprintf('%s%s%1.4E%1.4E%1.4E%1.4E%1.4E%1.4E',el0,ma0,G0.volx(i),...
        G0.ahtx(i),Kmodv(i),G0.x(i),G0.y(i),G0.z(i));
    end
    % write line
    fprintf(fid,'%s\n',line);
end
%% write the rest of the file after ELEME
% copy over from original MESH file
% find CONNE and copy line by line.

% we already know at what line CONNE will be, 
%so we don't have to look for it
fid0 = fopen(infile);
ELEMEsize=length(Kmodv); % this determines how many lines are already copied
% make sure to adjust - here I also copy the con00 line (29806'th element of MESH)
for ii=1:ELEMEsize
    s=fgetl(fid0);
end
s= fgetl(fid0); % to copy con00 line
fprintf(fid,'%s\n',s);

% print a blank line after end of ELEME
% and before CONNE
%fprintf(fid,'%s\n',''); 
s = fgetl(fid0);
% check we are at the right place in the file
if ~strfind(s,'CONNE')
     error('Incorrect format');
end
%%
% read and copy the rest of the file
% disp('Copying CONNE...');
fprintf(fid,'%s\n',s); 
s = fgetl(fid0);

while s>0
    fprintf(fid,'%s\n',s); 
    s = fgetl(fid0);
end

fclose(fid);
fclose(fid0);
disp('Modified file created');
% disp('Original and Modified file opening..');


% opentxt(infile);
% opentxt(outfile);


end

%% 
function strout=appstr(strin,nochar)

% function to append string as needed
% convert to string
sstrin=strin{1};
% now append
l=length(sstrin);
appl=nochar-l;
strin1=sstrin;
for s=1:appl
    strout1=[strin1,' '];
    strin1=strout1;
end
strout=strin1;

end


