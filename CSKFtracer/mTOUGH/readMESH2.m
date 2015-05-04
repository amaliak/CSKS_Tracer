function G = readMESH2(filename)
%G.elem: cell name
%G.ma: material name
%G.volx: element volume
%G.ahtx: interface area
%G.pmx: K modifier
%G.x, G.y, G.z

if isempty(filename)
    fid = fopen('MESH');
else
    %filename = [folder,filename];
    fid = fopen(filename);
end

v = fgetl(fid);
if ~strfind(v,'ELEME')
    error('Incorrect format');
end
% Primary variables
x = []; y = []; z = [];
volx = []; ahtx = []; pmx = [];
elem = []; ma = [];
s = fgetl(fid); % current line
s1 = fgetl(fid); % next line
  count = 1;
while 1==1
%   changed AK May 28
%   reading everything as a character, so that MESH files 
%   with empty columns can be read too
%   first do the characters ...

    if length(s)<80
        disp('End of ELEME reached');
        break;
    end
    tmp1 = sscanf(s,'%8c%7c%5c%10c%10c%10c%10c%10c%10c')';
    
    elem = [elem; char(tmp1(1:8)')];
    [~] = char(tmp1(9:15));
    ma = [ma; char(tmp1(16:20)')];
%   now the numbe'rs
%   numstart is the location in the string where the numbers start
    numstart = 20;

    volx = [volx ; str2num(char(tmp1(numstart+1:numstart+10)'))];
    ahtx = [ahtx ; str2num(char(tmp1(numstart+11:numstart+20)'))];
    pmx  = [pmx  ; str2num(char(tmp1(numstart+21:numstart+30)'))];
    x    = [x    ; str2num(char(tmp1(numstart+31:numstart+40)'))];
    y    = [y    ; str2num(char(tmp1(numstart+41:numstart+50)'))];
    z    = [z    ; str2num(char(tmp1(numstart+51:numstart+60)'))];
    

    s = s1;
    s1 = fgetl(fid);
    count=count+1;
    

end
fclose(fid);

%convert a character array into cell of string

G.elem = cellstr(elem); 
G.ma = cellstr(ma);
G.volx = volx;
G.ahtx = ahtx;
G.pmx = pmx;
G.x = x;
G.y = y;
G.z = z;

end