function eSol=readSAVEEOS1T(varargin)
if nargin ==2
    eSol = varargin{1};
    param = varargin{2};
end

fid = fopen('SAVE');
fgetl(fid); %ignore title line
% Primary variables
pressure=[]; T = []; X = [];
elem = []; por = [];

s = fgetl(fid); % read INCON.1
s1=fgetl(fid);  % read INCON.2
while 1==1
    tmp1 = sscanf(s,'%8c%f')';
    elem = [elem; char(tmp1(1:8))];
    por = [por;tmp1(9)];
    tmp = sscanf(s1,'%f %f %f',[3,1])';
    pressure=[pressure;tmp(1)];
    T = [T; tmp(2)];
    X = [X; tmp(3)];
    s = fgetl(fid);
    s1=fgetl(fid);
    if  (~isempty(strfind(s,'+++'))) || (isempty(s))
        %        disp('It is the end of this file');
        break;
    end
end
fclose(fid);
eSol.pressure = pressure;
eSol.T = T;
eSol.elem = cellstr(elem); %convert a character array into cell of stirng
eSol.por = por;
eSol.X = X; 
%if nargin == 2
%eSol.vec = eSol.transform(eSol,param);
%end

end
