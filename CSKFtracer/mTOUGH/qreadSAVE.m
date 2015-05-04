function eSol = qreadSAVE()
% For SAVE files with nvar columns of numeric data 
% Adjust eSol structure (end of file) accordingly

nvar = 4; % porosity, pressure, temperature, mole fraction (EOS1) TRAVER
nchar = 5; % number of charcters in element 

if nargin == 2
    eSol = varargin{1};
    param = varargin{2};
end

% read files
fid = fopen('SAVE');

fgetl(fid); % discard first line
eval(['data = fscanf(fid, ''%s %e\n',repmat('%e ',1,nvar-1),''');']);
eof=find(data==43,1,'first');
data=data(1:eof-1);
nrows = length(data)/(nvar+nchar); % transpose result 

fclose(fid);

% extract numbers
start = 6; %numbers start after 5th caracter
for ivar = 1:nvar
    num_idx(ivar,:) = start:(nchar+nvar):length(data); 
    start = start + 1;
end
    
num_idx = reshape(num_idx,(numel(num_idx)),1);

% extract characters (assumes total of 5 characters)
elem_idx = setxor((1:length(data)),num_idx); % complement of idx
elem = data(elem_idx);

elem = reshape(elem,[5, numel(elem)/5])'; % check numel(char_sensor)/3 is integer
elem = char(elem);

% for three columns of data in the SAVE file
num_data = data(num_idx);
num_data = reshape(num_data,[nvar,numel(num_data)/nvar])'; % check numel(obs)/3 is integer


% put to Sol structure for KF code
% retain standard order of eSol fields, P, T, elem, por, pmx, X, transform, vec
eSol.pressure = num_data(:,2);
eSol.T = num_data(:,3);
eSol.elem = cellstr(elem);
eSol.por = num_data(:,1);
eSol.pmx = [];
if nvar==4 % tracer test case
    eSol.X = num_data(:,4);
elseif nvar==5 % co2 case (judith please check)
    eSol.sat = num_data(:,4);
    eSol.sal = num_data(:,5);
end
