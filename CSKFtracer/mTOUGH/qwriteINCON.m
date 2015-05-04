function qwriteINCON(eSol,param,phase)

% this function converts the structure eSol to a TOUGH2 compatible INCON
% file

% param: needed for restart files to determine the time of current phase
% that needs to be written at the end of file (only param.dT needed)

% phase: which phase is this INCON for
% phase is read separately from param in the input because 
%         it is different 
%         between the TOUGH2obs and TOUGH2update

% this file works for all modules and supports both negative and positive
% numbers in the INCON
% Generates warning if there are negative (nonphysical) numbers

% if warning = 1, warning is generated for number of negative values in INCON
% takes some additional time (0.1 sec per file)

warning = 1 ; 

nofields = length(fieldnames(eSol));
numfields = nofields - 5; % excluding elem, transform, por, vec and pmx (not needed for INCON)
idfields = fieldnames(eSol);

eSol.elem = strtrim(eSol.elem); %this is likely not necessary and will 
% cause problems if the element name is longer than 5 characters (see fprintf, spaces are hardcoded)
% eSol.elem needs to be of lenth 5 - this is hard coding - should be
% corrected 
if length(eSol.elem{1})>5; keyboard; end;

%%
% firstpass is a flag set to 1 for user-generated INCON (as opposed to
% INCON generated from structure eSol)
firstpass = 0; 
 
if firstpass
    disp('Manual uniform initial conditions');
    % this is to create initial conditions manually
    % adjust accordingly for each case
    pressure_bc = 0.1 * 10^8; 
    temperature_bc = 90 ;
    pressure  = pressure_bc*ones(size(eSol.pressure));
    temperature = temperature_bc*ones(size(eSol.pressure));
end
%%
% File header
fid1 =fopen('INCON', 'w');
fprintf(fid1, 'INCON\n');

% in all modules
c1 = eSol.elem;
c2 = cellstr(num2str(eSol.por,'%16.8E\n'));
str = [c1';c2']; % holds element and porosity (first line of INCON)

% variable length of modules
for i = 1:numfields
    % numeric fields in eSol
    locnum = [1;2;6]; % denotes the location of the numeric fields in eSol, specific to tracer, correct hardcoding
    % Generate warning - takes an extra 0.1 sec 
    if warning
        eval(['negatives = length(','eSol.',char(idfields(locnum(i))),'(eSol.',char(idfields(locnum(i))),'<0)',');'])
        disp(['Warning!! There are ',num2str(negatives),' ',char(idfields(locnum(i))),' negative values in INCON'])
    end
    % Generate strings for fprintf should hold state variables in right order! for Tracer P T X
    c = eval(['cellstr(num2str(','eSol.',char(idfields(locnum(i))),',''%+19.13E\n''','));']);
    str = [str;c'];
end

fmtstring = repmat('%s', 1, numfields);
fprintf(fid1,['%s           %s\n',fmtstring,'\n'],str{:});

% File ending for restart conditions

time = cumsum(param.dT) * 86400; 
if phase == 1 
	fprintf(fid1,'\n'); % must include a blank record at the end, or '+++' with restart time information
else 
    % this part needs to be hardcoded depending on the simulation
    % 11: number of materials + SEED
    % last number: the initial time of the simulation 
	fprintf(fid1,'+++\n');
    fprintf(fid1,'             11 0.00000000E+00 %14.8e\n',time(phase-1));
end

fclose(fid1);
