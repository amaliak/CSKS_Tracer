function writeINCONEOS1X(eSol,param,phase)

% function to write the very first INCON for EOS1
% need to provide the matrix eSol corresponding to
% the same domain. 

% eSol should be created with 
% readSAVEEOS1 using any SAVE file that corresponds
% to the same domain (e.g. from test runs)
% the SAVE file is only used for the dimensions of 
% the problem and for element names

n=length(eSol.pressure);
eSol.elem = strtrim(eSol.elem);
% eSol.elem needs to be of lenth 5
if length(eSol.elem{1})>5; keyboard; end;

firstpass = 0; 
if firstpass == 1; disp('Manual uniform initial conditions'); end; 
if firstpass
    pressure_bc = 0.1 * 10^8; 
    temperature_bc = 90 ;
    pressure  = pressure_bc*ones(size(eSol.pressure));
    temperature = temperature_bc*ones(size(eSol.pressure));
else
    pressure = eSol.pressure;
    concentration = eSol.X;
    temperature = abs(eSol.T);
    %temperature=76*ones(size(temperature));
end

fid1 =fopen('INCON', 'w');
fprintf(fid1, 'INCON\n');
for i=1:1:n
    fprintf(fid1, '%5s           %14.8E\n',eSol.elem{i},eSol.por(i)); % cell index and porosity
    fprintf(fid1, ' %.13E ',pressure(i));% pressure
    fprintf(fid1, '%.13E ',temperature(i)); %temperature
    fprintf(fid1, '%.13E\n',concentration(i));
end
% phase is given in the input because it is different between the TOUGH2obs and TOUGH2update
%phase = input('Which phase is this (now finishing INCON)');
%phase = param.phase; 
time = cumsum(param.dT) * 86400; 
if phase == 1 
	fprintf(fid1,'\n'); % must include a blank record at the end, or '+++' with restart time information
else 
	fprintf(fid1,'+++\n');
        fprintf(fid1,'             11 0.00000000E+00 %14.8e\n',time(phase-1));
end

fclose(fid1);
end
