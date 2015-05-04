function writeINCON(eSol)
n=length(eSol.pressure);
fid1 =fopen('INCON', 'w');
fprintf(fid1, 'INCON\n');
for i=1:1:n
    fprintf(fid1, '%8s        %14.8E\n',eSol.elem{i},eSol.por(i)); % cell index and porosity
    fprintf(fid1, ' %.13E ',eSol.pressure(i));% pressure
    fprintf(fid1, '%.13E ',eSol.Xsalt(i)); % salinity
    fprintf(fid1, '%.13E ',eSol.s(i,1)+10); % CO2 saturation + 10
    fprintf(fid1, '%.13E\n',eSol.T(i));% Temprature
end
fprintf(fid1,'\n'); % must include a blank record at the end, or '+++' with restart time information
fclose(fid1);
end
