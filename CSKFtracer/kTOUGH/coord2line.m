function LNUM=coord2line(x,y,z)

%x,y,z provided as vectors of corresponding xyz points

%x = [281:2:359];
%y = [281:2:359];
%z = [1019:-2:1001];

% LNUM is the line number in the MESH FILE
% Note that if MESH file has been extracted to a vector the +1 for the header needs to be subtracted

% set up coordinates

flag = 1 ; % for non-regular grid

if ~flag
    % when x,y,z are meshgrid type vectors 
    no_b = length(x)*length(y)*length(z);
    coordinates=zeros(no_b,3);

    c=1;
    for i=1:length(x)
        for j=1:length(y)
            for k=1:length(z)
                coordinates(c,1) = x(i);
                coordinates(c,2) = y(j);
                coordinates(c,3) = z(k);
                c=c+1;
            end
        end
    end

else
    
    coordinates(:,1)=x;
    coordinates(:,2)=y;
    coordinates(:,3)=z;
    no_b = length(coordinates);
end

%%%%
meshpath = '/home/amaliak/Pumping/INV/tmp/state_update';
if ~exist([meshpath,'/MESH'])
    disp(['Please copy MESH file to directory RUNID'])
else
    % convert coordinates to ID's
    if exist([meshpath,'/MESHmap.mat'])
        load([meshpath,'/MESHmap.mat'])
    else
        G = readMESH([meshpath,'/MESH']);
        save([meshpath,'/MESHmap.mat'],G)
    end
    

    for i = 1:no_b % for each possible observation location
        locx=logical(G.x==coordinates(i,1));
        locy=logical(G.y==coordinates(i,2));
        locz=logical(G.z==coordinates(i,3));
        loc=find(locx.*locy.*locz ==1); % location in G matrix and MESHmap.mat, maybe useful later
        if length(loc)>1 loc=loc(1); end;
        %BHID(i,1)=G.elem(loc);
        LNUM(i,1)=loc+1;
    end
end
