fid = fopen('TIME','w');
%% ROCK
fprintf(fid,'ROCKS----1----*----2----*----3----*----4----*----5----*----6----*----7----*----8\n');
fprintf(fid,'\n'); % must add a blank record at the end of ROCK section

% PARAM
PARAM.TIMAX = 86400*90; %simulation time
PARAM.DELTEN = -1.0; %
PARAM.MOP = zeros(24,1);
PARAM.MOP(16) = 3; % ITER<=3, then double time step length

% INCON
% ignore disk file INCON
fprintf(fid,'INCON----1----*----2----*----3----*----4----*----5----*----6----*----7----*----8\n');
fprintf(fid,'\n');

% TIME 
fprintf(fid,'TIMES----1----*----2----*----3----*----4----*----5----*----6----*----7----*----8\n');
ITI = 7;
fprintf(fid,'%5d\n',ITI);
TIS = 86400.*[1 5 10 20 30 60 90];
fprintf(fid,'%10.4E%10.4E%10.4E%10.4E%10.4E%10.4E%10.4E\n',TIS);

%% GOFT
fprintf(fid,'GOFT----1----*----2----*----3----*----4----*----5----*----6----*----7----*----8\n');
for i = 1981:2025
    fprintf(fid,'%08d\n',i)
end
%%
fclose(fid);