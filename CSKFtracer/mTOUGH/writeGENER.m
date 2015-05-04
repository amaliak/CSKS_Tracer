fid = fopen('GENER','w');
fprintf(fid,'GENER----1----*----2----*----3----*----4----*----5----*----6----*----7----*----8\n');
table.index = [1:45 1981:1:2025];
table.no = 1:90;
% %% MOP(14) = 1, 8 characters
% % injection well with BHP = 220 bar at left boundary
% for i = 1:45
%     well.type = 'COM3';
%     well.prop = [0.05,0.0]; %BHP control well
%     fprintf(fid,'%08d  %3d       %10d% 9s %10.4E%10.4E\n',table.index(i),table.no(i),0,well.type,well.prop(1),well.prop(2));   
% end
% for i = 46:90
%     well.type = 'DELV';
%     well.prop = [6.73e-10,2.0e7]; %BHP control well
%     fprintf(fid,'%08d  %3d       %10d% 9s %10.4E%10.4E\n',table.index(i),table.no(i),0,well.type,well.prop(1),well.prop(2));
% end
% fprintf(fid,'\n');
% 
% % FOFT
% fprintf(fid,'FOFT ----1----*----2----*----3----*----4----*----5----*----6----*----7----*----8\n');
% for i = 1:90
%     fprintf(fid,'%08d\n',table.index(i));
% end
% fprintf(fid,'\n');
% 
% fprintf(fid,'GOFT ----1----*----2----*----3----*----4----*----5----*----6----*----7----*----8\n');
% % GOFT
% for i = 1:90
%     fprintf(fid,'%08d\n',table.index(i));
% end
% fprintf(fid,'\n');
% fclose(fid);
%% MOP(14) = 0, 5 characters
%Format (A3, I2, A3, I2, 4I5, 5X, A4, A1, 3E10.4)
%EL, NE, SL, NS, NSEQ, NADD, NADS, LTAB, TYPE, ITAB, GX, EX, HX
% injection well with BHP = 220 bar at left boundary

for i = 1:45
    well.type = 'COM3';
    well.prop = [0.05,0.0]; %BHP control well
    fprintf(fid,'A%04dinj%2d% 29s %10.4E%10.4E\n',table.index(i),table.no(i),well.type,well.prop(1),well.prop(2));   
end
for i = 46:90
    well.type = 'DELV';
    well.prop = [6.73e-10,2.0e7]; %BHP control well
    fprintf(fid,'A%04dpro%2d% 29s %10.4E%10.4E\n',table.index(i),table.no(i),well.type,well.prop(1),well.prop(2));   
end
fprintf(fid,'\n');

% FOFT
fprintf(fid,'FOFT ----1----*----2----*----3----*----4----*----5----*----6----*----7----*----8\n');
for i = 1:90
    fprintf(fid,'A%04d\n',table.index(i));
end
fprintf(fid,'\n');

fprintf(fid,'GOFT ----1----*----2----*----3----*----4----*----5----*----6----*----7----*----8\n');
% GOFT
for i = 1:90
    fprintf(fid,'A%04d\n',table.index(i));
end
fprintf(fid,'\n');
fclose(fid);
