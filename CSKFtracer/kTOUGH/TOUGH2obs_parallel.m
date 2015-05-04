function obs1 = TOUGH2obs_parallel(eSols,param,r)
obs1 = cell(r,1);
phase = param.phase; 
if exist('tmp')==7
tmpPATH = './tmp/';
else
    disp('No tmp folder');
    keyboard;
end

if exist('node')~=2
    disp('There is no node file!');
    keyboard;
end

fid = fopen('node');
% write eSol(initial condition:eSol.s(:,1) eSol.pressure) into INCON
for j = 1:1:r
    foldernamek=['obsfolder',num2str(j)];
    if exist(foldernamek)==7
        rmdir(foldernamek,'s');
    end
    mkdir(foldernamek);
    copyfile([tmpPATH,'obs_update/*'],foldernamek);
    cd(foldernamek);
    % write processer_name to nodek
    filename = [ 'node' num2str(j)];
    writeNodefile(filename,fid);
    % write INCON
    if exist('INCON')==2
        delete('INCON')
    end
    writeINCONEOS1X(eSols{j},param,param.phase);
    %qwriteINCON(eSols{j},param,param.phase);
    Gmod=writeMESHPMXfast('../tmp/obs_update/MESH','MESH',eSols{j}.pmx,param.G0);
    %Gmod = qwriteMESH('../tmp/obs_update/MESH','MESH',eSols{j}.pmx);
    % copies INFILE corresponding to run that ends at the time of the measurements being assimulated
    switchINFILE(phase);
    cd ../;
end
%keyboard;
tic;
filename = 'PBSParallel.sh';
workdir = pwd; % is it necessary?
writePBSfile(filename,workdir,r);
setenv('LD_LIBRARY_PATH','/share/apps/intel/intel-14/lib/intel64')
system(['sh ',filename]);
t1 = toc;
disp(['time for TOUGH2obs_parallel is',num2str(t1)]);
% system('sh PBSParallel.sh')
%% read TOUGH2 output file
% save each SAVE to eSol.pressure, eSol.s

for k = 1:r
    foldernamek = ['obsfolder',num2str(k)];
    cd(foldernamek);
    Sol = readSAVEEOS1X();
    %Sol = qreadSAVE();
    Sol.pmx = eSols{k}.pmx;% careful! the vec is not changed 
    obs = param.h(Sol);
    obs.vec = [obs.pressure;obs.X; log(obs.pmx) + log(param.pm0)];    
    obs1{k} = obs;
    cd ../
end
    function writePBSfile(filename,workdir,r)
        fid = fopen(filename,'w');
        fprintf(fid,'PBS_O_WORKDIR=''%s''\n',workdir);
        fprintf(fid,'         export PBS_O_WORKDIR\n');
        fprintf(fid,'\n');
        fprintf(fid,'for ((i=1;i<=%d;i++))\n',r);
        fprintf(fid,'do\n');
        fprintf(fid,'cd obsfolder$i/\n');
        fprintf(fid,'cmd="mpiexec -hostfile node$i -env MV2_ENABLE_AFFINITY=0 /home/amaliak/eos1/tough2-mp-eos1.debug"\n');
        fprintf(fid,'echo "Running MPIEXEC with: $cmd in directory "`pwd`\n');
        if strcmp(pwd,workdir) == 1
            mkdir('log');
        else
            Error('Not in the working directory');
        end
        fprintf(fid,'$cmd >& $PBS_O_WORKDIR/log/TOUGHobs$i.log &\n');
        fprintf(fid,'cd ..\n');
        fprintf(fid,'done\n');
        fprintf(fid,'wait\n');
        fclose(fid);
        disp(sprintf('%s is created successfully',filename));
    end

    function writeNodefile(filename,fid)
        fid1 =fopen(filename, 'w');
        for i=1:1:param.pps   % this is number of processor for each simulation
            t=fgetl(fid);
            tmp=sscanf(t,'%s');
            fprintf(fid1, '%s\n',tmp);
        end
        fclose(fid1);
    end
end
%% Run sequentially
% for k = 1:r
%     foldernamek = ['obsfolder',num2str(k)];
%     cd(foldernamek);
%     % write state varibles to INCON file
%     delete('INCON');
%     writeINCON(eSols{k});
%     callTOUGH2i(k);
%     obs = h();
%     obs1{k} = obs;
%     cd ../
% end
