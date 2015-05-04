function eSols = TOUGH2update_parallel(eSols,r,param)
% go the main folder ( best esimate of Sol, INFILE)

% INPUT: eSol is a structure containing the initial
% condition for pressure and saturation

% call simulator to make prediction
%%%%%%%%assert(max(s)<1+eps && min(s) > -eps);

% for k = 1:r
%     foldernamek = ['testfolder',num2str(k)];
%     cd(foldernamek);
%     % write state varibles to INCON file
%     writeINCON(eSols{k});
%     callTOUGH2s();
%     eSols{k} = readSAVE();
%     cd ../
% end
% param.pps = 4;  % change later

% % remove all the testfolder* folder
% files = dir(fullfile('./','testfolder*'));
% filenames = {files.name}';

% for i = 1:length(filenames)
%     if exist(filenames{i})==7 % a directory
%     rmdir(filenames{i},'s');
%     end
% end
phase = param.phase;
fclose all
 
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

% update INCON file
for k = 1:r
    foldernamek = ['testfolder',num2str(k)];
    if exist(foldernamek)==7
    	rmdir(foldernamek,'s');
    end
    mkdir(foldernamek);
    copyfile([tmpPATH,'state_update/*'],foldernamek);
    cd(foldernamek);
    % write processer_name to nodek
    filename = [ 'node' num2str(k)];
    writeNodefile(filename,fid);
    % write state varibles to INCON file
    if exist('INCON')==2
    	delete('INCON')
    end
    writeINCONEOS1X(eSols{k},param,param.phase);
    %qwriteINCON(eSols{k},param,param.phase);
    Gmod=writeMESHPMXfast('../tmp/state_update/MESH','./MESH',eSols{k}.pmx,param.G0);
    %Gmod = qwriteMESH('../tmp/state_update/MESH','./MESH',eSols{k}.pmx);
    pmx_tmp{k} = eSols{k}.pmx;
    switchINFILE(phase);
    cd ../
end

% call TOUGH2 in parallel
% if t1 is small, may indicate: 1. INCON contains NaN 2. node file is broken
% if t1 is very long, there would be convergence issue in specific folder
tic;
filename = 'PBSParallelPredict.sh';
workdir = pwd; % is it necessary?
writePBSfile(filename,workdir,r);
setenv('LD_LIBRARY_PATH','/share/apps/intel/intel-14/lib/intel64')
system(['sh ',filename]);
t1 = toc;
disp(['time for TOUGH2update_parallel is',num2str(t1)])

%% read TOUGH2 output file
% save each SAVE to eSol.pressure, eSol.s
for k = 1:r
    foldernamek = ['testfolder',num2str(k)];
    cd(foldernamek);
    eSols{k} = readSAVEEOS1X();
    %eSols{k} = qreadSAVE();
    eSols{k}.pmx = pmx_tmp{k};
    eSols{k}.transform = @transform_sol_tr;
    eSols{k}.vec = eSols{k}.transform(eSols{k},param);
    cd ../
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%functions%%%%%%%%%%%%%%%%%%%%%%%%%%%
	function writePBSfile(filename,workdir,r)
		fid = fopen(filename,'w');
		fprintf(fid,'PBS_O_WORKDIR=''%s''\n',workdir);
		fprintf(fid,'         export PBS_O_WORKDIR\n');
		fprintf(fid,'\n');
		fprintf(fid,'for ((i=1;i<=%d;i++))\n',r);
		fprintf(fid,'do\n');
		fprintf(fid,'cd testfolder$i/\n');
		%fprintf(fid,'cmd="mpiexec -hostfile node$i /home/yuel/CSKFparallel/FKF/state_update/tough2-mp-eco2n.debug"\n');
        	fprintf(fid,'cmd="mpiexec -hostfile node$i -env MV2_ENABLE_AFFINITY=0 /home/amaliak/eos1/tough2-mp-eos1.debug"\n');
		fprintf(fid,'echo "Running MPIEXEC with: $cmd in directory "`pwd`\n');
		if strcmp(pwd,workdir) == 1
		    mkdir('log');
		else
		    Error('Not in the working directory');
		end
		fprintf(fid,'$cmd >& $PBS_O_WORKDIR/log/TOUGHtest$i.log &\n');
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
