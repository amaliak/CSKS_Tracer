function creating_directories()
    % Judith 6/2/2014
    
% remove old folders
files = dir(fullfile('./','test*'));
filenames = {files.name}';

for i = 1:length(filenames)
    if exist(filenames{i})==7
    rmdir(filenames{i},'s');
    end
end
files = dir(fullfile('./','obs*'));
filenames = {files.name}';
for i = 1:length(filenames)
    if exist(filenames{i})==7
    rmdir(filenames{i},'s');
    end
end

if exist('state_update')==7
    rmdir('state_update','s');
end
if exist('obs_update')==7
    rmdir('obs_update','s');
end

% copy obs_update, state_update, PBSParallelPredict.sh , PBSParallel.sh to current directory
tmpPATH = './tmp/';
mkdir('state_update');
mkdir('obs_update');
copyfile([tmpPATH,'*'],'./');

% % establish folders
% for i=1:1:k    
%     % for parallel_update
%     filename = [ 'test' num2str(i)];
    
%     if exist(filename)==7
%     rmdir(filename,'s');
%     end

%     mkdir(filename);
%     copyfile([tmpPATH,'state_update/*'],filename);
%     %for parallel_observation
%     filename = [ 'obs' num2str(i)];
%     if exist(filename)==7
%     rmdir(filename,'s');
%     end
    
%     mkdir(filename);
%     copyfile([tmpPATH,'obs_update/*'],filename);
% end

disp(['removed all directories']);
disp('created state_update and obs_update');
end




