function writing_processors_name(param,r)
%nodes: number of computer nodes (24 ppn)    
%ppn: processor per node, default is 24
%pps: processor per simulation
%r : number of simulation

if param.nodes*param.ppn < param.pps*r
error('Not enough processor!')
end

% system(['qsub -I -lnodes=',num2str(nodes),':ppn=',num2str(ppn)]);
delete('node')
system('cat $PBS_NODEFILE >node');

fid = fopen('node');
k=1;
% while ~feof(fid)
    % assign processor to the kth node
    filename = [ 'node' num2str(k)];
    fid1 =fopen(filename, 'w');
    
    for i=1:1:param.pps   % this is number of processor for each simulation
        t=fgetl(fid);
        tmp=sscanf(t,'%s');
        fprintf(fid1, '%s\n',tmp);
    end
    fclose(fid1);
    
    % Foldername=sprintf('test%d',k);
    % Foldernameobs=sprintf('obs%d',k);
    % copyfile(filename,Foldername);
    % if k==1 % assign the same number of node
        copyfile(filename,'state_update');
        copyfile(filename,'obs_update');
    % end
    % movefile(filename,Foldernameobs);
    
    % k=k+1;
    % if k>r
        % break;
    % end
% end

fclose(fid);

disp(sprintf('the total # of processors available is %d \n',param.nodes*param.ppn));
disp(sprintf('%d of nodes are assigned to each simulation\n',param.pps));
% disp(sprintf('A maximum of %d simulations can be ran in parallel\n',r));

end


