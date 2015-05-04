function eSol = TOUGH2update(ieSol,param)
% go the main folder ( best esimate of Sol, INFILE)

% INPUT: ieSol is a structure containing the initial
% condition 
% inputfile is an id used to indicate which iINFILE to be used if multiple
% INFILES

% OUTPUT: eSol is a structure containing the updated state after the
% forward model has run

%%%%%%%%%%%%%
% write eSol(initial condition:eSol.s(:,1) eSol.pressure) into INCON
%keyboard;
if param.init == 1
	phase = param.phase ; 
else
	% at INV mode, when touGH updates sol, it runs the previous period since k corresponds to which data are being used, not which period is being estimated
	phase = param.phase% - 1;
end
cd state_update/
delete('INFILE','FOFT*','GENER','TABLE','MESHA','MESHB','INCON','OUTPUT_DATA','OUTPUT');
% write INCON file only after initial phase has been run
%writeINCONEOS1(eSol,param);

if param.init == 0
% only write INCON and MESH if in INV mode 
	writeINCONEOS1X(ieSol,param,param.phase);
	%qwriteINCON(ieSol,param,param.phase); 
	Gmod=writeMESHPMXfast('../tmp/state_update/MESH','./MESH',ieSol.pmx(1:end),param.G0);
	%Gmod = qwriteMESH('../tmp/state_update/MESH','./MESH',ieSol.pmx);
else
	% copy INCON from previous period
	if phase > 1; copyfile('SAVE','INCON'); end;
end
% choose appropriate INFILE for appropriate data assimilation period
copyfile(['INFILE',num2str(phase)],'INFILE');

% callTOUGH2()
param.f()
% read TOUGH2 output file for updated state variables

eSol = readSAVEEOS1X(ieSol,param);
%eSol = qreadSAVE();

eSol.pmx = ieSol.pmx;
eSol.transform = @transform_sol_tr;
eSol.vec = eSol.transform(eSol,param);

cd ../

end
