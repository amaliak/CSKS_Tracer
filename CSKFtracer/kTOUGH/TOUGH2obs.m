function obs = TOUGH2obs(eSol,param)
% go the main folder ( best esimate of Sol, INFILE)

% INPUT: eSol is a structure containing the initial
% condition for pressure and saturation
% eSol is the solution of the PREVIOUS step
% it needs to correspond to the whole TOUGH2 domain

% running TOUGH2 to propagate measurements to time k+1

%assert(max(eSol.s(:,1))<1+eps && min(eSol.s(:,1)) > -eps);
cd obs_update/
% if in INIT mode, do not run TOUGH2 again, if in INV mode, run TOUGH2
if param.init
	% Solution is at the current step
	copyfile ../state_update/SAVE ./
	Sol = eSol;
else 
	% Solution will be obtained at the next step
	delete('FOFT*','GENER','TABLE','MESHA','MESHB','INCON','OUTPUT_DATA','OUTPUT');
	% write initial conditionsc
	writeINCONEOS1X(eSol,param,param.phase+1);
	%qwriteINCON(eSol,param,param.phase+1);
	
	% create permeability field
	Gmod=writeMESHPMXfast('../tmp/obs_update/MESH','./MESH',eSol.pmx,param.G0);
	%Gmod=qwriteMESH('../tmp/obs_update/MESH','./MESH',eSol.pmx);
	% copy appropriate INFILE
	%TOUGH2obs in INV phase runs the next phase that's why we dont do phase-1 like on TOUGH2_update 

	switchINFILE(param.phase+1);
	% callTOUGH2()
	param.f();

	% obs = exp6_obs();
	% get obervations from TOUGH2 output
	Sol = readSAVEEOS1X();
	%Sol = qreadSAVE();
	% Sol only has P,T,elem, por, missiong pmxi
	% that's fine because Sol, is only used to extract pressure, X and pmx then discarded
	Sol.pmx = eSol.pmx ; 
end
% note: param.h = extract_obs_p(eSol)
obs = param.h(Sol); % function to extract observations from Sol structure
% convert pmx to logk
obs.vec = [obs.pressure; obs.X; log(obs.pmx) + log(param.pm0)];
cd ../          % uncomment it
end
