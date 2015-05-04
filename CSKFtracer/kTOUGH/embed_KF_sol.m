function big_new = embed_KF_sol(small,big_old,x,y,z,savearg)

% This function is used to embed the  solution from the update step of the 
% KF corresponding to the smaller
% domain, into the bigger domain. The big new vector generated can then be
% used as the initial conditions to rerun the forecast part of the KF

% find line numbers corresponding to small solutio

% big and small are vectors, not structures!
% if x,y,z are param.x,y,z then savearg=1 and coord2line does not have to be rerun every time
if ~savearg
lineno=coord2line(x,y,z);
else
load lineno.mat
end

big_new=big_old;

big_new(lineno-1) = small; 

% This needs to be done for eSol.pressure and then to be written in INCON
% and also needs to be done for eSol.pmx and then written in MESH

