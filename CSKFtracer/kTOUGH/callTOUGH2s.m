function callTOUGH2s()
delete('FOFT_P.*','MESHA','MESHB','GENER');
system('mpiexec /home/yuel/CSKFparallel/FKF/state_update/tough2-mp-eco2n.debug')
end