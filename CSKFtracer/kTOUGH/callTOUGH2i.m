function callTOUGH2i(k)
delete('FOFT_P.*','MESHA','MESHB','GENER');
system(sprintf('mpiexec -hostfile node%d /home/yuel/CSKFparallel/FKF/state_update/tough2-mp-eco2n.debug',k))
end