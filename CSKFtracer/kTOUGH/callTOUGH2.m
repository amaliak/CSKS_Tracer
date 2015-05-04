function callTOUGH2()
delete('FOFT_P.*','MESHA','MESHB','GENER');
setenv('LD_LIBRARY_PATH','/share/apps/intel/intel-14/lib/intel64')
system('mpiexec -hostfile node1 /home/amaliak/eos1heat/tough2-mp-eos1q.debug')
end
