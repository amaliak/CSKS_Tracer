PBS_O_WORKDIR='/home/amaliak/Tracer/INV'
         export PBS_O_WORKDIR

for ((i=1;i<=30;i++))
do
cd obsfolder$i/
cmd="mpiexec -hostfile node$i -env MV2_ENABLE_AFFINITY=0 /home/amaliak/eos1/tough2-mp-eos1.debug"
echo "Running MPIEXEC with: $cmd in directory "`pwd`
$cmd >& $PBS_O_WORKDIR/log/TOUGHobs$i.log &
cd ..
done
wait
