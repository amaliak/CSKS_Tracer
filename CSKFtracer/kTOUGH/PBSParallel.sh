PBS_O_WORKDIR='/home/yuel/CSKFparallel'
         export PBS_O_WORKDIR

for ((i=1;i<=20;i++))
do
cd obs$i/
cmd="mpiexec -hostfile node$i tough2-mp-eco2n.debug"
echo "Running MPIEXEC with: $cmd in directory "`pwd`
$cmd >& $PBS_O_WORKDIR/TOUGH2.log &
cd ..
done
wait
