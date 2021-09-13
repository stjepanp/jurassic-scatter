make clean -C ../../../src
make libutil.a -C ../../../src
##objcopy --prefix-symbols=foo_ ../../../src/libutil.a

sleep 2
make clean -C ../../src
make -C ../../src

rm out
#n=$(< aux/submission_index)
#m=$(( n + 1 ))
#echo $m > aux/submission_index
#echo last: $n, new: $m

#src=../../src
#w1=785
#w2=798
### run forward model
#srun --account=slmet --partition=gpus --gres=gpu:1 --nodes=1 --tasks-per-node=1 --cpus-per-task=1 $src/formod cloud-${w1}-${w2}.ctl obs33.tab atm.tab submissions/rad-${m}.tab AEROFILE aero.tab DIRLIST aux/first_001

read -n 1 -s -r -p "Press any key to continue"
sbatch jurun-ice-785.sh

watch -n 1 squeue -u pozgaj1

python3 diff.py
read -n 1 -s -r -p "Press any key to continue"
less out
