make clean -C ../../jurassic-gpu/src/
make libutil.a -C ../../jurassic-gpu/src/

sleep 2
make clean -C ../../src
make -C ../../src

rm out
n=$(< aux/submission_index)
m=$(( n + 1 ))
echo $m > aux/submission_index
echo last: $n, new: $m

read -n 1 -s -r -p "Press any key to continue"
sbatch jurun-ice-785.sh ../../src $m

watch -n 1 squeue -u pozgaj1

python3 diff.py
read -n 1 -s -r -p "Press any key to continue"
less out
