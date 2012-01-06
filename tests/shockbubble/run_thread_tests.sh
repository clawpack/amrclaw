#!/bin/sh

TIME_FILE="timing_${2}"
rm -f $TIME_FILE.txt

for num_threads in $(seq 1 ${1}) 
do
    echo "--------------------------------" | cat >> ${TIME_FILE}.txt
    echo "Threads = ${num_threads}" | cat >> ${TIME_FILE}.txt
    export OMP_NUM_THREADS=$num_threads
    /usr/bin/time --output=${TIME_FILE}.txt --append --verbose make output > log.txt
done 
