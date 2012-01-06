#!/bin/sh

TIME_FILE="timing_${2}"
rm -f $TIME_FILE.txt

for thread_method in GRID_THREADING SWEEP_THREADING; do
	export THREADING_METHOD=$thread_method

	for grid_max in 60 120 180 240; do
		export MAX1D=$grid_max
	
	    for num_threads in $(seq 1 ${1}) 
	    do
	        echo "--------------------------------" | cat >> ${TIME_FILE}.txt
			echo "thread_method = ${THREAD_METHOD}" | cat >> ${TIME_FILE}.txt
			echo "max1d = ${grid_max}" | cat >> ${TIME_FILE}.txt
	        echo "Threads = ${num_threads}" | cat >> ${TIME_FILE}.txt
	        export OMP_NUM_THREADS=$num_threads
	        /usr/bin/time --output=${TIME_FILE}.txt --append --verbose make output > log.txt
	    done 
	done
done