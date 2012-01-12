#!/bin/sh

function write_line {
	echo $1 | cat >> $2.txt
}

BUILD_FILE="build_log"
TIME_FILE="timing"
LOF_FILE="log"
rm -f ${BUILD_FILE}.txt
rm -f ${TIME_FILE}.txt
rm -f ${LOG_FILE}.txt

for thread_method in grid sweep; do
	export THREADING_METHOD=$thread_method

	for grid_max in 60 120; do
		export MAX1D=$grid_max
		# We have to remake every time we change THREADING_METHOD or MAX1D
		write_line "--------------------------------" ${BUILD_FILE}
		write_line "thread_method = ${THREAD_METHOD}" ${BUILD_FILE}
		write_line "max1d = ${grid_max}" ${BUILD_FILE}
		make new >> ${BUILD_FILE}.txt 
		
	    for num_threads in $(seq 1 ${1}) 
	    do
			write_line "--------------------------------" ${TIME_FILE}
			write_line "thread_method = ${THREAD_METHOD}" ${TIME_FILE}
			write_line "max1d = ${grid_max}" ${TIME_FILE}
			write_line "Threads = ${num_threads}" ${TIME_FILE}
			write_line "--------------------------------" ${LOG_FILE}
			write_line "thread_method = ${THREAD_METHOD}" ${LOG_FILE}
			write_line "max1d = ${grid_max}" ${LOG_FILE}
			write_line "Threads = ${num_threads}" ${LOG_FILE}
	        export OMP_NUM_THREADS=$num_threads
			if [ `uname` == 'Darwin' ]; then
				# Mac OS X does not include the extra info from time (not GNU time)
		        /usr/bin/time make output > ${LOG_FILE}.txt
			else
				/usr/bin/time --output=${TIME_FILE}.txt --append --verbose make output > ${LOG_FILE}.txt
			fi
	    done 
	done
done

