#!/usr/bin/env python
r"""Run AMRClaw thread tests"""

import sys
import os
import subprocess
import time

import setrun

if len(sys.argv) < 2:
    print "Must provide number of threads upper bound."
    print sys.argv
    sys.exit(2)

# Open and initialize log files
build_file_path = "./build_log.txt"
time_file_path = "./timings.txt"
log_file_path = "./log.txt"
build_file = open(build_file_path,'w')
time_file = open(time_file_path,"./timing.txt",'w')
log_file = open(log_file_path,"./log.txt",'w')
tm = time.localtime()
year = str(tm[0]).zfill(4)
month = str(tm[1]).zfill(2)
day = str(tm[2]).zfill(2)
hour = str(tm[3]).zfill(2)
minute = str(tm[4]).zfill(2)
second = str(tm[5]).zfill(2)
date = 'Started %s/%s/%s-%s:%s.%s\n' % (year,month,day,hour,minute,second)
build_file.write(date)
time_file.write(date)
log_file.write(date)

# Construct run command
RUNCLAW_CMD = "python $CLAWUTIL/src/python/runclaw.py xamr ./_output/ T F"
if os.uname()[0] == 'Darwin':
    time_file.write(" *** WARNING *** On Darwin timings are located in the log file (%s)." % log_file_path)
    run_cmd = "/usr/bin/time %s" % RUNCLAW_CMD
else:
    cmd = "/usr/bin/time --output=%s --append --verbose %s" % (time_file_path,RUNCLAW_CMD)
    

# Generate base data
amrdata = setrun.setrun().clawdata
amrdata.mx = 40
amrdata.my = 10
amrdata.verbosity = 0
amrdata.mxnest = 1

for thread_method in ["grid","sweep"]:
    os.environ["THREADING_METHOD"] = thread_method

    for grid_max in [60]:
        os.environ["MAX1D"] = str(grid_max)
        build_file.write("--------------------------------\n")
        build_file.write("thread_method = %s\n" % thread_method)
        build_file.write("max1d = %s\n" % grid_max)
        build_file.flush()
        subprocess.Popen("make new\n",shell=True,stdout=build_file).wait()

        for num_threads in range(0,int(sys.argv[1])):
            time_file.write("--------------------------------\n")
            time_file.write("thread_method = %s\n" % thread_method)
            time_file.write("max1d = %s\n" % grid_max)
            time_file.write("threads = %s\n" % (num_threads+1))
            log_file.write("--------------------------------\n")
            log_file.write("thread_method = %s\n" % thread_method)
            log_file.write("max1d = %s\n" % grid_max)
            log_file.write("threads = %s\n" % (num_threads+1))
            log_file.flush()
            time_file.flush()
            log_file.flush()

            os.environ["OMP_NUM_THREADS"] = str(num_threads + 1)
            subprocess.Popen(cmd,shell=True,stdout=log_file).wait()