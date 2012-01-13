#!/usr/bin/env python
# encoding: utf-8
r"""Run AMRClaw thread tests

p = num_threads

Single Grid Tests
=================
  Single grid tests of each approach
  
  Sweep Threading
  ---------------
    Vary (mx,my) and threads, mx = 40, 60, 80, 100, 120

  Grid Threading
  --------------
    Cut single grid up into even multiples forcing a specific number of grids
                mx
    |-----|-----|-----|-----|
    |     |     |     |     |  my
    |-----|-----|-----|-----|
    
    EWT = mx*my / p + 2 * my
    ECT = mx*my + 2*(p-1) * my
    
Adaptive Grid Tests
===================
  Tests for both sweep and grid threading and for all p
  
"""

import sys
import os
import subprocess
import time

import setrun

# Open and initialize log files
build_file_path = "./build_log.txt"
time_file_path = "./timings.txt"
log_file_path = "./log.txt"
build_file = open(build_file_path,'w')
time_file = open(time_file_path,'w')
log_file = open(log_file_path,'w')

# Construct run command
RUNCLAW_CMD = "python $CLAWUTIL/src/python/runclaw.py xamr ./_output/ T F"
if os.uname()[0] == 'Darwin':
    time_file.write(" *** WARNING *** On Darwin timings are located in the log file (%s)." % log_file_path)
    run_cmd = "/usr/bin/time %s" % RUNCLAW_CMD
else:
    cmd = "/usr/bin/time --output=%s --append --verbose %s" % (time_file_path,RUNCLAW_CMD)

class BaseThreadTest(object):
    
    def __init__(self,name,mx=40,mxnest=3,grid_max=60,thread_method="grid",max_threads=2):
        self.name = name
        self.grid_max = grid_max
        self.thread_method = thread_method
        self.max_threads = max_threads
        
        # Generate base data
        self.amrdata = setrun.setrun().clawdata
        self.amrdata.mx = mx
        self.amrdata.my = mx / 4
        self.amrdata.verbosity = 0
        self.amrdata.mxnest = mxnest

        dt = 0.0001
        tfinal = 0.75
        self.amrdata.dt_initial = dt
        self.amrdata.dt_variable = 0
        
        self.amrdata.outstyle = 3
        num_steps = int(tfinal / dt) + 1
        self.amrdata.max_steps = num_steps + 1
        self.amrdata.iout = [num_steps,num_steps]
        
    def __str__(self):
        output = "Test name: %s" % self.name
        output += "  mx = %s" % self.amrdata.mx
        output += "  mxnest = %s" % self.amrdata.mxnest
        output += "  max1d = %s" % self.grid_max
        output += "  thread_method = %s" % self.thread_method
        output += "  thread number = [1,%s]" % self.max_threads
        return output
        
    def write_info(self,out_file,num_threads):
        tm = time.localtime()
        year = str(tm[0]).zfill(4)
        month = str(tm[1]).zfill(2)
        day = str(tm[2]).zfill(2)
        hour = str(tm[3]).zfill(2)
        minute = str(tm[4]).zfill(2)
        second = str(tm[5]).zfill(2)
        date = 'Started %s/%s/%s-%s:%s.%s\n' % (year,month,day,hour,minute,second)
        out_file.write("--------------------------------\n")
        out_file.write(date)
        out_file.write(str(self))
        # out_file.write("thread_method = %s\n" % self.thread_method)
        # out_file.write("max1d = %s\n" % self.grid_max)
        # out_file.write("threads = %s" % num_threads)
        
    def flush_log_files(self):
        build_file.flush()
        log_file.flush()
        time_file.flush()
        
    def build_test(self):
        os.environ["THREADING_METHOD"] = self.thread_method
        os.environ["MAX1D"] = str(self.grid_max)
        self.flush_log_files()
        subprocess.Popen("make new -j %s" % self.max_threads,shell=True,stdout=build_file).wait()
        
    def run_tests(self):
        tm = time.localtime()
        year = str(tm[0]).zfill(4)
        month = str(tm[1]).zfill(2)
        day = str(tm[2]).zfill(2)
        hour = str(tm[3]).zfill(2)
        minute = str(tm[4]).zfill(2)
        second = str(tm[5]).zfill(2)
        date = 'Started %s/%s/%s-%s:%s.%s\n' % (year,month,day,hour,minute,second)
        build_file.write("--------------------------------\n")
        build_file.write(date)
        build_file.write("thread_method = %s\n" % self.thread_method)
        build_file.write("max1d = %s\n" % self.grid_max)
        log_file.write("Building...")
        self.build_test()
        log_file.write("Build completed.")

        # Write out data file
        self.amrdata.write()
        
        for num_threads in range(0,self.max_threads):
            os.environ["OMP_NUM_THREADS"] = str(num_threads + 1)
            self.write_info(log_file,num_threads + 1)
            self.write_info(time_file,num_threads + 1)
            
            # Finally, run test
            log_file.write("Running simulation with %s threads..." % int(num_threads + 1))
            self.flush_log_files()
            subprocess.Popen(RUNCLAW_CMD,shell=True,stdout=log_file).wait()
            log_file.write("Simulation completed.")
            

tests = []
max_threads = int(os.environ['OMP_NUM_THREADS'])

# Single grid sweep timings
for mx in [40,60,80,100,120]:
    tests.append(BaseThreadTest("Single Grid Sweep Threading, grid_max=%s" % mx,mxnest=1,mx=mx,grid_max=mx,thread_method="sweep",max_threads=max_threads))

if __name__ == "__main__":
    if len(sys.argv) > 1:
        if sys.argv[1].lower() == 'all':
            tests_to_be_run = tests
        else:
            tests_to_be_run = []
            for test in sys.argv[1:]:
                tests_to_be_run.append(tests[int(test)])
    else:
        for test in tests:
            print test
        sys.exit(0)
    
    
    # Execute tests
    for (i,test) in enumerate(tests_to_be_run):
        log_file.write("Running test %s, test %s of %s." % (test.name,int(i+1),int(len(tests_to_be_run))))
        test.run_tests()
