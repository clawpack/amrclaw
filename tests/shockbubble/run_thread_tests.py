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

def make_path(path):
    if not os.path.exists(path):
        os.makedirs(path)
    elif not os.path.isdir(path):
        os.remove(path)
        os.makedirs(path)


# Open and initialize log files
BUILD_PATH_BASE = "./logs"
TIME_PATH_BASE = "./logs"
LOG_PATH_BASE = "./logs"

for path in [BUILD_PATH_BASE,TIME_PATH_BASE,LOG_PATH_BASE]:
    make_path(path)

# Construct run command based on log
def construct_run_cmd(time_file_path):
    RUNCLAW_CMD = "python $CLAWUTIL/src/python/runclaw.py xamr ./_output/ T F"
    if os.uname()[0] == 'Darwin':
        # raise exceptions.Warning("On Darwin timings are located in the log file (%s).\n" % time_file_path)
        run_cmd = "/usr/bin/time %s" % RUNCLAW_CMD
    else:
        run_cmd = "/usr/bin/time --output=%s --append --verbose %s" % (time_file_path,RUNCLAW_CMD)
    return run_cmd

class BaseThreadTest(object):
    
    def __init__(self,name,mx=40,mxnest=3,grid_max=60,thread_method="grid",max_threads=2):
        # Basic test information
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
        
        # Open log files and write headers
        file_label = "_%s_m%s_gm%s_%s" % (self.name,mx,self.grid_max,self.thread_method)
        self._build_file_path = os.path.join(BUILD_PATH_BASE,"build%s.txt" % file_label)
        self._time_file_path = os.path.join(TIME_PATH_BASE,"time%s.txt" % file_label)
        self._log_file_path = os.path.join(LOG_PATH_BASE,"log%s.txt" % file_label)
        self.build_file = open(self._build_file_path,'w')
        self.time_file = open(self._time_file_path,'w')
        self.log_file = open(self._log_file_path,'w')
        self.write_info(self.build_file,self.max_threads)
        self.write_info(self.time_file,self.max_threads)
        self.write_info(self.log_file,self.max_threads)
        
    def __str__(self):
        output = "Test name: %s\n" % self.name
        output += "  mx = %s\n" % self.amrdata.mx
        output += "  mxnest = %s\n" % self.amrdata.mxnest
        output += "  max1d = %s\n" % self.grid_max
        output += "  thread_method = %s\n" % self.thread_method
        output += "  thread number = [1,%s]\n" % self.max_threads
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
        self.build_file.flush()
        self.log_file.flush()
        self.time_file.flush()
        
    def build_test(self):
        os.environ["THREADING_METHOD"] = self.thread_method
        os.environ["MAX1D"] = str(self.grid_max)
        self.flush_log_files()
        subprocess.Popen("make new -j %s" % self.max_threads,shell=True,stdout=self.build_file,stderr=self.build_file).wait()
        
    def run_tests(self):
        # tm = time.localtime()
        # year = str(tm[0]).zfill(4)
        # month = str(tm[1]).zfill(2)
        # day = str(tm[2]).zfill(2)
        # hour = str(tm[3]).zfill(2)
        # minute = str(tm[4]).zfill(2)
        # second = str(tm[5]).zfill(2)
        # date = 'Started %s/%s/%s-%s:%s.%s\n' % (year,month,day,hour,minute,second)
        # build_file.write("--------------------------------\n")
        # build_file.write(date)
        # build_file.write("thread_method = %s\n" % self.thread_method)
        # build_file.write("max1d = %s\n" % self.grid_max)
        self.log_file.write("Building...\n")
        self.build_test()
        self.log_file.write("Build completed.\n")

        # Write out data file
        self.amrdata.write()
        
        for num_threads in range(0,self.max_threads):
            os.environ["OMP_NUM_THREADS"] = str(num_threads + 1)
            self.time_file.write("\n *** OMP_NUM_THREADS = %s\n\n" % int(num_threads+1))
            self.log_file.write("\n *** OMP_NUM_THREADS = %s\n\n" % int(num_threads+1))

            run_cmd = construct_run_cmd(self._time_file_path)
            self.log_file.write(run_cmd + "\n")
            self.flush_log_files()
            self.time_file.close()
            subprocess.Popen(run_cmd,shell=True,stdout=self.log_file,stderr=self.log_file).wait()
            self.log_file.write("Simulation completed.\n")
            self.time_file = open(self._time_file_path,'aw')

tests = []
max_threads = int(os.environ['OMP_NUM_THREADS'])

# Single grid sweep timings
for mx in [40,60,80,100,120]:
    tests.append(BaseThreadTest("single_grid",mxnest=1,mx=mx,grid_max=mx,thread_method="sweep",max_threads=max_threads))

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
        test.run_tests()
