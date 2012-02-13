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
import math
import glob

import setrun

# =============================================================================
# Parameters for script
if os.environ.has_key('FC'):
    LOG_PATH_BASE = "./logs_%s" % os.environ['FC']
else:
    LOG_PATH_BASE = "./logs"
    
# End time for all simulations
TFINAL = 0.75

# =============================================================================
# Construct run command based on log
def construct_run_cmd(time_file_path):
    RUNCLAW_CMD = "python $CLAWUTIL/src/python/runclaw.py xamr ./_output/ T F"
    if os.uname()[0] == 'Darwin':
        # raise exceptions.Warning("On Darwin timings are located in the log file (%s).\n" % time_file_path)
        run_cmd = "/usr/bin/time %s" % RUNCLAW_CMD
    else:
        run_cmd = "/usr/bin/time --output=%s --append --verbose %s" % (time_file_path,RUNCLAW_CMD)
    return run_cmd

# =============================================================================
#  Base test class
class BaseThreadingTest(object):
    
    def __init__(self,name,threads,thread_method,mx,mxnest,grid_max):
        # Basic test information
        self.name = name
        self.threads = threads
        self.thread_method = thread_method
        self.grid_max = grid_max

        # Generate base data
        self.setrun = setrun.setrun()
        self.amrdata = self.setrun.clawdata
        self.probdata = self.setrun.probdata
        self.amrdata.verbosity = 0

        # Construct this test case's data
        self.amrdata.mxnest = mxnest
        self.amrdata.mx = mx
        self.amrdata.my = mx / 4
        
        # Output style
        dt = 0.001
        self.amrdata.dt_initial = dt
        self.amrdata.dt_variable = 1
        self.amrdata.outstyle = 2
        self.amrdata.nout = 1
        self.amrdata.tout = [TFINAL]
        
        # File log label
        self.file_label = "_%s_m%s_g%s_n%s" % (self.name,
                                               str(self.amrdata.mx).zfill(3),
                                               str(self.grid_max).zfill(3),
                                               str(self.amrdata.mxnest).zfill(2))
        
    def __str__(self):
        output = "Test name: %s - %s\n" % (str(test.__class__).strip("<class '__main__.").strip("'>"),self.name)
        output += "  file_label = %s\n" % self.file_label
        output += "  threads    = %s\n" % self.threads
        output += "  mx         = %s\n" % self.amrdata.mx
        output += "  mxnest     = %s\n" % self.amrdata.mxnest
        output += "  grid_max   = %s\n" % self.grid_max
        return output
        
        
    def write_info(self,out_file):
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
        
    def flush_log_files(self):
        self.build_file.flush()
        self.log_file.flush()
        self.time_file.flush()
        
    def open_log_files(self):
        # Open log files and write headers
        self._build_file_path = os.path.join(LOG_PATH_BASE,"build%s.txt" % self.file_label)
        self._time_file_path = os.path.join(LOG_PATH_BASE,"time%s.txt" % self.file_label)
        self._log_file_path = os.path.join(LOG_PATH_BASE,"log%s.txt" % self.file_label)
        self.build_file = open(self._build_file_path,'w')
        self.time_file = open(self._time_file_path,'w')
        self.log_file = open(self._log_file_path,'w')
        self.write_info(self.build_file)
        self.write_info(self.time_file)
        self.write_info(self.log_file)
        
    def close_log_files(self):
        # Close log files
        self.build_file.close()
        self.time_file.close()
        self.log_file.close()
    
    def write_data(self):
        # Write out data file
        self.amrdata.write()
        self.probdata.write()
        
    def run_tests(self):
        self.open_log_files()
        
        # Build binary
        self.log_file.write("Building...\n")
        os.environ["THREADING_METHOD"] = self.thread_method
        os.environ["MAX1D"] = str(self.grid_max)
        self.flush_log_files()
        subprocess.Popen("make new -j %s" % self.threads,shell=True,
                                stdout=self.build_file,stderr=self.build_file,
                                env=os.environ).wait()
        self.log_file.write("Build completed.\n")
        
        # Run tests
        for num_threads in self.threads:
            os.environ["OMP_NUM_THREADS"] = str(num_threads)
            self.time_file.write("\n *** OMP_NUM_THREADS = %s\n\n" % int(num_threads))
            self.log_file.write("\n *** OMP_NUM_THREADS = %s\n\n" % int(num_threads))

            self.write_data()

            run_cmd = construct_run_cmd(self._time_file_path)
            self.log_file.write(run_cmd + "\n")
            self.flush_log_files()
            self.time_file.close()
            # subprocess.Popen("make .data",shell=True,stdout=self.log_file,stderr=self.log_file,
            #                     env=os.environ).wait()
            subprocess.Popen(run_cmd,shell=True,stdout=self.log_file,stderr=self.log_file,
                                env=os.environ).wait()
            self.log_file.write("Simulation completed.\n")
            self.time_file = open(self._time_file_path,'aw')

        self.close_log_files()
        
        
class SingleGridThreadingTest(BaseThreadingTest):
        
    def __init__(self,name,threads,thread_method,mx):
        
        # Set base values
        super(SingleGridThreadingTest,self).__init__(name,threads,thread_method,mx,1,mx)

        # Setup non variable time stepping and output
        dt = 0.0001
        self.amrdata.dt_initial = dt
        self.amrdata.dt_variable = 0
        
        self.amrdata.outstyle = 3
        num_steps = int(TFINAL / dt) + 1
        self.amrdata.max_steps = num_steps + 1
        self.amrdata.iout = [num_steps,num_steps]

        
class WeakThreadingTest(BaseThreadingTest):

    def __init__(self,name,threads,thread_method,mx,grid_max):
        # Set base values
        super(WeakThreadingTest,self).__init__(name,threads,thread_method,mx,1,grid_max)

        # Setup non variable time stepping and output
        dt = 0.0001
        self.amrdata.dt_initial = dt
        self.amrdata.dt_variable = 0
        
        self.amrdata.outstyle = 3
        num_steps = int(TFINAL / dt) + 1
        self.amrdata.max_steps = num_steps + 1
        self.amrdata.iout = [num_steps,num_steps]
        
    def write_data(self):
        from math import sqrt

        # Set mx and my based on the number of threads
        self.amrdata.mx = self.grid_max * sqrt(num_threads)
        self.amrdata.my = self.amrdata.mx / 4
        self.amrdata.write()
        self.probdata.write()
    
        
# =============================================================================

# =============================================================================
# Create tests
# =============================================================================
tests = []
max_threads = int(os.environ['OMP_NUM_THREADS'])

# Test ranges
threads = [1,2,4,8,12,16]
sqrt_threads = [1,4,9,16]
for (i,count) in enumerate(threads):
    if count > max_threads:
        threads = threads[:i]
        break
for (i,count) in enumerate(sqrt_threads):
    if count > max_threads:
        sqrt_threads = sqrt_threads[:i]
        break
        
single_grid_mx = [N - 4 for N in [64,128,256,512]]
grid_max_tests = [64,128,256,512]

# Single Grid Tests
# =================
#   Single grid tests of each approach
#   
#   Sweep Threading
#   ---------------
#     Vary (mx,my) and threads
for mx in single_grid_mx:
    tests.append(SingleGridThreadingTest("single_sweep",threads,'sweep',mx))
# Weak scaling test
for mx in single_grid_mx:
    tests.append(WeakThreadingTest("weak_sweep",sqrt_threads,'sweep',mx,mx))
    
# Grid Threading
# --------------
#   Cut single grid up into even multiples forcing a specific number of grids
#               mx
#   |-----|-----|-----|-----|
#   |     |     |     |     |  my
#   |-----|-----|-----|-----|
#   
#   EWT = mx*my / p + 2 * my
#   ECT = mx*my + 2*(p-1) * my
for grid_max in grid_max_tests:
    tests.append(WeakThreadingTest('weak_grid',sqrt_threads,'grid',64,grid_max))

# Adaptive Grid Tests
# ===================
#   Tests for both sweep and grid threading and for all p
for grid_max in grid_max_tests:
    tests.append(BaseThreadingTest('amr_sweep',threads,'sweep',64,3,grid_max))
for grid_max in grid_max_tests:
    tests.append(BaseThreadingTest('amr_grid',threads,'grid',64,3,grid_max))
    

# =============================================================================
#  Command line support
if __name__ == "__main__":
    
    # Construct list of tests to be run
    if len(sys.argv) > 1:
        if sys.argv[1].lower() == 'all':
            tests_to_be_run = tests
        else:
            tests_to_be_run = []
            for test in sys.argv[1:]:
                tests_to_be_run.append(tests[int(test)])
    # or just print all the tests available
    else:
        for (i,test) in enumerate(tests):
            print "==== %s ==================" % i
            print test
        sys.exit(0)
    
    # Output tests to be performed
    for (i,test) in enumerate(tests_to_be_run):
        print "==== %s ==================" % i
        print test
    
    # Create log output directory
    if os.path.exists(LOG_PATH_BASE):
        for log_file in glob.glob(os.path.join(LOG_PATH_BASE,'*.txt')):
            os.remove(log_file)
        os.rmdir(LOG_PATH_BASE)
    os.makedirs(LOG_PATH_BASE)
    
    # Execute tests
    for (i,test) in enumerate(tests_to_be_run):
        test.run_tests()
