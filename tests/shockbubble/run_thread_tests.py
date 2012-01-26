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
    
    def __init__(self,name,threads,thread_method,grid_max):
        # Basic test information
        self.name = name
        self.threads = threads
        self.thread_method = thread_method
        self.grid_max = grid_max

        # Generate base data
        self.amrdata = setrun.setrun().clawdata
        self.amrdata.verbosity = 0
        
        # Set file base name
        self.file_label = "_%s" % (name)
        
    def __str__(self):
        output = "Test name: %s - %s\n" % (str(test.__class__).strip("<class '__main__.").strip("'>"),self.name)
        output += "  file_label = %s\n" % self.file_label
        output += "  threads    = %s\n" % self.threads
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
        # Create log file directory if not there already
        if not os.path.exists(LOG_PATH_BASE):
            os.makedirs(LOG_PATH_BASE)
        elif not os.path.isdir(LOG_PATH_BASE):
            os.remove(LOG_PATH_BASE)
            os.makedirs(LOG_PATH_BASE)
        
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

        # Write out data file
        self.amrdata.write()
        
        # Run tests
        for num_threads in self.threads:
            os.environ["OMP_NUM_THREADS"] = str(num_threads)
            self.time_file.write("\n *** OMP_NUM_THREADS = %s\n\n" % int(num_threads))
            self.log_file.write("\n *** OMP_NUM_THREADS = %s\n\n" % int(num_threads))

            run_cmd = construct_run_cmd(self._time_file_path)
            self.log_file.write(run_cmd + "\n")
            self.flush_log_files()
            self.time_file.close()
            subprocess.Popen("make .data",shell=True,stdout=self.log_file,stderr=self.log_file,
                                env=os.environ).wait()
            subprocess.Popen(run_cmd,shell=True,stdout=self.log_file,stderr=self.log_file,
                                env=os.environ).wait()
            self.log_file.write("Simulation completed.\n")
            self.time_file = open(self._time_file_path,'aw')

        self.close_log_files()
        
    
    
class SweepThreadingTest(BaseThreadingTest):
    
    def __init__(self,name,threads,mx=40,mxnest=3,grid_max=60):
        
        # Call super class to set name
        super(SweepThreadingTest,self).__init__(name,threads,"sweep",grid_max=grid_max)
        
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
        output = super(SweepThreadingTest,self).__str__()
        output += "  mx         = %s\n" % self.amrdata.mx
        output += "  mxnest     = %s\n" % self.amrdata.mxnest
        output += "  grid_max   = %s\n" % self.grid_max
            
        return output
        
        
class SingleGridSweepThreadingTest(SweepThreadingTest):
        
    def __init__(self,name,threads,mx=40):
        
        # Set base values
        super(SingleGridSweepThreadingTest,self).__init__(name,threads,mx=mx,mxnest=1,grid_max=mx)

        # Setup non variable time stepping and output
        dt = 0.0001
        self.amrdata.dt_initial = dt
        self.amrdata.dt_variable = 0
        
        self.amrdata.outstyle = 3
        num_steps = int(TFINAL / dt) + 1
        self.amrdata.max_steps = num_steps + 1
        self.amrdata.iout = [num_steps,num_steps]

        # File log label
        self.file_label = "_%s_m%s" % (self.name,str(self.amrdata.mx).zfill(3))
        
        
class GridThreadingTest(BaseThreadingTest):
    
    def __init__(self,name,threads,mx=40,mxnest=3,grid_max=60):
        
        # Call super class to set name
        super(GridThreadingTest,self).__init__(name,threads,"grid",grid_max=grid_max)
        
        # Construct this test case's data
        self.amrdata.mxnest = mxnest
        self.amrdata.mx = mx
        self.amrdata.my = mx / 4
        
        # Variable time stepping and output
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
        output = super(GridThreadingTest,self).__str__()
        output += "  mx         = %s\n" % self.amrdata.mx
        output += "  mxnest     = %s\n" % self.amrdata.mxnest
        output += "  grid_max   = %s\n" % self.grid_max
            
        return output
        
class StaticGridThreadingTest(GridThreadingTest):
    
    def __init__(self,name,threads,grid_max=40):
        
        # Call super class to set base parameters and load data
        super(StaticGridThreadingTest,self).__init__(name,threads,grid_max=grid_max,mxnest=1)
        
        # File log label
        self.file_label = "_%s_g%s" % (self.name,str(self.grid_max).zfill(3))

        # Setup non variable time stepping and output
        dt = 0.0001
        self.amrdata.dt_initial = dt
        self.amrdata.dt_variable = 0
        
        self.amrdata.outstyle = 3
        num_steps = int(TFINAL / dt) + 1
        self.amrdata.max_steps = num_steps + 1
        self.amrdata.iout = [num_steps,num_steps]
        

    def __str__(self):
        output = "Test name: %s - %s\n" % (str(test.__class__).strip("<class '__main__.").strip("'>"),self.name)
        output += "  file_label = %s\n" % self.file_label
        output += "  threads    = %s\n" % self.threads
        output += "  grid_max   = %s\n" % self.grid_max
            
        return output
        
        
    def run_tests(self):
        self.open_log_files()

        # Build binary
        os.environ["THREADING_METHOD"] = self.thread_method
        os.environ["MAX1D"] = str(self.grid_max)
        self.log_file.write("Building...\n")
        self.log_file.write("  max1d = %s" % self.grid_max)
        self.flush_log_files()
        subprocess.Popen("make new -j %s" % self.threads,shell=True,
                                stdout=self.build_file,stderr=self.build_file,
                                env=os.environ).wait()
        self.log_file.write("Build completed.\n")
        
        # Run tests
        for num_threads in self.threads:
            # Determine mx
            self.amrdata.mx = self.grid_max * num_threads
            self.amrdata.my = self.amrdata.mx / 4
            
            # Write out data file
            self.amrdata.write()

            # Set environment variables
            os.environ["OMP_NUM_THREADS"] = str(num_threads)
            self.time_file.write("\n *** OMP_NUM_THREADS = %s\n\n" % int(num_threads))
            self.log_file.write("\n *** OMP_NUM_THREADS = %s\n\n" % int(num_threads))

            run_cmd = construct_run_cmd(self._time_file_path)
            self.log_file.write(run_cmd + "\n")
            self.flush_log_files()
            self.time_file.close()
            subprocess.Popen("make .data",shell=True,stdout=self.log_file,stderr=self.log_file,
                                env=os.environ).wait()
            subprocess.Popen(run_cmd,shell=True,stdout=self.log_file,stderr=self.log_file,
                                env=os.environ).wait()
            self.log_file.write("Simulation completed.\n")
            self.time_file = open(self._time_file_path,'aw')

        self.close_log_files()
    
        
# =============================================================================

# =============================================================================
# Create tests
# =============================================================================
tests = []
max_threads = int(os.environ['OMP_NUM_THREADS'])

# Test ranges
threads = [1,2,4,8,12,16]
for (i,count) in enumerate(threads):
    if count > max_threads:
        threads = threads[:i]
        break
        
single_grid_mx = [40,60,80,100,120,140,160,180]
grid_max_tests = [40,60,80,100,120,140,160,180]

# Single Grid Tests
# =================
#   Single grid tests of each approach
#   
#   Sweep Threading
#   ---------------
#     Vary (mx,my) and threads
for mx in single_grid_mx:
    tests.append(SingleGridSweepThreadingTest("single_sweep",threads,mx=mx))
    
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
    tests.append(StaticGridThreadingTest("static_grid",threads,grid_max=grid_max))

# Adaptive Grid Tests
# ===================
#   Tests for both sweep and grid threading and for all p
 for grid_max in grid_max_tests:
     tests.append(SweepThreadingTest("amr_sweep",threads,mx=40,mxnest=3,grid_max=grid_max))
     tests.append(GridThreadingTest("amr_grid",threads,mx=40,mxnest=3,grid_max=grid_max))

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
    
    # Execute tests
    for (i,test) in enumerate(tests_to_be_run):
        test.run_tests()
