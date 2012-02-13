#!/usr/bin/env python
# encoding: utf-8
r"""Run AMRClaw thread tests

p = num_threads

Single Grid Tests
=================
  Single grid tests of each approach

  Sweep Threading
  ---------------
    Vary (mx,my) and threads, mx = ]

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
import getopt
import subprocess
import time
import math
import glob

import setrun
    
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
        output = "Test name: %s - %s\n" % (str(self.__class__).strip("<class '__main__.").strip("'>"),self.name)
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
        
    def open_log_files(self,log_base):
        # Open log files and write headers
        self._build_file_path = os.path.join(log_base,"build%s.txt" % self.file_label)
        self._time_file_path = os.path.join(log_base,"time%s.txt" % self.file_label)
        self._log_file_path = os.path.join(log_base,"log%s.txt" % self.file_label)
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
    
    def write_data(self,num_threads):
        # Write out data file
        self.amrdata.write()
        self.probdata.write()
        
    def run_tests(self,log_path='./logs',compiler='gcc'):
        self.open_log_files(log_path)
        
        # Build binary
        self.log_file.write("Building...\n")
        os.environ["THREADING_METHOD"] = self.thread_method
        os.environ["MAX1D"] = str(self.grid_max)
        self.flush_log_files()
        subprocess.Popen("set_devel %s opt; make new -j %s" % (compiler,
                                self.threads),shell=True,
                                stdout=self.build_file,stderr=self.build_file,
                                env=os.environ).wait()
        self.log_file.write("Build completed.\n")
        
        # Run tests
        for num_threads in self.threads:
            os.environ["OMP_NUM_THREADS"] = str(num_threads)
            self.time_file.write("\n *** OMP_NUM_THREADS = %s\n\n" % int(num_threads))
            self.log_file.write("\n *** OMP_NUM_THREADS = %s\n\n" % int(num_threads))

            self.write_data(num_threads)

            run_cmd = construct_run_cmd(self._time_file_path)
            self.log_file.write(run_cmd + "\n")
            self.flush_log_files()
            self.time_file.close()
            subprocess.Popen(('set_devel %s opt' % compiler, run_cmd),
                                shell=True,stdout=self.log_file,
                                stderr=self.log_file,
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
        
    def write_data(self,num_threads):
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
tests = {}
max_threads = int(os.environ['OMP_NUM_THREADS'])

# Test ranges
host_name = os.uname()[1]
if host_name == 'irene':
    threads = [1,2,3,4]
    sqrt_threads = [1,4,9,16]
elif host_name == 'tiberius':
    threads = [1,2,3,4]
    sqrt_threads = [1,4,9,16]
elif host_name == 'juniper':
    threads = [1,2,4,8,12,16]
    sqrt_threads = [1,4,9,16]
else:
    threads = [1,2,4,8,12,16]
    sqrt_threads = [1,4,9,16,25]
for (i,count) in enumerate(threads):
    if count > max_threads:
        threads = threads[:i]
        break
for (i,count) in enumerate(sqrt_threads):
    if count > max_threads:
        sqrt_threads = sqrt_threads[:i]
        break
        
single_grid_mx = [N - 4 for N in [128,256,512,1024]]
grid_max_tests = [32,64,128,256,512,1024]

# Single Grid Tests
# =================
#   Single grid tests of each approach
#   
#   Sweep Threading
#   ---------------
#     Vary (mx,my) and threads
tests['single_sweep'] = []
for mx in single_grid_mx:
    tests['single_sweep'].append(SingleGridThreadingTest("single_sweep",threads,'sweep',mx))
# Weak scaling test
tests['weak_sweep'] = []
for mx in single_grid_mx:
    tests['weak_sweep'].append(WeakThreadingTest("weak_sweep",sqrt_threads,'sweep',mx,mx))
    
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
tests['weak_grid'] = []
for grid_max in grid_max_tests:
    tests['weak_grid'].append(WeakThreadingTest('weak_grid',sqrt_threads,'grid',64,grid_max))

# Adaptive Grid Tests
# ===================
#   Tests for both sweep and grid threading and for all p
tests['amr_sweep'] = []
for grid_max in grid_max_tests:
    tests['amr_sweep'].append(BaseThreadingTest('amr_sweep',threads,'sweep',64,3,grid_max))
tests['amr_grid'] = []
for grid_max in grid_max_tests:
    tests['amr_grid'].append(BaseThreadingTest('amr_grid',threads,'grid',64,3,grid_max))
    

# =============================================================================
#  Command line support
available_test_groups = '\t' + ', '.join(tests.keys())
help_message = r"""Run threading tests

Usage:  run_thread_tests.py [options] test_group [test_rank], ...

test_group = Test group to run of which multiple can be given.  To run all 
             available tests use 'all'.  Valid test groups include:
             
%s

test_rank = Tests to run from within the test groups specified

Command line options:
    -v, --verbose - Verbose output (default = False)
    -c, --compiler - Which compiler suite to use (either gcc, intel, or all),
                     (default = gcc)
    -l, --logs - Location to put logs (default = ./logs_[compiler_suite])
    -h, --help - Display this help message
""" % available_test_groups

class Usage(Exception):
    def __init__(self,msg):
        self.msg = msg   

if __name__ == "__main__":
    # Parse command line arguments
    try:
        try:
            opts,args = getopt.getopt(sys.argv[1:],
                            "hvc:l:",
                            ['help','verbose','compiler','logs'])
                                
        except getopt.error, msg:
            raise Usage(msg)
            
        # Default arguments
        verbose = False
        compiler = 'gcc'
        log_base = './logs'
        
        # Option parsing
        for option,value in opts:
            if option in ('-v','--verbose'):
                verbose = True
            if option in ('-c','--compiler'):
                if value in ('gcc','intel','all'):
                    compiler = value
                else:
                    raise Exception("Invalid compiler choice, must be 'intel' or 'gcc'.")
            if option in ('-l','--logs'):
                log_base = os.path.expandvars(os.path.expanduser(value))
            if option in ('-h','--help'):
                raise Usage(help_message)
    except Usage, err:
        print >> sys.stderr, sys.argv[0].split('/')[-1] + ": " + str(err.msg)
        sys.exit(1)
        
    # Construct test set
    if len(args) == 0:
        # Print out available test groups and ranks
        for (group,group_tests) in tests.iteritems():
            print "===== %s ========================" % group
            for (i,test) in enumerate(group_tests):
                print "Test (%s)" % i
                print test
        sys.exit(0)
    elif len(args) > 0:
        # Parse argument list
        tests_to_run = {}
        for arg in args:
            try:
                test_index = int(arg)
                for group in current_groups:
                    tests_to_run[group].append(test_index)
            except ValueError, err:
                current_groups = []
                if arg == 'all':
                    current_groups = tests.keys()
                    for group in current_groups:
                        tests_to_run[group] = []
                else:
                    for group_name in tests.keys():
                        if arg in group_name:
                            current_groups.append(group_name)
                            if not tests_to_run.has_key(group_name):
                                tests_to_run[group_name] = []
    
    # Print out request tests, also add all test in a group at this point if 
    # none were specifically requested
    test_count = 0
    for group in tests_to_run.keys():
        if len(tests_to_run[group]) == 0:
            tests_to_run[group] = range(len(tests[group]))
        print "===== %s ========================" % group
        for test_index in tests_to_run[group]:
            print "Test (%s)" % test_index
            print tests[group][test_index]
            test_count += 1
    print "Total tests = %s" % test_count
    
    # Run requested tests
    if compiler == 'all':
        compilers = ['gcc','intel']
    else:
        compilers = [compiler]
    for FC in compilers:
        # Remove old log files if they exist in this location
        log_path = '_'.join((log_base,FC))
        if os.path.exists(log_path):
            if os.path.isdir(log_path):
                for log_file in glob.glob(os.path.join(log_path,'*.txt')):
                    os.remove(log_file)
                os.rmdir(log_path)
        os.makedirs(log_path)
        
        # Finally, run the tests
        for group in tests_to_run.iterkeys():
            for test_index in tests_to_run[group]:
                tests[group][test_index].run_tests(log_path,FC)
