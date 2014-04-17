"""
Run all examples and then compare to gallery versions, if available.
"""

from clawpack.clawutil import regression_tests, make_all
import os

my_env = os.environ
my_env['GIT_STATUS'] = 'True'
my_env['FFLAGS'] = '-O2 -fopenmp'
my_env['OMP_NUM_THREADS'] = '3'

make_all.make_all(make_clean_first=True, env=my_env)

print "\n-----------------------------------------------------------\n"

all_ok = regression_tests.test_subdirs()
if all_ok:
    print "===> All tests pass"
else:
    print "===> Some test(s) failed"
