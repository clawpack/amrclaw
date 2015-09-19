#!/usr/bin/env python
# encoding: utf-8

"""
Script to run a code or restart it from a checkpoint file automatically.

Needs to be cleaned up and made more general.

Assumes code is compiled already.

    python run_with_restart.py

Will do the following:

 - Check whether _output exists with data from a run that appeared to
   complete (by checking last line of fort.amr -- more robust way?).

 - If so, it quits,

 - If _output exists with at least one checkpoint file, it will use the 
   more recent one to restart the code.

 - If _output does not exist, it runs the code from scratch.

So you should be able to type the command above, hit Ctrl-C, and repeat this
process an arbitrary number of times and eventually the full output will be
generated.  

Notes:

 - This requires replacing `check.f` by `check_twofiles.f` in the Makefile
   to use the version that alternates writing two sets of checkpoint files
   rather than accumulating files.

 - The setrun.py file is used also for the restart.  The clawdata.restart
   value is set to True and written to claw.data explicitly from this script.

 - check_twofiles.f also flushes buffers for output to fort.gauge,
   fort.amr, fort.debug so output is not lost.

 - Before restarting, it moves any fort.gauge file to fort.gauge_<datetime>
   so that eventually these can be catenated together (after editing out
   repeated entries).

"""

import subprocess
import os, sys
from setrun import setrun

outdir = '_output'

# set any desired environment flags:

env = os.environ
#env['FFLAGS'] = '-O2 -fopenmp'  # currently assume code is already compiled.

# runtime environment variables:
env['OMP_NUM_THREADS'] = '3'

# The next line insures that stdout is not buffered so if the code dies
# the output sent to run_output.txt so the error message is visible:
env['GFORTRAN_UNBUFFERED_PRECONNECTED'] = 'y'  

def examine_outdir(outdir='_output'):
    """
    Check the outdir to see if the code has already run to completion
    (in which case nothing is done) or needs to be restarted.
    If outdir does not exist, run from scratch.
    """

    from numpy import Inf
    fortamr = os.path.join(outdir,'fort.amr')
    try:
        f = open(fortamr).readlines()
        finished = ('end of' in f[-1])  # examine last line for ending message
    except:
        finished = False

    try:
        cfile = os.path.join(outdir,'fort.tckaaaaa')
        f = open(cfile).readlines()
        ta = float(f[0][29:])
    except:
        ta = -Inf

    try:
        cfile = os.path.join(outdir,'fort.tckbbbbb')
        f = open(cfile).readlines()
        tb = float(f[0][29:])
    except:
        tb = -Inf

    if (ta == -Inf) and (tb == -Inf):
        print "Could not read fort.tckaaaaa or fort.tckbbbbb in outdir %s" \
            % outdir
        latest = None
        t_latest = None
    elif ta > tb:
        latest = 'aaaaa'
        t_latest = ta
    else:
        latest = 'bbbbb'
        t_latest = tb

    return finished, latest, t_latest


def run_code_or_restart():

    import time
    tm = time.localtime()
    year = str(tm[0]).zfill(4)
    month = str(tm[1]).zfill(2)
    day = str(tm[2]).zfill(2)
    hour = str(tm[3]).zfill(2)
    minute = str(tm[4]).zfill(2)
    second = str(tm[5]).zfill(2)
    timestamp = '%s-%s-%s-%s%s%s'  % (year,month,day,hour,minute,second)

    finished, latest, t_latest = examine_outdir(outdir)

    if finished:
        print "Code has finished running, remove %s to run again" % outdir
        return

    restart = (latest is not None)

    fname_output = 'run_output.txt'
    fname_errors = 'run_errors.txt'

    if restart:
        print "Will attempt to restart using checkpoint file %s at t = %s" \
                % (latest, t_latest)
        print "Appending output stream to %s" % fname_output
        access = 'a'
    else:
        print "Will run code -- no restart"
        print "Writing output stream to %s" % fname_output
        access = 'w'

    fout = open(fname_output, access)
    ferr = open(fname_errors, access)

    if restart:
        fout.flush()
        fout.write("\n=========== RESTART =============\n" + \
                "Local time: %s\n" % timestamp + \
                "Will attempt to restart using checkpoint file %s at t = %s\n" \
                % (latest, t_latest))
        fout.flush()
        make_args = ['make','output','RESTART=True']
    else:
        make_args = ['make','output']

    if restart:
        fortgauge = os.path.join(outdir,'fort.gauge')
        fortgauge2 = os.path.join(outdir,'fort.gauge_%s' % timestamp)
        os.system("mv %s %s" % (fortgauge,fortgauge2))
        fout.write("Moving %s to %s \n" % (fortgauge,fortgauge2))
        fout.flush()

    rundata = setrun('amrclaw')
    rundata.clawdata.restart = restart
    rundata.clawdata.restart_file = 'fort.chk' + str(latest)
    if restart:
        rundata.clawdata.output_t0 = False  # to avoid plotting at restart times
    rundata.write()
    
    job = subprocess.Popen(make_args, stdout=fout,stderr=ferr,env=env)
    return_code = job.wait()
    
    if return_code == 0:
        print "Successful run\n"
    else:
        print "Problem running code\n"
    print "See %s and %s" % (fname_output,fname_errors)
    
    fout.close()
    ferr.close()


if __name__ == "__main__":

    run_code_or_restart()
            
