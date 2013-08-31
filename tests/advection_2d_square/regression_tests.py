
from clawpack.visclaw import data
import os
from numpy import allclose
import subprocess

def test1():
    #os.system("make clean; make .output")
    job = subprocess.Popen(['make', 'clean'])
    return_code = job.wait()
    assert return_code == 0, "Problem with 'make clean'"
    job = subprocess.Popen(['make', '.output'])
    return_code = job.wait()
    assert return_code == 0, "Problem with 'make .output'"


def test2():
    plotdata = data.ClawPlotData()
    plotdata.outdir = open('.output').readline().strip()
    assert plotdata.outdir == '_output', "Unexpected contents of .output"

    g1 = plotdata.getgauge(1)
    tsum = g1.t.sum()
    qsum = g1.q[0,:].sum()
    tsum_expected = 15.844
    qsum_expected = 46.3254633
    tol = 1e-14
    assert allclose(tsum,tsum_expected,tol), \
        "gauge 1: tsum = %s, expected: %s"  % (tsum, tsum_expected)
    assert allclose(qsum,qsum_expected,tol), \
        "gauge 1: qsum = %s, expected: %s"  % (qsum, qsum_expected)
    print "Gauge 1 OK"
    

def test3():
    plotdata = data.ClawPlotData()
    plotdata.outdir = open('.output').readline().strip()
    assert plotdata.outdir == '_output', "Unexpected contents of .output"

    g2 = plotdata.getgauge(2)
    tsum = g2.t.sum()
    qsum = g2.q[0,:].sum()
    tsum_expected = 8.907
    qsum_expected = 4.9890462
    tol = 1e-14
    assert allclose(tsum,tsum_expected,tol), \
        "gauge 2: tsum = %s, expected: %s"  % (tsum, tsum_expected)
    assert allclose(qsum,qsum_expected,tol), \
        "gauge 2: qsum = %s, expected: %s"  % (qsum, qsum_expected)
    print "Gauge 2 OK"
    
if __name__=="__main__":
    test1()
    test2()
    test3()


