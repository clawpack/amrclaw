
from clawpack.visclaw import data
import os
from numpy import allclose

def test1():
    os.system("make clean; make .output")
    plotdata = data.ClawPlotData()
    plotdata.outdir = '_output'

    g0 = plotdata.getgauge(0)
    tsum = g0.t.sum()
    qsum = g0.q[0,:].sum()
    tsum_expected = 4.652
    qsum_expected = 20.9349558403
    tol = 1e-14
    assert allclose(tsum,tsum_expected,tol), \
        "gauge 0: tsum = %s, expected: %s"  % (tsum, tsum_expected)
    assert allclose(qsum,qsum_expected,tol), \
        "gauge 0: qsum = %s, expected: %s"  % (qsum, qsum_expected)
    print "Gauge 0 OK"

    g1 = plotdata.getgauge(1)
    tsum = g1.t.sum()
    qsum = g1.q[0,:].sum()
    tsum_expected = 4.7
    qsum_expected = 19.178640866
    tol = 1e-14
    assert allclose(tsum,tsum_expected,tol), \
        "gauge 1: tsum = %s, expected: %s"  % (tsum, tsum_expected)
    assert allclose(qsum,qsum_expected,tol), \
        "gauge 1: qsum = %s, expected: %s"  % (qsum, qsum_expected)
    print "Gauge 1 OK"
    
    g2 = plotdata.getgauge(2)
    tsum = g2.t.sum()
    qsum = g2.q[0,:].sum()
    tsum_expected = 1.1
    qsum_expected = 4.82442404
    tol = 1e-14
    assert allclose(tsum,tsum_expected,tol), \
        "gauge 2: tsum = %s, expected: %s"  % (tsum, tsum_expected)
    assert allclose(qsum,qsum_expected,tol), \
        "gauge 2: qsum = %s, expected: %s"  % (qsum, qsum_expected)
    print "Gauge 2 OK"
    os.system("rm -rf _output")
    
if __name__=="__main__":
    test1()


