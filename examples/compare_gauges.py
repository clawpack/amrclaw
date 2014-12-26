
from clawpack.visclaw.data import ClawPlotData
from clawpack.visclaw import gaugetools
from pylab import *
    
setgauges = gaugetools.read_setgauges('.')
gaugenos = setgauges.gauge_numbers
if len(gaugenos)==0:
    print "No gauges found in gauges.data"
else:
    plotdata1 = ClawPlotData()
    plotdata1.outdir = '_output'
    plotdata2 = ClawPlotData()
    plotdata2.outdir = '_output_old_gauges'
    for gaugeno in gaugenos:
        g1 = plotdata1.getgauge(gaugeno)
        q1 = g1.q
        g2 = plotdata2.getgauge(gaugeno)
        q2 = g2.q
        figure(gaugeno)
        clf()
        meqn = q1.shape[0]
        for m in range(meqn):
            subplot(meqn,1,m+1)
            dq = abs(q1[m,:]-q2[m,:]).max()
            print  "Max difference in q[%s] at gauge %s is %g" % (m,gaugeno,dq)
            plot(g1.t,g1.q[m,:],'b',label='1')
            plot(g2.t,g2.q[m,:],'r',label='2')
            legend()
        subplot(1,1,1)
        title('Comparison of gauge number %s' % gaugeno)
        



