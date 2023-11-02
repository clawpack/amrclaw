#!/usr/bin/env python

import sys

import numpy
import matplotlib.pyplot as plt

import clawpack.pyclaw.gauges as gauges

# Load all gauges

def check_old_gauge_data(path, gauge_id):

    # Load old gauge data
    data = numpy.loadtxt(path)
    old_ids = numpy.asarray(data[:, 0], dtype=int)
    gauge_indices = numpy.nonzero(old_ids == gauge_id)[0]
    q = data[gauge_indices, 3:]

    # Load new data
    gauge = gauges.GaugeSolution(gauge_id, "./regression_data/")

    print(numpy.linalg.norm(q - gauge.q.transpose(), ord=2))
    print(numpy.argmax(q - gauge.q.transpose()))

    fig = plt.figure()
    for i in range(gauge.q.shape[0]):
        axes = fig.add_subplot(1, gauge.q.shape[0], i + 1)
        axes.plot(q[:, i] - gauge.q[i, :])
        axes.set_title("q[%s, :] comparison" % i)

    return fig

if __name__ == "__main__":

    old_files = []
    gauge_ids = []

    if len(sys.argv) > 1:
        # Assume arguments are pairs of file locations and ids
        for i in range(1, len(sys.argv), 2):
            old_files.append(sys.argv[i])
            gauge_ids.append(int(sys.argv[i + 1]))
    else:
        raise ValueError("Need at least one pair to compare")

    print(old_files)
    print(gauge_ids)

    figures = []
    for i in range(len(old_files)):
        figures.append(check_old_gauge_data(old_files[i], gauge_ids[i]))

    plt.show()