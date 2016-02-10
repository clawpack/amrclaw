#!/usr/bin/env python

"""Module contains class definitions related to dealing with gauge data"""

import os

import numpy

class GaugeSolution(object):
    r"""Object representing a gauge station

    Gauges are measurement points in a simulation that can be set arbitrarily 
    relative to the computational grid.  They record all time step updates for 
    the finest grid they are located in so can contain a much more detailed 
    view on a particular point.  Each gauge file contains an identifying 
    ID number placed in the file name and by default recorded by the object.
    Additionally the time """


    def __init__(self, gauge_id=None, path=None):
        r"""GaugeSolution Initiatlization Routine
        
        See :class:`GaugeSolution` for more info.
        """

        # Gauge descriptors
        self.id = None
        self.location = None

        # Gauge data
        self.level = None
        self.t = None
        self.q = None

        # Read in gauge data from file
        if gauge_id is not None:
            if path is None:
                path = os.getcwd()
            self.read(gauge_id, path)


    def read(self, gauge_id, path=None):
        r"""Read the gauge file at path into this object

        """

        # Construct path to gauge
        if path is None:
            path = os.getcwd()
        gauge_file_name = "gauge%s.txt" % str(gauge_id).zfill(5)
        gauge_path = os.path.join(path, gauge_file_name)

        # Read header info
        with open(gauge_path, 'r') as gauge_file:
            # First line
            data = gauge_file.readline().split()
            self.id = int(data[2])
            self.location = (float(data[4]), float(data[5]))
            num_eqn = int(data[8])

        # Check to see if the gauge file name ID and that inside of the gauge
        # file are the same
        gauge_file_name = os.path.basename(path)
        if self.id != gauge_id:
            raise ValueError("Gauge ID requested does not match ID inside ",
                             "file!")

        # Read gauge data
        data = numpy.loadtxt(gauge_path, comments="#")
        self.level = data[:, 0].view('int64')
        self.t = data[:, 1]
        self.q = data[:, 2:].transpose()

        if num_eqn != self.q.shape[0]:
            raise ValueError("Number of fields in gauge file does not match",
                             "recorded number in header.")


    def write(self, path):
        r"""Write the data from this gauge to a file in `path`

        """

        if not self.is_valid():
            raise ValueError("Gauge is not initialized properly.")

        gauge_file_name = "gauge%s.txt" % str(self.id).zfill(5)
        with open(os.path.join(path, gauge_file_name), "w") as gauge_file:

            gauge_file.write("# gauge_id= %s location=( %s %s ) num_eqn= %s\n" %
                 (self.id, self.location[0], self.location[1], self.q.shape[0]))
            gauge_file.write("# Columns: level time q(1 ... num_eqn)\n")

            for i in xrange(self.q.shape[1]):
                gauge_file.write("%s %s %s\n" % (self.level[i], self.t[i], 
                              " ".join([str(value) for value in self.q[:, i]])))


    def is_valid(self):
        r"""Check to see if the gauge data has all been set."""

        if ((self.id is not None) and (self.location is not None) and 
            (self.level is not None) and (self.t is not None) and
            (self.q is not None)):

            return True

        return False


    def __repr__(self):

        if self.is_valid():
            output = "%4i" % self.id
            for j in range(len(self.location)):
                output = " ".join((output,"%17.10e" % self.location[j]))
            output = " ".join((output,"%13.6e" % self.t[0]))
            output = " ".join((output,"%13.6e\n" % self.t[-1]))
        else:
            output = None

        return output


    def __str__(self):
        return ("Gauge %s: location = %s, t = [%s, %s]" % 
                                    (self.id, self.location,
                                     self.t[0], self.t[-1]))


if __name__ == "__main__":

    pass