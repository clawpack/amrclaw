#!/usr/bin/env python

"""Base AMRClaw data class for writing out data parameter files."""

import os

import clawpack.clawutil.data

class AmrclawInputData(clawpack.clawutil.data.ClawData):
    r"""
    Data object containing AMRClaw input data.

    Extends ClawInputData adding necessary data for AMR.
    """
    def __init__(self, clawdata):

        # Some defaults are inherited from ClawInputData:
        super(AmrclawInputData,self).__init__()

        # Need to have a pointer to this so we can get num_dim and num_aux
        self._clawdata = clawdata
        
        # Refinement control
        self.add_attribute('amr_levels_max',1)
        self.add_attribute('refinement_ratios_x',[1])
        self.add_attribute('refinement_ratios_y',[1])
        if self._clawdata.num_dim == 3:
            self.add_attribute('refinement_ratios_z',[1])
        if self._clawdata.num_dim == 1:
            raise NotImplementedError("1d AMR not yet supported")
        self.add_attribute('variable_dt_refinement_ratios',False)

        self.add_attribute('refinement_ratios_t',[1])
        self.add_attribute('aux_type',[])

        # Flagging control
        self.add_attribute('flag_richardson',False)
        self.add_attribute('flag_richardson_tol',-1.0)
        self.add_attribute('flag2refine',True)
        self.add_attribute('flag2refine_tol',0.05)
        self.add_attribute('regrid_interval',2)
        self.add_attribute('regrid_buffer_width',3)
        self.add_attribute('clustering_cutoff',0.7)
        self.add_attribute('verbosity_regrid',3)

        
        # debugging flags:
        self.add_attribute('dprint',False)
        self.add_attribute('eprint',False)
        self.add_attribute('edebug',False)
        self.add_attribute('gprint',False)
        self.add_attribute('nprint',False)
        self.add_attribute('pprint',False)
        self.add_attribute('rprint',False)
        self.add_attribute('sprint',False)
        self.add_attribute('tprint',False)
        self.add_attribute('uprint',False)


    def write(self, out_file='amr.data', data_source='setrun.py'):

        self.open_data_file(out_file, data_source)
    
        self.data_write('amr_levels_max')

        num_ratios = max(abs(self.amr_levels_max)-1, 1)
        if len(self.refinement_ratios_x) < num_ratios:
            raise ValueError("*** Error in data parameter: " + \
                  "require len(refinement_ratios_x) >= %s " % num_ratios)
        if len(self.refinement_ratios_y) < num_ratios:
            raise ValueError("*** Error in data parameter: " + \
                  "require len(refinement_ratios_y) >= %s " % num_ratios)
        self.data_write('refinement_ratios_x')
        self.data_write('refinement_ratios_y')
        if self._clawdata.num_dim == 3:
            if len(self.refinement_ratios_z) < num_ratios:
                    raise ValueError("*** Error in data parameter: " + \
                      "require len(refinement_ratios_z) >= %s " % num_ratios)
            self.data_write('refinement_ratios_z')
        if len(self.refinement_ratios_t) < num_ratios:
            raise ValueError("*** Error in data parameter: " + \
                  "require len(refinement_ratios_t) >= %s " % num_ratios)
        self.data_write('refinement_ratios_t')
        self.data_write()  # writes blank line

        if self._clawdata.num_aux > 0:
            self.data_write(file, self.aux_type, 'aux_type')
        self.data_write()

        self.data_write('flag_richardson')
        if self.flag_richardson == False:
            # Still need to add flag_richardson to fortran, for now this works:
            self.flag_richardson_tol = -1.   
        self.data_write('flag_richardson_tol')
        self.data_write('flag2refine')
        self.data_write('flag2refine_tol')
        self.data_write('regrid_interval')
        self.data_write('regrid_buffer_width')
        self.data_write('clustering_cutoff')
        self.data_write('verbosity_regrid')
        self.data_write()

        self.data_write()

        self.data_write('dprint')
        self.data_write('eprint')
        self.data_write('edebug')
        self.data_write('gprint')
        self.data_write('nprint')
        self.data_write('pprint')
        self.data_write('rprint')
        self.data_write('sprint')
        self.data_write('tprint')
        self.data_write('uprint')
        self.data_write()

        self.close_data_file()


# ==============================================================================
#  Region data object
class RegionData(clawpack.clawutil.data.ClawData):
    r""""""

    def __init__(self,regions=None,num_dim=2):

        super(RegionData,self).__init__()

        if regions is None or not isinstance(regions,list):
            self.add_attribute('regions',[])
        else:
            self.add_attribute('regions',regions)
        self.add_attribute('num_dim',num_dim)


    def write(self,out_file='regions.data',data_source='setrun.py'):


        self.open_data_file(out_file,data_source)

        self.data_write(value=len(self.regions),alt_name='num_regions')
        if (self.num_dim == 3) and (len(self.regions) > 0):
            raise NotImplementedError("*** Regions not yet implemented in 3d")
        for regions in self.regions:
            self._out_file.write(8*"%g  " % tuple(regions) +"\n")
        self.close_data_file()


    def read(self, path, force=False):
        r"""Read in region data from file at path
        """

        data_file = open(os.path.abspath(path),'r')

        # Read past comments and blank lines
        header_lines = 0
        ignore_lines = True
        while ignore_lines:
            line = data_file.readline()
            if line[0] == "#" or len(line.strip()) == 0:
                header_lines += 1
            else:
                break

        # Read in number of regions        
        num_regions, tail = line.split("=:")
        num_regions = int(num_regions)
        varname = tail.split()[0]
        if varname != "num_regions":
            raise IOError("It does not appear that this file contains region data.")

        # Read in each region
        self.regions = []
        for n in xrange(num_regions):
            line = data_file.readline().split()
            self.regions.append([int(line[0]), int(line[1]),
                                 float(line[2]), float(line[3]), 
                                 float(line[4]), float(line[5]),
                                 float(line[6]), float(line[7])])

        data_file.close()


# ==============================================================================
        
# ==============================================================================
#  Gauge data object
class GaugeData(clawpack.clawutil.data.ClawData):
    r""""""

    @property
    def gauge_numbers(self):
        if len(self.gauges) == 1:
            return [self.gauges[0][0]]
        else:
            return [gauge[0] for gauge in self.gauges]

    def __init__(self, num_dim=2):
        super(GaugeData,self).__init__()

        self.add_attribute('num_dim',num_dim)
        self.add_attribute('gauges',[])

    def __str__(self):
        output = "Gauges: %s\n" % len(self.gauges)
        for gauge in self.gauges:
            output = "\t".join((output,"%4i:" % gauge[0]))
            output = " ".join((output,"%19.10e" % gauge[1]))
            output = " ".join((output,"%17.10e" % gauge[2]))
            output = " ".join((output,"%13.6e" % gauge[3]))
            output = " ".join((output,"%13.6e\n" % gauge[4]))
        return output

    def write(self,out_file='gauges.data',data_source='setrun.py'):
        r"""Write out gague information data file."""

        if (self.num_dim == 3) and (len(self.gauges) > 0):
            raise NotImplementedError("*** Gauges not yet implemented in 3d")

        # Check to make sure we have only unique gauge numebrs
        if len(self.gauges) > 0:
            if len(self.gauge_numbers) != len(set(self.gauge_numbers)):
                raise Exception("Non unique gauge numbers specified.")

        # Write out gauge data file
        self.open_data_file(out_file,data_source)
        self.data_write(name='ngauges',value=len(self.gauges))
        for gauge in self.gauges:
            self._out_file.write("%4i %19.10e  %17.10e  %13.6e  %13.6e\n" % tuple(gauge))
        self.close_data_file()

    def read(self,data_path="./",file_name='gauges.data'):
        r"""Read gauge data file"""
        path = os.path.join(data_path, file_name)
        gauge_file = open(path,'r')

        # Read past comments and blank lines
        header_lines = 0
        ignore_lines = True
        while ignore_lines:
            line = gauge_file.readline()
            if line[0] == "#" or len(line.strip()) == 0:
                header_lines += 1
            else:
                break

        # Read number of gauges, should be line that was last read in
        num_gauges = int(line.split()[0])

        # Read in each gauge line
        for n in xrange(num_gauges):
            line = gauge_file.readline().split()
            self.gauges.append([int(line[0]),float(line[1]),float(line[2]),
                                             float(line[3]),float(line[4])])

        gauge_file.close()

#  Gauge data objects
# ==============================================================================


if __name__ == '__main__':
    raise Exception("Not unit tests have been defined for this module.")

