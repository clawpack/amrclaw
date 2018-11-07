#!/usr/bin/env python

"""Base AMRClaw data class for writing out data parameter files."""

from __future__ import absolute_import
import os

import clawpack.clawutil.data
import six
from six.moves import range

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
        if self._clawdata.num_dim >= 2:
            self.add_attribute('refinement_ratios_y',[1])
        if self._clawdata.num_dim == 3:
            self.add_attribute('refinement_ratios_z',[1])
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
        self.data_write('refinement_ratios_x')
        if self._clawdata.num_dim >= 2:
            if len(self.refinement_ratios_y) < num_ratios:
                raise ValueError("*** Error in data parameter: " + \
                      "require len(refinement_ratios_y) >= %s " % num_ratios)
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

        if len(self.aux_type) < self._clawdata.num_aux:
            raise ValueError("*** length of aux_type array should be " +\
                  "same as num_aux = %i" % self._clawdata.num_aux)
        if self._clawdata.num_aux > 0:
            self.data_write('aux_type')
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

    def __init__(self,regions=None,num_dim=None):

        super(RegionData,self).__init__()

        if regions is None or not isinstance(regions,list):
            self.add_attribute('regions',[])
        else:
            self.add_attribute('regions',regions)
        self.add_attribute('num_dim',num_dim)


    def write(self,out_file='regions.data',data_source='setrun.py'):


        self.open_data_file(out_file,data_source)

        self.data_write(value=len(self.regions),alt_name='num_regions')
        for regions in self.regions:
            # write minlevel,maxlevel as integers, remaining floats in e format:
            self._out_file.write("%i  %i " % (regions[0], regions[1]))
            self._out_file.write((len(regions)-2)*"%20.14e  " % tuple(regions[2:]) +"\n")
            
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
        for n in range(num_regions):
            line = data_file.readline().split()
            self.regions.append([int(line[0]), int(line[1])] + [float(a) for a in line[2:]])

        data_file.close()


# ==============================================================================
        
# ==============================================================================
#  Gauge data object
class GaugeData(clawpack.clawutil.data.ClawData):
    r""""""

    defaults = {"file_format":"ascii", "display_format":"e15.7",
                "q_out_fields":"all", "aux_out_fields":"none",
                "min_time_increment":0.0}

    @property
    def gauge_numbers(self):
        if len(self.gauges) == 1:
            return [self.gauges[0][0]]
        else:
            return [gauge[0] for gauge in self.gauges]


    def __init__(self, num_dim=None):
        # num_dim is an argument for backward compatibility, but no longer used
        super(GaugeData,self).__init__()
        self.add_attribute('gauges',[])

        for (value, default) in six.iteritems(self.defaults):
            self.add_attribute(value, default)


    def __str__(self):
        output = "Gauges: %s\n" % len(self.gauges)
        for gauge in self.gauges:
            output = "\t".join((output,"%4i:" % gauge[0]))
            output = " ".join((output,"%19.10e" % gauge[1]))
            output = " ".join((output,"%17.10e" % gauge[2]))
            output = " ".join((output,"%13.6e" % gauge[3]))
            output = " ".join((output,"%13.6e\n" % gauge[4]))
        output = "\n\n".join((output, "Output Format: %s\n" % self.file_format))
        output = "\t".join((output, "display: %s" % self.display_format))
        output = "\n\t".join((output, "q fields: %s" % self.q_out_fields))
        output = "\n\t".join((output, "aux fields: %s" % self.aux_out_fields))
        output = "\n\t".join((output, "min. time increment: %s" 
                                                     % self.min_time_increment))
        return output


    def write(self, num_eqn, num_aux, out_file='gauges.data', 
                                      data_source='setrun.py'):
        r"""Write out gauge information data file."""

        # Check to make sure we have only unique gauge numbers
        if len(self.gauges) > 0:
            if len(self.gauge_numbers) != len(set(self.gauge_numbers)):
                raise Exception("Non unique gauge numbers specified.")

        # Write out gauge data file
        self.open_data_file(out_file,data_source)
        self.data_write(name='ngauges',value=len(self.gauges))
        for gauge in self.gauges:
            format = "%4i" + (len(gauge)-3) * "  %17.10e" + 2 * "  %13.6e" + "\n"
            self._out_file.write(format % tuple(gauge))
        self.data_write()
        
        # Expand all gauge format option dictionaries
        for key in self.defaults.keys():
            self.expand_gauge_format_option(key)

        # File format
        self._out_file.write("# File format\n")
        format_map = {'ascii':1, 'binary': 2}
        for gauge_num in self.gauge_numbers:
            try:
                file_format = format_map[self.file_format[gauge_num].lower()]
            except KeyError:
                raise ValueError("Invalid file format %s requested." % self.file_format[gauge_num])
            self._out_file.write("%s " % file_format)
        self._out_file.write("\n")
        self.data_write()

        # Display format for each gauge
        self._out_file.write("# Display format\n")
        for gauge_num in self.gauge_numbers:
            self._out_file.write("%s " % self.display_format[gauge_num])
        self._out_file.write("\n\n")

        # Minimum time increment output
        self._out_file.write("# Minimum output increment\n")
        for gauge_num in self.gauge_numbers:
            self._out_file.write("%s " % self.min_time_increment[gauge_num])
        self._out_file.write("\n\n")

        # Which q fields to output
        self._out_file.write("# q fields\n")
        for gauge_num in self.gauge_numbers:
            # Handle special values of "all" and "none"
            if isinstance(self.q_out_fields[gauge_num], six.string_types):
                if self.q_out_fields[gauge_num].lower() == 'all':
                    self._out_file.write("%s\n" % " ".join(['True'] * num_eqn))
                elif self.q_out_fields[gauge_num].lower() == 'none': 
                    self._out_file.write("%s\n" % " ".join(['False'] * num_eqn))
                else:
                    raise ValueError("Unknown q field string specified, '%s'" 
                                                 % self.q_out_fields[gauge_num])
            else:
                # Specify by field number
                if not isinstance(self.q_out_fields[gauge_num], list):
                    self.q_out_fields[gauge_num] = [self.q_out_fields[gauge_num]]
                bool_list = [n in self.q_out_fields[gauge_num] 
                                                       for n in range(num_eqn)]
                bool_list = [str(value) for value in bool_list]
                self._out_file.write("%s\n" % (" ".join(bool_list)))
        self.data_write()

        # Which aux fields to output
        if num_aux > 0:
            self._out_file.write("# aux fields\n")
            for gauge_num in self.gauge_numbers:
                # Handle special values of "all" and "none"
                if isinstance(self.aux_out_fields[gauge_num], six.string_types):
                    if self.aux_out_fields[gauge_num].lower() == 'all':
                        self._out_file.write("%s\n" % " ".join(['True'] * num_aux))
                    elif self.aux_out_fields[gauge_num].lower() == 'none': 
                        self._out_file.write("%s\n" % " ".join(['False'] * num_aux))
                    else:
                        raise ValueError("Unknown q field string specified, '%s'" 
                                                     % self.aux_out_fields[gauge_num])
                else:
                    # Specify by field number
                    if not isinstance(self.aux_out_fields[gauge_num], list):
                        self.aux_out_fields[gauge_num] = [self.aux_out_fields[gauge_num]]
                    bool_list = [n in self.aux_out_fields[gauge_num] 
                                                           for n in range(num_aux)]
                    bool_list = [str(value) for value in bool_list]
                    self._out_file.write("%s\n" % (" ".join(bool_list)))

        self.close_data_file()


    def expand_gauge_format_option(self, param_name):
        r"""Construct the full gauge output specification for *param_name*

        Also handles if each *param_name* is set to a single value (not a dict)
        and then assumes that all gauges should have this parameter.  For
        example if a user set

        > file_format = 'ascii'

        this would set all gauges ot have the 'ascii' format.

        """

        default_param = self.defaults[param_name]
        if not isinstance(getattr(self, param_name), dict):
            # Convert into dict for file
            default_param = getattr(self, param_name)
            setattr(self, param_name, {})

        # Check to make sure that gauges listed are actually gauges
        for gauge_num in getattr(self, param_name).keys():
            if gauge_num not in self.gauge_numbers:
                raise ValueError("Gauge number listed in format option not a ",
                                 "gauge.  gauge_id = %s" % gauge_num)

        for gauge_num in self.gauge_numbers:
            getattr(self, param_name).setdefault(gauge_num, default_param)


    def read(self, data_path="./", file_name='gauges.data'):
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
        for n in range(num_gauges):
            line = gauge_file.readline().split()
            self.gauges.append([int(line[0])] + [float(a) for a in line[1:]])

        # TODO:  Read in format data

        gauge_file.close()


#  Gauge data objects
# ==============================================================================


if __name__ == '__main__':
    raise Exception("Not unit tests have been defined for this module.")

