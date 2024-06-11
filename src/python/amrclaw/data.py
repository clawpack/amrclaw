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
        
        # set max1d, maximum size of patch in each dimension
        # and memsize, initial length of alloc array 
        #    use default values so that repeated doubling gets as close as
        #    possible to 2**30-1 = 2147483647, beyond which integer(kind=4)
        #    indexes overflow and it would be necessary to use kind=8.
        if self._clawdata.num_dim == 3:
            self.add_attribute('memsize',8388607)
            self.add_attribute('max1d',32)
        elif self._clawdata.num_dim == 2:
            self.add_attribute('memsize',4194303)
            self.add_attribute('max1d',60)
        else:
            self.add_attribute('memsize',1048575)
            self.add_attribute('max1d',500)

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
    
        self.data_write('memsize')
        self.data_write('max1d')
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
#  FlagRegion object (new in v5.7.0)

class FlagRegion(clawpack.clawutil.data.ClawData):
    
    def __init__(self, num_dim, region=None):
        r"""
        Note: for backward compatibility, can pass in a old-style region 
        in the form (in 2d):
            [minlevel,maxlevel,t1,t2,x1,x2,y1,y2]
        and it will be converted to the new form.
        
        """
        
        super(FlagRegion,self).__init__()
        self.add_attribute('name','')
        self.add_attribute('num_dim',num_dim)
        self.add_attribute('minlevel',None)
        self.add_attribute('maxlevel',None)
        self.add_attribute('t1',None)
        self.add_attribute('t2',None)
        self.add_attribute('spatial_region_type',None)
        self.add_attribute('spatial_region',None)
        self.add_attribute('spatial_region_file',None)
        
        if region is not None:
            self.convert_old_region(region)

        
    def convert_old_region(self, region):
        """
        Take a list region = [minlevel, maxlevel, t1, t2, x1, x2, y1, y2]
        in the old style and convert to a new flagregion.
        """
        self.minlevel = region[0]
        self.maxlevel = region[1]
        self.t1 = region[2]
        self.t2 = region[3]
        self.spatial_region_type = 1  # rectangle
        self.spatial_region = region[4:]
        
    def read_spatial_region(self, fname=None):
        """
        Reads a ruled rectangle and assigns to self.spatial_region
        But this isn't really supported yet, since cannot write to 
        flagregions.data.
        """
        from clawpack.amrclaw import region_tools
        
        if fname is None:
            fname = self.spatial_region_file
        if not os.path.isfile(fname):
            print('*** Could not find file ',fname)
        else:
            self.spatial_region = region_tools.RuledRectangle(fname)
            
        
# ==============================================================================
#  FlagRegionData object  (new in v5.7.0)

class FlagRegionData(clawpack.clawutil.data.ClawData):
    r"""To replace RegionData, allowing more flexibility."""

    def __init__(self,flagregions=None,num_dim=2):

        super(FlagRegionData,self).__init__()

        if flagregions is None or not isinstance(flagregions,list):
            self.add_attribute('flagregions',[])
        else:
            self.add_attribute('flagregions',flagregions)
        self.add_attribute('num_dim',num_dim)


    def write(self,out_file='flagregions.data',data_source='setrun.py'):
        """
        Write out the flagregion.data file that will be read by the Fortran
        code.  
        """

        self.open_data_file(out_file,data_source)

        self.data_write(value=len(self.flagregions),alt_name='num_flagregions')
        for flagregion in self.flagregions:
            flagregion._out_file = self._out_file
            # write minlevel,maxlevel as integers, t1,t2 as floats:
            flagregion.data_write()
            flagregion.data_write('name')
            flagregion.data_write('minlevel')
            flagregion.data_write('maxlevel')
            flagregion.data_write('t1')
            flagregion.data_write('t2')
            flagregion.data_write('spatial_region_type')
            
            if flagregion.spatial_region_type == 1:
                # spatial_region is an extent [x1,x2,y1,y2]
                assert len(flagregion.spatial_region) == 4, \
                    "*** FlagRegionData.write requires len(spatial_region) = 4"
                flagregion.data_write('spatial_region')
                
            elif flagregion.spatial_region_type == 2:
                if type(flagregion.spatial_region_file) is str:
                    # assumed to be path to RuledRectangle file
                    flagregion.data_write('spatial_region_file')
                else:
                    raise ValueError('*** FlagRegionData.write requires path to RuledRectangle')
                    
        self.close_data_file()


    def read(self, path='flagregions.data', force=False):
        r"""Read in flagregion data from file at path
        """

        with open(os.path.abspath(path),'r') as data_file:
            # Read past comments and blank lines
            header_lines = 0
            ignore_lines = True
            while ignore_lines:
                line = data_file.readline()
                if line[0] == "#" or len(line.strip()) == 0:
                    header_lines += 1
                else:
                    break
    
            # Read in number of flagregions        
            num_flagregions, tail = line.split("=:")
            num_flagregions = int(num_flagregions)
            varname = tail.split()[0]
            if varname != "num_flagregions":
                raise IOError("It does not appear that this file contains flagregion data.")
    
            print('Will read %i flagregions' % num_flagregions)
            
            # Read in each flagregion
            self.flagregions = []
            for n in range(num_flagregions):
                flagregion = FlagRegion(num_dim=2)
                line = data_file.readline() # blank
                line = data_file.readline().split()
                flagregion.name = str(line[0])[1:-1] # strip extra quotes
                line = data_file.readline().split()
                flagregion.minlevel = int(line[0])
                line = data_file.readline().split()
                flagregion.maxlevel = int(line[0])
                line = data_file.readline().split()
                flagregion.t1 = float(line[0])
                line = data_file.readline().split()
                flagregion.t2 = float(line[0])
                line = data_file.readline().split()
                flagregion.spatial_region_type = int(line[0])
                
                line = data_file.readline().split()
                if flagregion.spatial_region_type == 1:
                    flagregion.spatial_region = [float(val) for val in line[:4]]
                elif flagregion.spatial_region_type == 2:
                    # the next line is assumed to be path to a
                    # RuledRectangle file, strip quotes
                    flagregion.spatial_region_file = line[0][1:-1] 
                    # read the file: 
                    flagregion.read_spatial_region()
                else:
                    raise ValueError('*** Unrecognized spatial_region_type')
                
                self.flagregions.append(flagregion)
    
# ==============================================================================
#  Gauge data object
class GaugeData(clawpack.clawutil.data.ClawData):
    r""""""

    defaults = {"file_format":"ascii", "display_format":"e15.7",
                "q_out_fields":"all", "aux_out_fields":"none",
                "min_time_increment":0.0, "gtype":"stationary"}

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

        for (value, default) in self.defaults.items():
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
        output = "\n\n".join((output, "Gauge type: %s\n" % self.gtype))
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

        for gaugeno in self.gauge_numbers:
            if abs(gaugeno) >= 2**31:
                raise Exception("Gauge number %i is too large, must be < 2**31"\
                    % gaugeno)


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
        self._out_file.write("# File format (1=ascii, 2=binary64, 3=binary32)\n")
        format_map = {'ascii':1, 'binary':2, 'binary64':2, 'binary32':3}
        for gauge_num in self.gauge_numbers:
            try:
                file_format = format_map[self.file_format[gauge_num].lower()]
            except KeyError:
                raise ValueError("Invalid file format %s requested." % self.file_format[gauge_num])
            self._out_file.write("%s " % file_format)
        self._out_file.write("\n")
        self.data_write()

        # Display format for each gauge
        self._out_file.write("# Display format (for ascii files)\n")
        for gauge_num in self.gauge_numbers:
            self._out_file.write("%s " % self.display_format[gauge_num])
        self._out_file.write("\n\n")

        # Minimum time increment output
        self._out_file.write("# Minimum output increment\n")
        for gauge_num in self.gauge_numbers:
            self._out_file.write("%s " % self.min_time_increment[gauge_num])
        self._out_file.write("\n\n")

        # Gauge type
        self._out_file.write("# Gauge type (1=stationary, 2=Lagrangian)\n")
        gtype_map = {'stationary':1, 'lagrangian': 2}
        for gauge_num in self.gauge_numbers:
            try:
                gtype = gtype_map[self.gtype[gauge_num].lower()]
            except KeyError:
                raise ValueError("Invalid gauge type %s requested." \
                        % self.gtype[gauge_num])
            self._out_file.write("%s " % gtype)
        self._out_file.write("\n")
        self.data_write()

        # Which q fields to output
        self._out_file.write("# q fields\n")
        for gauge_num in self.gauge_numbers:
            # Handle special values of "all" and "none"
            if isinstance(self.q_out_fields[gauge_num], str):
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
                if isinstance(self.aux_out_fields[gauge_num], str):
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

        this would set all gauges to have the 'ascii' format.

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


# ==============================================================================
#  Adjoint data object
class AdjointData(clawpack.clawutil.data.ClawData):
    r""""""

    def __init__(self, use_adjoint=False, num_dim=None):

        super(AdjointData,self).__init__()

        self.add_attribute('use_adjoint',use_adjoint)
        self.add_attribute('adjoint_outdir','')
        self.add_attribute('adjoint_files',[])
        self.add_attribute('numadjoints',0)
        self.add_attribute('t1',None)
        self.add_attribute('t2',None)
        self.add_attribute('innerprod_index', None)


    def write(self,out_file='adjoint.data',data_source='setrun.py'):

        self.open_data_file(out_file,data_source)

        self.data_write('use_adjoint')
        self.data_write('adjoint_outdir')
        self.data_write('t1')
        self.data_write('t2')
        self.data_write('innerprod_index')

        self.set_adjoint_files()
        self.data_write('numadjoints')
            
        for file in self.adjoint_files:
            # if path is relative in setrun, assume it's relative to the
            # same directory that out_file comes from
            fname = os.path.abspath(os.path.join(os.path.dirname(out_file),file))
            self._out_file.write("'%s'\n" % fname)
        self.close_data_file()

    def set_adjoint_files(self):
        import glob
        from clawpack.pyclaw.fileio.ascii import read_t

        self.adjoint_files = []
        if self.use_adjoint:
            files = glob.glob(os.path.join(self.adjoint_outdir,"fort.b*"))
            files.sort()
            for file in files:
                frameno = int(file[-4:])
                #print('+++ file, frameno: ',file, frameno)
                [t,num_eqn,nstates,num_aux,num_dim,num_ghost,file_format] \
                    = read_t(frameno, self.adjoint_outdir)
                if t <= self.t2:
                    #print('+++   using this file')
                    self.adjoint_files.append(file)
            self.numadjoints = len(self.adjoint_files)
            if (len(self.adjoint_files) == 0):
                print("*** WARNING: No binary files found for adjoint output!")


    def read_adjoint_files(self):
        self.set_adjoint_files()
        # read them in for Python processing...  add if we need it



if __name__ == '__main__':
    raise Exception("Not unit tests have been defined for this module.")

