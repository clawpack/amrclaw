"""
Runs "make .plots" and then compares the .png files in _plots with the
archived results in the Clawpack gallery. 

Function test returns True if all .png files agree in the two directories.

Image differences can be viewed by opening a browser to 
    _image_diff/_ImageDiffIndex.html

For this to work you need to first clone and/or fetch the latest gallery
results from
     git://github.com/clawpack/clawpack.github.com

"""

import os, subprocess
from clawpack.clawutil import imagediff

verbose = True
relocatable = False  # True ==> copies all images to subdirs of _image_diff
                     # so this directory can be moved, e.g. posted on web

def test(stdout=None,stderr=None):
    """
    Redirect output and errors if stdout or stderr passed in.
    """
    
    import sys
    if stdout is not None:
        sys.stdout = stdout
    if stderr is not None:
        sys.stderr = stderr
        
    try:
        CLAW = os.environ['CLAW']
    except:
        raise Exception("Environment variable CLAW not set")

    # Compare plots to what is in the Clawpack gallery:
    gallery_dir = CLAW + "/clawpack.github.com/doc/_static/"

    # For testing with gallery built locally, instead use:
    # gallery_dir = CLAW + "/doc/doc/gallery/"  

    this_example = os.path.split(os.getcwd())[-1]
    
    gallery_plots = gallery_dir + "amrclaw/examples/" + this_example + "/_plots"
    if not os.path.isdir(gallery_plots):
        error_msg = "Missing directory %s\n Need to clone clawpack.github.com\n"\
                     % gallery_plots
        #raise Exception(error_msg)
        sys.stderr.write(error_msg)
        gallery_found = False
    else:
        gallery_found = True
        

    # Run the code and create _plots directory:
    cmd = "make clean; make .output"
    status=subprocess.Popen(cmd,shell=True).wait()
    if status != 0:
        #raise Exception("Problem running the code, status = %s" % status)
        error_msg = "Problem running the code, status = %s" % status
        sys.stderr.write(error_msg)
    else:
        cmd = "make .plots"
        status=subprocess.Popen(cmd,shell=True).wait()
        if status != 0:
            error_msg = "Problem running the code, status = %s" % status
            sys.stderr.write(error_msg)

    # Compare the resulting plots to the gallery version:
    if (status==0) and gallery_found:
        try:
            regression_ok = imagediff.imagediff_dir('_plots',gallery_plots, \
                                    relocatable=relocatable,overwrite=True, \
                                    verbose=verbose)
            if not regression_ok:
                sys.stderr.write("***Regression files are not all identical")
        except:
            error_msg = "Error running imagediff with directories \n  %s\n  %s\n" \
                        % ('_plots',gallery_plots)
            #raise Exception(error_msg)
            sys.stderr.write(error_msg)
            regression_ok = False
    else:
        regression_ok = False
            
    return regression_ok
    
if __name__=="__main__":
    regression_ok = test()
    print "regression_ok = ",regression_ok


