
# Script used to get rid of maxmx and maxmy dependencies when converting
# from 4.x to 5.0 form.  Executed in library and application directories.

#
# Fix a set of target files in directory tree rootdir by replacing
# oldpat with newpat. 
#
# Now supports wildcards in list of targetfiles.
#

from __future__ import absolute_import
from __future__ import print_function
import os,sys,glob
from six.moves import zip

rootdir = '.'
targetfiles = ['*.f*']

oldpat_list = ["1-mbc:maxmx", "1-mbc:maxmy", "maxmx,maxmy,"]
newpat_list = ["1-mbc:mx",    "1-mbc:my",    ""]

for oldpat,newpat in zip(oldpat_list, newpat_list):
    print("============================================")
    print('Replacing "%s" with "%s"' % (oldpat,newpat))
    print("============================================")
    for (dirpath, subdirs, files) in os.walk(rootdir):
        currentdir = os.path.abspath(os.getcwd())
        os.chdir(os.path.abspath(dirpath))
        tfiles = []
        for fpat in targetfiles:
            for f in glob.glob(fpat):
                tfiles.append(f)
        for file in tfiles:

            infile = open(file,'r')
            lines = infile.read()
            infile.close()

            if lines.find(oldpat) > -1:
                lines = lines.replace(oldpat, newpat)
                print("Fixed file   ",dirpath + '/' + file)
            else:
                print("No change to ",dirpath + '/' + file)

            outfile = open(file,'w')
            outfile.write(lines)
            outfile.close()

        os.chdir(currentdir)

