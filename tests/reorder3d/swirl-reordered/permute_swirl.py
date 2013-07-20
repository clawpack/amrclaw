from permute import *
import sys

if __name__ == '__main__':
    fname = sys.argv[1]
    if (fname == 'driver.f'):
        plist = [('q', [4,1,2,3]),
                 ('aux', [4,1,2,3])]
    elif (fname == 'qinit.f'):
        plist = [('q', [4,1,2,3]),
                 ('aux', [4,1,2,3])]
    elif (fname == 'setaux.f'):
        plist = [('aux', [4,1,2,3])]

    permute_file(sys.argv[1], sys.argv[2], plist, False)
