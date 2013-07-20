from permute import *
import sys

if __name__ == '__main__':
    fname = sys.argv[1]
    if (fname == 'rpn3adv.f'):
        plist = [('wave', [2, 3, 1]),
                 ('ql', [2, 1]),
                 ('qr', [2, 1]),
                 ('amdq', [2, 1]),
                 ('apdq', [2, 1]),
                 ('auxl', [2, 1]),
                 ('auxr', [2, 1])]
    elif (fname == 'rpt3adv.f'):
        plist = [('ql', [2, 1]),
                 ('qr', [2, 1]),
                 ('asdq', [2, 1]),
                 ('bmasdq', [2, 1]),
                 ('bpasdq', [2, 1]),
                 ('aux1', [2, 1, 3]),
                 ('aux2', [2, 1, 3]),
                 ('aux3', [2, 1, 3])]
    elif (fname == 'rptt3adv.f'):
        plist = [('ql', [2, 1]),
                 ('qr', [2, 1]),
                 ('bsasdq', [2, 1]),
                 ('cmbsasdq', [2, 1]),
                 ('cpbsasdq', [2, 1]),
                 ('aux1', [2, 1, 3]),
                 ('aux2', [2, 1, 3]),
                 ('aux3', [2, 1, 3])]

    permute_file(sys.argv[1], sys.argv[2], plist, False)
