# run from the example dir
import numpy as np
from clawpack.amrclaw.util import load_bin

for g in ["_output/gauge00000.bin", "_output/gauge00001.bin"]:
    lvl, t, q = load_bin(g, num_eqn=1, has_header=True)
    print(g, len(t), "rows; first3:", list(zip(lvl[:3], t[:3], q[:3])))
