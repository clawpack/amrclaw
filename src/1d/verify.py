# run from the example dir
import numpy as np

def load_bin(path, num_eqn=1, has_header=True, endian="<"):
    with open(path, "rb") as f:
        if has_header: f.readline()
        raw = f.read()
    dt = np.dtype([("level",  endian+"i4"),
                   ("time",   endian+"f4"),
                   ("q",      (endian+"f4", (num_eqn,)))])
    arr = np.frombuffer(raw, dtype=dt)
    return arr["level"], arr["time"], arr["q"][:,0] if num_eqn==1 else arr["q"]

for g in ["_output/gauge00000.bin", "_output/gauge00001.bin"]:
    lvl,t,q = load_bin(g, num_eqn=1, has_header=True)
    print(g, len(t), "rows; first3:", list(zip(lvl[:3], t[:3], q[:3])))
