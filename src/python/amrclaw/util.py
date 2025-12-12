"""
Utility functions for AMRClaw.

This module contains utility functions for reading and processing AMRClaw output.
"""

import numpy as np


def load_bin(path, num_eqn=1, has_header=True, endian="<"):
    """Load binary gauge file data.
    
    Reads a binary gauge file and returns level, time, and q arrays.
    
    :Input:
     - *path* (str) - Path to the binary gauge file (.bin file)
     - *num_eqn* (int) - Number of equations/components in q. Default is 1.
     - *has_header* (bool) - Whether the file has a header line to skip. 
       Default is True.
     - *endian* (str) - Byte order, '<' for little-endian, '>' for big-endian.
       Default is '<'.
    
    :Output:
     - *level* (ndarray) - Array of refinement levels
     - *time* (ndarray) - Array of time values
     - *q* (ndarray) - Array of q values. If num_eqn==1, returns 1D array,
       otherwise returns 2D array with shape (npoints, num_eqn)
    
    :Example:
    
        >>> level, time, q = load_bin("_output/gauge00000.bin", num_eqn=1, has_header=True)
        >>> print(f"Loaded {len(time)} time points")
    """
    with open(path, "rb") as f:
        if has_header:
            f.readline()
        raw = f.read()
    
    dt = np.dtype([("level", endian + "i4"),
                   ("time", endian + "f4"),
                   ("q", (endian + "f4", (num_eqn,)))])
    arr = np.frombuffer(raw, dtype=dt)
    
    if num_eqn == 1:
        return arr["level"], arr["time"], arr["q"][:, 0]
    else:
        return arr["level"], arr["time"], arr["q"]

