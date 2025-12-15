#!/usr/bin/env python
"""
Test script to verify load_bin function works for different num_eqn values.

This tests the load_bin function in util.py which should work regardless
of spatial dimensions (1D, 2D, 3D) since it only depends on num_eqn.
"""

import sys
import os
import tempfile
import struct
import numpy as np

# Add path to import load_bin
# We're already in src/python/amrclaw, so we can import directly
from amrclaw.util import load_bin


def create_test_binary_file(path, num_eqn=1, has_header=True, num_points=10):
    """Create a test binary gauge file"""
    with open(path, 'wb') as f:
        if has_header:
            f.write(b"# Test header\n")
        
        # Write binary data: for each point, write (level, time, q[0], q[1], ...)
        for i in range(num_points):
            # level (int32)
            f.write(struct.pack('<i', 1))
            # time (float32)
            f.write(struct.pack('<f', float(i) * 0.1))
            # q values (float32 each)
            for j in range(num_eqn):
                f.write(struct.pack('<f', float(i * num_eqn + j) * 0.01))
    
    return path


def test_load_bin_1eqn():
    """Test load_bin with num_eqn=1 (1D case)"""
    print("Testing load_bin with num_eqn=1...")
    
    temp_dir = tempfile.mkdtemp()
    test_file = os.path.join(temp_dir, 'test_gauge.bin')
    
    try:
        create_test_binary_file(test_file, num_eqn=1, num_points=5)
        level, time, q = load_bin(test_file, num_eqn=1, has_header=True)
        
        assert len(level) == 5, f"Expected 5 points, got {len(level)}"
        assert len(time) == 5, f"Expected 5 time points, got {len(time)}"
        assert q.ndim == 1, f"Expected 1D array for num_eqn=1, got {q.ndim}D"
        assert len(q) == 5, f"Expected 5 q values, got {len(q)}"
        
        print(f"  ✓ Loaded {len(time)} points")
        print(f"  ✓ Level shape: {level.shape}, Time shape: {time.shape}, Q shape: {q.shape}")
        print("  ✓ num_eqn=1 test passed\n")
    finally:
        import shutil
        shutil.rmtree(temp_dir)


def test_load_bin_multieqn():
    """Test load_bin with num_eqn>1 (e.g., 2D or 3D systems)"""
    print("Testing load_bin with num_eqn=3...")
    
    temp_dir = tempfile.mkdtemp()
    test_file = os.path.join(temp_dir, 'test_gauge.bin')
    
    try:
        create_test_binary_file(test_file, num_eqn=3, num_points=5)
        level, time, q = load_bin(test_file, num_eqn=3, has_header=True)
        
        assert len(level) == 5, f"Expected 5 points, got {len(level)}"
        assert len(time) == 5, f"Expected 5 time points, got {len(time)}"
        assert q.ndim == 2, f"Expected 2D array for num_eqn>1, got {q.ndim}D"
        assert q.shape == (5, 3), f"Expected shape (5, 3), got {q.shape}"
        
        print(f"  ✓ Loaded {len(time)} points")
        print(f"  ✓ Level shape: {level.shape}, Time shape: {time.shape}, Q shape: {q.shape}")
        print("  ✓ num_eqn=3 test passed\n")
    finally:
        import shutil
        shutil.rmtree(temp_dir)


def test_load_bin_no_header():
    """Test load_bin without header"""
    print("Testing load_bin without header...")
    
    temp_dir = tempfile.mkdtemp()
    test_file = os.path.join(temp_dir, 'test_gauge.bin')
    
    try:
        create_test_binary_file(test_file, num_eqn=1, has_header=False, num_points=3)
        level, time, q = load_bin(test_file, num_eqn=1, has_header=False)
        
        assert len(level) == 3, f"Expected 3 points, got {len(level)}"
        print("  ✓ No header test passed\n")
    finally:
        import shutil
        shutil.rmtree(temp_dir)


if __name__ == '__main__':
    print("=" * 60)
    print("Testing load_bin function (dimension-independent)")
    print("=" * 60)
    print()
    
    results = []
    results.append(("load_bin with num_eqn=1", test_load_bin_1eqn()))
    results.append(("load_bin with num_eqn=3", test_load_bin_multieqn()))
    results.append(("load_bin without header", test_load_bin_no_header()))
    
    print("=" * 60)
    print("Summary:")
    print("=" * 60)
    all_passed = True
    for name, passed in results:
        status = "✓ PASSED" if passed else "✗ FAILED"
        print(f"  {name}: {status}")
        if not passed:
            all_passed = False
    
    print()
    if all_passed:
        print("All tests passed! ✓")
        print()
        print("The load_bin function works correctly and is dimension-independent.")
        print("It only depends on num_eqn, not on spatial dimensions (1D/2D/3D).")
    else:
        print("Some tests failed. ✗")
    
    exit(0 if all_passed else 1)

