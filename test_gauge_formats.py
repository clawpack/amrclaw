#!/usr/bin/env python
"""
Test script to verify gauge format writing works for 1D, 2D, and 3D.

This tests the format string calculation in data.py line 478.
"""

import sys
import os
import tempfile
import shutil

# Add the source to the path
amrclaw_dir = os.path.dirname(__file__)
clawpack_dir = os.path.dirname(os.path.dirname(amrclaw_dir))
sys.path.insert(0, os.path.join(clawpack_dir, 'clawutil', 'src', 'python'))
sys.path.insert(0, os.path.join(amrclaw_dir, 'src', 'python'))

from clawpack.clawutil import data


def test_gauge_format_1d():
    """Test gauge format writing for 1D"""
    print("Testing 1D gauge format...")
    
    rundata = data.ClawRunData('amrclaw', 1)
    rundata.gaugedata.gauges = []
    # 1D format: [gaugeno, x, t1, t2]
    rundata.gaugedata.gauges.append([0, 0.2, 0.0, 1e9])
    rundata.gaugedata.gauges.append([1, 0.5, 0.0, 1e9])
    
    # Test different formats
    rundata.gaugedata.file_format = ['ascii', 'binary32']
    
    temp_dir = tempfile.mkdtemp()
    old_cwd = os.getcwd()
    try:
        os.chdir(temp_dir)
        rundata.write()
        
        # Check that gauges.data was created
        gauges_file = os.path.join(temp_dir, 'gauges.data')
        assert os.path.exists(gauges_file), "gauges.data file not created"
        
        # Read and verify format
        with open(gauges_file, 'r') as f:
            lines = f.readlines()
            # Skip header lines (those with '=:' are data_write format), find gauge data
            gauge_lines = []
            for line in lines:
                if line.strip() and not line.startswith('#') and '=:' not in line:
                    parts = line.split()
                    # Gauge data lines should have the format: gaugeno x t1 t2 (4 fields for 1D)
                    if len(parts) >= 4:  # At least 4 fields for gauge data
                        gauge_lines.append(parts)
            
            assert len(gauge_lines) == 2, f"Expected 2 gauge lines, got {len(gauge_lines)}"
            for parts in gauge_lines:
                assert len(parts) == 4, f"1D gauge should have 4 fields, got {len(parts)}: {parts}"
                print(f"  ✓ 1D gauge line: {' '.join(parts)}")
        
        print("  ✓ 1D test passed\n")
    finally:
        os.chdir(old_cwd)
        shutil.rmtree(temp_dir)


def test_gauge_format_2d():
    """Test gauge format writing for 2D"""
    print("Testing 2D gauge format...")
    
    rundata = data.ClawRunData('amrclaw', 2)
    rundata.gaugedata.gauges = []
    # 2D format: [gaugeno, x, y, t1, t2]
    rundata.gaugedata.gauges.append([1, 0.6, 0.4, 0.0, 10.0])
    rundata.gaugedata.gauges.append([2, 0.3, 0.7, 0.0, 10.0])
    
    rundata.gaugedata.file_format = ['binary32', 'binary64']
    
    temp_dir = tempfile.mkdtemp()
    old_cwd = os.getcwd()
    try:
        os.chdir(temp_dir)
        rundata.write()
        
        gauges_file = os.path.join(temp_dir, 'gauges.data')
        assert os.path.exists(gauges_file), "gauges.data file not created"
        
        with open(gauges_file, 'r') as f:
            lines = f.readlines()
            # Skip header lines (those with '=:' are data_write format), find gauge data
            gauge_lines = []
            for line in lines:
                if line.strip() and not line.startswith('#') and '=:' not in line:
                    parts = line.split()
                    # Gauge data lines should have the format: gaugeno x y t1 t2 (5 fields for 2D)
                    if len(parts) >= 5:  # At least 5 fields for gauge data
                        gauge_lines.append(parts)
            
            assert len(gauge_lines) == 2, f"Expected 2 gauge lines, got {len(gauge_lines)}"
            for parts in gauge_lines:
                assert len(parts) == 5, f"2D gauge should have 5 fields, got {len(parts)}: {parts}"
                print(f"  ✓ 2D gauge line: {' '.join(parts)}")
        
        print("  ✓ 2D test passed\n")
    finally:
        os.chdir(old_cwd)
        shutil.rmtree(temp_dir)


def test_gauge_format_3d():
    """Test gauge format writing for 3D"""
    print("Testing 3D gauge format...")
    
    rundata = data.ClawRunData('amrclaw', 3)
    rundata.gaugedata.gauges = []
    # 3D format: [gaugeno, x, y, z, t1, t2]
    rundata.gaugedata.gauges.append([1, 1.0, 0.1, 0.1, 0.0, 1e9])
    rundata.gaugedata.gauges.append([2, 0.1, 1.0, 0.1, 0.0, 1e9])
    
    rundata.gaugedata.file_format = ['ascii', 'binary32']
    
    temp_dir = tempfile.mkdtemp()
    old_cwd = os.getcwd()
    try:
        os.chdir(temp_dir)
        rundata.write()
        
        gauges_file = os.path.join(temp_dir, 'gauges.data')
        assert os.path.exists(gauges_file), "gauges.data file not created"
        
        with open(gauges_file, 'r') as f:
            lines = f.readlines()
            # Skip header lines (those with '=:' are data_write format), find gauge data
            gauge_lines = []
            for line in lines:
                if line.strip() and not line.startswith('#') and '=:' not in line:
                    parts = line.split()
                    # Gauge data lines should have the format: gaugeno x y z t1 t2 (6 fields for 3D)
                    if len(parts) >= 6:  # At least 6 fields for gauge data
                        gauge_lines.append(parts)
            
            assert len(gauge_lines) == 2, f"Expected 2 gauge lines, got {len(gauge_lines)}"
            for parts in gauge_lines:
                assert len(parts) == 6, f"3D gauge should have 6 fields, got {len(parts)}: {parts}"
                print(f"  ✓ 3D gauge line: {' '.join(parts)}")
        
        print("  ✓ 3D test passed\n")
    finally:
        os.chdir(old_cwd)
        shutil.rmtree(temp_dir)


def test_format_string_calculation():
    """Test the format string calculation directly"""
    print("Testing format string calculation...")
    
    # Test cases
    gauge_1d = [0, 0.2, 0, 1e9]
    gauge_2d = [1, 0.6, 0.4, 0.0, 10.0]
    gauge_3d = [1, 1.0, 0.1, 0.1, 0.0, 1e9]
    
    for name, gauge in [('1D', gauge_1d), ('2D', gauge_2d), ('3D', gauge_3d)]:
        format_str = "%4i" + (len(gauge)-3) * "  %17.10e" + 2 * "  %13.6e" + "\n"
        result = format_str % tuple(gauge)
        expected_fields = len(gauge)
        actual_fields = len(result.split())
        assert actual_fields == expected_fields, \
            f"{name}: Expected {expected_fields} fields, got {actual_fields}"
        print(f"  ✓ {name}: {len(gauge)} fields formatted correctly")
    
    print("  ✓ Format string calculation test passed\n")


if __name__ == '__main__':
    print("=" * 60)
    print("Testing gauge format writing for 1D, 2D, and 3D")
    print("=" * 60)
    print()
    
    results = []
    results.append(("Format string calculation", test_format_string_calculation()))
    results.append(("1D gauge format", test_gauge_format_1d()))
    results.append(("2D gauge format", test_gauge_format_2d()))
    results.append(("3D gauge format", test_gauge_format_3d()))
    
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
        sys.exit(0)
    else:
        print("Some tests failed. ✗")
        sys.exit(1)

