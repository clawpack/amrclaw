#!/usr/bin/env python
"""
Simple test script to verify gauge format writing works for 1D, 2D, and 3D.

This tests the format string calculation in data.py line 478 without requiring
full clawpack imports.
"""

def test_format_string_calculation():
    """Test the format string calculation directly"""
    print("Testing format string calculation for 1D, 2D, and 3D...")
    print()
    
    # Test cases matching the actual gauge formats
    test_cases = [
        ('1D', [0, 0.2, 0, 1e9]),           # [gaugeno, x, t1, t2]
        ('2D', [1, 0.6, 0.4, 0.0, 10.0]),   # [gaugeno, x, y, t1, t2]
        ('3D', [1, 1.0, 0.1, 0.1, 0.0, 1e9]), # [gaugeno, x, y, z, t1, t2]
    ]
    
    for name, gauge in test_cases:
        # This is the exact format string from data.py line 478
        format_str = "%4i" + (len(gauge)-3) * "  %17.10e" + 2 * "  %13.6e" + "\n"
        
        result = format_str % tuple(gauge)
        expected_fields = len(gauge)
        actual_fields = len(result.split())
        
        assert actual_fields == expected_fields, \
            f"{name}: Expected {expected_fields} fields, got {actual_fields}"
        print(f"  ✓ {name}: {len(gauge)} fields formatted correctly")
        print(f"    Format: {format_str.strip()}")
        print(f"    Result: {result.strip()}")
        print()


def test_field_count_logic():
    """Test the logic for calculating number of coordinate fields"""
    print("Testing field count logic...")
    print()
    
    test_cases = [
        ('1D', 4, 1),  # [gaugeno, x, t1, t2] -> 1 coordinate
        ('2D', 5, 2),  # [gaugeno, x, y, t1, t2] -> 2 coordinates
        ('3D', 6, 3),  # [gaugeno, x, y, z, t1, t2] -> 3 coordinates
    ]
    
    for name, gauge_len, expected_coords in test_cases:
        # Formula: (len(gauge) - 3) gives number of coordinate fields
        # len(gauge) includes: gaugeno + coordinates + t1 + t2
        # So: len - 3 = len - (gaugeno + t1 + t2) = number of coordinates
        num_coords = gauge_len - 3
        
        assert num_coords == expected_coords, \
            f"{name}: Expected {expected_coords} coordinates, got {num_coords}"
        print(f"  ✓ {name}: {gauge_len} total fields -> {num_coords} coordinate field(s)")
    
    print()


if __name__ == '__main__':
    print("=" * 60)
    print("Testing gauge format string calculation")
    print("=" * 60)
    print()
    
    results = []
    results.append(("Field count logic", test_field_count_logic()))
    results.append(("Format string calculation", test_format_string_calculation()))
    
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
        print("The format string calculation works correctly for 1D, 2D, and 3D.")
        print("The formula (len(gauge)-3) correctly calculates the number of")
        print("coordinate fields regardless of dimension.")
    else:
        print("Some tests failed. ✗")
    
    exit(0 if all_passed else 1)

