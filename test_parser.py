"""
Basic smoke test for the .all parser using a synthetic file.
Run with: python3 test_parser.py
"""

import sys
import os
import math

sys.path.insert(0, os.path.dirname(__file__))

from generate_test_all import generate
from all_parser import parse_all_file


def test_basic_parse():
    test_file = "/tmp/test_mbe.all"
    generate(output_path=test_file, n_pings=10, n_beams=64)

    result = parse_all_file(test_file)

    assert result.n_pings == 10, f"Expected 10 pings, got {result.n_pings}"
    assert len(result.soundings) > 0, "No soundings extracted"
    assert result.n_position_datagrams == 10, f"Expected 10 positions, got {result.n_position_datagrams}"

    # Check depths are sensible
    depths = [s.depth for s in result.soundings]
    assert all(50 < d < 150 for d in depths), f"Depths out of range: {min(depths):.1f} – {max(depths):.1f}"

    # Check lat/lon have been assigned
    lats = [s.latitude for s in result.soundings]
    lons = [s.longitude for s in result.soundings]
    assert all(56 < la < 58 for la in lats), f"Latitudes unexpected: {min(lats):.4f} – {max(lats):.4f}"

    print(f"PASS: {len(result.soundings)} soundings, "
          f"depth {min(depths):.1f}–{max(depths):.1f} m, "
          f"lat {min(lats):.4f}–{max(lats):.4f}")

    if result.errors:
        print(f"  Parser warnings ({len(result.errors)}):")
        for e in result.errors[:5]:
            print(f"    {e}")

    os.unlink(test_file)


if __name__ == "__main__":
    test_basic_parse()
    print("All tests passed.")
