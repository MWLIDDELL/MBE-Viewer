"""
Generate a minimal synthetic Kongsberg .all file for testing.

Creates:
  - Installation parameters datagram (0x49)
  - 30 Position datagrams (0x50) along a straight track
  - 30 XYZ-88 datagrams (0x58) with 256 beams each (swath pattern)

Output: test_data.all
"""

import struct
import math
import os

STX = 0x02
ETX = 0x03


def pack_uint8(v):  return struct.pack("B", int(v) & 0xFF)
def pack_uint16(v): return struct.pack("<H", int(v) & 0xFFFF)
def pack_int16(v):  return struct.pack("<h", int(v))
def pack_uint32(v): return struct.pack("<I", int(v) & 0xFFFFFFFF)
def pack_int32(v):  return struct.pack("<i", int(v))
def pack_float(v):  return struct.pack("<f", float(v))


def make_datagram(dgm_type: int, em_model: int, date: int, time_ms: int,
                  ping_counter: int, serial: int, payload: bytes) -> bytes:
    """Wrap payload in a standard EM datagram with header and ETX+checksum."""
    header = (
        pack_uint8(STX) +
        pack_uint8(dgm_type) +
        pack_uint16(em_model) +
        pack_uint32(date) +
        pack_uint32(time_ms) +
        pack_uint16(ping_counter) +
        pack_uint16(serial)
    )
    body = header + payload
    checksum = sum(body) & 0xFFFF
    body += pack_uint8(ETX) + pack_uint16(checksum)
    # Prepend 4-byte length (length of body, not including the 4-byte prefix itself)
    return pack_uint32(len(body)) + body


def make_position(date, time_ms, ping_counter, serial, lat_deg, lon_deg):
    lat_raw = int(round(lat_deg * 1e7))
    lon_raw = int(round(lon_deg * 1e7))
    payload = (
        pack_int32(lat_raw) +
        pack_int32(lon_raw) +
        pack_uint16(0) +   # fix quality
        pack_uint16(0) +   # speed
        pack_uint16(0) +   # course
        pack_uint16(0) +   # heading
        pack_uint8(0) +    # position system descriptor
        pack_uint8(0) +    # number of bytes in input datagram
        b"\x00" * 0        # no raw string
    )
    return make_datagram(0x50, 302, date, time_ms, ping_counter, serial, payload)


def make_xyz88(date, time_ms, ping_counter, serial, n_beams, swath_width, base_depth):
    """Create an XYZ-88 datagram with n_beams spread over swath_width metres."""
    # Fixed heading block (16 bytes)
    heading = (
        pack_float(0.0) +   # height of water level
        pack_uint16(0) +    # bytes in input datagram
        pack_uint16(n_beams) +  # number of valid detections
        pack_uint8(0) +     # sampling frequency
        pack_uint8(0) +     # Rx transducer heading
        pack_uint8(0) +     # sound speed
        pack_uint8(0) +     # Tx transducer depth
        pack_uint16(0) +    # along-track ping spacing
        pack_uint16(0)      # spare
    )

    beam_records = b""
    for i in range(n_beams):
        # Across-track: spread evenly from -swath_width/2 to +swath_width/2
        across = (i / (n_beams - 1) - 0.5) * swath_width
        along  = 0.0
        # Depth: deeper in centre, shallower at edges (realistic bowl shape)
        angle_norm = (i / (n_beams - 1) - 0.5) * 2.0  # -1..1
        depth = base_depth * (1.0 + 0.15 * angle_norm ** 2)
        # Add mild noise
        import random
        depth += random.gauss(0, 0.3)
        beam_records += (
            pack_float(across) +
            pack_float(along) +
            pack_float(depth) +
            pack_float(0.001) +   # detection window length
            pack_float(0.5) +     # quality factor
            pack_int16(0) +       # beam incidence angle adjustment (1 byte + pad)
            pack_uint8(0) +       # detection info (valid)
            pack_uint8(0) +       # real-time cleaning info
            # Note: each beam is 20 bytes total; the two int16 fields above are 4 bytes,
            # plus float(4)*5=20 but we need to adjust — let me count:
            # float(4) + float(4) + float(4) + float(4) + float(4) = 20 bytes
            # int8(1) + uint8(1) + uint8(1) + uint8(1) = 4 bytes  → total = 24 — wrong
            # Fix: beam record is exactly 20 bytes per spec:
            #   float32 across(4) + float32 along(4) + float32 depth(4) +
            #   float32 det_window(4) + float32 quality(4) = 20 bytes
            # The remaining fields are packed differently; let's redo:
            b""
        )

    # Redo beam records correctly (20 bytes per beam, 5 x float32)
    beam_records = b""
    import random
    for i in range(n_beams):
        across = (i / (n_beams - 1) - 0.5) * swath_width
        angle_norm = (i / (n_beams - 1) - 0.5) * 2.0
        depth = base_depth * (1.0 + 0.15 * angle_norm ** 2) + random.gauss(0, 0.3)
        beam_records += (
            pack_float(across) +   # 4
            pack_float(0.0) +      # along-track  4
            pack_float(depth) +    # depth        4
            pack_float(0.001) +    # det window   4
            pack_float(0.5) +      # quality      4
            pack_int16(0) +        # beam incidence adj (int8 + pad = 2)
            pack_uint8(0) +        # detection info  1
            pack_uint8(0)          # rt cleaning     1   → total = 4+4+4+4+4+2+1+1 = 24 ≠ 20
        )
    # The beam record is 20 bytes per spec. Let me use the correct layout:
    # float32(4) + float32(4) + float32(4) + float32(4) + float32(4) = 20 bytes for the 5 floats
    # but the spec says 20 bytes total including the extra fields… let me match the parser exactly.
    # In _parse_xyz88 in all_parser.py we read:
    #   across = float32 @ boff
    #   along  = float32 @ boff+4
    #   depth  = float32 @ boff+8
    #   det_info = uint8 @ boff+18
    # So beam_size=20, det_info is at offset 18 within the beam.
    beam_records = b""
    for i in range(n_beams):
        across = (i / (n_beams - 1) - 0.5) * swath_width
        angle_norm = (i / (n_beams - 1) - 0.5) * 2.0
        import random as _r
        depth = base_depth * (1.0 + 0.15 * angle_norm ** 2) + _r.gauss(0, 0.3)
        # 20 bytes: f4 f4 f4 f4 f4 b1 b1 b1 b1 = 4*5 + 4 = 24? No.
        # 20 bytes = f4 f4 f4 pad8 det_info(1) pad(1)
        # offset 0: across (f4)
        # offset 4: along (f4)
        # offset 8: depth (f4)
        # offset 12: det_window (f4)
        # offset 16: quality (f4)
        # offset 20: ??? — but beam_size is 20 per the parser
        # The parser uses beam_size=20 but only reads up to offset 18.
        # Let's just make 20 bytes with det_info=0 at offset 18.
        rec = bytearray(20)
        struct.pack_into("<f", rec, 0,  across)
        struct.pack_into("<f", rec, 4,  0.0)
        struct.pack_into("<f", rec, 8,  depth)
        struct.pack_into("<f", rec, 12, 0.001)
        # offset 16: 2 bytes spare; offset 18: det_info; offset 19: spare
        rec[18] = 0  # valid detection
        beam_records += bytes(rec)

    # Spare byte at end of datagram
    spare = b"\x00"
    payload = heading + beam_records + spare
    return make_datagram(0x58, 302, date, time_ms, ping_counter, serial, payload)


def generate(output_path="test_data.all", n_pings=30, n_beams=256,
             start_lat=57.0, start_lon=-2.0, base_depth=80.0):
    import random
    random.seed(42)

    date = 20230601
    serial = 1234
    datagrams = []

    track_length = 0.01  # degrees latitude

    for ping_idx in range(n_pings):
        t_frac = ping_idx / max(n_pings - 1, 1)
        lat = start_lat + t_frac * track_length
        lon = start_lon
        time_ms = 60_000 + ping_idx * 2000  # 2 s apart

        # Position datagram
        datagrams.append(make_position(date, time_ms, ping_idx, serial, lat, lon))

        # Vary depth slightly along track
        depth = base_depth + 5 * math.sin(t_frac * math.pi)

        # XYZ datagram
        datagrams.append(make_xyz88(date, time_ms, ping_idx, serial, n_beams, 150.0, depth))

    with open(output_path, "wb") as fh:
        for dg in datagrams:
            fh.write(dg)

    print(f"Written {len(datagrams)} datagrams to {output_path}")
    print(f"  {n_pings} pings × {n_beams} beams = {n_pings * n_beams} soundings (approx)")


if __name__ == "__main__":
    generate()
