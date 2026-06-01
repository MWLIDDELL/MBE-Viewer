"""
Generate a minimal synthetic Kongsberg .all file for testing.

Creates:
  - Installation parameters datagram (0x49)
  - N Position datagrams (0x50) along a straight track
  - N XYZ-88 datagrams (0x58) with n_beams each (swath pattern)

Beam record layout matches real EM2040 files (verified against hardware data):
  bytes  0-3:   float32 Z  depth      (m, positive down)
  bytes  4-7:   float32 Y  across-track (m, +starboard)
  bytes  8-11:  float32 X  along-track  (m, +forward)
  bytes 12-15:  uint16 detection_window + uint8 quality + uint8 spare
  bytes 16-19:  uint8 beam_angle_adj + uint8 detection_info + uint8 rt_clean + uint8 reflectivity

XYZ header block (20 bytes):
  bytes  0-3:   float32 height of water level (~0)
  bytes  4-7:   spare
  bytes  8-9:   uint16 n_valid (number of valid beams)
  bytes 10-19:  spare

Output: test_data.all
"""

import struct
import math
import random
import os


STX = 0x02
ETX = 0x03


def _u8(v):   return struct.pack("B", int(v) & 0xFF)
def _u16(v):  return struct.pack("<H", int(v) & 0xFFFF)
def _i16(v):  return struct.pack("<h", int(v))
def _u32(v):  return struct.pack("<I", int(v) & 0xFFFFFFFF)
def _i32(v):  return struct.pack("<i", int(v))
def _f32(v):  return struct.pack("<f", float(v))


def make_datagram(dgm_type, em_model, date, time_ms, ping_counter, serial, payload):
    header = (
        _u8(STX) + _u8(dgm_type) + _u16(em_model) +
        _u32(date) + _u32(time_ms) + _u16(ping_counter) + _u16(serial)
    )
    body = header + payload
    checksum = sum(body) & 0xFFFF
    body += _u8(ETX) + _u16(checksum)
    return _u32(len(body)) + body


def make_position(date, time_ms, ping_counter, serial, lat_deg, lon_deg):
    payload = (
        _i32(int(round(lat_deg * 1e7))) +
        _i32(int(round(lon_deg * 1e7))) +
        _u16(0) + _u16(0) + _u16(0) + _u16(0) +
        _u8(0) + _u8(0)
    )
    return make_datagram(0x50, 2040, date, time_ms, ping_counter, serial, payload)


def make_xyz88(date, time_ms, ping_counter, serial, n_beams, swath_half_width, base_depth):
    """Build an XYZ-88 datagram matching the real EM2040 field layout."""
    # XYZ header block (20 bytes)
    xyz_header = bytearray(20)
    struct.pack_into("<f", xyz_header, 0, 0.0)          # height of water level
    struct.pack_into("<H", xyz_header, 8, n_beams)      # n_valid at offset +8

    beam_records = bytearray()
    for i in range(n_beams):
        # across-track: spread from -swath_half_width to +swath_half_width
        t = i / (n_beams - 1) - 0.5          # -0.5 to +0.5
        across = t * 2 * swath_half_width     # meters, + = starboard
        along  = 3.5                          # typical inter-ping forward movement
        # depth: slightly deeper at outer beams (realistic bowl)
        depth = base_depth * (1.0 + 0.08 * t ** 2) + random.gauss(0, 0.2)

        rec = bytearray(20)
        struct.pack_into("<f", rec, 0,  depth)   # Z
        struct.pack_into("<f", rec, 4,  across)  # Y
        struct.pack_into("<f", rec, 8,  along)   # X
        # bytes 12-15: detection window (uint16=0) + quality (uint8=1) + spare
        struct.pack_into("<H", rec, 12, 0)
        rec[14] = 1   # quality
        # bytes 16-19: beam_angle_adj, detection_info(0=valid), rt_clean, reflectivity
        rec[16] = 0
        rec[17] = 0   # detection_info = 0 → valid
        rec[18] = 0
        rec[19] = 128  # arbitrary reflectivity

        beam_records += rec

    spare = b"\x00"
    payload = bytes(xyz_header) + bytes(beam_records) + spare
    return make_datagram(0x58, 2040, date, time_ms, ping_counter, serial, payload)


def generate(output_path="test_data.all", n_pings=30, n_beams=256,
             start_lat=10.536, start_lon=0.289, base_depth=35.0):
    random.seed(42)
    date   = 20121205
    serial = 201

    datagrams = []
    track_len = 0.002  # degrees latitude

    for ping_idx in range(n_pings):
        t_frac   = ping_idx / max(n_pings - 1, 1)
        lat      = start_lat + t_frac * track_len
        time_ms  = 60_000 + ping_idx * 2000

        datagrams.append(make_position(date, time_ms, ping_idx, serial, lat, start_lon))

        depth = base_depth + 5 * math.sin(t_frac * math.pi)
        datagrams.append(make_xyz88(date, time_ms, ping_idx, serial, n_beams, 100.0, depth))

    with open(output_path, "wb") as fh:
        for dg in datagrams:
            fh.write(dg)

    print(f"Written {len(datagrams)} datagrams to {output_path}")
    print(f"  {n_pings} pings × {n_beams} beams = {n_pings * n_beams} soundings")


if __name__ == "__main__":
    generate()
