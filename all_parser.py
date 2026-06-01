"""
Kongsberg .all file parser.

Parses EM-series multibeam echosounder binary files and extracts
georeferenced XYZ sounding data from XYZ datagrams (datagram ID 88 / 0x58)
and depth datagrams (datagram ID 68 / 0x44).

References:
  Kongsberg EM Series Multibeam Echo Sounder - EM datagram formats
  (document 850-160692 / Rev.U)
"""

import struct
import math
from dataclasses import dataclass, field
from typing import List, Optional
import os


# ---------------------------------------------------------------------------
# Data structures
# ---------------------------------------------------------------------------

@dataclass
class Sounding:
    """A single georeferenced depth sounding."""
    longitude: float   # degrees
    latitude: float    # degrees
    depth: float       # metres, positive down
    across_track: float = 0.0   # metres
    along_track: float  = 0.0   # metres
    beam_number: int    = 0
    ping_counter: int   = 0


@dataclass
class ParseResult:
    soundings: List[Sounding] = field(default_factory=list)
    n_pings: int = 0
    n_position_datagrams: int = 0
    file_name: str = ""
    errors: List[str] = field(default_factory=list)


# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

DATAGRAM_TYPES = {
    0x41: "PU Status",
    0x43: "Clock",
    0x44: "Depth (beam depths)",
    0x45: "Single beam depth",
    0x46: "Raw range and angle (F)",
    0x47: "Surface sound speed",
    0x48: "Heading",
    0x49: "Installation parameters",
    0x4E: "Raw range and angle (N)",
    0x50: "Position",
    0x52: "Runtime parameters",
    0x53: "Seabed image",
    0x55: "Sound speed profile",
    0x57: "Operator station parameters",
    0x58: "XYZ 88",
    0x59: "Seabed image 89",
    0x6B: "Watercolumn",
    0x6E: "Network attitude",
    0x70: "Position (old)",
    0x73: "Scanning fish detector",
    0x78: "XYZ 120",
    0x79: "Seabed image 120",
    0x7A: "Bathy and snippet 120",
}

STX = 0x02
ETX = 0x03


# ---------------------------------------------------------------------------
# Low-level reader
# ---------------------------------------------------------------------------

def _read_uint8(buf, offset):
    return struct.unpack_from("B", buf, offset)[0]


def _read_uint16(buf, offset):
    return struct.unpack_from("<H", buf, offset)[0]


def _read_int16(buf, offset):
    return struct.unpack_from("<h", buf, offset)[0]


def _read_uint32(buf, offset):
    return struct.unpack_from("<I", buf, offset)[0]


def _read_int32(buf, offset):
    return struct.unpack_from("<i", buf, offset)[0]


# ---------------------------------------------------------------------------
# Datagram header
# ---------------------------------------------------------------------------

def _parse_header(buf, offset):
    """Returns (num_bytes, stx, datagram_type, em_model, date, time_ms, ping_counter, serial) or None."""
    if offset + 16 > len(buf):
        return None
    num_bytes   = _read_uint32(buf, offset)      # number of bytes in datagram (excl. 4-byte length field)
    stx         = _read_uint8(buf, offset + 4)
    dgm_type    = _read_uint8(buf, offset + 5)
    em_model    = _read_uint16(buf, offset + 6)
    date        = _read_uint32(buf, offset + 8)  # YYYYMMDD
    time_ms     = _read_uint32(buf, offset + 12) # milliseconds since midnight
    ping_cnt    = _read_uint16(buf, offset + 16)
    serial      = _read_uint16(buf, offset + 18)
    return num_bytes, stx, dgm_type, em_model, date, time_ms, ping_cnt, serial


# ---------------------------------------------------------------------------
# Position datagram (0x50)
# ---------------------------------------------------------------------------

def _parse_position(buf, offset, num_bytes):
    """Extract lat/lon from a Position datagram. Returns (lat_deg, lon_deg)."""
    # After the 20-byte header:
    #   int32  latitude   (1e-7 deg)
    #   int32  longitude  (1e-7 deg)
    base = offset + 20
    if base + 8 > len(buf):
        return None
    lat_raw = _read_int32(buf, base)
    lon_raw = _read_int32(buf, base + 4)
    return lat_raw * 1e-7, lon_raw * 1e-7


# ---------------------------------------------------------------------------
# XYZ 88 datagram (0x58)
# ---------------------------------------------------------------------------

def _parse_xyz88(buf, offset, num_bytes, ping_counter):
    """
    Parse XYZ 88 datagram and return list of Sounding objects.

    Datagram layout (after 20-byte common header):
      XYZ header block (20 bytes):
        float32  height of water level re vessel (m)         [+0]
        uint16   sampling / spare                            [+4]
        uint16   spare                                       [+6]
        uint16   number of valid detections                  [+8]
        uint16   sampling frequency                          [+10]
        uint32   spare                                       [+12]
        uint32   spare                                       [+16]
      -- then N beam records (20 bytes each) --
    Each beam record (20 bytes):
        float32  Z depth (m, positive down)                  [+0]
        float32  Y across-track (m, +starboard)              [+4]
        float32  X along-track  (m, +forward)                [+8]
        uint16   detection window                            [+12]
        uint8    quality                                     [+14]
        uint8    spare                                       [+15]
        uint8    beam incidence angle adj                    [+16]
        uint8    detection info (bit 7: 0=valid, 1=invalid)  [+17]
        uint8    real-time cleaning info                     [+18]
        uint8    reflectivity                                [+19]
    """
    base = offset + 20
    heading_size = 20
    if base + heading_size > len(buf):
        return []

    # Number of valid detections is at byte offset +8 within the XYZ header
    n_valid = _read_uint16(buf, base + 8)

    beam_base = base + heading_size
    beam_size = 20
    soundings = []

    end_of_datagram = offset + 4 + num_bytes
    max_beams = min(n_valid, (end_of_datagram - beam_base - 2) // beam_size)

    for i in range(max_beams):
        boff = beam_base + i * beam_size
        if boff + beam_size > len(buf):
            break
        depth  = struct.unpack_from("<f", buf, boff)[0]
        across = struct.unpack_from("<f", buf, boff + 4)[0]
        along  = struct.unpack_from("<f", buf, boff + 8)[0]
        det_info = _read_uint8(buf, boff + 17)
        # bit 7 of detection_info: 0 = valid detection
        if det_info & 0x80:
            continue
        if math.isnan(depth) or math.isinf(depth) or depth <= 0 or depth > 12000:
            continue
        soundings.append(Sounding(
            longitude=0.0, latitude=0.0,
            depth=depth,
            across_track=across,
            along_track=along,
            beam_number=i,
            ping_counter=ping_counter,
        ))
    return soundings


# ---------------------------------------------------------------------------
# Depth datagram (0x44)
# ---------------------------------------------------------------------------

def _parse_depth(buf, offset, num_bytes, ping_counter):
    """
    Parse legacy Depth datagram (0x44) and return list of Sounding objects.

    After 20-byte common header:
      uint8   heading (1/100 deg)
      uint8   sound speed (dm/s, 14.5 offset)
      uint8   transducer depth (cm)
      uint8   max beams
      uint8   valid beams
      uint8   z-resolution (cm)
      uint8   xy-resolution (cm)
      uint8   sampling rate (100 Hz units)
      -- then N beam records (6 bytes each) --
    Each beam:
      int16  depth (in z-resolution units)
      int16  across-track (in xy-resolution units)
      int16  along-track  (in xy-resolution units)
    """
    base = offset + 20
    if base + 8 > len(buf):
        return []

    # heading_raw  = _read_uint8(buf, base)
    # ss_raw       = _read_uint8(buf, base + 1)
    tx_depth_cm  = _read_uint8(buf, base + 2)
    # max_beams    = _read_uint8(buf, base + 3)
    valid_beams  = _read_uint8(buf, base + 4)
    z_res_cm     = _read_uint8(buf, base + 5)
    xy_res_cm    = _read_uint8(buf, base + 6)

    z_scale  = z_res_cm  / 100.0
    xy_scale = xy_res_cm / 100.0
    tx_depth = tx_depth_cm / 100.0

    beam_base = base + 8
    beam_size = 6
    soundings = []

    for i in range(valid_beams):
        boff = beam_base + i * beam_size
        if boff + beam_size > len(buf):
            break
        depth_raw  = _read_int16(buf, boff)
        across_raw = _read_int16(buf, boff + 2)
        along_raw  = _read_int16(buf, boff + 4)
        depth  = depth_raw * z_scale + tx_depth
        across = across_raw * xy_scale
        along  = along_raw  * xy_scale
        if abs(depth) > 12000:
            continue
        soundings.append(Sounding(
            longitude=0.0, latitude=0.0,
            depth=depth,
            across_track=across,
            along_track=along,
            beam_number=i,
            ping_counter=ping_counter,
        ))
    return soundings


# ---------------------------------------------------------------------------
# Main parser
# ---------------------------------------------------------------------------

def parse_all_file(filepath: str, max_soundings: int = 2_000_000) -> ParseResult:
    """
    Parse a Kongsberg .all file and return a ParseResult containing soundings
    with approximate lat/lon positions interpolated from Position datagrams.
    """
    result = ParseResult(file_name=os.path.basename(filepath))

    with open(filepath, "rb") as fh:
        buf = fh.read()

    n = len(buf)
    offset = 0

    # We accumulate (ping_counter -> [soundings]) then marry with positions
    ping_soundings: dict = {}
    # Position interpolation: list of (byte_offset, lat, lon)
    pos_records: list = []

    while offset < n - 4:
        # Find STX
        if offset + 5 >= n:
            break

        # Read the 4-byte length prefix
        try:
            num_bytes = _read_uint32(buf, offset)
        except Exception:
            offset += 1
            continue

        if num_bytes < 6 or offset + 4 + num_bytes > n:
            offset += 1
            continue

        stx = _read_uint8(buf, offset + 4)
        if stx != STX:
            offset += 1
            continue

        # Parse header
        hdr = _parse_header(buf, offset)
        if hdr is None:
            offset += 1
            continue

        num_bytes, stx, dgm_type, em_model, date, time_ms, ping_cnt, serial = hdr

        datagram_end = offset + 4 + num_bytes

        try:
            if dgm_type == 0x50:  # Position
                pos = _parse_position(buf, offset, num_bytes)
                if pos is not None:
                    pos_records.append((offset, pos[0], pos[1]))
                    result.n_position_datagrams += 1

            elif dgm_type == 0x58:  # XYZ 88
                snds = _parse_xyz88(buf, offset, num_bytes, ping_cnt)
                if snds:
                    ping_soundings.setdefault(ping_cnt, []).extend(snds)
                    result.n_pings += 1

            elif dgm_type == 0x44:  # Depth
                snds = _parse_depth(buf, offset, num_bytes, ping_cnt)
                if snds:
                    ping_soundings.setdefault(ping_cnt, []).extend(snds)
                    result.n_pings += 1

        except Exception as e:
            result.errors.append(f"offset={offset} type=0x{dgm_type:02X}: {e}")

        offset = datagram_end

    # ------------------------------------------------------------------
    # Assign positions to soundings.
    # Strategy: each ping is assigned the nearest position record by
    # file-offset proximity (a crude but effective proxy for time when
    # the file is written in chronological order).
    # ------------------------------------------------------------------
    if not pos_records:
        # No position data — assign dummy 0,0 and return as-is (relative XY only)
        result.errors.append(
            "No Position datagrams found; soundings will have lat=0 lon=0. "
            "Visualization will use relative across/along-track coordinates."
        )
        for snds in ping_soundings.values():
            result.soundings.extend(snds)
            if len(result.soundings) >= max_soundings:
                break
        return result

    pos_offsets = [p[0] for p in pos_records]

    def nearest_pos(target_offset):
        # Binary search for nearest position record
        lo, hi = 0, len(pos_offsets) - 1
        while lo < hi:
            mid = (lo + hi) // 2
            if pos_offsets[mid] < target_offset:
                lo = mid + 1
            else:
                hi = mid
        # lo is insertion point; check lo-1 and lo
        best = lo
        if lo > 0 and abs(pos_offsets[lo - 1] - target_offset) < abs(pos_offsets[lo] - target_offset):
            best = lo - 1
        return pos_records[best][1], pos_records[best][2]

    # Build a mapping from ping_counter -> approximate byte offset.
    # We re-scan the buffer lightly to record XYZ88/Depth datagram positions.
    ping_offsets: dict = {}
    offset2 = 0
    while offset2 < n - 4:
        try:
            nb = _read_uint32(buf, offset2)
        except Exception:
            break
        if nb < 6 or offset2 + 4 + nb > n:
            offset2 += 1
            continue
        stx2 = _read_uint8(buf, offset2 + 4)
        if stx2 != STX:
            offset2 += 1
            continue
        dgm2 = _read_uint8(buf, offset2 + 5)
        if dgm2 in (0x58, 0x44):
            pc = _read_uint16(buf, offset2 + 16)
            if pc not in ping_offsets:
                ping_offsets[pc] = offset2
        offset2 += 4 + nb

    total = 0
    for ping_cnt, snds in ping_soundings.items():
        poff = ping_offsets.get(ping_cnt, 0)
        lat, lon = nearest_pos(poff)
        # Convert across/along-track offsets to approximate lat/lon deltas.
        # 1 degree latitude ≈ 111_111 m; 1 degree longitude ≈ 111_111 * cos(lat) m
        cos_lat = math.cos(math.radians(lat)) or 1e-9
        for s in snds:
            s.latitude  = lat  + s.along_track  / 111_111.0
            s.longitude = lon  + s.across_track / (111_111.0 * cos_lat)
            result.soundings.append(s)
            total += 1
            if total >= max_soundings:
                return result

    return result
