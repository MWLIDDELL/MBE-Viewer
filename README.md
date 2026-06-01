# MBE Viewer

Interactive visualiser for Kongsberg EM-series multibeam echosounder `.all` files.

## Features

- **Spatial depth map** — georeferenced scatter plot colour-coded by depth
- **Depth histogram** — colour-coded distribution of all soundings
- **Cross-track swath profile** — across-track distance vs depth for all pings
- Dark-themed Tkinter GUI with embedded Matplotlib toolbar (pan, zoom, save)
- Controls for colourmap, point size and decimation (for large files)
- Background file loading with progress indicator

## Supported datagrams

| ID | Type |
|----|------|
| `0x50` | Position (used for georeferencing) |
| `0x58` | XYZ-88 (primary depth source) |
| `0x44` | Depth legacy datagram |

## Requirements

```
python >= 3.11
numpy
matplotlib
scipy        # optional, not imported at runtime
tkinter      # standard library (ensure python3-tk is installed)
```

Install dependencies:

```bash
pip install numpy matplotlib scipy
```

On Debian/Ubuntu you may also need:

```bash
sudo apt install python3-tk
```

## Usage

```bash
# Interactive — opens file dialog
python3 viewer.py

# Direct load
python3 viewer.py path/to/survey.all
```

## Generating a synthetic test file

```bash
python3 generate_test_all.py   # writes test_data.all (30 pings × 256 beams)
python3 viewer.py test_data.all
```

## Running tests

```bash
python3 test_parser.py
```

## Project layout

```
MBE-Viewer/
├── all_parser.py          # Binary .all file parser
├── viewer.py              # Tkinter + Matplotlib GUI
├── generate_test_all.py   # Synthetic .all file generator (for testing)
├── test_parser.py         # Smoke tests
└── README.md
```

## Notes on georeferencing

Position datagrams (`0x50`) are matched to each ping by nearest-neighbour
interpolation in file-offset space (a proxy for time in chronologically-ordered
files). Beam across-track and along-track offsets are converted to
latitude/longitude using a flat-Earth approximation accurate to ~1 m for
typical swath widths.

If no Position datagrams are present the tool still displays relative XY data
with `latitude = 0, longitude = 0`.
