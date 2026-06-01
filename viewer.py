"""
MBE Viewer — Kongsberg .all file visualiser.

Usage:
    python3 viewer.py [path/to/file.all]

Opens a Tkinter window with:
  • A spatial depth map (scatter plot, colour-coded by depth)
  • A depth histogram
  • A cross-track swath profile for the selected ping range
  • File info panel
  • Colourmap / decimation controls
"""

import sys
import os
import math
import threading
import tkinter as tk
from tkinter import ttk, filedialog, messagebox

import numpy as np
import matplotlib
matplotlib.use("TkAgg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from matplotlib.colors import Normalize
import matplotlib.cm as cm

from all_parser import parse_all_file, ParseResult


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

COLORMAPS = ["viridis_r", "plasma_r", "cividis_r", "jet", "coolwarm_r", "Blues", "deep_r"]


def _depth_stats(depths):
    if len(depths) == 0:
        return 0, 0, 0
    return float(np.min(depths)), float(np.max(depths)), float(np.mean(depths))


# ---------------------------------------------------------------------------
# Main application
# ---------------------------------------------------------------------------

class MBEViewerApp:
    def __init__(self, root: tk.Tk):
        self.root = root
        self.root.title("MBE Viewer — Kongsberg .all")
        self.root.geometry("1280x800")
        self.root.configure(bg="#1e1e2e")

        self.result: ParseResult | None = None
        self._xs = np.empty(0)
        self._ys = np.empty(0)
        self._zs = np.empty(0)

        self._build_ui()

        # If a file was passed on the command line, load it straight away
        if len(sys.argv) > 1 and os.path.isfile(sys.argv[1]):
            self.root.after(200, lambda: self._load_file(sys.argv[1]))

    # ------------------------------------------------------------------
    # UI construction
    # ------------------------------------------------------------------

    def _build_ui(self):
        style = ttk.Style()
        style.theme_use("clam")
        style.configure("TFrame", background="#1e1e2e")
        style.configure("TLabel", background="#1e1e2e", foreground="#cdd6f4")
        style.configure("TButton", background="#313244", foreground="#cdd6f4", relief="flat", padding=4)
        style.map("TButton", background=[("active", "#45475a")])
        style.configure("TCombobox", fieldbackground="#313244", background="#313244", foreground="#cdd6f4")
        style.configure("TScale", background="#1e1e2e", troughcolor="#313244")
        style.configure("Horizontal.TProgressbar", troughcolor="#313244", background="#89b4fa")

        # ---- Top toolbar ------------------------------------------------
        toolbar = ttk.Frame(self.root, padding=4)
        toolbar.pack(side=tk.TOP, fill=tk.X)

        ttk.Button(toolbar, text="Open .all file…", command=self._open_file).pack(side=tk.LEFT, padx=4)

        ttk.Label(toolbar, text="Colourmap:").pack(side=tk.LEFT, padx=(12, 2))
        self._cmap_var = tk.StringVar(value="viridis_r")
        cmap_cb = ttk.Combobox(toolbar, textvariable=self._cmap_var, values=COLORMAPS,
                               width=14, state="readonly")
        cmap_cb.pack(side=tk.LEFT)
        cmap_cb.bind("<<ComboboxSelected>>", lambda e: self._redraw())

        ttk.Label(toolbar, text="Max points:").pack(side=tk.LEFT, padx=(12, 2))
        self._max_pts_var = tk.StringVar(value="200000")
        pts_entry = ttk.Combobox(toolbar, textvariable=self._max_pts_var,
                                 values=["10000", "50000", "100000", "200000", "500000", "all"],
                                 width=8, state="readonly")
        pts_entry.pack(side=tk.LEFT)
        pts_entry.bind("<<ComboboxSelected>>", lambda e: self._redraw())

        ttk.Label(toolbar, text="Point size:").pack(side=tk.LEFT, padx=(12, 2))
        self._pt_size = tk.DoubleVar(value=1.5)
        size_spin = ttk.Spinbox(toolbar, from_=0.5, to=10.0, increment=0.5,
                                textvariable=self._pt_size, width=5,
                                command=self._redraw)
        size_spin.pack(side=tk.LEFT)

        self._status_var = tk.StringVar(value="No file loaded.")
        ttk.Label(toolbar, textvariable=self._status_var, foreground="#a6e3a1").pack(
            side=tk.RIGHT, padx=8)

        # ---- Progress bar (hidden until loading) ------------------------
        self._progress = ttk.Progressbar(self.root, mode="indeterminate", length=200)

        # ---- Info panel -------------------------------------------------
        info_frame = ttk.Frame(self.root, padding=(8, 2))
        info_frame.pack(side=tk.TOP, fill=tk.X)
        self._info_var = tk.StringVar(value="")
        ttk.Label(info_frame, textvariable=self._info_var, font=("Courier", 9),
                  foreground="#89dceb").pack(side=tk.LEFT)

        # ---- Figure area ------------------------------------------------
        self._fig = plt.Figure(figsize=(14, 7), facecolor="#1e1e2e")
        self._fig.subplots_adjust(left=0.07, right=0.97, top=0.93, bottom=0.1,
                                  wspace=0.35, hspace=0.45)

        # Axes
        self._ax_map    = self._fig.add_subplot(1, 2, 1)
        self._ax_hist   = self._fig.add_subplot(2, 2, 2)
        self._ax_swath  = self._fig.add_subplot(2, 2, 4)

        for ax in (self._ax_map, self._ax_hist, self._ax_swath):
            ax.set_facecolor("#181825")
            for spine in ax.spines.values():
                spine.set_edgecolor("#45475a")
            ax.tick_params(colors="#cdd6f4", labelsize=8)
            ax.xaxis.label.set_color("#cdd6f4")
            ax.yaxis.label.set_color("#cdd6f4")
            ax.title.set_color("#cdd6f4")

        canvas = FigureCanvasTkAgg(self._fig, master=self.root)
        canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)

        nav_frame = tk.Frame(self.root, bg="#1e1e2e")
        nav_frame.pack(side=tk.BOTTOM, fill=tk.X)
        toolbar2 = NavigationToolbar2Tk(canvas, nav_frame)
        toolbar2.config(background="#1e1e2e")
        toolbar2.update()

        self._canvas = canvas
        self._draw_empty_plots()

    # ------------------------------------------------------------------
    # File I/O
    # ------------------------------------------------------------------

    def _open_file(self):
        path = filedialog.askopenfilename(
            title="Open Kongsberg .all file",
            filetypes=[("Kongsberg ALL files", "*.all *.ALL"), ("All files", "*.*")]
        )
        if path:
            self._load_file(path)

    def _load_file(self, path: str):
        self._status_var.set(f"Loading {os.path.basename(path)}…")
        self._progress.pack(side=tk.TOP, fill=tk.X)
        self._progress.start(10)

        def worker():
            try:
                result = parse_all_file(path)
                self.root.after(0, lambda: self._on_load_complete(result))
            except Exception as exc:
                self.root.after(0, lambda: self._on_load_error(str(exc)))

        threading.Thread(target=worker, daemon=True).start()

    def _on_load_complete(self, result: ParseResult):
        self._progress.stop()
        self._progress.pack_forget()
        self.result = result

        n = len(result.soundings)
        if n == 0:
            self._status_var.set("No soundings found.")
            messagebox.showwarning("MBE Viewer", "No depth soundings could be extracted from the file.")
            return

        self._xs = np.array([s.longitude   for s in result.soundings], dtype=np.float64)
        self._ys = np.array([s.latitude    for s in result.soundings], dtype=np.float64)
        self._zs = np.array([s.depth       for s in result.soundings], dtype=np.float64)

        dmin, dmax, dmean = _depth_stats(self._zs)
        self._info_var.set(
            f"File: {result.file_name}  |  "
            f"Soundings: {n:,}  |  Pings: {result.n_pings:,}  |  "
            f"Depth: {dmin:.1f} – {dmax:.1f} m  (mean {dmean:.1f} m)  |  "
            f"Positions: {result.n_position_datagrams:,}"
        )
        self._status_var.set(f"Loaded  {n:,} soundings from {result.file_name}")
        if result.errors:
            print(f"[Parser warnings] {len(result.errors)} issue(s):")
            for e in result.errors[:10]:
                print(" ", e)

        self._redraw()

    def _on_load_error(self, msg: str):
        self._progress.stop()
        self._progress.pack_forget()
        self._status_var.set("Error loading file.")
        messagebox.showerror("MBE Viewer", f"Failed to parse file:\n{msg}")

    # ------------------------------------------------------------------
    # Drawing
    # ------------------------------------------------------------------

    def _get_decimated(self):
        max_str = self._max_pts_var.get()
        if max_str == "all":
            return self._xs, self._ys, self._zs
        max_pts = int(max_str)
        n = len(self._xs)
        if n <= max_pts:
            return self._xs, self._ys, self._zs
        idx = np.random.choice(n, max_pts, replace=False)
        return self._xs[idx], self._ys[idx], self._zs[idx]

    def _draw_empty_plots(self):
        for ax, title in [(self._ax_map, "Spatial Depth Map"),
                          (self._ax_hist, "Depth Histogram"),
                          (self._ax_swath, "Cross-track Swath Profile")]:
            ax.cla()
            ax.set_facecolor("#181825")
            ax.set_title(title, fontsize=9, color="#cdd6f4")
            ax.text(0.5, 0.5, "No data", ha="center", va="center",
                    transform=ax.transAxes, color="#585b70", fontsize=11)
        self._canvas.draw_idle()

    def _redraw(self):
        if self.result is None or len(self._xs) == 0:
            self._draw_empty_plots()
            return

        xs, ys, zs = self._get_decimated()
        cmap_name = self._cmap_var.get()
        pt_size   = float(self._pt_size.get())

        zmin, zmax = float(np.nanmin(zs)), float(np.nanmax(zs))
        norm = Normalize(vmin=zmin, vmax=zmax)
        cmap = plt.get_cmap(cmap_name)

        # ------ Map -------------------------------------------------------
        ax = self._ax_map
        ax.cla()
        ax.set_facecolor("#181825")
        sc = ax.scatter(xs, ys, c=zs, cmap=cmap, norm=norm,
                        s=pt_size, linewidths=0, rasterized=True, zorder=2)

        # Colourbar
        try:
            self._cbar.remove()
        except Exception:
            pass
        self._cbar = self._fig.colorbar(sc, ax=ax, fraction=0.03, pad=0.01)
        self._cbar.set_label("Depth (m)", color="#cdd6f4", fontsize=8)
        self._cbar.ax.yaxis.set_tick_params(color="#cdd6f4", labelsize=7)
        plt.setp(self._cbar.ax.yaxis.get_ticklabels(), color="#cdd6f4")

        ax.set_title("Spatial Depth Map", fontsize=9, color="#cdd6f4")
        ax.set_xlabel("Longitude (°)", fontsize=8)
        ax.set_ylabel("Latitude (°)", fontsize=8)
        ax.tick_params(colors="#cdd6f4", labelsize=7)
        ax.set_facecolor("#181825")
        for sp in ax.spines.values():
            sp.set_edgecolor("#45475a")

        # Grid
        ax.grid(True, color="#313244", linewidth=0.4, zorder=1)
        ax.set_axisbelow(True)

        # Format tick labels as degrees
        ax.xaxis.set_major_formatter(plt.FuncFormatter(lambda v, _: f"{v:.4f}°"))
        ax.yaxis.set_major_formatter(plt.FuncFormatter(lambda v, _: f"{v:.4f}°"))
        ax.tick_params(axis="x", rotation=30)

        # ------ Histogram -------------------------------------------------
        ax2 = self._ax_hist
        ax2.cla()
        ax2.set_facecolor("#181825")
        n_bins = min(80, max(20, len(np.unique(np.round(zs, 0)))))
        counts, edges, patches = ax2.hist(zs, bins=n_bins, orientation="horizontal",
                                          color="#89b4fa", edgecolor="none", alpha=0.85)
        # Colour patches by depth
        bin_centers = 0.5 * (edges[:-1] + edges[1:])
        for patch, d in zip(patches, bin_centers):
            patch.set_facecolor(cmap(norm(d)))

        ax2.set_xlabel("Count", fontsize=8, color="#cdd6f4")
        ax2.set_ylabel("Depth (m)", fontsize=8, color="#cdd6f4")
        ax2.set_title("Depth Histogram", fontsize=9, color="#cdd6f4")
        ax2.tick_params(colors="#cdd6f4", labelsize=7)
        ax2.invert_yaxis()
        for sp in ax2.spines.values():
            sp.set_edgecolor("#45475a")

        # ------ Swath profile (all soundings, across-track vs depth) ------
        ax3 = self._ax_swath
        ax3.cla()
        ax3.set_facecolor("#181825")

        ats = np.array([s.across_track for s in self.result.soundings], dtype=np.float32)
        deps = self._zs

        # Decimate for swath plot too
        n_swath = min(len(ats), 100_000)
        if len(ats) > n_swath:
            idx2 = np.random.choice(len(ats), n_swath, replace=False)
            ats2, deps2 = ats[idx2], deps[idx2]
        else:
            ats2, deps2 = ats, deps

        ax3.scatter(ats2, deps2, c=deps2, cmap=cmap, norm=norm,
                    s=0.5, linewidths=0, rasterized=True, alpha=0.6)
        ax3.set_xlabel("Across-track (m)", fontsize=8, color="#cdd6f4")
        ax3.set_ylabel("Depth (m)", fontsize=8, color="#cdd6f4")
        ax3.set_title("Cross-track Swath Profile", fontsize=9, color="#cdd6f4")
        ax3.tick_params(colors="#cdd6f4", labelsize=7)
        ax3.invert_yaxis()
        for sp in ax3.spines.values():
            sp.set_edgecolor("#45475a")

        self._canvas.draw_idle()


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

def main():
    root = tk.Tk()
    app = MBEViewerApp(root)
    root.mainloop()


if __name__ == "__main__":
    main()
