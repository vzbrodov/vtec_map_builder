#!/usr/bin/env python3
"""Plot VTEC maps and optionally build a daily animation from Fortran CSV output.

The Fortran program writes a long CSV table:
time_ut, latitude, longitude, vtec.
This script reshapes each epoch to a 2-D lon/lat image and saves a PNG/GIF.
"""

from __future__ import annotations

import argparse
import csv
import os
from pathlib import Path

# In some sandboxed/editor environments the default matplotlib config directory is
# not writable. Point it to /tmp so plotting works without extra setup.
os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib-vtec")

import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter

try:
    import cartopy.crs as ccrs
    import cartopy.feature as cfeature
except Exception:  # Cartopy is optional; a plain lon/lat plot still works.
    ccrs = None
    cfeature = None


def read_grid(path: Path) -> tuple[dict[float, dict[tuple[float, float], float]], list[float], list[float], list[float]]:
    """Read the long grid CSV into dictionaries indexed by time and grid point."""
    data: dict[float, dict[tuple[float, float], float]] = {}
    latitudes: set[float] = set()
    longitudes: set[float] = set()

    with path.open(newline="") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            time_ut = float(row["time_ut"])
            lat = float(row["latitude"])
            lon = float(row["longitude"])
            vtec = float(row["vtec"])
            data.setdefault(time_ut, {})[(lat, lon)] = vtec
            latitudes.add(lat)
            longitudes.add(lon)

    return data, sorted(data), sorted(latitudes), sorted(longitudes)


def frame_to_array(data, time_ut: float, latitudes: list[float], longitudes: list[float]):
    """Convert one epoch from sparse dictionary form to a 2-D image array."""
    frame = data[time_ut]
    return [[frame.get((lat, lon), float("nan")) for lon in longitudes] for lat in latitudes]


def make_axes(use_cartopy: bool):
    """Create map axes. Cartopy is used when installed; otherwise use plain lon/lat."""
    if use_cartopy:
        fig = plt.figure(figsize=(12, 6))
        ax = plt.axes(projection=ccrs.PlateCarree())
        ax.set_global()
        ax.coastlines(linewidth=0.8)
        ax.add_feature(cfeature.BORDERS, linewidth=0.3)
        ax.gridlines(draw_labels=True, linewidth=0.3, alpha=0.5)
        return fig, ax
    fig, ax = plt.subplots(figsize=(12, 6))
    ax.set_xlabel("Longitude, deg")
    ax.set_ylabel("Latitude, deg")
    ax.grid(True, linewidth=0.3, alpha=0.5)
    return fig, ax


def plot_one(data, time_ut, latitudes, longitudes, output: Path, use_cartopy: bool, vmin: float | None, vmax: float | None, interpolation: str):
    """Save one PNG map for the requested UT epoch."""
    values = frame_to_array(data, time_ut, latitudes, longitudes)
    fig, ax = make_axes(use_cartopy)
    image_args = dict(
        extent=[min(longitudes), max(longitudes), min(latitudes), max(latitudes)],
        origin="lower",
        cmap="turbo",
        aspect="auto",
        vmin=vmin,
        vmax=vmax,
        interpolation=interpolation,
    )
    if use_cartopy:
        mesh = ax.imshow(values, transform=ccrs.PlateCarree(), **image_args)
    else:
        mesh = ax.imshow(values, **image_args)
    ax.set_title(f"VTEC, UT={time_ut:.2f} h")
    fig.colorbar(mesh, ax=ax, label="TECU", shrink=0.82)
    fig.tight_layout()
    fig.savefig(output, dpi=160)
    plt.close(fig)


def animate(data, times, latitudes, longitudes, output: Path, use_cartopy: bool, vmin: float | None, vmax: float | None, interpolation: str):
    """Save a GIF animation over all UT epochs."""
    fig, ax = make_axes(use_cartopy)
    first = frame_to_array(data, times[0], latitudes, longitudes)
    image_args = dict(
        extent=[min(longitudes), max(longitudes), min(latitudes), max(latitudes)],
        origin="lower",
        cmap="turbo",
        aspect="auto",
        vmin=vmin,
        vmax=vmax,
        interpolation=interpolation,
    )
    if use_cartopy:
        image = ax.imshow(first, transform=ccrs.PlateCarree(), **image_args)
    else:
        image = ax.imshow(first, **image_args)
    title = ax.set_title("")
    fig.colorbar(image, ax=ax, label="TECU", shrink=0.82)

    def update(i):
        values = frame_to_array(data, times[i], latitudes, longitudes)
        image.set_data(values)
        title.set_text(f"VTEC, UT={times[i]:.2f} h")
        return image, title

    anim = FuncAnimation(fig, update, frames=len(times), interval=250, blit=False)
    anim.save(output, writer=PillowWriter(fps=4))
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--grid", default="output/vtec_grid.csv", type=Path)
    parser.add_argument("--out-dir", default="output/figures", type=Path)
    parser.add_argument("--time", type=float, default=None, help="UT hour for a single PNG")
    parser.add_argument("--vmin", type=float, default=0.0, help="lower color scale limit; use nan for autoscale")
    parser.add_argument("--vmax", type=float, default=100.0, help="upper color scale limit; use nan for autoscale")
    parser.add_argument("--interpolation", default="bilinear", help="matplotlib image interpolation mode")
    parser.add_argument("--animate", action="store_true", help="write GIF animation")
    args = parser.parse_args()

    args.out_dir.mkdir(parents=True, exist_ok=True)
    data, times, latitudes, longitudes = read_grid(args.grid)
    use_cartopy = ccrs is not None
    # argparse cannot parse None directly, so "nan" is used as a small escape hatch
    # for matplotlib autoscaling.
    vmin = None if args.vmin != args.vmin else args.vmin
    vmax = None if args.vmax != args.vmax else args.vmax

    selected_time = times[0] if args.time is None else min(times, key=lambda t: abs(t - args.time))
    plot_one(
        data,
        selected_time,
        latitudes,
        longitudes,
        args.out_dir / f"vtec_{selected_time:05.2f}ut.png",
        use_cartopy,
        vmin,
        vmax,
        args.interpolation,
    )
    if args.animate:
        animate(data, times, latitudes, longitudes, args.out_dir / "vtec_day.gif", use_cartopy, vmin, vmax, args.interpolation)


if __name__ == "__main__":
    main()
