#!/usr/bin/env python3
"""Blend a detailed VTEC grid with a smooth background using station coverage.

This is an experimental post-processor. Near stations it keeps the detailed map;
far from stations it gradually falls back to a low-degree background map. It can
help suppress unphysical oscillations in oceans, but it is not part of the core
spherical-harmonic least-squares solver.
"""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path


def read_station_coords(path: Path) -> list[tuple[float, float]]:
    """Read station latitude/longitude from station_coordinates_from_tec.csv."""
    stations: list[tuple[float, float]] = []
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            stations.append((float(row["latitude"]), float(row["longitude"])))
    return stations


def angular_distance_rad(lat1: float, lon1: float, lat2: float, lon2: float) -> float:
    """Great-circle angular distance between two geographic points."""
    lat1r = math.radians(lat1)
    lat2r = math.radians(lat2)
    dlon = math.radians(lon1 - lon2)
    c = math.sin(lat1r) * math.sin(lat2r) + math.cos(lat1r) * math.cos(lat2r) * math.cos(dlon)
    return math.acos(max(-1.0, min(1.0, c)))


def coverage_weight(lat: float, lon: float, stations: list[tuple[float, float]], radius_deg: float) -> float:
    """Return a 0..1 weight based on distance to the nearest station."""
    radius = math.radians(radius_deg)
    nearest = min(angular_distance_rad(lat, lon, slat, slon) for slat, slon in stations)
    return math.exp(-0.5 * (nearest / radius) ** 2)


def read_grid(path: Path) -> dict[tuple[float, float, float], float]:
    """Read a VTEC grid CSV into a dictionary keyed by (time, lat, lon)."""
    values: dict[tuple[float, float, float], float] = {}
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            key = (float(row["time_ut"]), float(row["latitude"]), float(row["longitude"]))
            values[key] = float(row["vtec"])
    return values


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--detail-grid", required=True, type=Path)
    parser.add_argument("--background-grid", required=True, type=Path)
    parser.add_argument("--coords", default="station_coordinates_from_tec.csv", type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--radius", type=float, default=12.0)
    parser.add_argument("--floor-zero", action="store_true")
    args = parser.parse_args()

    stations = read_station_coords(args.coords)
    detail = read_grid(args.detail_grid)
    background = read_grid(args.background_grid)
    args.output.parent.mkdir(parents=True, exist_ok=True)

    # The coverage weight depends only on grid position, not on time.
    weights: dict[tuple[float, float], float] = {}
    grid_points = sorted({(lat, lon) for _, lat, lon in detail})
    for lat, lon in grid_points:
        weights[(lat, lon)] = coverage_weight(lat, lon, stations, args.radius)

    with args.output.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["time_ut", "latitude", "longitude", "vtec"])
        for key in sorted(detail):
            time_ut, lat, lon = key
            weight = weights[(lat, lon)]
            # Convex combination: detail near stations, background in sparse regions.
            value = weight * detail[key] + (1.0 - weight) * background[key]
            if args.floor_zero:
                value = max(0.0, value)
            writer.writerow([f"{time_ut:.3f}", f"{lat:.3f}", f"{lon:.3f}", f"{value:.8E}"])


if __name__ == "__main__":
    main()
