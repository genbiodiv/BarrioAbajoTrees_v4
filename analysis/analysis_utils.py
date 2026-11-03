#!/usr/bin/env python3
"""Shared data, spatial, and plotting utilities for the thematic workflow."""

from __future__ import annotations

from dataclasses import dataclass
import json
import math
import os
from pathlib import Path
import shutil
import tempfile
from typing import Iterable, Optional, Sequence, Tuple

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.collections import PatchCollection
from matplotlib.colors import Normalize
from matplotlib.patches import Polygon as MplPolygon
import numpy as np
import pandas as pd
from pyproj import Transformer
from scipy.special import gammaln
import shapefile
from shapely.geometry import Point, box, shape
from shapely.ops import nearest_points, transform as shapely_transform, unary_union


ROOT = Path(__file__).resolve().parents[1]
DATABASE = ROOT / "data" / "Data_S1_Barrio_Abajo_tree_inventory.xlsx"
DATA_DIR = ROOT / "data"
SPATIAL_DIR = DATA_DIR / "spatial"

WGS84 = 4326
UTM18N = 32618
PROJECT_TO_UTM = Transformer.from_crs(WGS84, UTM18N, always_xy=True)
PROJECT_TO_WGS84 = Transformer.from_crs(UTM18N, WGS84, always_xy=True)

PERIMETER_FILE = SPATIAL_DIR / "Perimetro del Barrio abajo.shp"
SECTOR_FILE = SPATIAL_DIR / "sectores.shp"
BLOCK_FILE = SPATIAL_DIR / "Manzanas Barrio Abajo.shp"
ROAD_CENTERLINES_FILE = DATA_DIR / "barrio_abajo_roads.geojson"


@dataclass
class FigureData:
    """Canonical census and projected geometries used by every figure."""

    trees: pd.DataFrame
    traits: pd.DataFrame
    blocks: pd.DataFrame
    perimeter: object
    sectors: list
    sector_outline: object
    grid: pd.DataFrame
    bounds: Tuple[float, float, float, float]


def configure_style() -> None:
    """Set a compact journal-style Matplotlib theme."""
    mpl.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 8.5,
            "axes.titlesize": 9.5,
            "axes.labelsize": 8.5,
            "axes.linewidth": 0.7,
            "xtick.labelsize": 7.5,
            "ytick.labelsize": 7.5,
            "legend.fontsize": 7.5,
            "legend.title_fontsize": 8,
            "figure.facecolor": "white",
            "savefig.facecolor": "white",
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )


def expected_rarefied_richness(counts: pd.Series, sample_size: int) -> float:
    """Return expected richness for an individual-based rarefied sample."""
    counts = counts.astype(int)
    total = int(counts.sum())
    if total < sample_size or sample_size <= 0:
        return float("nan")

    def log_choose(n: int, k: int) -> float:
        if k < 0 or k > n:
            return float("-inf")
        return float(gammaln(n + 1) - gammaln(k + 1) - gammaln(n - k + 1))

    denominator = log_choose(total, sample_size)
    richness = 0.0
    for abundance in counts:
        remaining = total - int(abundance)
        absence_probability = 0.0
        if remaining >= sample_size:
            absence_probability = math.exp(
                log_choose(remaining, sample_size) - denominator
            )
        richness += 1.0 - absence_probability
    return richness


def hill_numbers(species: pd.Series) -> Tuple[int, float, float]:
    """Return q0, q1, and q2 as effective numbers of species."""
    proportions = species.value_counts(normalize=True)
    q0 = int(len(proportions))
    q1 = float(np.exp(-(proportions * np.log(proportions)).sum()))
    q2 = float(1.0 / np.square(proportions).sum())
    return q0, q1, q2


def load_data() -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Read and validate the two source sheets used by the analysis."""
    trees = pd.read_excel(DATABASE, sheet_name="ArbolesIndividuales_OK")
    traits = pd.read_excel(DATABASE, sheet_name="CaracteristicasEspecies")

    required = {
        "ID",
        "Species_name",
        "Genus",
        "Family",
        "X",
        "Y",
        "Sector",
        "Manzana",
        "TRUNK",
        "AVGDBH",
        "ESTH",
    }
    missing = required - set(trees.columns)
    if missing:
        raise ValueError(f"Missing required individual columns: {sorted(missing)}")
    if trees["ID"].duplicated().any():
        raise ValueError("The individual-tree sheet contains duplicate IDs")

    trait_columns = ["Species_name", "Origin", "Form"]
    if traits["Species_name"].duplicated().any():
        raise ValueError("The species-trait sheet contains duplicate species names")
    trees = trees.merge(
        traits[trait_columns], on="Species_name", how="left", validate="many_to_one"
    )
    unmatched = trees[trait_columns[1:]].isna().all(axis=1)
    if unmatched.any():
        names = sorted(trees.loc[unmatched, "Species_name"].unique())
        raise ValueError(f"Species without trait records: {names}")

    trees["Sector"] = trees["Sector"].astype(str).str.strip()
    trees["Manzana"] = pd.to_numeric(trees["Manzana"], errors="raise").astype(int)
    trees["DBH_cm"] = pd.to_numeric(trees["AVGDBH"], errors="coerce")
    trees["stem_count"] = pd.to_numeric(trees["TRUNK"], errors="coerce")
    trees["height_m"] = pd.to_numeric(trees["ESTH"], errors="coerce")
    trees.loc[trees["height_m"] <= 0, "height_m"] = np.nan
    trees["basal_area_m2"] = (
        trees["stem_count"] * np.pi * (trees["DBH_cm"] / 200.0) ** 2
    )
    return trees, traits


def read_shapes(
    path: Path, source_crs: int, target_crs: Optional[int] = None
) -> Tuple[list, list]:
    """Read a shapefile and optionally transform all non-null geometries."""
    reader = shapefile.Reader(str(path), encoding="utf-8")
    fields = [field[0] for field in reader.fields[1:]]
    converter = None
    if target_crs is not None and source_crs != target_crs:
        converter = Transformer.from_crs(
            source_crs, target_crs, always_xy=True
        ).transform

    geometries = []
    records = []
    for item in reader.iterShapeRecords():
        if item.shape.shapeType == shapefile.NULL:
            continue
        geometry = shape(item.shape.__geo_interface__)
        if converter is not None:
            geometry = shapely_transform(converter, geometry)
        geometries.append(geometry)
        records.append(dict(zip(fields, list(item.record))))
    return geometries, records


def load_all_blocks() -> pd.DataFrame:
    """Read cadastral blocks and correct the known G-40 geometry label."""
    geometries, records = read_shapes(BLOCK_FILE, WGS84, UTM18N)
    rows = []
    for geometry, record in zip(geometries, records):
        sector = str(record["Sector"]).strip()
        block_number = int(record["Manzana"])
        cadastral_code = str(record["COD_CATAST"])
        if sector == "G" and cadastral_code.endswith("0040"):
            block_number = 40
        rows.append(
            {
                "Sector": sector,
                "Manzana": block_number,
                "geometry": geometry,
                "COD_CATAST": cadastral_code,
            }
        )
    blocks = pd.DataFrame(rows)
    duplicate = blocks.duplicated(["Sector", "Manzana"], keep=False)
    if duplicate.any():
        pairs = blocks.loc[duplicate, ["Sector", "Manzana"]].to_dict("records")
        raise ValueError(f"Duplicate corrected block labels: {pairs}")
    return blocks


def build_grid(
    bounds_geometry: object,
    clip_geometry: object,
    n_x: int = 20,
    n_y: int = 20,
) -> pd.DataFrame:
    """Create the retained regular 20 by 20 grid and clip it to the barrio."""
    min_x, min_y, max_x, max_y = bounds_geometry.bounds
    width = (max_x - min_x) / n_x
    height = (max_y - min_y) / n_y
    rows = []
    for row in range(n_y):
        for column in range(n_x):
            raw = box(
                min_x + column * width,
                min_y + row * height,
                min_x + (column + 1) * width,
                min_y + (row + 1) * height,
            )
            clipped = raw.intersection(clip_geometry)
            if not clipped.is_empty and clipped.area > 0:
                rows.append(
                    {
                        "grid_id": row * n_x + column,
                        "row": row,
                        "column": column,
                        "geometry": clipped,
                        "raw_geometry": raw,
                    }
                )
    return pd.DataFrame(rows)


def assign_grid_ids(
    trees: pd.DataFrame,
    bounds_geometry: object,
    n_x: int = 20,
    n_y: int = 20,
) -> np.ndarray:
    """Assign points to the regular grid without a spatial join."""
    min_x, min_y, max_x, max_y = bounds_geometry.bounds
    width = (max_x - min_x) / n_x
    height = (max_y - min_y) / n_y
    columns = np.floor((trees["x_m"].to_numpy() - min_x) / width).astype(int)
    rows = np.floor((trees["y_m"].to_numpy() - min_y) / height).astype(int)
    columns = np.clip(columns, 0, n_x - 1)
    rows = np.clip(rows, 0, n_y - 1)
    return rows * n_x + columns


def prepare_data() -> FigureData:
    """Load the canonical census and construct projected analysis geometries."""
    trees, traits = load_data()
    x_m, y_m = PROJECT_TO_UTM.transform(
        trees["X"].to_numpy(dtype=float), trees["Y"].to_numpy(dtype=float)
    )
    trees = trees.assign(x_m=x_m, y_m=y_m)

    blocks = load_all_blocks()
    perimeter_geometries, _ = read_shapes(PERIMETER_FILE, WGS84, UTM18N)
    sector_geometries, _ = read_shapes(SECTOR_FILE, WGS84, UTM18N)
    perimeter = unary_union(perimeter_geometries)
    sector_outline = unary_union(sector_geometries)
    block_outline = unary_union(blocks["geometry"].tolist())
    grid = build_grid(block_outline, perimeter)
    trees["grid_id"] = assign_grid_ids(trees, block_outline)

    missing = set(trees["grid_id"]) - set(grid["grid_id"])
    if missing:
        raise ValueError(f"Tree points assigned outside clipped grid: {sorted(missing)}")

    return FigureData(
        trees=trees,
        traits=traits,
        blocks=blocks,
        perimeter=perimeter,
        sectors=sector_geometries,
        sector_outline=sector_outline,
        grid=grid,
        bounds=perimeter.bounds,
    )


def geometry_patches(geometry: object) -> list:
    """Convert Polygon/MultiPolygon exteriors to Matplotlib patches."""
    if geometry.is_empty:
        return []
    if geometry.geom_type == "Polygon":
        polygons = [geometry]
    elif geometry.geom_type == "MultiPolygon":
        polygons = list(geometry.geoms)
    elif hasattr(geometry, "geoms"):
        polygons = [part for part in geometry.geoms if part.geom_type == "Polygon"]
    else:
        return []
    return [
        MplPolygon(np.asarray(polygon.exterior.coords), closed=True)
        for polygon in polygons
    ]


def add_geometry_collection(
    ax: plt.Axes,
    geometries: Iterable,
    *,
    facecolor="none",
    edgecolor: str = "0.35",
    linewidth: float = 0.45,
    alpha: float = 1.0,
    zorder: float = 1,
) -> PatchCollection:
    """Draw one or more shapely polygon geometries."""
    patches = []
    for geometry in geometries:
        patches.extend(geometry_patches(geometry))
    collection = PatchCollection(
        patches,
        facecolor=facecolor,
        edgecolor=edgecolor,
        linewidth=linewidth,
        alpha=alpha,
        zorder=zorder,
    )
    ax.add_collection(collection)
    return collection


def set_map_extent(ax: plt.Axes, geometry: object, padding: float = 0.025) -> None:
    min_x, min_y, max_x, max_y = geometry.bounds
    dx = max_x - min_x
    dy = max_y - min_y
    ax.set_xlim(min_x - padding * dx, max_x + padding * dx)
    ax.set_ylim(min_y - padding * dy, max_y + padding * dy)
    ax.set_aspect("equal")


def draw_choropleth(
    ax: plt.Axes,
    geometries: Sequence,
    values: Sequence[float],
    *,
    cmap,
    colorbar_label: str,
    boundary: object,
    vmin: Optional[float] = None,
    vmax: Optional[float] = None,
) -> None:
    """Draw a polygon choropleth with a compact color bar."""
    patches = []
    expanded = []
    for geometry, value in zip(geometries, values):
        parts = geometry_patches(geometry)
        patches.extend(parts)
        expanded.extend([float(value)] * len(parts))
    array = np.asarray(expanded)
    if vmin is None:
        vmin = float(np.nanmin(array))
    if vmax is None:
        vmax = float(np.nanmax(array))
    if np.isclose(vmin, vmax):
        vmax = vmin + 1.0
    collection = PatchCollection(
        patches,
        cmap=cmap,
        norm=Normalize(vmin=vmin, vmax=vmax),
        edgecolor="white",
        linewidth=0.25,
        zorder=1,
    )
    collection.set_array(array)
    ax.add_collection(collection)
    add_geometry_collection(
        ax,
        [boundary],
        facecolor="none",
        edgecolor="0.25",
        linewidth=0.65,
        zorder=4,
    )
    set_map_extent(ax, boundary, padding=0.015)
    ax.set_axis_off()
    colorbar = ax.figure.colorbar(collection, ax=ax, fraction=0.038, pad=0.012)
    colorbar.set_label(colorbar_label, fontsize=7.5)
    colorbar.ax.tick_params(labelsize=6.8, length=2)


def save_figure(figure: plt.Figure, stem: str, output: Path, dpi: int) -> None:
    """Save matching PDF and PNG figures atomically."""
    output.mkdir(parents=True, exist_ok=True)
    pdf_path = output / f"{stem}.pdf"
    png_path = output / f"{stem}.png"

    def complete_pdf(path: Path) -> bool:
        with path.open("rb") as handle:
            handle.seek(max(0, path.stat().st_size - 2048))
            return b"%%EOF" in handle.read()

    def complete_png(path: Path) -> bool:
        with path.open("rb") as handle:
            return handle.read(8) == b"\x89PNG\r\n\x1a\n"

    def copy_atomically(source: Path, destination: Path, validator) -> None:
        for attempt in range(1, 4):
            staging = output / f".{stem}.{attempt}.staging{destination.suffix}"
            staging.unlink(missing_ok=True)
            with source.open("rb") as source_handle, staging.open("wb") as target_handle:
                shutil.copyfileobj(source_handle, target_handle)
                target_handle.flush()
                os.fsync(target_handle.fileno())
            if validator(staging):
                os.replace(staging, destination)
                return
            staging.unlink(missing_ok=True)
        raise OSError(f"Could not copy a complete figure: {destination}")

    try:
        with tempfile.TemporaryDirectory(prefix=f"barrio-{stem}-") as temp_dir:
            temp_root = Path(temp_dir)
            for attempt in range(1, 4):
                temporary_pdf = temp_root / f"{stem}.{attempt}.pdf"
                figure.savefig(temporary_pdf, bbox_inches="tight")
                if complete_pdf(temporary_pdf):
                    break
            else:
                raise OSError(f"Could not write a complete PDF: {pdf_path}")

            temporary_png = temp_root / f"{stem}.png"
            figure.savefig(temporary_png, dpi=dpi, bbox_inches="tight")
            if not complete_png(temporary_png):
                raise OSError(f"Invalid PNG signature: {temporary_png}")
            copy_atomically(temporary_pdf, pdf_path, complete_pdf)
            copy_atomically(temporary_png, png_path, complete_png)
    finally:
        for attempt in range(1, 4):
            (output / f".{stem}.{attempt}.staging.pdf").unlink(missing_ok=True)
            (output / f".{stem}.{attempt}.staging.png").unlink(missing_ok=True)
        plt.close(figure)


def load_road_centerlines() -> pd.DataFrame:
    """Load the frozen road extract and project it to UTM zone 18N."""
    collection = json.loads(ROAD_CENTERLINES_FILE.read_text(encoding="utf-8"))
    rows = []
    for feature in collection["features"]:
        properties = dict(feature["properties"])
        geometry = shapely_transform(
            PROJECT_TO_UTM.transform, shape(feature["geometry"])
        )
        if not geometry.is_empty:
            rows.append({**properties, "geometry": geometry})
    roads = pd.DataFrame(rows)
    expected = {"Higher-order roads", "Local/secondary streets"}
    observed = set(roads["analysis_class"])
    if observed != expected:
        raise ValueError(f"Unexpected road classes: {sorted(observed)}")
    return roads


def road_unions(roads: pd.DataFrame) -> Tuple[object, object]:
    """Return merged higher-order and local/secondary centerlines."""
    higher = unary_union(
        roads.loc[
            roads["analysis_class"].eq("Higher-order roads"), "geometry"
        ].tolist()
    )
    local = unary_union(
        roads.loc[
            roads["analysis_class"].eq("Local/secondary streets"), "geometry"
        ].tolist()
    )
    if higher.is_empty or local.is_empty:
        raise ValueError("Both road classes require at least one centerline")
    return higher, local


def assign_tree_road_context(
    data: FigureData,
    trees: pd.DataFrame,
    higher_roads: object,
    local_roads: object,
    distance_m: float,
) -> pd.DataFrame:
    """Assign nearby records to the class nearest their block frontage."""
    result = trees.copy()
    points = [Point(x, y) for x, y in zip(result["x_m"], result["y_m"])]
    distance_higher = np.array([point.distance(higher_roads) for point in points])
    distance_local = np.array([point.distance(local_roads) for point in points])
    result["nearest_road_distance_m"] = np.minimum(distance_higher, distance_local)
    result["road_associated"] = result["nearest_road_distance_m"].le(distance_m)

    block_boundaries = {
        (str(row.Sector), int(row.Manzana)): row.geometry.boundary
        for row in data.blocks.itertuples(index=False)
    }
    road_classes = []
    for row, point in zip(result.itertuples(index=False), points):
        key = (str(row.Sector), int(row.Manzana))
        if key not in block_boundaries:
            raise KeyError(f"No cadastral geometry for tree block {key}")
        frontage_point = nearest_points(point, block_boundaries[key])[1]
        road_classes.append(
            "Higher-order roads"
            if frontage_point.distance(higher_roads)
            <= frontage_point.distance(local_roads)
            else "Local/secondary streets"
        )
    result["road_class"] = road_classes
    return result


def block_frontage_lengths(
    data: FigureData,
    higher_roads: object,
    local_roads: object,
    distance_m: float,
    sample_interval_m: float = 1.0,
) -> pd.DataFrame:
    """Estimate road-associated frontage by sampling block boundaries."""
    observed_blocks = {
        (str(sector), int(block))
        for sector, block in data.trees[["Sector", "Manzana"]]
        .drop_duplicates()
        .itertuples(index=False, name=None)
    }
    rows = []
    for block in data.blocks.itertuples(index=False):
        key = (str(block.Sector), int(block.Manzana))
        if key not in observed_blocks:
            continue
        lengths = {"Higher-order roads": 0.0, "Local/secondary streets": 0.0}
        boundary = block.geometry.boundary
        parts = list(boundary.geoms) if hasattr(boundary, "geoms") else [boundary]
        for line in parts:
            count = max(1, int(math.ceil(line.length / sample_interval_m)))
            represented_length = line.length / count
            for index in range(count):
                point = line.interpolate((index + 0.5) / count * line.length)
                distance_higher = point.distance(higher_roads)
                distance_local = point.distance(local_roads)
                if min(distance_higher, distance_local) > distance_m:
                    continue
                road_class = (
                    "Higher-order roads"
                    if distance_higher <= distance_local
                    else "Local/secondary streets"
                )
                lengths[road_class] += represented_length
        for road_class, frontage_m in lengths.items():
            if frontage_m > 0:
                rows.append(
                    {
                        "Sector": str(block.Sector),
                        "Manzana": int(block.Manzana),
                        "block_id": f"{block.Sector}-{int(block.Manzana)}",
                        "road_class": road_class,
                        "frontage_m": frontage_m,
                        "frontage_km": frontage_m / 1_000.0,
                    }
                )
    return pd.DataFrame(rows)


def draw_line_geometry(
    ax: plt.Axes,
    geometry: object,
    *,
    color: str,
    linewidth: float,
    alpha: float = 1.0,
    zorder: float = 2,
) -> None:
    """Draw a LineString or multipart line geometry."""
    if geometry.is_empty:
        return
    if geometry.geom_type == "LineString":
        lines = [geometry]
    elif hasattr(geometry, "geoms"):
        lines = [part for part in geometry.geoms if part.geom_type == "LineString"]
    else:
        lines = []
    for line in lines:
        coordinates = np.asarray(line.coords)
        ax.plot(
            coordinates[:, 0],
            coordinates[:, 1],
            color=color,
            linewidth=linewidth,
            alpha=alpha,
            zorder=zorder,
        )
