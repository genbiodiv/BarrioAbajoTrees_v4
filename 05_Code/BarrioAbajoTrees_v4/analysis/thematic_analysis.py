#!/usr/bin/env python3
"""Thematically organized analysis of the Barrio Abajo tree census.

The workflow starts only from ``ArbolesIndividuales_OK`` and the documented
spatial layers.  It preserves the recorded cadastral labels, applies the 25
approved block corrections to separate analysis fields, and then treats each
spatial framework independently:

* whole neighborhood;
* above and below Avenida Murillo;
* higher-order and local/secondary roads;
* residential and tourist-cultural reference areas;
* sectors, cadastral blocks, and a clipped 20 x 20 grid.

The output is primarily descriptive.  It calculates abundance, Hill diversity,
dominance, origin composition, Jaccard/Bray-Curtis beta diversity, basal area,
and standing above-ground carbon stock.  Targeted inference is limited to the
replicated block comparisons: a negative-binomial rate model for the Murillo
contrast and a paired, sector-clustered bootstrap for road-frontage density.
The old easting/northing regression suite is not fitted.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass, replace
import math
from pathlib import Path
from typing import Dict, Iterable, List, Mapping, Optional, Sequence, Tuple
import warnings

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from matplotlib.ticker import MaxNLocator
import numpy as np
import pandas as pd
from scipy.optimize import minimize
from scipy.special import gammaln
from scipy.stats import binomtest, gaussian_kde, norm
from shapely.geometry import LineString, Point
from shapely.ops import split, transform as shapely_transform, unary_union

import analysis_utils as base
import carbon_stock
from analysis_utils import expected_rarefied_richness, hill_numbers


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_OUTPUT = ROOT / "outputs"
CONTEXT_DIR = ROOT / "data" / "context"
COUNTRY_FILE = CONTEXT_DIR / "colombia_sinSA.shp"
DEPARTMENT_FILE = CONTEXT_DIR / "4326_atpio.shp"
CITY_FILE = CONTEXT_DIR / "quilla_mpio.shp"
SEED = 23_052_025
RAREFACTION_DRAWS = 4_000
BOOTSTRAP_ITERATIONS = 20_000

BLOCK_CORRECTION_DISTANCE_M = 25.0
ROAD_ASSOCIATION_DISTANCE_M = 25.0
CONTEXT_RADIUS_M = 150.0
MURILLO_NAME = "Avenida Murillo"

CONTEXT_ANCHORS = {
    "Residential": [("La Castellana Park", -74.7871, 10.9895)],
    "Tourist-cultural": [
        ("Old Customs Square", -74.7786, 10.9890),
        ("Caribbean Cultural Park", -74.7780, 10.9860),
    ],
}

ZONE_COLORS = {
    "Above Murillo": "#7b3294",
    "Below Murillo": "#008837",
    "Higher-order roads": "#d95f02",
    "Local/secondary streets": "#4c78a8",
    "Not assigned (>25 m)": "#9e9e9e",
    "Residential": "#984ea3",
    "Tourist-cultural": "#e66101",
    "Other/unclassified": "#bdbdbd",
    "Native": "#1b9e77",
    "Exotic": "#d95f02",
    "Unknown": "#9e9e9e",
}


@dataclass
class AnalysisBundle:
    spatial: base.FigureData
    trees: pd.DataFrame
    roads: pd.DataFrame
    murillo_roads: pd.DataFrame
    murillo_divider: object
    murillo_zones: Dict[str, object]
    context_buffers: Dict[str, object]
    context_anchors: pd.DataFrame
    sector_geometries: Dict[str, object]
    road_frontage: pd.DataFrame
    carbon_inventory: pd.DataFrame
    carbon_uncertainty: pd.DataFrame
    carbon_scenarios: pd.DataFrame
    height_diagnostics: pd.DataFrame
    height_information: Dict[str, float]


def iter_line_parts(geometry: object) -> Iterable[object]:
    if geometry.geom_type == "LineString":
        yield geometry
    elif hasattr(geometry, "geoms"):
        for part in geometry.geoms:
            if part.geom_type == "LineString":
                yield part


def fit_murillo_divider(
    data: base.FigureData, roads: pd.DataFrame
) -> Tuple[pd.DataFrame, np.ndarray, float, object, Dict[str, object]]:
    """Fit the approved curved Murillo divider and name its two sides."""
    murillo = roads.loc[roads["name"].fillna("").eq(MURILLO_NAME)].copy()
    coordinates: List[Tuple[float, float]] = []
    for geometry in murillo["geometry"]:
        for part in iter_line_parts(geometry):
            coordinates.extend(part.coords)
    values = np.asarray(coordinates, dtype=float)
    center_y = float(values[:, 1].mean())
    scaled_y = (values[:, 1] - center_y) / 1_000.0
    coefficients = np.polyfit(scaled_y, values[:, 0], 2)
    y_values = np.linspace(data.bounds[1] - 500, data.bounds[3] + 500, 500)
    x_values = np.polyval(coefficients, (y_values - center_y) / 1_000.0)
    divider = LineString(np.column_stack([x_values, y_values]))

    pieces = list(split(data.perimeter, divider).geoms)
    zones: Dict[str, List[object]] = {"Above Murillo": [], "Below Murillo": []}
    for piece in pieces:
        point = piece.representative_point()
        predicted_x = np.polyval(coefficients, (point.y - center_y) / 1_000.0)
        label = "Above Murillo" if point.x < predicted_x else "Below Murillo"
        zones[label].append(piece)
    merged = {label: unary_union(parts) for label, parts in zones.items()}
    return murillo, coefficients, center_y, divider, merged


def assign_murillo(
    trees: pd.DataFrame, coefficients: np.ndarray, center_y: float
) -> pd.DataFrame:
    result = trees.copy()
    scaled_y = (result["y_m"].to_numpy(dtype=float) - center_y) / 1_000.0
    predicted_x = np.polyval(coefficients, scaled_y)
    residual = result["x_m"].to_numpy(dtype=float) - predicted_x
    derivative = np.polyval(np.polyder(coefficients), scaled_y) / 1_000.0
    result["murillo_zone"] = np.where(
        residual < 0, "Above Murillo", "Below Murillo"
    )
    result["murillo_distance_m"] = np.abs(residual) / np.sqrt(1 + derivative**2)
    return result


def apply_block_corrections(
    trees: pd.DataFrame, blocks: pd.DataFrame
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Preserve source labels and correct only the 25 approved implausible labels."""
    result = trees.copy()
    result["recorded_sector"] = result["Sector"].astype(str)
    result["recorded_manzana"] = result["Manzana"].astype(int)
    result["recorded_block_id"] = (
        result["recorded_sector"]
        + "-"
        + result["recorded_manzana"].astype(str)
    )
    block_lookup = {
        (str(row.Sector), int(row.Manzana)): row.geometry
        for row in blocks.itertuples(index=False)
    }
    correction_rows: List[Dict[str, object]] = []
    analysis_sector: List[str] = []
    analysis_manzana: List[int] = []

    block_geometries = blocks["geometry"].tolist()
    for row in result.itertuples(index=False):
        point = Point(float(row.x_m), float(row.y_m))
        recorded_key = (str(row.recorded_sector), int(row.recorded_manzana))
        distances = np.asarray(
            [point.distance(geometry) for geometry in block_geometries], dtype=float
        )
        nearest_index = int(np.argmin(distances))
        nearest_row = blocks.iloc[nearest_index]
        nearest_key = (str(nearest_row["Sector"]), int(nearest_row["Manzana"]))
        recorded_distance = float(point.distance(block_lookup[recorded_key]))
        should_correct = (
            recorded_key != nearest_key
            and recorded_distance > BLOCK_CORRECTION_DISTANCE_M
        )
        chosen_key = nearest_key if should_correct else recorded_key
        analysis_sector.append(chosen_key[0])
        analysis_manzana.append(chosen_key[1])
        correction_rows.append(
            {
                "ID": row.ID,
                "recorded_sector": recorded_key[0],
                "recorded_manzana": recorded_key[1],
                "recorded_block_id": f"{recorded_key[0]}-{recorded_key[1]}",
                "nearest_sector": nearest_key[0],
                "nearest_manzana": nearest_key[1],
                "nearest_block_id": f"{nearest_key[0]}-{nearest_key[1]}",
                "nearest_block_distance_m": float(distances[nearest_index]),
                "recorded_block_distance_m": recorded_distance,
                "block_assignment_corrected": bool(should_correct),
                "analysis_sector": chosen_key[0],
                "analysis_manzana": chosen_key[1],
                "analysis_block_id": f"{chosen_key[0]}-{chosen_key[1]}",
            }
        )

    result["Sector"] = analysis_sector
    result["Manzana"] = np.asarray(analysis_manzana, dtype=int)
    result["analysis_block_id"] = (
        result["Sector"].astype(str) + "-" + result["Manzana"].astype(str)
    )
    audit = pd.DataFrame(correction_rows)
    result = result.merge(
        audit[
            [
                "ID",
                "nearest_sector",
                "nearest_manzana",
                "nearest_block_id",
                "nearest_block_distance_m",
                "recorded_block_distance_m",
                "block_assignment_corrected",
            ]
        ],
        on="ID",
        how="left",
        validate="one_to_one",
    )
    corrected_count = int(result["block_assignment_corrected"].sum())
    if corrected_count != 25:
        raise ValueError(
            f"Expected 25 approved block corrections, obtained {corrected_count}"
        )
    return result, audit


def load_sector_geometries() -> Dict[str, object]:
    geometries, records = base.read_shapes(base.SECTOR_FILE, base.WGS84, base.UTM18N)
    grouped: Dict[str, List[object]] = {}
    for geometry, record in zip(geometries, records):
        label = str(record["layer"]).replace("sector", "").upper()
        grouped.setdefault(label, []).append(geometry)
    return {
        label: unary_union(parts) for label, parts in sorted(grouped.items())
    }


def build_contexts(
    data: base.FigureData,
) -> Tuple[Dict[str, object], pd.DataFrame]:
    buffers: Dict[str, object] = {}
    anchors: List[Dict[str, object]] = []
    for context, sites in CONTEXT_ANCHORS.items():
        parts = []
        for site, longitude, latitude in sites:
            x_m, y_m = base.PROJECT_TO_UTM.transform(longitude, latitude)
            point = Point(x_m, y_m)
            parts.append(point.buffer(CONTEXT_RADIUS_M))
            anchors.append(
                {
                    "context": context,
                    "site": site,
                    "longitude": longitude,
                    "latitude": latitude,
                    "x_m": x_m,
                    "y_m": y_m,
                }
            )
        buffers[context] = unary_union(parts).intersection(data.perimeter)
    return buffers, pd.DataFrame(anchors)


def assign_contexts(
    trees: pd.DataFrame, buffers: Mapping[str, object]
) -> pd.DataFrame:
    result = trees.copy()
    points = [Point(x, y) for x, y in zip(result["x_m"], result["y_m"])]
    residential = np.asarray(
        [buffers["Residential"].covers(point) for point in points], dtype=bool
    )
    tourist = np.asarray(
        [buffers["Tourist-cultural"].covers(point) for point in points], dtype=bool
    )
    if np.any(residential & tourist):
        raise ValueError("Approved contextual buffers unexpectedly overlap")
    result["context_area"] = np.select(
        [residential, tourist],
        ["Residential", "Tourist-cultural"],
        default="Other/unclassified",
    )
    return result


def attach_carbon(
    trees: pd.DataFrame,
) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame, Dict[str, float], pd.DataFrame]:
    parameters, adjustments, height_information, height_diagnostics = (
        carbon_stock.fit_height_system(trees)
    )
    inventory = carbon_stock.calculate_inventory_estimates(
        trees, parameters, adjustments, height_information
    )
    uncertainty = carbon_stock.monte_carlo_uncertainty(
        inventory, height_information["cv_log_rse"]
    )
    scenarios = carbon_stock.scenario_table(inventory)
    carbon = inventory[
        [
            "ID",
            "height_model_m",
            "height_used_m",
            "height_source",
            "wood_density_g_cm3",
            "wood_density_sd_g_cm3",
            "estimation_method",
            "agb_kg",
            "agc_kg",
            "co2e_kg",
        ]
    ].copy()
    carbon["tree_form_carbon_kg"] = carbon["agc_kg"].where(
        carbon["estimation_method"].eq("Chave 2014 D-H-WD"), 0.0
    )
    result = trees.merge(carbon, on="ID", how="left", validate="one_to_one")
    result["tree_form_carbon_kg"] = result["tree_form_carbon_kg"].fillna(0.0)
    return (
        result,
        inventory,
        uncertainty,
        scenarios,
        height_information,
        height_diagnostics,
    )


def prepare_bundle() -> Tuple[AnalysisBundle, pd.DataFrame]:
    data = base.prepare_data()
    trees, correction_audit = apply_block_corrections(data.trees, data.blocks)

    roads = base.load_road_centerlines()
    murillo, coefficients, center_y, divider, zones = fit_murillo_divider(data, roads)
    trees = assign_murillo(trees, coefficients, center_y)

    buffers, anchors = build_contexts(data)
    trees = assign_contexts(trees, buffers)

    corrected_data = replace(data, trees=trees)
    higher_roads, local_roads = base.road_unions(roads)

    recorded_road_assignment = base.assign_tree_road_context(
        data,
        data.trees,
        higher_roads,
        local_roads,
        ROAD_ASSOCIATION_DISTANCE_M,
    )
    recorded_road_assignment["road_group_using_recorded_block"] = (
        recorded_road_assignment["road_class"].where(
            recorded_road_assignment["road_associated"], "Not assigned (>25 m)"
        )
    )
    trees = base.assign_tree_road_context(
        corrected_data,
        trees,
        higher_roads,
        local_roads,
        ROAD_ASSOCIATION_DISTANCE_M,
    )
    trees["road_group"] = trees["road_class"].where(
        trees["road_associated"], "Not assigned (>25 m)"
    )
    trees = trees.merge(
        recorded_road_assignment[["ID", "road_group_using_recorded_block"]],
        on="ID",
        how="left",
        validate="one_to_one",
    )
    trees["road_class_changed_by_block_correction"] = (
        trees["road_group"] != trees["road_group_using_recorded_block"]
    )
    corrected_data = replace(corrected_data, trees=trees)
    frontage = base.block_frontage_lengths(
        corrected_data,
        higher_roads,
        local_roads,
        ROAD_ASSOCIATION_DISTANCE_M,
    )

    (
        trees,
        carbon_inventory,
        carbon_uncertainty,
        carbon_scenarios,
        height_information,
        height_diagnostics,
    ) = attach_carbon(trees)
    corrected_data = replace(corrected_data, trees=trees)

    bundle = AnalysisBundle(
        spatial=corrected_data,
        trees=trees,
        roads=roads,
        murillo_roads=murillo,
        murillo_divider=divider,
        murillo_zones=zones,
        context_buffers=buffers,
        context_anchors=anchors,
        sector_geometries=load_sector_geometries(),
        road_frontage=frontage,
        carbon_inventory=carbon_inventory,
        carbon_uncertainty=carbon_uncertainty,
        carbon_scenarios=carbon_scenarios,
        height_diagnostics=height_diagnostics,
        height_information=height_information,
    )
    return bundle, correction_audit


def ordered_species_counts(values: pd.Series) -> pd.Series:
    """Return descending counts with alphabetical tie-breaking."""
    return (
        values.value_counts(sort=False)
        .sort_index(kind="mergesort")
        .sort_values(ascending=False, kind="mergesort")
    )


def summarize_group(
    group: pd.DataFrame,
    top_five: set,
    rarefaction_n: Optional[int],
) -> Dict[str, object]:
    n_records = len(group)
    if n_records == 0:
        return {
            "n_records": 0,
            "richness_q0": 0,
            "diversity_q1": np.nan,
            "diversity_q2": np.nan,
            "rarefaction_n": rarefaction_n,
            "rarefied_richness": np.nan,
            "top_five_fraction": np.nan,
            "dominant_species": "",
            "dominant_species_fraction": np.nan,
            "native_records": 0,
            "exotic_records": 0,
            "unknown_origin_records": 0,
            "native_fraction_known": np.nan,
            "basal_area_m2": 0.0,
            "tree_form_carbon_Mg": 0.0,
            "tree_form_co2e_Mg": 0.0,
        }
    q0, q1, q2 = hill_numbers(group["Species_name"])
    counts = ordered_species_counts(group["Species_name"])
    known = group["Origin"].isin(["Native", "Exotic"])
    native = int(group["Origin"].eq("Native").sum())
    exotic = int(group["Origin"].eq("Exotic").sum())
    rarefied = np.nan
    if rarefaction_n is not None and n_records >= rarefaction_n:
        rarefied = expected_rarefied_richness(counts, rarefaction_n)
    carbon_mg = float(group["tree_form_carbon_kg"].sum()) / 1_000.0
    return {
        "n_records": n_records,
        "richness_q0": q0,
        "diversity_q1": q1,
        "diversity_q2": q2,
        "rarefaction_n": rarefaction_n,
        "rarefied_richness": rarefied,
        "top_five_fraction": float(group["Species_name"].isin(top_five).mean()),
        "dominant_species": str(counts.index[0]),
        "dominant_species_fraction": float(counts.iloc[0] / n_records),
        "native_records": native,
        "exotic_records": exotic,
        "unknown_origin_records": int(n_records - native - exotic),
        "native_fraction_known": float(native / int(known.sum()))
        if known.any()
        else np.nan,
        "basal_area_m2": float(group["basal_area_m2"].fillna(0).sum()),
        "tree_form_carbon_Mg": carbon_mg,
        "tree_form_co2e_Mg": carbon_mg * 44.0 / 12.0,
    }


def framework_metrics(
    trees: pd.DataFrame,
    framework: str,
    group_column: Optional[str],
    group_order: Sequence[object],
    exposures: Optional[Mapping[object, Tuple[float, str]]] = None,
    rarefaction_n: Optional[int] = None,
) -> pd.DataFrame:
    top_five = set(ordered_species_counts(trees["Species_name"]).head(5).index)
    rows: List[Dict[str, object]] = []
    for group_name in group_order:
        group = trees if group_column is None else trees.loc[trees[group_column].eq(group_name)]
        values = summarize_group(group, top_five, rarefaction_n)
        exposure_value = np.nan
        exposure_unit = ""
        if exposures is not None and group_name in exposures:
            exposure_value, exposure_unit = exposures[group_name]
        row = {
            "framework": framework,
            "group": str(group_name),
            "exposure_value": exposure_value,
            "exposure_unit": exposure_unit,
            **values,
        }
        if np.isfinite(exposure_value) and exposure_value > 0:
            row["records_per_exposure"] = values["n_records"] / exposure_value
            row["carbon_Mg_per_exposure"] = (
                values["tree_form_carbon_Mg"] / exposure_value
            )
            row["basal_area_m2_per_exposure"] = (
                values["basal_area_m2"] / exposure_value
            )
        else:
            row["records_per_exposure"] = np.nan
            row["carbon_Mg_per_exposure"] = np.nan
            row["basal_area_m2_per_exposure"] = np.nan
        rows.append(row)
    return pd.DataFrame(rows)


def murillo_exposures(bundle: AnalysisBundle) -> Dict[str, Tuple[float, str]]:
    result: Dict[str, Tuple[float, str]] = {}
    boundaries = unary_union(
        [geometry.boundary for geometry in bundle.spatial.blocks["geometry"]]
    )
    for label, zone in bundle.murillo_zones.items():
        result[label] = (float(boundaries.intersection(zone).length) / 1_000.0, "km block frontage")
    return result


def sector_exposures(bundle: AnalysisBundle) -> Dict[str, Tuple[float, str]]:
    rows = bundle.spatial.blocks.copy()
    rows["perimeter_km"] = rows["geometry"].map(lambda geometry: geometry.length / 1_000.0)
    totals = rows.groupby("Sector")["perimeter_km"].sum()
    return {str(index): (float(value), "km block perimeter") for index, value in totals.items()}


def block_exposures(bundle: AnalysisBundle) -> Dict[str, Tuple[float, str]]:
    return {
        f"{row.Sector}-{int(row.Manzana)}": (float(row.geometry.length) / 1_000.0, "km block perimeter")
        for row in bundle.spatial.blocks.itertuples(index=False)
    }


def grid_exposures(bundle: AnalysisBundle) -> Dict[int, Tuple[float, str]]:
    return {
        int(row.grid_id): (float(row.geometry.area) / 10_000.0, "ha clipped cell")
        for row in bundle.spatial.grid.itertuples(index=False)
    }


def context_exposures(bundle: AnalysisBundle) -> Dict[str, Tuple[float, str]]:
    union = unary_union(list(bundle.context_buffers.values()))
    other = bundle.spatial.perimeter.difference(union)
    return {
        "Residential": (bundle.context_buffers["Residential"].area / 10_000.0, "ha"),
        "Tourist-cultural": (bundle.context_buffers["Tourist-cultural"].area / 10_000.0, "ha"),
        "Other/unclassified": (other.area / 10_000.0, "ha"),
    }


def road_exposures(bundle: AnalysisBundle) -> Dict[str, Tuple[float, str]]:
    totals = bundle.road_frontage.groupby("road_class")["frontage_km"].sum()
    return {
        str(index): (float(value), "km mapped frontage")
        for index, value in totals.items()
    }


def pairwise_beta(community: pd.DataFrame) -> Dict[str, float]:
    matrix = community.to_numpy(dtype=float)
    binary = matrix > 0
    rows: List[Tuple[float, float, float, float]] = []
    for first in range(len(matrix) - 1):
        for second in range(first + 1, len(matrix)):
            shared = int(np.sum(binary[first] & binary[second]))
            only_first = int(np.sum(binary[first] & ~binary[second]))
            only_second = int(np.sum(~binary[first] & binary[second]))
            denominator = shared + only_first + only_second
            if denominator == 0:
                continue
            jaccard = (only_first + only_second) / denominator
            minimum_unique = min(only_first, only_second)
            turnover_denominator = shared + 2 * minimum_unique
            turnover = 2 * minimum_unique / turnover_denominator if turnover_denominator else 0.0
            bray_denominator = float(np.sum(matrix[first] + matrix[second]))
            bray = float(np.sum(np.abs(matrix[first] - matrix[second])) / bray_denominator) if bray_denominator else np.nan
            rows.append((jaccard, turnover, jaccard - turnover, bray))
    values = np.asarray(rows, dtype=float)
    return {
        "units": int(len(matrix)),
        "pairwise_comparisons": int(len(rows)),
        "mean_jaccard_dissimilarity": float(np.nanmean(values[:, 0])),
        "mean_turnover": float(np.nanmean(values[:, 1])),
        "mean_nestedness": float(np.nanmean(values[:, 2])),
        "mean_bray_curtis_dissimilarity": float(np.nanmean(values[:, 3])),
    }


def beta_diversity_table(trees: pd.DataFrame) -> pd.DataFrame:
    definitions = [
        ("Murillo sides", "murillo_zone", ["Above Murillo", "Below Murillo"]),
        ("Road classes", "road_group", ["Higher-order roads", "Local/secondary streets"]),
        ("Residential vs tourist-cultural", "context_area", ["Residential", "Tourist-cultural"]),
        ("Sectors", "Sector", None),
        ("Blocks", "analysis_block_id", None),
        ("Occupied grid cells", "grid_id", None),
    ]
    rows = []
    for framework, column, allowed in definitions:
        subset = trees if allowed is None else trees.loc[trees[column].isin(allowed)]
        community = pd.crosstab(subset[column], subset["Species_name"])
        rows.append({"framework": framework, **pairwise_beta(community)})
    return pd.DataFrame(rows)


def road_frontage_unit_table(bundle: AnalysisBundle) -> pd.DataFrame:
    trees = bundle.trees.loc[bundle.trees["road_associated"]].copy()
    rows: List[Dict[str, object]] = []
    for unit in bundle.road_frontage.itertuples(index=False):
        group = trees.loc[
            trees["Sector"].eq(unit.Sector)
            & trees["Manzana"].astype(int).eq(int(unit.Manzana))
            & trees["road_class"].eq(unit.road_class)
        ]
        known = group["Origin"].isin(["Native", "Exotic"])
        native = int(group["Origin"].eq("Native").sum())
        carbon_mg = float(group["tree_form_carbon_kg"].sum()) / 1_000.0
        rows.append(
            {
                "block_id": unit.block_id,
                "Sector": unit.Sector,
                "Manzana": int(unit.Manzana),
                "road_class": unit.road_class,
                "frontage_km": float(unit.frontage_km),
                "n_records": len(group),
                "records_per_frontage_km": len(group) / float(unit.frontage_km),
                "richness_q0": int(group["Species_name"].nunique()),
                "diversity_q1": hill_numbers(group["Species_name"])[1] if len(group) else np.nan,
                "native_fraction_known": native / int(known.sum()) if known.any() else np.nan,
                "tree_form_carbon_Mg": carbon_mg,
                "tree_form_carbon_Mg_per_frontage_km": carbon_mg / float(unit.frontage_km),
            }
        )
    return pd.DataFrame(rows)


def road_paired_summary(units: pd.DataFrame) -> pd.DataFrame:
    paired_ids = (
        units.groupby("block_id")["road_class"].nunique().loc[lambda value: value.eq(2)].index
    )
    paired = units.loc[units["block_id"].isin(paired_ids)].copy()
    rows = []
    metrics = [
        "records_per_frontage_km",
        "diversity_q1",
        "native_fraction_known",
        "tree_form_carbon_Mg_per_frontage_km",
    ]
    for metric in metrics:
        pivot = paired.pivot(index="block_id", columns="road_class", values=metric).dropna()
        higher = pivot["Higher-order roads"]
        local = pivot["Local/secondary streets"]
        difference = higher - local
        rows.append(
            {
                "metric": metric,
                "paired_blocks_with_values": len(pivot),
                "median_higher_order": float(higher.median()),
                "median_local_secondary": float(local.median()),
                "median_paired_difference_higher_minus_local": float(difference.median()),
                "blocks_higher_greater": int((difference > 0).sum()),
                "blocks_local_greater": int((difference < 0).sum()),
                "ties": int((difference == 0).sum()),
            }
        )
    return pd.DataFrame(rows)


def murillo_block_table(bundle: AnalysisBundle) -> pd.DataFrame:
    """Return all cadastral blocks, including zeros, for block-level inference."""
    tree_groups = {
        str(block_id): group
        for block_id, group in bundle.trees.groupby("analysis_block_id")
    }
    rows: List[Dict[str, object]] = []
    for block in bundle.spatial.blocks.itertuples(index=False):
        block_id = f"{block.Sector}-{int(block.Manzana)}"
        overlap = {
            label: float(block.geometry.intersection(zone).area)
            for label, zone in bundle.murillo_zones.items()
        }
        zone = max(overlap, key=overlap.get)
        group = tree_groups.get(block_id, bundle.trees.iloc[0:0])
        native = int(group["Origin"].eq("Native").sum())
        exotic = int(group["Origin"].eq("Exotic").sum())
        unknown = int(group["Origin"].eq("Unknown").sum())
        perimeter_km = float(block.geometry.length) / 1_000.0
        rows.append(
            {
                "block_id": block_id,
                "Sector": str(block.Sector),
                "Manzana": int(block.Manzana),
                "murillo_zone": zone,
                "zone_area_fraction": overlap[zone] / float(block.geometry.area),
                "block_perimeter_km": perimeter_km,
                "n_records": int(len(group)),
                "records_per_block_perimeter_km": len(group) / perimeter_km,
                "native_records": native,
                "exotic_records": exotic,
                "unknown_origin_records": unknown,
                "native_fraction_known": (
                    native / (native + exotic) if native + exotic else np.nan
                ),
                "tree_form_carbon_Mg": float(group["tree_form_carbon_kg"].sum())
                / 1_000.0,
            }
        )
    return pd.DataFrame(rows).sort_values(["Sector", "Manzana"]).reset_index(drop=True)


def murillo_negative_binomial(blocks: pd.DataFrame) -> pd.DataFrame:
    """Fit one NB2 count model with log block perimeter as an offset.

    The model estimates the rate ratio for above versus below Murillo.  HC3
    sandwich covariance is calculated for the regression coefficients, with
    the maximum-likelihood dispersion estimate treated as a nuisance value.
    """
    y = blocks["n_records"].to_numpy(dtype=float)
    above = blocks["murillo_zone"].eq("Above Murillo").to_numpy(dtype=float)
    design = np.column_stack([np.ones(len(blocks)), above])
    offset = np.log(blocks["block_perimeter_km"].to_numpy(dtype=float))

    grouped = blocks.groupby("murillo_zone").agg(
        records=("n_records", "sum"), exposure=("block_perimeter_km", "sum")
    )
    below_rate = float(
        grouped.loc["Below Murillo", "records"]
        / grouped.loc["Below Murillo", "exposure"]
    )
    above_rate = float(
        grouped.loc["Above Murillo", "records"]
        / grouped.loc["Above Murillo", "exposure"]
    )
    beta_start = np.array([math.log(below_rate), math.log(above_rate / below_rate)])
    mu_start = np.exp(design @ beta_start + offset)
    alpha_start = max(
        float(np.sum(np.square(y - mu_start) - y) / np.sum(np.square(mu_start))),
        0.05,
    )

    def negative_log_likelihood(parameters: np.ndarray) -> float:
        beta = parameters[:2]
        alpha = math.exp(float(parameters[2]))
        size = 1.0 / alpha
        mu = np.exp(np.clip(design @ beta + offset, -30.0, 30.0))
        values = (
            gammaln(y + size)
            - gammaln(size)
            - gammaln(y + 1.0)
            + size * (math.log(size) - np.log(size + mu))
            + y * (np.log(mu) - np.log(size + mu))
        )
        return -float(np.sum(values))

    fit = minimize(
        negative_log_likelihood,
        np.r_[beta_start, math.log(alpha_start)],
        method="L-BFGS-B",
        bounds=[(-20.0, 20.0), (-20.0, 20.0), (-12.0, 6.0)],
    )
    if not fit.success:
        raise RuntimeError(f"Murillo negative-binomial model failed: {fit.message}")

    beta = fit.x[:2]
    alpha = math.exp(float(fit.x[2]))
    mu = np.exp(design @ beta + offset)
    weights = mu / (1.0 + alpha * mu)
    bread = np.linalg.pinv(design.T @ (weights[:, None] * design))
    leverage = weights * np.sum((design @ bread) * design, axis=1)
    score_scalar = (y - mu) / (1.0 + alpha * mu)
    adjusted_score = score_scalar / np.clip(1.0 - leverage, 0.05, None)
    meat = design.T @ (np.square(adjusted_score)[:, None] * design)
    covariance = bread @ meat @ bread
    standard_error = float(math.sqrt(max(covariance[1, 1], 0.0)))
    z_value = float(beta[1] / standard_error)
    critical = float(norm.ppf(0.975))
    lower = float(math.exp(beta[1] - critical * standard_error))
    upper = float(math.exp(beta[1] + critical * standard_error))

    return pd.DataFrame(
        [
            {
                "model": "Negative binomial count model with log(block perimeter km) offset",
                "comparison": "Above Murillo vs Below Murillo",
                "blocks": int(len(blocks)),
                "blocks_above": int((above == 1).sum()),
                "blocks_below": int((above == 0).sum()),
                "coefficient_log_rate_ratio": float(beta[1]),
                "HC3_standard_error": standard_error,
                "rate_ratio": float(math.exp(beta[1])),
                "rate_ratio_95pct_lower": lower,
                "rate_ratio_95pct_upper": upper,
                "wald_z": z_value,
                "wald_p": float(2.0 * norm.sf(abs(z_value))),
                "dispersion_alpha": alpha,
                "log_likelihood": -float(fit.fun),
                "converged": bool(fit.success),
            }
        ]
    )


def road_density_paired_inference(units: pd.DataFrame) -> pd.DataFrame:
    """Estimate the paired road-density effect with sector-clustered resampling."""
    pivot = units.pivot(
        index=["Sector", "block_id"],
        columns="road_class",
        values="records_per_frontage_km",
    ).dropna()
    difference = (
        pivot["Higher-order roads"] - pivot["Local/secondary streets"]
    ).rename("difference").reset_index()
    sector_values = {
        sector: group["difference"].to_numpy(dtype=float)
        for sector, group in difference.groupby("Sector")
    }
    sectors = np.asarray(sorted(sector_values), dtype=object)
    rng = np.random.default_rng(SEED + 101)
    bootstrap = np.empty(BOOTSTRAP_ITERATIONS)
    for index in range(BOOTSTRAP_ITERATIONS):
        selected = rng.choice(sectors, size=len(sectors), replace=True)
        sample = np.concatenate([sector_values[str(sector)] for sector in selected])
        bootstrap[index] = float(np.median(sample))

    values = difference["difference"].to_numpy(dtype=float)
    positive = int((values > 0).sum())
    negative = int((values < 0).sum())
    test = binomtest(positive, positive + negative, p=0.5, alternative="two-sided")
    proportion_interval = test.proportion_ci(confidence_level=0.95, method="exact")
    return pd.DataFrame(
        [
            {
                "comparison": "Higher-order minus local/secondary frontage",
                "paired_blocks": int(len(values)),
                "sectors_resampled": int(len(sectors)),
                "median_paired_difference_records_per_km": float(np.median(values)),
                "cluster_bootstrap_iterations": BOOTSTRAP_ITERATIONS,
                "cluster_bootstrap_95pct_lower": float(np.quantile(bootstrap, 0.025)),
                "cluster_bootstrap_95pct_upper": float(np.quantile(bootstrap, 0.975)),
                "blocks_higher_greater": positive,
                "blocks_local_greater": negative,
                "ties": int((values == 0).sum()),
                "proportion_higher_greater": positive / len(values),
                "exact_proportion_95pct_lower": float(proportion_interval.low),
                "exact_proportion_95pct_upper": float(proportion_interval.high),
                "exact_sign_test_p": float(test.pvalue),
            }
        ]
    )


def rarefied_sample_coverage(counts: np.ndarray, sample_size: int) -> float:
    """Estimate coverage for an individual-based rarefied sample (m < n)."""
    counts = np.asarray(counts, dtype=int)
    total = int(counts.sum())
    if sample_size <= 0 or sample_size >= total:
        raise ValueError("Coverage-standardized rarefaction requires 0 < m < n")
    log_denominator = float(
        gammaln(total) - gammaln(sample_size + 1) - gammaln(total - sample_size)
    )
    missing_probability = 0.0
    for abundance in counts:
        remaining = total - int(abundance)
        if remaining < sample_size:
            continue
        log_numerator = float(
            gammaln(remaining + 1)
            - gammaln(sample_size + 1)
            - gammaln(remaining - sample_size + 1)
        )
        missing_probability += (abundance / total) * math.exp(
            log_numerator - log_denominator
        )
    return float(1.0 - missing_probability)


def coverage_standardized_diversity(bundle: AnalysisBundle) -> pd.DataFrame:
    """Rarefy paired assemblages to common empirical sample coverage."""
    definitions = [
        (
            "Murillo sides",
            bundle.trees,
            "murillo_zone",
            ["Above Murillo", "Below Murillo"],
        ),
        (
            "Road classes",
            bundle.trees.loc[bundle.trees["road_associated"]],
            "road_class",
            ["Higher-order roads", "Local/secondary streets"],
        ),
        (
            "Context areas",
            bundle.trees,
            "context_area",
            ["Residential", "Tourist-cultural"],
        ),
    ]
    rng = np.random.default_rng(SEED + 202)
    rows: List[Dict[str, object]] = []
    for framework, data, column, groups in definitions:
        abundance = {
            group: data.loc[data[column].eq(group), "Species_name"]
            .value_counts()
            .to_numpy(dtype=int)
            for group in groups
        }
        full_coverage = {
            group: rarefied_sample_coverage(values, int(values.sum()) - 1)
            for group, values in abundance.items()
        }
        target = float(min(full_coverage.values()))
        for group in groups:
            counts = abundance[group]
            total = int(counts.sum())
            candidates = np.arange(1, total, dtype=int)
            coverage_values = np.asarray(
                [rarefied_sample_coverage(counts, int(size)) for size in candidates]
            )
            sample_size = int(candidates[np.argmin(np.abs(coverage_values - target))])
            achieved = float(rarefied_sample_coverage(counts, sample_size))
            samples = rng.multivariate_hypergeometric(
                counts, sample_size, size=RAREFACTION_DRAWS
            )
            q0 = np.sum(samples > 0, axis=1).astype(float)
            proportions = samples / sample_size
            with np.errstate(divide="ignore", invalid="ignore"):
                entropy_terms = np.where(
                    proportions > 0, proportions * np.log(proportions), 0.0
                )
            q1 = np.exp(-np.sum(entropy_terms, axis=1))
            q2 = 1.0 / np.sum(np.square(proportions), axis=1)
            for order, values in ((0, q0), (1, q1), (2, q2)):
                rows.append(
                    {
                        "framework": framework,
                        "group": group,
                        "n_records": total,
                        "sample_coverage_at_n_minus_1": full_coverage[group],
                        "common_target_coverage": target,
                        "standardized_sample_size": sample_size,
                        "achieved_coverage": achieved,
                        "hill_order_q": order,
                        "effective_species_mean": float(np.mean(values)),
                        "rarefaction_2_5pct": float(np.quantile(values, 0.025)),
                        "rarefaction_97_5pct": float(np.quantile(values, 0.975)),
                        "rarefaction_draws": RAREFACTION_DRAWS,
                    }
                )
    return pd.DataFrame(rows)


def road_threshold_sensitivity(bundle: AnalysisBundle) -> pd.DataFrame:
    """Repeat the road comparison at 15, 25, and 40 m thresholds."""
    higher_roads, local_roads = base.road_unions(bundle.roads)
    rows: List[Dict[str, object]] = []
    for threshold in (15.0, 25.0, 40.0):
        trees = base.assign_tree_road_context(
            bundle.spatial,
            bundle.trees,
            higher_roads,
            local_roads,
            threshold,
        )
        spatial = replace(bundle.spatial, trees=trees)
        frontage = base.block_frontage_lengths(
            spatial, higher_roads, local_roads, threshold
        )
        temporary = replace(bundle, spatial=spatial, trees=trees, road_frontage=frontage)
        units = road_frontage_unit_table(temporary)
        paired = road_paired_summary(units).set_index("metric").loc[
            "records_per_frontage_km"
        ]
        associated = trees.loc[trees["road_associated"]]
        counts = associated["road_class"].value_counts()
        exposure = frontage.groupby("road_class")["frontage_km"].sum()
        higher_count = int(counts.get("Higher-order roads", 0))
        local_count = int(counts.get("Local/secondary streets", 0))
        higher_frontage = float(exposure.get("Higher-order roads", np.nan))
        local_frontage = float(exposure.get("Local/secondary streets", np.nan))
        higher_density = higher_count / higher_frontage
        local_density = local_count / local_frontage
        rows.append(
            {
                "threshold_m": threshold,
                "higher_order_records": higher_count,
                "local_secondary_records": local_count,
                "unassigned_records": int(len(trees) - len(associated)),
                "higher_order_frontage_km": higher_frontage,
                "local_secondary_frontage_km": local_frontage,
                "higher_order_records_per_km": higher_density,
                "local_secondary_records_per_km": local_density,
                "aggregate_density_ratio_higher_vs_local": higher_density
                / local_density,
                "paired_blocks": int(paired["paired_blocks_with_values"]),
                "median_paired_difference_records_per_km": float(
                    paired["median_paired_difference_higher_minus_local"]
                ),
                "blocks_higher_greater": int(paired["blocks_higher_greater"]),
                "blocks_local_greater": int(paired["blocks_local_greater"]),
            }
        )
    return pd.DataFrame(rows)


def grid_resolution_sensitivity(bundle: AnalysisBundle) -> pd.DataFrame:
    """Repeat the grid beta-diversity summary at coarser and finer grains."""
    rows: List[Dict[str, object]] = []
    block_outline = unary_union(bundle.spatial.blocks["geometry"].tolist())
    min_x, min_y, max_x, max_y = block_outline.bounds
    for dimension in (15, 20, 25):
        grid = base.build_grid(
            block_outline,
            bundle.spatial.perimeter,
            n_x=dimension,
            n_y=dimension,
        )
        grid_id = base.assign_grid_ids(
            bundle.trees,
            block_outline,
            n_x=dimension,
            n_y=dimension,
        )
        temporary = bundle.trees.assign(sensitivity_grid_id=grid_id)
        community = pd.crosstab(
            temporary["sensitivity_grid_id"], temporary["Species_name"]
        )
        beta = pairwise_beta(community)
        abundance = temporary.groupby("sensitivity_grid_id").size()
        richness = temporary.groupby("sensitivity_grid_id")["Species_name"].nunique()
        rows.append(
            {
                "grid_columns": dimension,
                "grid_rows": dimension,
                "nominal_cell_width_m": (max_x - min_x) / dimension,
                "nominal_cell_height_m": (max_y - min_y) / dimension,
                "clipped_cells": int(len(grid)),
                "occupied_cells": int(len(community)),
                "median_records_per_occupied_cell": float(abundance.median()),
                "median_richness_per_occupied_cell": float(richness.median()),
                **beta,
            }
        )
    return pd.DataFrame(rows)


def build_metric_tables(bundle: AnalysisBundle) -> Tuple[pd.DataFrame, Dict[str, pd.DataFrame]]:
    trees = bundle.trees
    tables: Dict[str, pd.DataFrame] = {}
    tables["whole"] = framework_metrics(
        trees, "Whole neighborhood", None, ["Barrio Abajo"], rarefaction_n=None
    )
    tables["murillo"] = framework_metrics(
        trees,
        "Murillo sides",
        "murillo_zone",
        ["Above Murillo", "Below Murillo"],
        murillo_exposures(bundle),
        int(trees.groupby("murillo_zone").size().min()),
    )
    road_trees = trees.loc[trees["road_associated"]].copy()
    tables["roads"] = framework_metrics(
        road_trees,
        "Road classes",
        "road_class",
        ["Higher-order roads", "Local/secondary streets"],
        road_exposures(bundle),
        int(road_trees.groupby("road_class").size().min()),
    )
    tables["contexts"] = framework_metrics(
        trees,
        "Context areas",
        "context_area",
        ["Residential", "Tourist-cultural", "Other/unclassified"],
        context_exposures(bundle),
        int(trees.groupby("context_area").size().min()),
    )
    sector_order = sorted(bundle.sector_geometries)
    sector_minimum = int(trees.groupby("Sector").size().min())
    tables["sectors"] = framework_metrics(
        trees,
        "Sectors",
        "Sector",
        sector_order,
        sector_exposures(bundle),
        sector_minimum,
    )
    block_order = [
        f"{row.Sector}-{int(row.Manzana)}"
        for row in bundle.spatial.blocks.sort_values(["Sector", "Manzana"]).itertuples(index=False)
    ]
    tables["blocks"] = framework_metrics(
        trees,
        "Cadastral blocks",
        "analysis_block_id",
        block_order,
        block_exposures(bundle),
        10,
    )
    grid_order = sorted(bundle.spatial.grid["grid_id"].astype(int).tolist())
    tables["grid"] = framework_metrics(
        trees,
        "20 x 20 grid",
        "grid_id",
        grid_order,
        grid_exposures(bundle),
        10,
    )
    combined = pd.concat(tables.values(), ignore_index=True)
    return combined, tables


def draw_base(ax: plt.Axes, bundle: AnalysisBundle, blocks: bool = True) -> None:
    base.add_geometry_collection(
        ax,
        [bundle.spatial.perimeter],
        facecolor="white",
        edgecolor="0.2",
        linewidth=0.8,
        zorder=0,
    )
    if blocks:
        base.add_geometry_collection(
            ax,
            bundle.spatial.blocks["geometry"],
            facecolor="none",
            edgecolor="0.82",
            linewidth=0.3,
            zorder=1,
        )
    base.set_map_extent(ax, bundle.spatial.perimeter, padding=0.02)
    ax.set_axis_off()


def plot_line_geometries(
    ax: plt.Axes,
    geometries: Iterable[object],
    color: str,
    linewidth: float,
    alpha: float = 1.0,
    zorder: float = 4,
) -> None:
    for geometry in geometries:
        base.draw_line_geometry(
            ax, geometry, color=color, linewidth=linewidth, alpha=alpha, zorder=zorder
        )


def add_lonlat_scale_bar(ax: plt.Axes, kilometres: float) -> None:
    """Draw an approximate scale bar on a small geographic inset."""
    x_min, x_max = ax.get_xlim()
    y_min, y_max = ax.get_ylim()
    latitude = 0.5 * (y_min + y_max)
    length_degrees = kilometres / (111.32 * np.cos(np.deg2rad(latitude)))
    start_x = x_min + 0.08 * (x_max - x_min)
    start_y = y_min + 0.06 * (y_max - y_min)
    ax.plot(
        [start_x, start_x + length_degrees],
        [start_y, start_y],
        color="black",
        linewidth=1.0,
        zorder=20,
    )
    ax.text(
        start_x + length_degrees / 2,
        start_y + 0.025 * (y_max - y_min),
        f"{kilometres:g} km",
        ha="center",
        va="bottom",
        fontsize=6.5,
        zorder=20,
    )


def plot_context_map(
    ax: plt.Axes,
    base_geometries: Sequence[object],
    highlighted: Sequence[object],
    title: str,
    scale_km: float,
) -> None:
    """Draw one geographic context inset for Figure 1."""
    extent = unary_union(base_geometries)
    base.add_geometry_collection(
        ax,
        base_geometries,
        facecolor="#e6e6e6",
        edgecolor="white",
        linewidth=0.25,
    )
    base.add_geometry_collection(
        ax,
        highlighted,
        facecolor="#b22222",
        edgecolor="#8b1a1a",
        linewidth=0.4,
        zorder=3,
    )
    base.set_map_extent(ax, extent, padding=0.02)
    ax.set_title(title, loc="left", fontweight="bold")
    ax.set_axis_off()
    add_lonlat_scale_bar(ax, scale_km)


def add_metric_scale_bar(ax: plt.Axes, length_m: float = 200) -> None:
    """Draw a metric scale bar on a projected neighborhood map."""
    x_min, x_max = ax.get_xlim()
    y_min, y_max = ax.get_ylim()
    start_x = x_min + 0.08 * (x_max - x_min)
    start_y = y_min + 0.06 * (y_max - y_min)
    ax.plot(
        [start_x, start_x + length_m],
        [start_y, start_y],
        color="black",
        linewidth=1.5,
        solid_capstyle="butt",
        zorder=30,
    )
    ax.text(
        start_x + length_m / 2,
        start_y + 0.018 * (y_max - y_min),
        f"{int(length_m)} m",
        ha="center",
        va="bottom",
        fontsize=7,
        zorder=30,
    )


def add_north_arrow(ax: plt.Axes, x: float = 0.08, y: float = 0.17) -> None:
    ax.annotate(
        "N",
        xy=(x, y + 0.08),
        xytext=(x, y),
        xycoords="axes fraction",
        textcoords="axes fraction",
        ha="center",
        va="center",
        fontsize=8,
        fontweight="bold",
        arrowprops={"arrowstyle": "-|>", "color": "black", "lw": 1.0},
        zorder=30,
    )


def figure_neighborhood_overview(
    bundle: AnalysisBundle, output: Path, dpi: int
) -> None:
    """Show geographic context, growth forms, and the complete census."""
    country, country_records = base.read_shapes(
        COUNTRY_FILE, base.WGS84, base.WGS84
    )
    municipalities, _ = base.read_shapes(
        DEPARTMENT_FILE, base.WGS84, base.WGS84
    )
    city, _ = base.read_shapes(CITY_FILE, 3857, base.WGS84)
    department = [
        geometry
        for geometry, record in zip(country, country_records)
        if str(record.get("DeCodigo", "")) == "08"
    ]
    if not department:
        department = [unary_union(municipalities)]
    barrio_wgs84 = [
        shapely_transform(base.PROJECT_TO_WGS84.transform, bundle.spatial.perimeter)
    ]

    figure = plt.figure(figsize=(11.0, 8.3), constrained_layout=True)
    layout = figure.add_gridspec(
        2,
        4,
        height_ratios=[0.34, 0.66],
        width_ratios=[1, 1, 1, 1.15],
    )
    top_axes = [figure.add_subplot(layout[0, index]) for index in range(4)]
    map_ax = figure.add_subplot(layout[1, :])

    plot_context_map(top_axes[0], country, department, "A  Colombia", 200)
    plot_context_map(top_axes[1], municipalities, city, "B  Atlántico", 25)
    plot_context_map(top_axes[2], city, barrio_wgs84, "C  Barranquilla", 5)

    form_order = ["Tree", "Palm", "Shrub", "Grass"]
    form_labels = ["Tree", "Palm", "Shrub", "Arborescent\nherb"]
    form_colors = {
        "Tree": "#087f5b",
        "Palm": "#d9a21b",
        "Shrub": "#7b3294",
        "Grass": "#b2182b",
    }
    form_counts = bundle.trees["Form"].value_counts().reindex(form_order, fill_value=0)
    bars = top_axes[3].bar(
        form_labels,
        form_counts.to_numpy(),
        color=[form_colors[value] for value in form_order],
    )
    top_axes[3].set_ylabel("Records")
    top_axes[3].set_title("D  Growth forms", loc="left", fontweight="bold")
    top_axes[3].spines[["top", "right"]].set_visible(False)
    top_axes[3].grid(axis="y", color="0.9", linewidth=0.5, zorder=0)
    top_axes[3].tick_params(axis="x", labelsize=7)
    for bar, count in zip(bars, form_counts):
        top_axes[3].text(
            bar.get_x() + bar.get_width() / 2,
            count + 18,
            f"{int(count):,}",
            ha="center",
            va="bottom",
            fontsize=7,
        )

    draw_base(map_ax, bundle)
    higher_roads = bundle.roads.loc[
        bundle.roads["analysis_class"].eq("Higher-order roads"), "geometry"
    ]
    plot_line_geometries(map_ax, higher_roads, "0.65", 1.0, 0.8, 2)
    for form in form_order:
        subset = bundle.trees.loc[bundle.trees["Form"].eq(form)]
        map_ax.scatter(
            subset["x_m"],
            subset["y_m"],
            s=5 if form == "Tree" else 8,
            color=form_colors[form],
            alpha=0.70,
            linewidths=0,
            label="Arborescent herb" if form == "Grass" else form,
            zorder=4,
        )
    road_labels = [
        ("Calle 53", -74.7896, 10.9900, 80),
        ("Carrera 54", -74.7845, 10.9945, 20),
        ("Vía 40", -74.7795, 10.9915, -80),
        ("Avenida Murillo", -74.7843, 10.9883, 75),
        ("Carrera 50", -74.7830, 10.98905, -11),
        ("Carrera 46", -74.7840, 10.9867, -11),
    ]
    for name, longitude, latitude, angle in road_labels:
        x_value, y_value = base.PROJECT_TO_UTM.transform(longitude, latitude)
        map_ax.text(
            x_value,
            y_value,
            name,
            rotation=angle,
            fontsize=7,
            fontweight="bold",
            ha="center",
            va="center",
            color="0.15",
            zorder=8,
        )
    map_ax.set_title(
        "E  Publicly accessible trees and arborescent plants (n = 1,256)",
        loc="left",
        fontweight="bold",
    )
    map_ax.legend(frameon=False, loc="lower right", ncol=4, markerscale=1.4)
    add_metric_scale_bar(map_ax, 200)
    add_north_arrow(map_ax)
    figure.suptitle("Study area and complete inventory", fontweight="bold")
    base.save_figure(figure, "Figure1_neighborhood_inventory", output, dpi)


def figure_top_species(
    bundle: AnalysisBundle, output: Path, dpi: int
) -> None:
    """Map every record of the ten most abundant species."""
    top_species = ordered_species_counts(bundle.trees["Species_name"]).head(10)
    figure, axes = plt.subplots(
        2,
        5,
        figsize=(13.0, 5.8),
        constrained_layout=True,
        sharex=True,
        sharey=True,
    )
    for panel_index, (ax, (species_name, count)) in enumerate(
        zip(axes.flat, top_species.items()), start=1
    ):
        draw_base(ax, bundle)
        ax.scatter(
            bundle.trees["x_m"],
            bundle.trees["y_m"],
            s=1.3,
            color="0.82",
            alpha=0.35,
            linewidths=0,
            zorder=2,
        )
        subset = bundle.trees.loc[bundle.trees["Species_name"].eq(species_name)]
        ax.scatter(
            subset["x_m"],
            subset["y_m"],
            s=8,
            color="#087f5b",
            alpha=0.82,
            linewidths=0,
            zorder=4,
        )
        panel = chr(64 + panel_index)
        ax.set_title(
            f"{panel}  {species_name}",
            loc="left",
            fontstyle="italic",
            fontsize=8.5,
        )
        ax.text(
            0.99,
            0.99,
            f"n = {int(count)}",
            transform=ax.transAxes,
            ha="right",
            va="top",
            fontsize=7.3,
            bbox={"facecolor": "white", "edgecolor": "none", "pad": 0.8, "alpha": 0.8},
            zorder=10,
        )
    figure.suptitle(
        "Spatial distribution of the ten most abundant species",
        fontweight="bold",
    )
    base.save_figure(figure, "Figure2_top_species_distribution", output, dpi)


def figure_distribution_murillo(
    bundle: AnalysisBundle, metrics: pd.DataFrame, output: Path, dpi: int
) -> None:
    rows = metrics.set_index("group").loc[["Above Murillo", "Below Murillo"]]
    colors = [ZONE_COLORS[label] for label in rows.index]
    figure = plt.figure(figsize=(11, 4.7), constrained_layout=True)
    layout = figure.add_gridspec(1, 3, width_ratios=[1.8, 0.85, 0.95])
    axes = [figure.add_subplot(layout[0, index]) for index in range(3)]

    draw_base(axes[0], bundle)
    for label, geometry in bundle.murillo_zones.items():
        base.add_geometry_collection(
            axes[0], [geometry], facecolor=ZONE_COLORS[label], edgecolor="none", alpha=0.15, zorder=1.5
        )
    for label in rows.index:
        subset = bundle.trees.loc[bundle.trees["murillo_zone"].eq(label)]
        axes[0].scatter(subset["x_m"], subset["y_m"], s=5, color=ZONE_COLORS[label], alpha=0.66, linewidths=0, zorder=4)
    plot_line_geometries(axes[0], bundle.murillo_roads["geometry"], "#e66101", 2.0, 0.95, 6)
    base.draw_line_geometry(axes[0], bundle.murillo_divider, color="black", linewidth=0.9, alpha=0.9, zorder=7)
    axes[0].set_title("A  Avenida Murillo and census records", loc="left", fontweight="bold")
    axes[0].legend(
        handles=[Patch(facecolor=ZONE_COLORS[label], alpha=0.6, label=label) for label in rows.index],
        frameon=False,
        loc="lower left",
    )

    axes[1].bar(rows.index, rows["n_records"], color=colors)
    axes[1].set_ylabel("Tree records")
    axes[1].set_title("B  Recorded abundance", loc="left", fontweight="bold")
    axes[1].tick_params(axis="x", rotation=25)
    for index, value in enumerate(rows["n_records"]):
        axes[1].text(index, value + 12, f"{int(value):,}", ha="center", fontsize=8)

    axes[2].bar(rows.index, rows["records_per_exposure"], color=colors)
    axes[2].set_ylabel("Records per km of block frontage")
    axes[2].set_title("C  Frontage-standardized abundance", loc="left", fontweight="bold")
    axes[2].tick_params(axis="x", rotation=25)
    for index, value in enumerate(rows["records_per_exposure"]):
        axes[2].text(index, value * 1.02, f"{value:.1f}", ha="center", fontsize=8)
    for ax in axes[1:]:
        ax.spines[["top", "right"]].set_visible(False)
        ax.grid(axis="y", color="0.9", linewidth=0.5, zorder=0)
    figure.suptitle("Tree abundance on the two sides of Avenida Murillo", fontweight="bold")
    base.save_figure(figure, "Figure4_Murillo_comparison", output, dpi)


def geometry_metric_table(
    geometries: pd.DataFrame,
    metric: pd.DataFrame,
    left_keys: Sequence[str],
    right_key: str,
) -> pd.DataFrame:
    left = geometries.copy()
    if len(left_keys) == 1:
        left["join_id"] = left[left_keys[0]].astype(str)
    else:
        left["join_id"] = left[left_keys[0]].astype(str) + "-" + left[left_keys[1]].astype(int).astype(str)
    right = metric.copy()
    right["join_id"] = right[right_key].astype(str)
    return left.merge(right.drop(columns=[right_key]), on="join_id", how="left", validate="one_to_one")


def figure_multiscale(
    bundle: AnalysisBundle,
    tables: Dict[str, pd.DataFrame],
    beta: pd.DataFrame,
    output: Path,
    dpi: int,
) -> None:
    figure, axes = plt.subplots(2, 2, figsize=(10.5, 8.5), constrained_layout=True)

    sector_frame = pd.DataFrame(
        {"sector": list(bundle.sector_geometries), "geometry": list(bundle.sector_geometries.values())}
    ).merge(
        tables["sectors"][["group", "richness_q0"]],
        left_on="sector",
        right_on="group",
        how="left",
        validate="one_to_one",
    )
    base.draw_choropleth(
        axes[0, 0], sector_frame["geometry"], sector_frame["richness_q0"], cmap="YlGn", colorbar_label="Observed species", boundary=bundle.spatial.perimeter, vmin=0
    )
    axes[0, 0].set_title("A  Sectors", loc="left", fontweight="bold")

    block_frame = geometry_metric_table(
        bundle.spatial.blocks,
        tables["blocks"],
        ["Sector", "Manzana"],
        "group",
    )
    base.draw_choropleth(
        axes[0, 1], block_frame["geometry"], block_frame["richness_q0"].fillna(0), cmap="YlGn", colorbar_label="Observed species", boundary=bundle.spatial.perimeter, vmin=0
    )
    axes[0, 1].set_title("B  Cadastral blocks", loc="left", fontweight="bold")

    grid_frame = bundle.spatial.grid.merge(
        tables["grid"][["group", "richness_q0"]].assign(group=lambda frame: frame["group"].astype(int)),
        left_on="grid_id",
        right_on="group",
        how="left",
        validate="one_to_one",
    )
    base.draw_choropleth(
        axes[1, 0], grid_frame["geometry"], grid_frame["richness_q0"].fillna(0), cmap="YlGn", colorbar_label="Observed species", boundary=bundle.spatial.perimeter, vmin=0
    )
    axes[1, 0].set_title("C  Clipped 20 x 20 grid", loc="left", fontweight="bold")

    selected = beta.set_index("framework").loc[["Sectors", "Blocks", "Occupied grid cells"]]
    positions = np.arange(len(selected))
    axes[1, 1].bar(positions, selected["mean_nestedness"], color="#fdae61", label="Nestedness")
    axes[1, 1].bar(positions, selected["mean_turnover"], bottom=selected["mean_nestedness"], color="#2b83ba", label="Turnover")
    axes[1, 1].set_xticks(positions, ["Sectors", "Blocks", "Grid cells"])
    axes[1, 1].set_ylabel("Mean Jaccard dissimilarity")
    axes[1, 1].set_title("D  Beta diversity across grains", loc="left", fontweight="bold")
    axes[1, 1].spines[["top", "right"]].set_visible(False)
    axes[1, 1].legend(frameon=False)
    for position, value in zip(positions, selected["mean_jaccard_dissimilarity"]):
        axes[1, 1].text(position, value + 0.015, f"{value:.2f}", ha="center", fontsize=8)
    figure.suptitle("Species richness and turnover across spatial grains", fontweight="bold")
    base.save_figure(figure, "Figure3_diversity_multiscale", output, dpi)


def figure_origin(
    bundle: AnalysisBundle,
    tables: Dict[str, pd.DataFrame],
    output: Path,
    dpi: int,
) -> None:
    figure = plt.figure(figsize=(11, 5.2), constrained_layout=True)
    layout = figure.add_gridspec(1, 3, width_ratios=[1.55, 1.25, 1.15])
    axes = [figure.add_subplot(layout[0, index]) for index in range(3)]

    draw_base(axes[0], bundle)
    for origin in ["Exotic", "Native"]:
        subset = bundle.trees.loc[bundle.trees["Origin"].eq(origin)]
        axes[0].scatter(subset["x_m"], subset["y_m"], s=5, color=ZONE_COLORS[origin], alpha=0.6, linewidths=0, label=origin, zorder=4)
    axes[0].set_title("A  Individual records", loc="left", fontweight="bold")
    axes[0].legend(frameon=False, loc="lower left")

    combined = pd.concat(
        [
            tables["whole"].assign(display="Whole neighborhood"),
            tables["murillo"].assign(display=lambda frame: frame["group"]),
            tables["roads"].assign(display=lambda frame: frame["group"]),
            tables["contexts"].loc[tables["contexts"]["group"].isin(["Residential", "Tourist-cultural"])].assign(display=lambda frame: frame["group"]),
        ],
        ignore_index=True,
    )
    labels = combined["display"].tolist()
    native_fraction = combined["native_records"] / combined["n_records"]
    exotic_fraction = combined["exotic_records"] / combined["n_records"]
    unknown_fraction = combined["unknown_origin_records"] / combined["n_records"]
    positions = np.arange(len(combined))
    axes[1].barh(positions, native_fraction, color=ZONE_COLORS["Native"], label="Native")
    axes[1].barh(positions, exotic_fraction, left=native_fraction, color=ZONE_COLORS["Exotic"], label="Exotic")
    axes[1].barh(positions, unknown_fraction, left=native_fraction + exotic_fraction, color=ZONE_COLORS["Unknown"], label="Unknown")
    axes[1].set_yticks(positions, labels)
    axes[1].invert_yaxis()
    axes[1].set_xlim(0, 1)
    axes[1].set_xlabel("Fraction of records")
    axes[1].set_title("B  Origin composition", loc="left", fontweight="bold")
    axes[1].spines[["top", "right"]].set_visible(False)
    axes[1].legend(frameon=False, ncol=3, loc="lower center", bbox_to_anchor=(0.5, -0.2))

    sector_frame = pd.DataFrame(
        {"sector": list(bundle.sector_geometries), "geometry": list(bundle.sector_geometries.values())}
    ).merge(
        tables["sectors"][["group", "native_fraction_known"]],
        left_on="sector",
        right_on="group",
        how="left",
        validate="one_to_one",
    )
    base.draw_choropleth(
        axes[2], sector_frame["geometry"], sector_frame["native_fraction_known"], cmap="YlGn", colorbar_label="Native fraction", boundary=bundle.spatial.perimeter, vmin=0, vmax=1
    )
    axes[2].set_title("C  Native fraction by sector", loc="left", fontweight="bold")
    figure.suptitle("Native and exotic trees across spatial frameworks", fontweight="bold")
    base.save_figure(figure, "Figure5_origin_composition", output, dpi)


def kde_density_on_boundary(
    coordinates: np.ndarray,
    boundary: object,
    resolution: int = 130,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    min_x, min_y, max_x, max_y = boundary.bounds
    x_axis = np.linspace(min_x, max_x, resolution)
    y_axis = np.linspace(min_y, max_y, resolution)
    xx, yy = np.meshgrid(x_axis, y_axis)
    positions = np.vstack([xx.ravel(), yy.ravel()])
    density = gaussian_kde(coordinates.T, bw_method=0.25)(positions).reshape(xx.shape)
    inside = np.fromiter(
        (boundary.covers(Point(float(x), float(y))) for x, y in zip(xx.ravel(), yy.ravel())),
        dtype=bool,
        count=xx.size,
    ).reshape(xx.shape)
    density[~inside] = np.nan
    return xx, yy, density


def draw_kde_panel(
    ax: plt.Axes,
    bundle: AnalysisBundle,
    subset: pd.DataFrame,
    boundary: object,
    color_map: str,
    title: str,
) -> None:
    base.add_geometry_collection(ax, [boundary], facecolor="#fafafa", edgecolor="0.3", linewidth=0.8, zorder=0)
    xx, yy, density = kde_density_on_boundary(subset[["x_m", "y_m"]].to_numpy(), boundary)
    ax.pcolormesh(xx, yy, density, cmap=color_map, shading="auto", alpha=0.82, rasterized=True, zorder=1)
    ax.scatter(subset["x_m"], subset["y_m"], s=5, color="black", alpha=0.35, linewidths=0, zorder=3)
    base.set_map_extent(ax, boundary, padding=0.05)
    ax.set_axis_off()
    ax.set_title(title, loc="left", fontweight="bold")


def figure_context_kde(
    bundle: AnalysisBundle,
    tables: Dict[str, pd.DataFrame],
    output: Path,
    dpi: int,
) -> None:
    figure, axes = plt.subplots(1, 3, figsize=(11, 4.3), constrained_layout=True)
    draw_base(axes[0], bundle)
    coordinates = bundle.trees[["x_m", "y_m"]].to_numpy()
    xx, yy, density = kde_density_on_boundary(coordinates, bundle.spatial.perimeter)
    valid = density[np.isfinite(density)]
    threshold = float(np.quantile(valid, 0.80))
    axes[0].contourf(xx, yy, np.where(density >= threshold, density, np.nan), levels=8, cmap="Greens", alpha=0.65, zorder=2)
    for context in ["Residential", "Tourist-cultural"]:
        base.add_geometry_collection(axes[0], [bundle.context_buffers[context]], facecolor="none", edgecolor=ZONE_COLORS[context], linewidth=1.5, zorder=5)
    for row in bundle.context_anchors.itertuples(index=False):
        axes[0].scatter(row.x_m, row.y_m, marker="*", s=45, color=ZONE_COLORS[row.context], edgecolor="white", linewidth=0.5, zorder=7)
    axes[0].set_title("A  Upper 20% of neighborhood KDE", loc="left", fontweight="bold")
    axes[0].legend(
        handles=[Patch(facecolor="none", edgecolor=ZONE_COLORS[label], label=label) for label in ["Residential", "Tourist-cultural"]],
        frameon=False,
        loc="lower left",
    )

    residential = bundle.trees.loc[bundle.trees["context_area"].eq("Residential")]
    tourist = bundle.trees.loc[bundle.trees["context_area"].eq("Tourist-cultural")]
    draw_kde_panel(
        axes[1], bundle, residential, bundle.context_buffers["Residential"], "Purples", f"B  Residential KDE (n={len(residential)})"
    )
    draw_kde_panel(
        axes[2], bundle, tourist, bundle.context_buffers["Tourist-cultural"], "Oranges", f"C  Tourist-cultural KDE (n={len(tourist)})"
    )
    figure.suptitle("Descriptive tree-density surfaces in provisional contextual areas", fontweight="bold")
    base.save_figure(figure, "Figure6_context_KDE", output, dpi)


def figure_roads(
    bundle: AnalysisBundle,
    road_metrics: pd.DataFrame,
    road_units: pd.DataFrame,
    output: Path,
    dpi: int,
) -> None:
    order = ["Higher-order roads", "Local/secondary streets"]
    table = road_metrics.set_index("group").loc[order]
    colors = [ZONE_COLORS[label] for label in order]
    figure = plt.figure(figsize=(11, 7.2), constrained_layout=True)
    layout = figure.add_gridspec(2, 3, height_ratios=[1.65, 1])
    map_ax = figure.add_subplot(layout[0, :])
    axes = [figure.add_subplot(layout[1, index]) for index in range(3)]

    draw_base(map_ax, bundle)
    for road_class in ["Local/secondary streets", "Higher-order roads"]:
        subset = bundle.roads.loc[bundle.roads["analysis_class"].eq(road_class), "geometry"]
        plot_line_geometries(
            map_ax,
            subset,
            ZONE_COLORS[road_class],
            0.55 if road_class.startswith("Local") else 2.0,
            0.55 if road_class.startswith("Local") else 0.95,
            4,
        )
    for road_class in order:
        subset = bundle.trees.loc[bundle.trees["road_associated"] & bundle.trees["road_class"].eq(road_class)]
        map_ax.scatter(subset["x_m"], subset["y_m"], s=5, color=ZONE_COLORS[road_class], alpha=0.62, linewidths=0, zorder=6)
    map_ax.set_title("A  Trees associated with mapped road frontage", loc="left", fontweight="bold")
    map_ax.legend(
        handles=[Line2D([0], [0], color=ZONE_COLORS[label], lw=2, label=label) for label in order],
        frameon=False,
        loc="lower left",
    )

    values = [
        ("records_per_exposure", "Records per frontage km", "B  Abundance"),
        ("native_fraction_known", "Native fraction", "C  Origin"),
        ("carbon_Mg_per_exposure", "Tree-form carbon (Mg C/km)", "D  Carbon stock"),
    ]
    for ax, (column, ylabel, title) in zip(axes, values):
        ax.bar(order, table[column], color=colors)
        ax.set_ylabel(ylabel)
        ax.set_title(title, loc="left", fontweight="bold")
        ax.tick_params(axis="x", rotation=20)
        ax.spines[["top", "right"]].set_visible(False)
        ax.grid(axis="y", color="0.9", linewidth=0.5)
        for index, value in enumerate(table[column]):
            ax.text(index, value * 1.02 if value else 0.01, f"{value:.2f}", ha="center", fontsize=8)
    figure.suptitle("Higher-order and local-road tree assemblages", fontweight="bold")
    base.save_figure(figure, "Figure7_road_hierarchy", output, dpi)


def figure_carbon(
    bundle: AnalysisBundle,
    tables: Dict[str, pd.DataFrame],
    output: Path,
    dpi: int,
) -> None:
    figure, axes = plt.subplots(2, 2, figsize=(10, 8), constrained_layout=True)
    tree_inventory = bundle.carbon_inventory.loc[
        bundle.carbon_inventory["estimation_method"].eq("Chave 2014 D-H-WD")
    ].copy()
    species = (
        tree_inventory.groupby("Species_name")["agc_kg"].sum().div(1_000).sort_values(ascending=False).head(12).sort_values()
    )
    axes[0, 0].barh(species.index, species.values, color="#4c78a8")
    axes[0, 0].set_xlabel("Standing above-ground carbon (Mg C)")
    axes[0, 0].set_title("A  Leading species contributions", loc="left", fontweight="bold")
    axes[0, 0].spines[["top", "right"]].set_visible(False)

    block_frame = geometry_metric_table(
        bundle.spatial.blocks,
        tables["blocks"],
        ["Sector", "Manzana"],
        "group",
    )
    base.draw_choropleth(
        axes[0, 1],
        block_frame["geometry"],
        block_frame["carbon_Mg_per_exposure"].fillna(0),
        cmap="YlOrBr",
        colorbar_label="Mg C per km block perimeter",
        boundary=bundle.spatial.perimeter,
        vmin=0,
    )
    axes[0, 1].set_title("B  Carbon density by block", loc="left", fontweight="bold")

    murillo = tables["murillo"].set_index("group").loc[["Above Murillo", "Below Murillo"]]
    axes[1, 0].bar(murillo.index, murillo["carbon_Mg_per_exposure"], color=[ZONE_COLORS[index] for index in murillo.index])
    axes[1, 0].set_ylabel("Tree-form carbon (Mg C/km block frontage)")
    axes[1, 0].set_title("C  Murillo sides", loc="left", fontweight="bold")
    axes[1, 0].tick_params(axis="x", rotation=20)
    axes[1, 0].spines[["top", "right"]].set_visible(False)

    uncertainty = bundle.carbon_uncertainty.set_index("component").loc["Tree-form records"]
    mean = float(uncertainty["mean_agc_Mg"])
    lower = float(uncertainty["p2_5_agc_Mg"])
    upper = float(uncertainty["p97_5_agc_Mg"])
    axes[1, 1].errorbar([0], [mean], yerr=[[mean - lower], [upper - mean]], fmt="o", color="#4c78a8", capsize=5, markersize=7)
    axes[1, 1].set_xlim(-0.7, 0.7)
    axes[1, 1].set_xticks([0], ["Tree-form records"])
    axes[1, 1].set_ylabel("Standing above-ground carbon (Mg C)")
    axes[1, 1].set_title("D  Conditional uncertainty interval", loc="left", fontweight="bold")
    axes[1, 1].text(
        0.03,
        0.95,
        f"{mean:.1f} Mg C\n95% interval: {lower:.1f}-{upper:.1f}",
        transform=axes[1, 1].transAxes,
        ha="left",
        va="top",
        fontsize=8,
    )
    axes[1, 1].spines[["top", "right"]].set_visible(False)
    figure.suptitle("Standing above-ground carbon stock", fontweight="bold")
    base.save_figure(figure, "Figure8_carbon_stock", output, dpi)


def figure_coverage_standardized_diversity(
    coverage: pd.DataFrame, output: Path, dpi: int
) -> None:
    """Plot coverage-standardized Hill diversity as a supplementary figure."""
    frameworks = ["Murillo sides", "Road classes", "Context areas"]
    colors = {
        "Above Murillo": ZONE_COLORS["Above Murillo"],
        "Below Murillo": ZONE_COLORS["Below Murillo"],
        "Higher-order roads": ZONE_COLORS["Higher-order roads"],
        "Local/secondary streets": ZONE_COLORS["Local/secondary streets"],
        "Residential": ZONE_COLORS["Residential"],
        "Tourist-cultural": ZONE_COLORS["Tourist-cultural"],
    }
    figure, axes = plt.subplots(1, 3, figsize=(10.6, 3.5), constrained_layout=True)
    for axis, framework in zip(axes, frameworks):
        subset = coverage.loc[coverage["framework"].eq(framework)]
        groups = list(dict.fromkeys(subset["group"].tolist()))
        offsets = np.linspace(-0.11, 0.11, len(groups))
        for offset, group in zip(offsets, groups):
            values = subset.loc[subset["group"].eq(group)].sort_values("hill_order_q")
            x = values["hill_order_q"].to_numpy(dtype=float) + offset
            mean = values["effective_species_mean"].to_numpy(dtype=float)
            lower = mean - values["rarefaction_2_5pct"].to_numpy(dtype=float)
            upper = values["rarefaction_97_5pct"].to_numpy(dtype=float) - mean
            axis.errorbar(
                x,
                mean,
                yerr=np.vstack([lower, upper]),
                fmt="o",
                ms=5,
                capsize=2.5,
                color=colors[group],
                label=group,
            )
        target = float(subset["common_target_coverage"].iloc[0])
        axis.set_title(f"{framework}\ncommon coverage = {target:.3f}")
        axis.set_xticks([0, 1, 2], ["q0", "q1", "q2"])
        axis.set_ylabel("Effective number of species")
        axis.set_ylim(bottom=0)
        axis.spines[["top", "right"]].set_visible(False)
        axis.legend(frameon=False, fontsize=7, loc="upper right")
    figure.suptitle(
        "Coverage-standardized diversity contrasts", fontweight="bold"
    )
    base.save_figure(
        figure, "FigureS1_coverage_standardized_diversity", output, dpi
    )


def write_report(
    bundle: AnalysisBundle,
    tables: Dict[str, pd.DataFrame],
    beta: pd.DataFrame,
    road_paired: pd.DataFrame,
    inference: Mapping[str, pd.DataFrame],
    correction_audit: pd.DataFrame,
    output: Path,
) -> None:
    whole = tables["whole"].iloc[0]
    murillo = tables["murillo"].set_index("group")
    roads = tables["roads"].set_index("group")
    contexts = tables["contexts"].set_index("group")
    beta_lookup = beta.set_index("framework")
    uncertainty = bundle.carbon_uncertainty.set_index("component").loc["Tree-form records"]
    road_density_pair = road_paired.set_index("metric").loc["records_per_frontage_km"]
    murillo_nb = inference["murillo_negative_binomial"].iloc[0]
    road_test = inference["road_paired_inference"].iloc[0]

    content = f"""# Thematically organized results

## Decisions applied

- Above Murillo is the west/northwest side of the mapped avenue, shown on the left of the standard map orientation.
- Below Murillo is the east/southeast side.
- The {int(correction_audit['block_assignment_corrected'].sum())} approved block corrections are applied to analysis fields. Original block labels remain in the record-level export.
- Road association uses the approved {ROAD_ASSOCIATION_DISTANCE_M:.0f} m threshold.
- Residential and tourist-cultural areas use the approved provisional {CONTEXT_RADIUS_M:.0f} m buffers.
- The regular clipped 20 x 20 grid is retained.
- No easting/northing trend model is used in this analysis.

## Whole neighborhood

The census contains {int(whole['n_records']):,} records representing {int(whole['richness_q0'])} species. Hill diversity is q1 = {whole['diversity_q1']:.2f} and q2 = {whole['diversity_q2']:.2f}. The five most abundant species account for {100 * whole['top_five_fraction']:.1f}% of records. Among records of known origin, {100 * whole['native_fraction_known']:.1f}% are native.

## Avenida Murillo

Above Murillo contains {int(murillo.loc['Above Murillo', 'n_records']):,} records and below Murillo contains {int(murillo.loc['Below Murillo', 'n_records']):,}. Frontage-standardized abundance is {murillo.loc['Above Murillo', 'records_per_exposure']:.1f} and {murillo.loc['Below Murillo', 'records_per_exposure']:.1f} records per kilometre, respectively. Their pooled species lists have Jaccard dissimilarity {beta_lookup.loc['Murillo sides', 'mean_jaccard_dissimilarity']:.3f} and Bray-Curtis dissimilarity {beta_lookup.loc['Murillo sides', 'mean_bray_curtis_dissimilarity']:.3f}.

Across all 81 mapped blocks, the negative-binomial model with log block perimeter as an offset estimated an above-versus-below rate ratio of {murillo_nb['rate_ratio']:.2f} (HC3 95% CI {murillo_nb['rate_ratio_95pct_lower']:.2f}-{murillo_nb['rate_ratio_95pct_upper']:.2f}). This is an association among blocks, not a causal effect of the avenue.

## Road hierarchy

At 25 m, {int(roads.loc['Higher-order roads', 'n_records']):,} records are associated with higher-order roads and {int(roads.loc['Local/secondary streets', 'n_records']):,} with local or secondary streets. Aggregate abundance is {roads.loc['Higher-order roads', 'records_per_exposure']:.1f} versus {roads.loc['Local/secondary streets', 'records_per_exposure']:.1f} records per frontage kilometre. Native fractions among known-origin records are {roads.loc['Higher-order roads', 'native_fraction_known']:.3f} and {roads.loc['Local/secondary streets', 'native_fraction_known']:.3f}. In the {int(road_density_pair['paired_blocks_with_values'])} blocks with both frontage types, higher-order density exceeds local-road density in {int(road_density_pair['blocks_higher_greater'])} blocks, local-road density is greater in {int(road_density_pair['blocks_local_greater'])}, and {int(road_density_pair['ties'])} are tied. The median paired density difference is {road_test['median_paired_difference_records_per_km']:.1f} records per kilometre (sector-clustered bootstrap 95% CI {road_test['cluster_bootstrap_95pct_lower']:.1f}-{road_test['cluster_bootstrap_95pct_upper']:.1f}).

Eight of the 25 corrected records change from local/secondary to higher-order frontage when road class is recalculated from the corrected cadastral block. Their distance-based 25 m road-association status does not change.

## Contextual reference areas

The residential area contains {int(contexts.loc['Residential', 'n_records']):,} records and the tourist-cultural areas contain {int(contexts.loc['Tourist-cultural', 'n_records']):,}. Their q1 diversities are {contexts.loc['Residential', 'diversity_q1']:.2f} and {contexts.loc['Tourist-cultural', 'diversity_q1']:.2f}, and their known-origin native fractions are {contexts.loc['Residential', 'native_fraction_known']:.3f} and {contexts.loc['Tourist-cultural', 'native_fraction_known']:.3f}. The KDE panels are descriptive density surfaces within the approved buffers; they do not establish a land-use effect.

## Spatial grains and beta diversity

All 17 sectors contain records. Of 81 mapped blocks, {int((tables['blocks']['n_records'] > 0).sum())} contain records; of {len(tables['grid'])} clipped grid cells, {int((tables['grid']['n_records'] > 0).sum())} contain records. Mean Jaccard dissimilarity is {beta_lookup.loc['Sectors', 'mean_jaccard_dissimilarity']:.3f} among sectors, {beta_lookup.loc['Blocks', 'mean_jaccard_dissimilarity']:.3f} among occupied blocks, and {beta_lookup.loc['Occupied grid cells', 'mean_jaccard_dissimilarity']:.3f} among occupied grid cells. Turnover is reported separately from nestedness in `beta_diversity.csv`.

## Standing carbon stock

The primary tree-form estimate is {uncertainty['mean_agc_Mg']:.1f} Mg C, with a conditional 95% uncertainty interval of {uncertainty['p2_5_agc_Mg']:.1f}-{uncertainty['p97_5_agc_Mg']:.1f} Mg C. This equals approximately {uncertainty['mean_co2e_Mg']:.0f} Mg CO2e. It is standing above-ground stock, not annual capture or sequestration. Palm estimates remain a separate sensitivity scenario because the available measurements do not fully satisfy the palm equation's calibration conditions.

## Interpretation limits

- No observations are simulated. Monte Carlo draws are used only to propagate stated carbon-estimation uncertainty.
- Empty blocks and cells mean that the census contains no mapped records there; they are not independent confirmation that no unrecorded tree existed.
- Road hierarchy is a proxy for infrastructure and management context. The database does not identify the caretaker of each tree.
- Context buffers are approved working definitions, not independently mapped land-use polygons.
- Formal inference is limited to the block-replicated Murillo and paired-road density contrasts. Other framework comparisons remain descriptive. The same record can legitimately belong to one Murillo side, one road class, one context, one sector, one block, and one grid cell.
"""
    (output / "thematic_analysis_report.md").write_text(content, encoding="utf-8")


def export_outputs(
    bundle: AnalysisBundle,
    correction_audit: pd.DataFrame,
    combined_metrics: pd.DataFrame,
    tables: Dict[str, pd.DataFrame],
    beta: pd.DataFrame,
    road_units: pd.DataFrame,
    road_paired: pd.DataFrame,
    inference: Mapping[str, pd.DataFrame],
    output: Path,
) -> None:
    assignments = bundle.trees[
        [
            "ID",
            "X",
            "Y",
            "x_m",
            "y_m",
            "recorded_sector",
            "recorded_manzana",
            "recorded_block_id",
            "Sector",
            "Manzana",
            "analysis_block_id",
            "block_assignment_corrected",
            "nearest_block_id",
            "murillo_zone",
            "murillo_distance_m",
            "road_group",
            "road_group_using_recorded_block",
            "road_class_changed_by_block_correction",
            "nearest_road_distance_m",
            "context_area",
            "grid_id",
        ]
    ].copy()
    assignments = assignments.rename(
        columns={"Sector": "analysis_sector", "Manzana": "analysis_manzana"}
    )
    assignments.to_csv(output / "tree_spatial_assignments_corrected.csv", index=False)
    correction_audit.to_csv(output / "block_corrections_applied.csv", index=False)
    combined_metrics.to_csv(output / "framework_metrics.csv", index=False)
    for name, table in tables.items():
        table.to_csv(output / f"metrics_{name}.csv", index=False)
    beta.to_csv(output / "beta_diversity.csv", index=False)
    road_units.to_csv(output / "road_frontage_unit_metrics.csv", index=False)
    road_paired.to_csv(output / "road_paired_summary.csv", index=False)
    for name, table in inference.items():
        table.to_csv(output / f"{name}.csv", index=False)
    bundle.road_frontage.to_csv(output / "road_frontage_exposure.csv", index=False)
    bundle.carbon_uncertainty.to_csv(output / "carbon_uncertainty.csv", index=False)
    bundle.carbon_scenarios.to_csv(output / "carbon_scenarios.csv", index=False)
    bundle.height_diagnostics.to_csv(output / "height_model_cross_validation.csv", index=False)
    pd.DataFrame(
        [{"parameter": key, "value": value} for key, value in bundle.height_information.items()]
    ).to_csv(output / "height_model_parameters.csv", index=False)
    carbon_columns = [
        "ID",
        "Species_name",
        "Family",
        "Origin",
        "Form",
        "DBH_cm",
        "height_m",
        "height_model_m",
        "height_used_m",
        "height_source",
        "wood_density_g_cm3",
        "estimation_method",
        "agb_kg",
        "agc_kg",
        "co2e_kg",
    ]
    bundle.carbon_inventory[carbon_columns].to_csv(
        output / "carbon_individual_estimates.csv", index=False
    )


def validate_analysis(
    bundle: AnalysisBundle,
    correction_audit: pd.DataFrame,
    tables: Dict[str, pd.DataFrame],
    inference: Mapping[str, pd.DataFrame],
) -> pd.DataFrame:
    """Reconcile every framework to the same 1,256 individual records."""
    whole = tables["whole"].iloc[0]
    checks: List[Dict[str, object]] = []

    def add(check: str, observed: float, expected: float, tolerance: float = 0.0) -> None:
        passed = bool(abs(float(observed) - float(expected)) <= tolerance)
        checks.append(
            {
                "check": check,
                "observed": observed,
                "expected": expected,
                "tolerance": tolerance,
                "passed": passed,
            }
        )

    add("Unique individual IDs", bundle.trees["ID"].nunique(), 1_256)
    add("Approved block corrections", correction_audit["block_assignment_corrected"].sum(), 25)
    add("Above Murillo records", (bundle.trees["murillo_zone"] == "Above Murillo").sum(), 772)
    add("Below Murillo records", (bundle.trees["murillo_zone"] == "Below Murillo").sum(), 484)
    add("Higher-order road records after correction", (bundle.trees["road_group"] == "Higher-order roads").sum(), 373)
    add("Local/secondary road records after correction", (bundle.trees["road_group"] == "Local/secondary streets").sum(), 722)
    add("Road-unassigned records", (bundle.trees["road_group"] == "Not assigned (>25 m)").sum(), 161)
    add("Road classes changed by block correction", bundle.trees["road_class_changed_by_block_correction"].sum(), 8)
    add("Residential records", (bundle.trees["context_area"] == "Residential").sum(), 148)
    add("Tourist-cultural records", (bundle.trees["context_area"] == "Tourist-cultural").sum(), 147)
    add("Other context records", (bundle.trees["context_area"] == "Other/unclassified").sum(), 961)
    add("Occupied sectors", (tables["sectors"]["n_records"] > 0).sum(), 17)
    add("Occupied blocks", (tables["blocks"]["n_records"] > 0).sum(), 72)
    add("Occupied grid cells", (tables["grid"]["n_records"] > 0).sum(), 188)
    for name in ["murillo", "contexts", "sectors", "blocks", "grid"]:
        add(f"{name} record total", tables[name]["n_records"].sum(), 1_256)
        add(
            f"{name} carbon total",
            tables[name]["tree_form_carbon_Mg"].sum(),
            whole["tree_form_carbon_Mg"],
            1e-9,
        )
    add(
        "Whole origin total",
        whole["native_records"] + whole["exotic_records"] + whole["unknown_origin_records"],
        1_256,
    )
    add("Whole species richness", whole["richness_q0"], 81)
    add("Tree-form point estimate (Mg C)", whole["tree_form_carbon_Mg"], 117.03374966482598, 1e-6)
    murillo_blocks = inference["murillo_block_metrics"]
    add("Murillo block model units", len(murillo_blocks), 81)
    add("Murillo block model record total", murillo_blocks["n_records"].sum(), 1_256)
    add(
        "Murillo negative-binomial convergence",
        int(inference["murillo_negative_binomial"].iloc[0]["converged"]),
        1,
    )
    add(
        "Paired-road inference blocks",
        inference["road_paired_inference"].iloc[0]["paired_blocks"],
        47,
    )
    threshold_25 = inference["road_threshold_sensitivity"].set_index("threshold_m").loc[25.0]
    add("Road sensitivity 25 m higher-order records", threshold_25["higher_order_records"], 373)
    add("Road sensitivity 25 m local records", threshold_25["local_secondary_records"], 722)
    grid_20 = inference["grid_resolution_sensitivity"].set_index("grid_columns").loc[20]
    add("Grid sensitivity 20 x 20 clipped cells", grid_20["clipped_cells"], 258)
    add("Grid sensitivity 20 x 20 occupied cells", grid_20["occupied_cells"], 188)
    add("Coverage-standardized diversity rows", len(inference["coverage_standardized_diversity"]), 18)
    result = pd.DataFrame(checks)
    if not result["passed"].all():
        failed = result.loc[~result["passed"]]
        raise ValueError(f"Analysis validation failed: {failed.to_dict('records')}")
    return result


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=DEFAULT_OUTPUT,
        help="Destination directory (default: outputs/ at the repository root)",
    )
    parser.add_argument("--dpi", type=int, default=300)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    output = args.output_dir.resolve()
    output.mkdir(parents=True, exist_ok=True)
    base.configure_style()
    mpl.rcParams.update({"figure.titlesize": 12, "axes.titlesize": 9.5})
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", message=".*invalid value encountered in cast.*")
        bundle, correction_audit = prepare_bundle()
        combined_metrics, tables = build_metric_tables(bundle)
        beta = beta_diversity_table(bundle.trees)
        road_units = road_frontage_unit_table(bundle)
        road_paired = road_paired_summary(road_units)

        murillo_blocks = murillo_block_table(bundle)
        inference = {
            "murillo_block_metrics": murillo_blocks,
            "murillo_negative_binomial": murillo_negative_binomial(murillo_blocks),
            "road_paired_inference": road_density_paired_inference(road_units),
            "coverage_standardized_diversity": coverage_standardized_diversity(bundle),
            "road_threshold_sensitivity": road_threshold_sensitivity(bundle),
            "grid_resolution_sensitivity": grid_resolution_sensitivity(bundle),
        }

        validation = validate_analysis(bundle, correction_audit, tables, inference)

        export_outputs(
            bundle,
            correction_audit,
            combined_metrics,
            tables,
            beta,
            road_units,
            road_paired,
            inference,
            output,
        )
        validation.to_csv(output / "validation_checks.csv", index=False)
        figure_neighborhood_overview(bundle, output, args.dpi)
        figure_top_species(bundle, output, args.dpi)
        figure_multiscale(bundle, tables, beta, output, args.dpi)
        figure_distribution_murillo(bundle, tables["murillo"], output, args.dpi)
        figure_origin(bundle, tables, output, args.dpi)
        figure_context_kde(bundle, tables, output, args.dpi)
        figure_roads(bundle, tables["roads"], road_units, output, args.dpi)
        figure_carbon(bundle, tables, output, args.dpi)
        figure_coverage_standardized_diversity(
            inference["coverage_standardized_diversity"], output, args.dpi
        )
        write_report(
            bundle,
            tables,
            beta,
            road_paired,
            inference,
            correction_audit,
            output,
        )
    print(f"Thematic analysis complete: {output}")


if __name__ == "__main__":
    main()
