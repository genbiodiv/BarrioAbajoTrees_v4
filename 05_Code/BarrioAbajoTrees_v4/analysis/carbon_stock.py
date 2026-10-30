#!/usr/bin/env python3
"""Standing above-ground carbon calculations for the Barrio Abajo census.

The workflow deliberately estimates standing stock, not annual sequestration.
Trees and palms use different allometries, and shrubs/herbaceous records are
excluded because the inventory lacks measurements and calibrated equations for
those growth forms.
"""

from __future__ import annotations

import math
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import least_squares


ROOT = Path(__file__).resolve().parents[1]
WOOD_DENSITY = ROOT / "data" / "wood_density_lookup.csv"

SEED = 23_052_025
HEIGHT_FOLDS = 10
SPECIES_SHRINKAGE_K = 20.0
N_MONTE_CARLO = 4_000

# Chave et al. (2014), equation 4. D is cm, H is m, WD is g cm-3,
# and the result is kg of oven-dry above-ground biomass.
CHAVE_COEFFICIENT = 0.0673
CHAVE_EXPONENT = 0.976
CHAVE_LOG_RSE = 0.3577

# Goodman et al. (2013), family-level diameter-only palm model. The
# published calibration has Hstem > 3 m and 6 <= D < 40 cm.
PALM_INTERCEPT = -3.3488
PALM_DBH_SLOPE = 2.7483
PALM_LOG_RSE = 0.588

# Tropical angiosperm stem carbon fraction summarized by Thomas & Martin
# (2012), as implemented by the BIOMASS uncertainty workflow.
CARBON_FRACTION = 0.4713
CARBON_FRACTION_SD = 0.0206


def weibull_height(dbh_cm: np.ndarray, parameters: np.ndarray) -> np.ndarray:
    """Three-parameter asymptotic height model with a 1.3 m lower anchor."""
    asymptote, rate, shape = parameters
    dbh_cm = np.asarray(dbh_cm, dtype=float)
    return 1.3 + asymptote * (1.0 - np.exp(-rate * dbh_cm**shape))


def fit_weibull(dbh_cm: np.ndarray, height_m: np.ndarray) -> np.ndarray:
    """Robustly fit the local height-diameter relationship on the log scale."""
    dbh_cm = np.asarray(dbh_cm, dtype=float)
    height_m = np.asarray(height_m, dtype=float)

    def residual(parameters: np.ndarray) -> np.ndarray:
        prediction = weibull_height(dbh_cm, parameters)
        return np.log(prediction) - np.log(height_m)

    fit = least_squares(
        residual,
        x0=np.array([10.0, 0.03, 1.0]),
        bounds=(np.array([0.5, 1e-5, 0.1]), np.array([80.0, 2.0, 3.0])),
        loss="soft_l1",
        max_nfev=10_000,
    )
    if not fit.success:
        raise RuntimeError(f"Height model did not converge: {fit.message}")
    return fit.x


def species_height_adjustments(
    species: pd.Series,
    observed_height: np.ndarray,
    base_prediction: np.ndarray,
) -> pd.Series:
    """Shrink species mean log-residuals toward the neighborhood curve."""
    residuals = np.log(observed_height) - np.log(base_prediction)
    grouped = (
        pd.DataFrame({"Species_name": species.to_numpy(), "residual": residuals})
        .groupby("Species_name")["residual"]
        .agg(["mean", "size"])
    )
    grouped["adjustment"] = (
        grouped["mean"] * grouped["size"] / (grouped["size"] + SPECIES_SHRINKAGE_K)
    )
    return grouped["adjustment"]


def height_metrics(observed: np.ndarray, predicted: np.ndarray) -> dict[str, float]:
    residual = observed - predicted
    return {
        "rmse_m": float(np.sqrt(np.mean(residual**2))),
        "r_squared": float(1.0 - np.sum(residual**2) / np.sum((observed - observed.mean()) ** 2)),
        "log_rmse": float(np.sqrt(np.mean((np.log(observed) - np.log(predicted)) ** 2))),
        "median_absolute_relative_error": float(
            np.median(np.abs(observed - predicted) / observed)
        ),
        "mean_bias_m": float(np.mean(predicted - observed)),
    }


def cross_validate_height(calibration: pd.DataFrame) -> tuple[pd.DataFrame, float]:
    """Ten-fold random CV for the global and species-adjusted height models."""
    rng = np.random.default_rng(SEED)
    fold = np.empty(len(calibration), dtype=int)
    fold[rng.permutation(len(calibration))] = np.arange(len(calibration)) % HEIGHT_FOLDS
    observed = calibration["height_m"].to_numpy(dtype=float)
    predictions = {
        "Local Weibull": np.empty(len(calibration)),
        "Local Weibull + shrunken species residual": np.empty(len(calibration)),
    }

    for held_out in range(HEIGHT_FOLDS):
        train = calibration.loc[fold != held_out]
        test = calibration.loc[fold == held_out]
        parameters = fit_weibull(
            train["DBH_cm"].to_numpy(), train["height_m"].to_numpy()
        )
        train_base = weibull_height(train["DBH_cm"].to_numpy(), parameters)
        test_base = weibull_height(test["DBH_cm"].to_numpy(), parameters)

        global_log_rse = np.sqrt(
            np.mean((np.log(train["height_m"].to_numpy()) - np.log(train_base)) ** 2)
        )
        predictions["Local Weibull"][fold == held_out] = (
            test_base * np.exp(0.5 * global_log_rse**2)
        )

        adjustments = species_height_adjustments(
            train["Species_name"], train["height_m"].to_numpy(), train_base
        )
        train_adjusted = train_base * np.exp(
            train["Species_name"].map(adjustments).fillna(0.0).to_numpy()
        )
        adjusted_log_rse = np.sqrt(
            np.mean(
                (
                    np.log(train["height_m"].to_numpy())
                    - np.log(train_adjusted)
                )
                ** 2
            )
        )
        predictions["Local Weibull + shrunken species residual"][fold == held_out] = (
            test_base
            * np.exp(test["Species_name"].map(adjustments).fillna(0.0).to_numpy())
            * np.exp(0.5 * adjusted_log_rse**2)
        )

    rows = []
    for model_name, predicted in predictions.items():
        rows.append(
            {
                "model": model_name,
                "folds": HEIGHT_FOLDS,
                "n_observations": len(calibration),
                **height_metrics(observed, predicted),
            }
        )
    diagnostics = pd.DataFrame(rows)
    selected_log_rmse = float(
        diagnostics.loc[
            diagnostics["model"].eq(
                "Local Weibull + shrunken species residual"
            ),
            "log_rmse",
        ].iloc[0]
    )
    return diagnostics, selected_log_rmse


def fit_height_system(
    trees: pd.DataFrame,
) -> tuple[np.ndarray, pd.Series, dict[str, float], pd.DataFrame]:
    calibration = trees.loc[
        trees["Form"].eq("Tree")
        & trees["DBH_cm"].ge(5.0)
        & trees["height_m"].ge(1.3),
        ["Species_name", "DBH_cm", "height_m"],
    ].copy()
    if len(calibration) < 100:
        raise ValueError("Too few tree height records for a local model")

    diagnostics, cv_log_rse = cross_validate_height(calibration)
    parameters = fit_weibull(
        calibration["DBH_cm"].to_numpy(), calibration["height_m"].to_numpy()
    )
    base = weibull_height(calibration["DBH_cm"].to_numpy(), parameters)
    adjustments = species_height_adjustments(
        calibration["Species_name"], calibration["height_m"].to_numpy(), base
    )
    adjusted_median = base * np.exp(
        calibration["Species_name"].map(adjustments).fillna(0.0).to_numpy()
    )
    in_sample_log_rse = float(
        np.sqrt(
            np.mean(
                (
                    np.log(calibration["height_m"].to_numpy())
                    - np.log(adjusted_median)
                )
                ** 2
            )
        )
    )
    information = {
        "n_calibration": len(calibration),
        "asymptote_parameter": float(parameters[0]),
        "rate_parameter": float(parameters[1]),
        "shape_parameter": float(parameters[2]),
        "in_sample_log_rse": in_sample_log_rse,
        "smearing_factor": float(np.exp(0.5 * in_sample_log_rse**2)),
        "cv_log_rse": cv_log_rse,
    }
    return parameters, adjustments, information, diagnostics


def attach_wood_density(trees: pd.DataFrame) -> pd.DataFrame:
    lookup = pd.read_csv(WOOD_DENSITY)
    if lookup["Species_name"].duplicated().any():
        raise ValueError("Wood-density lookup contains duplicate species")

    tree_species = set(trees.loc[trees["Form"].eq("Tree"), "Species_name"])
    lookup_species = set(lookup["Species_name"])
    missing = tree_species - lookup_species
    extra = lookup_species - tree_species
    if missing or extra:
        raise ValueError(
            f"Wood-density lookup mismatch; missing={sorted(missing)}, extra={sorted(extra)}"
        )

    return trees.merge(
        lookup.drop(columns="Genus"),
        on="Species_name",
        how="left",
        validate="many_to_one",
    )


def calculate_inventory_estimates(
    trees: pd.DataFrame,
    parameters: np.ndarray,
    adjustments: pd.Series,
    height_information: dict[str, float],
) -> pd.DataFrame:
    inventory = attach_wood_density(trees)
    valid_dbh = inventory["DBH_cm"].gt(0)
    tree_with_dbh = inventory["Form"].eq("Tree") & valid_dbh
    base = np.full(len(inventory), np.nan)
    base[tree_with_dbh] = weibull_height(
        inventory.loc[tree_with_dbh, "DBH_cm"].to_numpy(), parameters
    )
    adjustment = (
        inventory["Species_name"].map(adjustments).fillna(0.0).to_numpy()
    )
    inventory["height_model_m"] = (
        base * np.exp(adjustment) * height_information["smearing_factor"]
    )
    observed_height = inventory["height_m"].ge(1.3)
    inventory["height_used_m"] = inventory["height_m"].where(
        observed_height, inventory["height_model_m"]
    )
    inventory.loc[~inventory["Form"].eq("Tree"), "height_used_m"] = np.nan
    inventory["height_source"] = "not used"
    inventory.loc[
        inventory["Form"].eq("Tree") & observed_height, "height_source"
    ] = "observed"
    inventory.loc[
        inventory["Form"].eq("Tree") & ~observed_height & valid_dbh,
        "height_source",
    ] = "local model"
    inventory.loc[
        inventory["Form"].eq("Tree") & ~valid_dbh, "height_source"
    ] = "unavailable: missing DBH"

    inventory["tree_in_chave_diameter_range"] = (
        inventory["Form"].eq("Tree") & inventory["DBH_cm"].ge(5.0)
    )
    inventory["palm_in_goodman_diameter_range"] = (
        inventory["Form"].eq("Palm")
        & inventory["DBH_cm"].ge(6.0)
        & inventory["DBH_cm"].lt(40.0)
    )
    inventory["palm_total_height_gt_3m"] = (
        inventory["Form"].eq("Palm") & inventory["height_m"].gt(3.0)
    )

    valid_stems = inventory["stem_count"].gt(0)
    tree_mask = (
        inventory["Form"].eq("Tree")
        & valid_dbh
        & valid_stems
        & inventory["wood_density_g_cm3"].gt(0)
        & inventory["height_used_m"].ge(1.3)
    )
    palm_mask = inventory["Form"].eq("Palm") & valid_dbh & valid_stems

    inventory["agb_kg"] = np.nan
    tree_single_stem = CHAVE_COEFFICIENT * (
        inventory.loc[tree_mask, "wood_density_g_cm3"]
        * inventory.loc[tree_mask, "height_used_m"]
        * inventory.loc[tree_mask, "DBH_cm"] ** 2
    ) ** CHAVE_EXPONENT
    inventory.loc[tree_mask, "agb_kg"] = (
        inventory.loc[tree_mask, "stem_count"] * tree_single_stem
    )

    palm_single_stem = np.exp(
        PALM_INTERCEPT
        + PALM_DBH_SLOPE * np.log(inventory.loc[palm_mask, "DBH_cm"])
        + 0.5 * PALM_LOG_RSE**2
    )
    inventory.loc[palm_mask, "agb_kg"] = (
        inventory.loc[palm_mask, "stem_count"] * palm_single_stem
    )
    inventory["agc_kg"] = inventory["agb_kg"] * CARBON_FRACTION
    inventory["co2e_kg"] = inventory["agc_kg"] * 44.0 / 12.0

    inventory["agb_model_height_for_all_kg"] = np.nan
    model_height_single_stem = CHAVE_COEFFICIENT * (
        inventory.loc[tree_mask, "wood_density_g_cm3"]
        * inventory.loc[tree_mask, "height_model_m"]
        * inventory.loc[tree_mask, "DBH_cm"] ** 2
    ) ** CHAVE_EXPONENT
    inventory.loc[tree_mask, "agb_model_height_for_all_kg"] = (
        inventory.loc[tree_mask, "stem_count"] * model_height_single_stem
    )
    inventory["agb_one_mean_stem_kg"] = np.nan
    inventory.loc[tree_mask, "agb_one_mean_stem_kg"] = tree_single_stem

    inventory["estimation_method"] = "excluded: no calibrated growth-form model"
    inventory.loc[
        inventory["Form"].eq("Tree") & ~valid_dbh,
        "estimation_method",
    ] = "excluded: missing DBH"
    inventory.loc[tree_mask, "estimation_method"] = "Chave 2014 D-H-WD"
    inventory.loc[
        palm_mask & inventory["palm_in_goodman_diameter_range"],
        "estimation_method",
    ] = "Goodman 2013 diameter-only; diameter in calibration range"
    inventory.loc[
        palm_mask & ~inventory["palm_in_goodman_diameter_range"],
        "estimation_method",
    ] = "Goodman 2013 diameter-only; extrapolated diameter"
    return inventory


def scenario_table(inventory: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, object]] = []

    def add(
        group: str,
        scenario: str,
        mask: pd.Series,
        agb: pd.Series,
        interpretation: str,
    ) -> None:
        values = agb.loc[mask].dropna()
        rows.append(
            {
                "group": group,
                "scenario": scenario,
                "records": int(mask.sum()),
                "recorded_stems": int(inventory.loc[mask, "stem_count"].sum()),
                "agb_Mg": float(values.sum() / 1_000.0),
                "agc_Mg": float(values.sum() * CARBON_FRACTION / 1_000.0),
                "co2e_Mg": float(
                    values.sum() * CARBON_FRACTION * 44.0 / 12.0 / 1_000.0
                ),
                "interpretation": interpretation,
            }
        )

    tree = inventory["Form"].eq("Tree") & inventory["agb_kg"].notna()
    palm = inventory["Form"].eq("Palm") & inventory["agb_kg"].notna()
    palm_range = palm & inventory["palm_in_goodman_diameter_range"]
    giant_ceiba = (
        tree
        & inventory["Species_name"].eq("Ceiba pentandra")
        & inventory["DBH_cm"].ge(100.0)
    )
    observed_tree = tree & inventory["height_source"].eq("observed")

    add(
        "Trees",
        "Inventory estimate: observed height where available, otherwise local model",
        tree,
        inventory["agb_kg"],
        "Recommended headline stock for tree-form records; includes DBH < 5 cm extrapolation.",
    )
    add(
        "Trees",
        "Local height model used for every tree",
        tree,
        inventory["agb_model_height_for_all_kg"],
        "Checks sensitivity to field height estimates.",
    )
    add(
        "Trees",
        "Only DBH >= 5 cm",
        tree & inventory["DBH_cm"].ge(5.0),
        inventory["agb_kg"],
        "Restricts Chave allometry to its reported diameter calibration threshold.",
    )
    add(
        "Trees",
        "One mean stem per record",
        tree,
        inventory["agb_one_mean_stem_kg"],
        "Lower structural scenario that ignores additional recorded stems.",
    )
    add(
        "Trees",
        "Excluding Ceiba pentandra records with DBH >= 100 cm",
        tree & ~giant_ceiba,
        inventory["agb_kg"],
        "Influence check for four giant stems.",
    )
    add(
        "Trees",
        "Observed heights only",
        observed_tree,
        inventory["agb_kg"],
        "Incomplete-case subtotal; not an estimate for the full inventory.",
    )
    add(
        "Palms",
        "All palm diameters",
        palm,
        inventory["agb_kg"],
        "Exploratory; 24 records fall outside the Goodman diameter range.",
    )
    add(
        "Palms",
        "Only 6 <= DBH < 40 cm",
        palm_range,
        inventory["agb_kg"],
        "Diameter-supported subset, although stem height > 3 m cannot be verified from total height.",
    )
    add(
        "Combined",
        "Trees plus all palms",
        tree | palm,
        inventory["agb_kg"],
        "Illustrative total; palm component is substantially extrapolative.",
    )
    add(
        "Combined",
        "Trees plus diameter-range palms",
        tree | palm_range,
        inventory["agb_kg"],
        "More conservative combined subtotal; still cannot verify palm stem-height eligibility.",
    )
    return pd.DataFrame(rows)


def monte_carlo_uncertainty(
    inventory: pd.DataFrame, height_log_rse: float
) -> pd.DataFrame:
    """Propagate selected inputs; this is not a full urban-transfer interval."""
    rng = np.random.default_rng(SEED)
    tree = inventory.loc[
        inventory["Form"].eq("Tree") & inventory["agb_kg"].notna()
    ].reset_index(drop=True)
    palm = inventory.loc[
        inventory["Form"].eq("Palm") & inventory["agb_kg"].notna()
    ].reset_index(drop=True)

    species_table = (
        tree[
            ["Species_name", "wood_density_g_cm3", "wood_density_sd_g_cm3"]
        ]
        .drop_duplicates("Species_name")
        .reset_index(drop=True)
    )
    species_code = pd.Categorical(
        tree["Species_name"], categories=species_table["Species_name"]
    ).codes
    tree_total = np.empty(N_MONTE_CARLO)
    palm_total = np.empty(N_MONTE_CARLO)
    batch_size = 250

    for start in range(0, N_MONTE_CARLO, batch_size):
        stop = min(start + batch_size, N_MONTE_CARLO)
        size = stop - start

        species_wd = np.clip(
            rng.normal(
                species_table["wood_density_g_cm3"].to_numpy()[:, None],
                species_table["wood_density_sd_g_cm3"].to_numpy()[:, None],
                size=(len(species_table), size),
            ),
            0.08,
            1.39,
        )
        wd = species_wd[species_code, :]
        height_mean = tree["height_used_m"].to_numpy()[:, None]
        height = np.repeat(height_mean, size, axis=1)
        imputed = tree["height_source"].eq("local model").to_numpy()
        if imputed.any():
            height[imputed, :] = np.exp(
                rng.normal(
                    np.log(height_mean[imputed]) - 0.5 * height_log_rse**2,
                    height_log_rse,
                    size=(int(imputed.sum()), size),
                )
            )
        allometry_factor = np.exp(
            rng.normal(
                -0.5 * CHAVE_LOG_RSE**2,
                CHAVE_LOG_RSE,
                size=(len(tree), size),
            )
        )
        tree_agb = (
            tree["stem_count"].to_numpy()[:, None]
            * CHAVE_COEFFICIENT
            * (
                wd
                * height
                * tree["DBH_cm"].to_numpy()[:, None] ** 2
            )
            ** CHAVE_EXPONENT
            * allometry_factor
        )
        tree_carbon_fraction = np.clip(
            rng.normal(
                CARBON_FRACTION,
                CARBON_FRACTION_SD,
                size=(len(tree), size),
            ),
            0.0,
            1.0,
        )
        tree_total[start:stop] = np.sum(
            tree_agb * tree_carbon_fraction, axis=0
        ) / 1_000.0

        palm_agb = (
            palm["stem_count"].to_numpy()[:, None]
            * np.exp(
                PALM_INTERCEPT
                + PALM_DBH_SLOPE
                * np.log(palm["DBH_cm"].to_numpy()[:, None])
                + rng.normal(
                    0.0, PALM_LOG_RSE, size=(len(palm), size)
                )
            )
        )
        palm_carbon_fraction = np.clip(
            rng.normal(
                CARBON_FRACTION,
                CARBON_FRACTION_SD,
                size=(len(palm), size),
            ),
            0.0,
            1.0,
        )
        palm_total[start:stop] = np.sum(
            palm_agb * palm_carbon_fraction, axis=0
        ) / 1_000.0

    rows = []
    note = (
        "Propagates generic wood-density uncertainty by species, imputed-height "
        "residuals, equation residuals, and carbon-fraction variation. It excludes "
        "DBH/field-height measurement error, height-model parameter uncertainty, "
        "and systematic transfer error from forest to urban tree architecture."
    )
    for component, values in (
        ("Tree-form records", tree_total),
        ("Palms (all diameters)", palm_total),
        ("Combined illustrative total", tree_total + palm_total),
    ):
        rows.append(
            {
                "component": component,
                "iterations": N_MONTE_CARLO,
                "mean_agc_Mg": float(values.mean()),
                "sd_agc_Mg": float(values.std(ddof=1)),
                "p2_5_agc_Mg": float(np.quantile(values, 0.025)),
                "median_agc_Mg": float(np.median(values)),
                "p97_5_agc_Mg": float(np.quantile(values, 0.975)),
                "mean_co2e_Mg": float(values.mean() * 44.0 / 12.0),
                "p2_5_co2e_Mg": float(np.quantile(values, 0.025) * 44.0 / 12.0),
                "p97_5_co2e_Mg": float(np.quantile(values, 0.975) * 44.0 / 12.0),
                "scope_note": note,
            }
        )
    return pd.DataFrame(rows)


