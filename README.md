# Barrio Abajo urban-tree census and reproducible analysis

This repository contains the curated data, frozen spatial inputs, Python workflow, derived tables, publication figures, and manuscript associated with the census of publicly accessible trees and other arborescent plants in Barrio Abajo, Barranquilla, Colombia.

> **Release v4.** The current workflow is organized around the complete neighborhood, dominant species, multiscale diversity, Avenida Murillo, road hierarchy, provisional contextual areas, and standing above-ground carbon. The superseded easting/northing regression suite is not part of this release.

The repository is designed to support three uses:

- reproduce every table and figure from the frozen inputs;
- inspect the record-level spatial assignments, corrections, and carbon estimates; and
- reuse the curated census with a clear account of scope, units, analytical definitions, and limitations.

## Contents

- [Study at a glance](#study-at-a-glance)
- [Quick start](#quick-start)
- [Repository map](#repository-map)
- [Source data](#source-data)
- [Operational definitions](#operational-definitions)
- [Analytical design](#analytical-design)
- [Workflow](#workflow)
- [Figures](#figures)
- [Derived result files](#derived-result-files)
- [Verified release checkpoints](#verified-release-checkpoints)
- [Interpretation boundaries](#interpretation-boundaries)
- [Citation](#citation)
- [License and attribution](#license-and-attribution)

## Study at a glance

| Item | Description |
|---|---|
| Study area | Barrio Abajo, Barranquilla, Colombia; approximately 86 ha |
| Field period | 21 non-consecutive dates from 23 March to 13 October 2023 |
| Sampling target | Trees, palms, shrubs, and arborescent herbs taller than 1.3 m and visible from streets or other publicly accessible places |
| Census size | 1,256 unique records |
| Taxonomic coverage | 81 species, 68 genera, and 30 families |
| Growth forms | 1,037 trees, 193 palms, 23 shrubs, and 3 arborescent herbs |
| Spatial units | Whole neighborhood; Murillo sides; road classes; contextual areas; 17 sectors; 81 cadastral blocks; and a clipped 20-column × 20-row grid |
| Analytical CRS | WGS 84 / UTM zone 18N, EPSG:32618 |
| Associated manuscript | *Diversity, spatial structure, and standing carbon stock of the urban forest in Barrio Abajo, Barranquilla* |

The biological inventory is a census of the defined publicly accessible sampling frame at one period in time. It is not a complete census of private backyards and it is not a repeated-measures study.

## Quick start

Python 3.8 or newer is required. The release workflow has been verified with Python 3.12 and package versions satisfying `requirements.txt`.

From the repository root on macOS or Linux:

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install -r requirements.txt
python analysis/generate_all_figures.py
```

On Windows PowerShell, use:

```powershell
python -m venv .venv
.venv\Scripts\Activate.ps1
python -m pip install --upgrade pip
python -m pip install -r requirements.txt
python analysis/generate_all_figures.py
```

By default, the workflow writes to `outputs/`. The directory is intentionally ignored by Git so local reruns do not modify the versioned release files.

Useful options:

```bash
# Write to another directory
python analysis/generate_all_figures.py --output-dir test_outputs

# Change raster output resolution
python analysis/generate_all_figures.py --dpi 600

# Display all command-line options
python analysis/generate_all_figures.py --help
```

The analysis uses only files included in this repository after the Python packages are installed. It does not download observations, maps, or trait data while running.

## What a successful run produces

A complete run creates:

- eight main figures and one supplementary figure, each in PDF and PNG;
- 26 derived CSV files;
- a concise Markdown results report; and
- `validation_checks.csv`, which must show that all 36 checks passed.

The script stops with an error if a reconciliation or inferential-design check fails. The archived `results/` and `figures/` directories are the verified release snapshots. Locally regenerated files may differ at the final decimal place because of numerical-library versions, and PDFs or PNGs may differ in metadata or font rendering. These differences should not change the reported conclusions or cause validation failure.

## Repository map

| Path | Role |
|---|---|
| `data/Data_S1_Barrio_Abajo_tree_inventory.xlsx` | Curated biological inventory and species traits; the manuscript Data S1 file |
| `data/spatial/` | Neighborhood perimeter, sector, and cadastral-block geometries |
| `data/context/` | Colombia, Atlántico, and Barranquilla context layers used in Figure 1 |
| `data/barrio_abajo_roads.geojson` | Frozen OpenStreetMap motor-road centerlines with analysis classes and source metadata |
| `data/wood_density_lookup.csv` | Species- or genus-level wood-density values and uncertainty fields |
| `analysis/generate_all_figures.py` | Public one-command entry point |
| `analysis/thematic_analysis.py` | Spatial assignments, metrics, inferential analyses, sensitivity checks, validation, and figures |
| `analysis/carbon_stock.py` | Height model, tree and palm allometry, carbon scenarios, and Monte Carlo uncertainty |
| `analysis/analysis_utils.py` | Input validation, projections, grid construction, diversity functions, and shared plotting utilities |
| `results/` | Versioned, machine-readable outputs from the release run |
| `figures/` | Versioned manuscript figures in PDF and PNG |
| `manuscript/` | Manuscript v4 in Word and Markdown |
| `requirements.txt` | Supported Python dependency ranges |
| `LICENSE` | MIT license for repository code and documentation |

## Source data

### Data S1 workbook

All biological calculations begin with `data/Data_S1_Barrio_Abajo_tree_inventory.xlsx`.

| Sheet | Dimensions | Role in the workflow |
|---|---:|---|
| `README` | Documentation sheet | Survey dates, scope, record count, and source notes |
| `ArbolesIndividuales_OK` | 1,256 rows × 17 columns | One row per inventoried individual; primary biological input |
| `CaracteristicasEspecies` | 81 rows × 15 columns | One row per species; `Origin` and `Form` are joined to the individual records |

Derived workbook sheets are not used. No biological records are generated or simulated.

#### Individual-record fields

| Field or field group | Meaning and unit |
|---|---|
| `Status` | Quality-control status; all analytical rows are marked `OK` |
| `ID` | Unique individual identifier; uniqueness is checked before analysis |
| `Species_name`, `Genus`, `Family` | Standardized taxonomic fields |
| `Y`, `X` | Latitude and longitude in WGS 84 decimal degrees |
| `common_nam` | Recorded common name, when available |
| `Sector`, `Manzana`, `Predio`, `num_arbol` | Recorded cadastral sector, block, property, and within-property plant identifier |
| `TRUNK` | Number of recorded stems |
| `AVGCIR` | Mean circumference at breast height, cm |
| `AVGDBH` | Mean diameter at breast height, cm |
| `ESTH` | Estimated total height, m; missing values are retained and modeled only for the carbon analysis |
| `DESC` | Original record description and field notes |

`TRUNK`, `AVGCIR`, and `AVGDBH` each contain eight missing values; `ESTH` contains 404 missing values. Missing values are not silently converted to zero.

#### Species-trait fields

| Field or field group | Meaning |
|---|---|
| `Abundance`, `Species_name` | Census abundance and standardized species name |
| `Origin` | Native, exotic, or unknown origin |
| `Form` | Tree, palm, shrub, or arborescent herb |
| `Charismatic`, `Usos` | Charismatic-species indicator and documented use category |
| `Dispersion`, `Polinator` | Dispersal and pollinator categories |
| `Leaf`, `GrowthRate` | Leaf habit and growth-rate category |
| `Longevity_cat`, `Longevity_range` | Longevity categories |
| `MatHeight_cm` | Reported mature height, cm |
| `UICN_long`, `UICN_short` | Conservation category labels |

The present workflow uses only `Species_name`, `Origin`, and `Form` from the trait sheet. The other trait fields are retained for documentation and future analyses.

### Spatial and auxiliary inputs

| Input | Contents and treatment |
|---|---|
| `data/spatial/Perimetro del Barrio abajo.*` | One neighborhood perimeter polygon in WGS 84 |
| `data/spatial/Manzanas Barrio Abajo.*` | 81 cadastral block polygons in WGS 84 |
| `data/spatial/sectores.*` | Sector polygons in WGS 84; multipart labels are merged for analysis |
| `data/context/colombia_sinSA.*` | Colombia departmental context polygons in WGS 84 |
| `data/context/4326_atpio.*` | Atlántico municipal context polygons in WGS 84 |
| `data/context/quilla_mpio.*` | Barranquilla municipality context polygon supplied in Web Mercator |
| `data/barrio_abajo_roads.geojson` | 248 frozen motor-road features from an OpenStreetMap API extract downloaded 29 October 2026 |
| `data/wood_density_lookup.csv` | Wood density in g cm⁻³, uncertainty, match level, source taxon, and provenance for every tree-form species |

The analysis projects biological points and working geometries to EPSG:32618 before calculating distances, areas, frontage, or grid dimensions.

## Operational definitions

| Term | Definition used here |
|---|---|
| Census record | One inventoried plant taller than 1.3 m and visible from a street or public space; “record” is used because palms, shrubs, and arborescent herbs are included with trees |
| Above Murillo | West–northwest side of Avenida Murillo, shown on the left in the standard map orientation |
| Below Murillo | East–southeast side of Avenida Murillo |
| Corrected block | One of 25 records whose original block label was more than 25 m from the recorded point and less plausible than the nearest mapped block; original labels remain in the audit export |
| Road-associated record | A record within the selected distance of a frozen mapped motor-road centerline; 25 m is the primary threshold |
| Higher-order road | Vía 40, Avenida Murillo, Avenida Carrera 46, Carrera 54, or Carrera 50, following the Category I–II corridors crossing or delimiting the study area |
| Local/secondary street | Any other mapped motor road in the frozen extract |
| Road class assignment | Class nearest the boundary frontage of the corrected cadastral block; frontage is estimated from 1 m samples along block boundaries |
| Residential reference area | The neighborhood-clipped 150 m buffer around La Castellana Park |
| Tourist-cultural reference area | The neighborhood-clipped union of 150 m buffers around Old Customs Square and the Caribbean Cultural Park |
| Retained grid | A regular 20-column × 20-row grid clipped to the neighborhood, yielding 258 cells, of which 188 contain records |
| Known-origin fraction | Native or exotic records only; the two unknown-origin records are excluded from the denominator |
| Standing carbon | Estimated above-ground carbon present at the census date; not annual sequestration or avoided emissions |

> **Grid warning:** “20 × 20 grid” describes the number of columns and rows, not 20 m × 20 m cells. Before clipping, the retained cells are approximately 72.77 m wide × 53.42 m high.

## Analytical design

| Framework | Groups or units | Exposure or standardization | Inferential status |
|---|---|---|---|
| Whole neighborhood | All 1,256 records | None | Descriptive census summary |
| Murillo sides | Above versus below Murillo | Block-boundary frontage, km | Negative-binomial count model across all 81 blocks, including zero-record blocks; log frontage offset; HC3 covariance |
| Road hierarchy | Higher-order versus local/secondary frontage | Mapped block frontage, km | Paired comparison in 47 blocks containing both frontage types; sector-clustered bootstrap and exact sign test |
| Context areas | Residential, tourist-cultural, and other/unclassified | Area, ha | Descriptive; reference areas are provisional and unreplicated |
| Sectors | 17 occupied sector labels | None | Descriptive multiscale summary and pairwise beta diversity |
| Blocks | 81 mapped blocks, 72 occupied | None | Descriptive multiscale summary; blocks are also the replicated inferential units above |
| Grid | 258 clipped cells, 188 occupied | None | Descriptive multiscale summary and grain sensitivity |

Formal inference is intentionally limited to the replicated Murillo and road-frontage comparisons. Pooled species lists, context-area contrasts, kernel-density maps, carbon totals, and multiscale summaries are descriptive.

### Diversity, composition, and dominance

- `q0` is observed species richness.
- `q1` is the exponential of Shannon entropy and is interpreted as the effective number of common species.
- `q2` is the inverse Simpson concentration and emphasizes abundant species.
- Unequal pooled groups are compared with individual-based rarefaction and coverage-standardized Hill numbers.
- Jaccard dissimilarity describes presence–absence composition and is partitioned into turnover and nestedness.
- Bray–Curtis dissimilarity describes abundance composition.
- Dominance is reported as the share of the leading species and the combined share of the five most abundant species.

Coverage standardization uses 4,000 multivariate-hypergeometric rarefactions for each Hill order and paired framework. The intervals describe rarefaction variation conditional on the observed assemblages; they are not estimates of unseen species.

### Targeted inference and sensitivity

- The Murillo model is a negative-binomial count model with side as the predictor, `log(block perimeter km)` as the offset, maximum-likelihood dispersion, and HC3 sandwich covariance.
- The road model reports the median within-block difference in records per frontage kilometer. Its 95% interval uses 20,000 resamples clustered by the 17 sectors; an exact paired sign test is also reported.
- Road association is repeated at 15, 25, and 40 m.
- Grid composition is repeated at 15 × 15, 20 × 20, and 25 × 25 columns and rows.

### Basal area and carbon

Approximate basal area is calculated as:

```text
stem count × π × (DBH / 200)²
```

The database stores a mean DBH for multi-stem records rather than the diameter of each stem, so multi-stem basal area is approximate.

The primary carbon estimate includes 1,029 tree-form records with DBH. Above-ground biomass follows Chave et al. (2014), equation 4:

```text
AGB (kg) = 0.0673 × (wood density × height × DBH²)^0.976
```

Here wood density is in g cm⁻³, height is in m, and DBH is in cm. Observed height is used for 685 eligible records; height for 344 records is predicted with a local three-parameter Weibull model plus shrunken species residuals. The selected height model was evaluated with 10-fold cross-validation. Biomass is converted to carbon with a fraction of 0.4713 and to CO₂ equivalent with 44/12.

The tree-form uncertainty analysis uses 4,000 Monte Carlo iterations and propagates the stated wood-density, imputed-height, allometric-equation, and carbon-fraction uncertainty. It does not include DBH or field-height measurement error, height-model parameter uncertainty, or systematic transfer error from forest allometry to managed urban trees.

Palms are estimated separately with Goodman et al. (2013) and remain a sensitivity scenario because stem-height eligibility cannot be verified and some diameters fall outside the published calibration range. Shrubs and arborescent herbs are not assigned carbon estimates.

## Workflow

The public entry point performs these steps in order:

1. Read and validate the individual and species-trait sheets.
2. Project coordinates and geometries to EPSG:32618.
3. Preserve the recorded cadastral labels and apply the 25 approved analysis-only block corrections.
4. Assign every record independently to Murillo side, road context, provisional context area, sector, block, and grid cell.
5. Calculate abundance, richness, Hill diversity, rarefaction, dominance, origin, basal area, carbon, and beta diversity.
6. Fit the block-replicated Murillo and paired-road analyses.
7. Run road-threshold, grid-grain, coverage-standardization, height-model, and carbon-uncertainty analyses.
8. Reconcile every exhaustive framework to the same 1,256 records and tree-form carbon total.
9. Export all derived tables and render the figures.

Random procedures use the fixed seed `23052025`.

## Figures

| Figure | Files | Contents |
|---|---|---|
| 1 | `Figure1_neighborhood_inventory.*` | Study context, neighborhood boundary, blocks, and all records |
| 2 | `Figure2_top_species_distribution.*` | Distributions of the ten most abundant species |
| 3 | `Figure3_diversity_multiscale.*` | Richness, abundance, and turnover across sectors, blocks, and grid cells |
| 4 | `Figure4_Murillo_comparison.*` | Abundance and diversity above and below Murillo |
| 5 | `Figure5_origin_composition.*` | Native and exotic composition across frameworks |
| 6 | `Figure6_context_KDE.*` | Descriptive kernel-density surfaces for provisional contextual areas |
| 7 | `Figure7_road_hierarchy.*` | Road classes, associated records, frontage density, origin, and carbon |
| 8 | `Figure8_carbon_stock.*` | Species contributions and spatial distribution of standing carbon |
| S1 | `FigureS1_coverage_standardized_diversity.*` | Coverage-standardized `q0`, `q1`, and `q2` for the three paired frameworks |

Each wildcard represents both `.pdf` and `.png`.

## Derived result files

| File or group | Rows | Contents |
|---|---:|---|
| `framework_metrics.csv` | 364 | Combined metrics for every framework and group |
| `metrics_whole.csv` | 1 | Whole-neighborhood metrics |
| `metrics_murillo.csv` | 2 | Pooled Murillo-side metrics and frontage-standardized values |
| `metrics_roads.csv` | 2 | Pooled road-class metrics for records within 25 m |
| `metrics_contexts.csv` | 3 | Residential, tourist-cultural, and other-context metrics |
| `metrics_sectors.csv` | 17 | Sector-level metrics |
| `metrics_blocks.csv` | 81 | Block-level metrics, including zero-record blocks |
| `metrics_grid.csv` | 258 | Clipped grid-cell metrics, including empty cells |
| `tree_spatial_assignments_corrected.csv` | 1,256 | Original coordinates and labels plus every analysis assignment |
| `block_corrections_applied.csv` | 1,256 | Correction audit; 25 rows are flagged as corrected |
| `beta_diversity.csv` | 6 | Jaccard, turnover, nestedness, and Bray–Curtis summaries |
| `coverage_standardized_diversity.csv` | 18 | Coverage targets, standardized sample sizes, Hill numbers, and rarefaction intervals |
| `murillo_block_metrics.csv` | 81 | Block counts, frontage, origin, and carbon used for the Murillo model |
| `murillo_negative_binomial.csv` | 1 | Murillo rate ratio, HC3 interval, dispersion, likelihood, and convergence |
| `road_frontage_exposure.csv` | 119 | Block-by-road-class frontage exposure |
| `road_frontage_unit_metrics.csv` | 119 | Block-by-road-class abundance, diversity, origin, and carbon density |
| `road_paired_summary.csv` | 4 | Within-block paired summaries for abundance, diversity, origin, and carbon |
| `road_paired_inference.csv` | 1 | Clustered bootstrap interval and exact sign test |
| `road_threshold_sensitivity.csv` | 3 | Results for 15, 25, and 40 m road thresholds |
| `grid_resolution_sensitivity.csv` | 3 | Results for 15 × 15, 20 × 20, and 25 × 25 grids |
| `carbon_individual_estimates.csv` | 1,256 | Record-level height source, wood density, biomass, carbon, and CO₂e fields |
| `carbon_scenarios.csv` | 10 | Primary, influence, palm, and combined carbon scenarios |
| `carbon_uncertainty.csv` | 3 | Monte Carlo summaries for trees, palms, and an illustrative combined total |
| `height_model_parameters.csv` | 7 | Fitted local height-model parameters and residual diagnostics |
| `height_model_cross_validation.csv` | 2 | Ten-fold diagnostics for the global and species-adjusted height models |
| `validation_checks.csv` | 36 | Observed values, expected values, tolerances, and pass/fail status |
| `thematic_analysis_report.md` | — | Human-readable summary generated from the analytical outputs |

## Verified release checkpoints

These values provide a compact check that the intended data and definitions were used. The CSV files remain the authoritative machine-readable outputs.

| Checkpoint | Release value |
|---|---:|
| Unique records | 1,256 |
| Species | 81 |
| Ten most abundant species | 903 records, 71.9% |
| Native / exotic / unknown origin | 474 / 780 / 2 |
| Approximate basal area | 93.14 m² |
| Murillo records, above / below | 772 / 484 |
| Murillo block rate ratio, above versus below | 2.02; HC3 95% CI 1.30–3.14; `p = 0.0018` |
| Road records, higher-order / local-secondary / beyond 25 m | 373 / 722 / 161 |
| Aggregate road density | 59.23 versus 33.11 records km⁻¹; ratio 1.79 |
| Paired-road median difference | 19.87 records km⁻¹; sector-clustered 95% CI 11.06–32.84 |
| Occupied sectors / blocks / grid cells | 17 / 72 / 188 |
| Mean Jaccard dissimilarity, sectors / blocks / grid | 0.717 / 0.837 / 0.890 |
| Residential / tourist-cultural records | 148 / 147 |
| Primary tree-form standing carbon | 117.03 Mg C, or 429.12 Mg CO₂e |
| Conditional tree-form 95% carbon interval | 104.32–129.23 Mg C |
| Validation checks | 36 of 36 passed |

## Interpretation boundaries

- The census covers publicly accessible and publicly visible vegetation, not all vegetation on private property.
- A zero in a block or grid cell means no census record was assigned there; it does not independently prove that no unrecorded plant existed.
- Spatial frameworks overlap by design. The same record belongs to one Murillo side, one road group, one context, one sector, one block, and one grid cell.
- The Murillo rate ratio and paired-road interval describe heterogeneity among mapped blocks. They do not establish a causal effect of the avenue or road class.
- Road hierarchy is a proxy for infrastructure and possible management setting. The inventory does not identify the caretaker or planting authority for each record.
- Aggregate native fractions along road classes are heterogeneous among paired blocks and should not be generalized as a uniform road effect.
- The residential and tourist-cultural buffers are provisional working definitions, not independently mapped land-use polygons. Their comparisons and KDE panels are descriptive.
- Coverage-standardized intervals reflect conditional rarefaction variation, not field-sampling uncertainty and not estimates of unseen species.
- Carbon values are standing above-ground stock. They exclude roots, soil carbon, annual flux, avoided emissions, and future sequestration.
- Forest allometries may not fully represent open-grown, pruned, or otherwise managed urban trees.
- The palm estimate is a sensitivity result and is not part of the primary 117.03 Mg C headline.

## Reusing or extending the workflow

- Treat files in `data/` as source inputs and files in `results/` and `figures/` as derived outputs.
- Do not hand-edit a derived CSV. Change the documented source or code and rerun the public entry point.
- Preserve `ID`, recorded cadastral fields, and the correction audit when adding new spatial assignments.
- Use a new `--output-dir` when testing code or alternative definitions.
- If changing the road threshold, context radius, grid grain, allometric assumptions, or random seed, report the new definition and do not present the output as the v4 release result.
- When reporting a reproducibility issue, include the command, operating system, Python version, installed package versions, and `validation_checks.csv`.

## Key method references

- Baselga, A. (2010). Partitioning the turnover and nestedness components of beta diversity. *Global Ecology and Biogeography*, 19, 134–143. <https://doi.org/10.1111/j.1466-8238.2009.00490.x>
- Chao, A., and Jost, L. (2012). Coverage-based rarefaction and extrapolation: standardizing samples by completeness rather than size. *Ecology*, 93, 2533–2547. <https://doi.org/10.1890/11-1952.1>
- Chao, A., et al. (2014). Rarefaction and extrapolation with Hill numbers. *Ecological Monographs*, 84, 45–67. <https://doi.org/10.1890/13-0133.1>
- Chave, J., et al. (2014). Improved allometric models to estimate the aboveground biomass of tropical trees. *Global Change Biology*, 20, 3177–3190. <https://doi.org/10.1111/gcb.12629>
- Fischer, F. J., et al. (2026). Global Wood Density Database v.2. Zenodo. <https://doi.org/10.5281/zenodo.18262736>
- Goodman, R. C., et al. (2013). Amazon palm biomass and allometry. *Forest Ecology and Management*, 310, 994–1004. <https://doi.org/10.1016/j.foreco.2013.09.045>
- Jost, L. (2006). Entropy and diversity. *Oikos*, 113, 363–375. <https://doi.org/10.1111/j.2006.0030-1299.14714.x>
- Thomas, S. C., and Martin, A. R. (2012). Carbon content of tree tissues: a synthesis. *Forests*, 3, 332–352. <https://doi.org/10.3390/f3020332>

## Citation

Please cite the associated manuscript when its final bibliographic details are available. For prepublication reuse, cite this repository, identify release v4 or the exact commit hash, and cite Data S1 explicitly.

Associated manuscript authors: Rafik Neme, Damaris Grandas, Daniel Posada, Valeria Machacón, Leyla Falla, Jorge Moreno-Tilano, Yalmar Vargas, and María Cristina Martínez-Habibe.

## License and attribution

The source code and repository documentation are distributed under the MIT License in `LICENSE`.

The frozen road extract is derived from OpenStreetMap and remains subject to the Open Database License 1.0. Attribute it as `© OpenStreetMap contributors`. Its embedded metadata records the API source, bounding box, extraction date, and road-class definition.

The wood-density lookup records its source taxon, matching level, uncertainty fields, and source DOI (`10.5281/zenodo.18262736`). Institutional spatial layers and the biological inventory retain any terms imposed by their data owners. The MIT software license should not be interpreted as relicensing third-party or institutional data.

Before assigning a separate open-data license to Data S1 or the institutional spatial layers, confirm that all data owners have authorized that license.

## Contact

For scientific questions, data corrections, or reuse requests, contact:

**Rafik Neme**, Department of Natural Sciences, Universidad del Norte; Laboratorio de Genómica y Biodiversidad del Caribe Colombiano; rneme@uninorte.edu.co
