# Barrio Abajo urban trees

Reproducible data, analysis, figures, and manuscript for the census of publicly
accessible trees and other arborescent plants in Barrio Abajo, Barranquilla,
Colombia.

The analysis begins with the complete neighborhood and the distributions of
its ten most abundant species. It then examines richness and turnover across
sectors, cadastral blocks, and a clipped 20 × 20 grid before comparing the two
sides of Avenida Murillo, higher-order and local/secondary roads, and
provisional residential and tourist-cultural areas. It describes abundance,
dominance, Hill diversity, beta diversity, species origin, basal area, and
standing above-ground carbon stock. The former easting/northing model suite is
not part of this repository.

Inference is deliberately limited to comparisons with block-level replication.
The Murillo contrast uses a negative-binomial count model with log block
perimeter as an offset and HC3 covariance. Road-frontage density is compared
within the 47 blocks containing both road classes, with a sector-clustered
bootstrap interval. Coverage-standardized Hill diversity and threshold/grid
sensitivity analyses are reported as supporting results.

## Run everything

Python 3.8 or newer is recommended. From the repository root:

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install -r requirements.txt
python analysis/generate_all_figures.py
```

On Windows PowerShell, activate with:

```powershell
.venv\Scripts\Activate.ps1
```

The command creates `outputs/` with eight main figures and one supplementary
figure in PDF and PNG, all derived
tables, record-level spatial assignments, carbon estimates, and
`validation_checks.csv`. To use another destination:

```bash
python analysis/generate_all_figures.py --output-dir test_outputs
```

The workflow uses a fixed random seed for carbon uncertainty and should
reproduce the versioned results apart from harmless platform-level differences
in figure rendering.

## Source data and approved definitions

Every biological calculation starts from the `ArbolesIndividuales_OK` sheet in
`data/Data_S1_Barrio_Abajo_tree_inventory.xlsx`. Derived workbook sheets are
not analytical inputs, and no tree observations are simulated.

- Above Murillo is the west/northwest side, shown on the left of the standard
  map; below Murillo is the east/southeast side.
- Twenty-five implausible recorded block labels are corrected in analysis-only
  fields. Original labels are retained in the audit output.
- Road association uses a 25 m distance to a frozen mapped motor-road
  centerline. Road class follows the frontage of the corrected block.
- Higher-order roads follow Barranquilla District Decree 0317 of 2023; all
  other mapped motor roads are grouped as local/secondary.
- Provisional contextual areas are 150 m buffers around La Castellana Park,
  Old Customs Square, and the Caribbean Cultural Park.
- The retained clipped grid has 20 columns by 20 rows.

The KDE panels are descriptive. Road hierarchy is a proxy for infrastructure
and possible management setting; the census does not identify the caretaker of
each tree.

## Verified headline values

A successful run reports:

- 1,256 unique individuals and 81 species;
- 903 records (71.9%) in the ten most abundant species;
- 772 records above and 484 below Murillo;
- 373 higher-order-road records, 722 local/secondary-road records, and 161
  records beyond the 25 m road threshold;
- 148 residential, 147 tourist-cultural, and 961 other-context records;
- 17 occupied sectors, 72 occupied blocks, and 188 occupied grid cells;
- 117.03 Mg C in the primary tree-form standing-stock estimate; and
- a conditional 95% uncertainty interval of approximately 104.3–129.2 Mg C.

All reconciliation and inferential-design checks must pass before the run
completes.

## Repository contents

| Path | Purpose |
|---|---|
| `data/Data_S1_Barrio_Abajo_tree_inventory.xlsx` | Curated census and species traits; submitted as Data S1 |
| `data/context/` | Geographic context boundaries used in Figure 1 |
| `data/spatial/` | Neighborhood, sector, and cadastral-block geometries |
| `data/barrio_abajo_roads.geojson` | Frozen road centerlines and analysis hierarchy |
| `data/wood_density_lookup.csv` | Documented wood-density inputs |
| `analysis/generate_all_figures.py` | One-command public entry point |
| `analysis/thematic_analysis.py` | Spatial assignments, summaries, validation, and figures |
| `analysis/carbon_stock.py` | Height model, tree/palm allometry, and uncertainty |
| `analysis/analysis_utils.py` | Shared input, spatial, diversity, and plotting utilities |
| `results/` | Verified derived tables from the manuscript run |
| `figures/` | Verified manuscript figures in PDF and PNG |
| `manuscript/` | Manuscript v4 in Word and Markdown |
| `requirements.txt` | Python dependencies |

## Figure sequence

1. Study area and complete inventory
2. Individual distributions of the ten most abundant species
3. Species richness and turnover across spatial grains
4. Abundance above and below Avenida Murillo
5. Native and exotic composition
6. Residential and tourist-cultural KDEs
7. Higher-order and local-road assemblages
8. Standing above-ground carbon stock

Supplementary Figure S1 reports coverage-standardized Hill diversity for the
Murillo, road-class, and contextual contrasts. Supporting CSV files also report
the Murillo block model, paired-road inference, and road-threshold and grid-grain
sensitivity analyses.

## Interpretation limits

- Carbon values are standing above-ground stock, not annual capture or
  sequestration.
- Monte Carlo draws propagate stated estimation uncertainty; they do not add
  or replace observations.
- Palm carbon remains a sensitivity result because the census does not fully
  meet the palm equation's calibration conditions.
- Empty blocks and grid cells contain no census records; they are not an
  independent confirmation that no unrecorded individual was present.

Public repository: <https://github.com/genbiodiv/BarrioAbajoTrees_v3>

## License

Repository code is distributed under `LICENSE`. The frozen road extract is
derived from OpenStreetMap and remains subject to the Open Database License.
