# Thematically organized results

## Decisions applied

- Above Murillo is the west/northwest side of the mapped avenue, shown on the left of the standard map orientation.
- Below Murillo is the east/southeast side.
- The 25 approved block corrections are applied to analysis fields. Original block labels remain in the record-level export.
- Road association uses the approved 25 m threshold.
- Residential and tourist-cultural areas use the approved provisional 150 m buffers.
- The regular clipped 20 x 20 grid is retained.
- No easting/northing trend model is used in this analysis.

## Whole neighborhood

The census contains 1,256 records representing 81 species. Hill diversity is q1 = 24.07 and q2 = 13.24. The five most abundant species account for 55.3% of records. Among records of known origin, 37.8% are native.

## Avenida Murillo

Above Murillo contains 772 records and below Murillo contains 484. Frontage-standardized abundance is 53.9 and 28.4 records per kilometre, respectively. Their pooled species lists have Jaccard dissimilarity 0.531 and Bray-Curtis dissimilarity 0.387.

Across all 81 mapped blocks, the negative-binomial model with log block perimeter as an offset estimated an above-versus-below rate ratio of 2.02 (HC3 95% CI 1.30-3.14). This is an association among blocks, not a causal effect of the avenue.

## Road hierarchy

At 25 m, 373 records are associated with higher-order roads and 722 with local or secondary streets. Aggregate abundance is 59.2 versus 33.1 records per frontage kilometre. Native fractions among known-origin records are 0.421 and 0.289. In the 47 blocks with both frontage types, higher-order density exceeds local-road density in 35 blocks, local-road density is greater in 12, and 0 are tied. The median paired density difference is 19.9 records per kilometre (sector-clustered bootstrap 95% CI 11.1-32.8).

Eight of the 25 corrected records change from local/secondary to higher-order frontage when road class is recalculated from the corrected cadastral block. Their distance-based 25 m road-association status does not change.

## Contextual reference areas

The residential area contains 148 records and the tourist-cultural areas contain 147. Their q1 diversities are 15.61 and 8.50, and their known-origin native fractions are 0.378 and 0.592. The KDE panels are descriptive density surfaces within the approved buffers; they do not establish a land-use effect.

## Spatial grains and beta diversity

All 17 sectors contain records. Of 81 mapped blocks, 72 contain records; of 258 clipped grid cells, 188 contain records. Mean Jaccard dissimilarity is 0.717 among sectors, 0.837 among occupied blocks, and 0.890 among occupied grid cells. Turnover is reported separately from nestedness in `beta_diversity.csv`.

## Standing carbon stock

The primary tree-form estimate is 116.9 Mg C, with a conditional 95% uncertainty interval of 104.3-129.2 Mg C. This equals approximately 429 Mg CO2e. It is standing above-ground stock, not annual capture or sequestration. Palm estimates remain a separate sensitivity scenario because the available measurements do not fully satisfy the palm equation's calibration conditions.

## Interpretation limits

- No observations are simulated. Monte Carlo draws are used only to propagate stated carbon-estimation uncertainty.
- Empty blocks and cells mean that the census contains no mapped records there; they are not independent confirmation that no unrecorded tree existed.
- Road hierarchy is a proxy for infrastructure and management context. The database does not identify the caretaker of each tree.
- Context buffers are approved working definitions, not independently mapped land-use polygons.
- Formal inference is limited to the block-replicated Murillo and paired-road density contrasts. Other framework comparisons remain descriptive. The same record can legitimately belong to one Murillo side, one road class, one context, one sector, one block, and one grid cell.
