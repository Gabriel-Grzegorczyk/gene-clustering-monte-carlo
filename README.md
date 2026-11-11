# Gene proximity and genomic clustering randomization

A standalone Python workflow for testing whether curated gene groups are physically closer, or contain a larger fraction of clustered genes, than scaffold-stratified random samples of annotated protein-coding loci.

The workflow performs **two one-sided Monte Carlo randomization tests using the same sampled gene sets**:

1. **Proximity:** is the median nearest-neighbor distance smaller than expected?
2. **Clustering:** is the fraction of genes in genomic clusters larger than expected at a specified threshold (35,000 bp by default)?

No BLAST search, sequence alignment, phylogenetic inference, or new homology assignment is performed. Gene classifications are supplied by the user.

## What is included

- `gene_clustering.py`: the complete analysis, in one executable script.
- `requirements.txt`: minimum package requirements.
- `requirements-tested.txt`: exact package versions used for the supplied validation.
- `examples/genes.csv`: a small normalized input table with **synthetic** gene groups.
- `examples/a_mellifera_proteome_genomic.gff`: a **synthetic** GFF, not an Apis annotation; its filename demonstrates compatibility with the study convention.
- `tests/test_analysis.py`: numerical and end-to-end checks.
- `column_dictionary.csv`: machine-readable definitions of every summary column.

The study workbook and real genome annotations are not included. The example must not be used as biological evidence.

## Installation

Use Python 3.9 or newer. From the repository directory:

```bash
python3 -m venv .venv
source .venv/bin/activate
python3 -m pip install -r requirements.txt
```

On Windows PowerShell, activate with `.venv\Scripts\Activate.ps1`. The analysis itself is platform-independent. On WSL, use the Linux commands above.

Validation was performed with Python 3.12.14, NumPy 2.3.5, and openpyxl 3.1.5. For a matching package environment, install `requirements-tested.txt`. Other compatible versions are not guaranteed to generate byte-identical outputs. Python and dependency versions are recorded in each run.

## Quick start: synthetic example

```bash
python3 gene_clustering.py \
  --genes examples/genes.csv \
  --gff-dir examples \
  --iterations 99 \
  --output results/demo \
  --xlsx
```

The example has a close group, a dispersed group, and a scaffold-singleton group. The output directory must not already exist. Use a new output path when repeating the command, or omit `--output` to create a timestamped directory.

For a study run, use 9,999 draws (the default), not the 99-draw demonstration.

## Inputs

Choose **one** input mode: `--genes` or `--workbook`. Both modes use the same GFF reader, locus matching, calculations, and randomization.

### Mode A: normalized CSV (recommended for other datasets)

A comma-delimited UTF-8 CSV with these columns:

| Column | Required | Meaning |
|---|---|---|
| `species` | Yes | Species key used in the GFF filename, e.g. `a_mellifera`. Custom keys are allowed. |
| `superfamily` | Yes | Superfamily label; custom values are allowed. |
| `level` | Yes | `family`, `subfamily`, `clade`, or `class`. |
| `group` | Yes | Curated classification label. Empty labels are reported and excluded. |
| `gene` | Yes | Gene-level identifier or alias, not an individual protein isoform. |
| `scaffold` | Yes | Exact scaffold accession from the GFF. |
| `start` | Yes | Positive integer, 1-based inclusive start. |
| `end` | Yes | Positive integer, inclusive end, not smaller than start. |
| `alias` | No | One additional gene alias to assist matching. |

Each row assigns one locus to one group at one classification level. Repeat a gene across levels if you wish to analyze its subfamily and its family separately; coordinates must agree. Do not duplicate a gene within a level or assign multiple groups to the same gene at the same level. Different gene labels resolving to the same GFF locus within a group are rejected.

Use one locus per gene, not one row per transcript or protein isoform. Existing cluster IDs are not required and are not read.

```bash
python3 gene_clustering.py \
  --genes curated_genes.csv \
  --gff-dir annotations \
  --iterations 9999 \
  --output results/full_analysis \
  --xlsx
```

### Mode B: study workbook input

This input mode accepts the study workbook schema documented below. Input headers must match exactly; output headers are English.

Every selected sheet needs `Species_Code`, `Gene_Name`, `Scaffold`, `Start`, and `End`. `Alias_GFF` is optional. Required classification columns are:

| Superfamily | Sheet | Output level | Input classification column |
|---|---|---|---|
| GST | `S2-GSTs` | family | `familia` |
| P450 | `S3-P450` | subfamily | `Sub_familia` |
| P450 | `S3-P450` | family | `Familia` |
| P450 | `S3-P450` | class | `Classe` |
| CCE | `S5-CCEs` | clade | `clado` |
| CCE | `S5-CCEs` | class | `Classe` |
| UGT | `S6-UGTs` | family | `Familias` |
| ABC | `S7-ABCs` | subfamily | `Sub_familia` |

Species mappings:

| Code | File key |
|---|---|
| Ag | a_gambiae |
| Am | a_mellifera |
| Bt | b_terrestris |
| Dm | d_melanogaster |
| Fv | f_varia |
| Mb | m_bicolor |
| Mq | m_quadrifasciata |
| Mr | m_rotundata |
| Ta | t_angustula |

The workbook must contain each selected sheet. Use plain values; formula cells are read using cached values, and Python does not recalculate Excel formulas.

Example for a study workbook (run from the directory containing this script):

```bash
python3 gene_clustering.py \
  --workbook Dados_Clusters.xlsx \
  --gff-dir annotations \
  --iterations 9999 \
  --output results/full_analysis \
  --xlsx
```

Single-group pilot:

```bash
python3 gene_clustering.py \
  --workbook Dados_Clusters.xlsx \
  --gff-dir annotations \
  --species a_mellifera \
  --superfamily P450 \
  --level subfamily \
  --group CYP6AS \
  --iterations 999 \
  --output results/CYP6AS_pilot
```

### GFF requirements and eligible universe

For each selected species, the file must be named `SPECIES_proteome_genomic.gff` inside `--gff-dir`. Files must be uncompressed GFF3 with nine tab-separated columns.

The parser expects **RefSeq-style attributes**, not arbitrary GFF dialects:

- Only `gene` records with `gene_biotype=protein_coding` are eligible.
- Records marked `pseudo=true` are excluded.
- Organelle scaffolds are identified from `region` records whose `genome` attribute contains `mitochond`, `chloroplast`, or `plastid`; their genes are excluded after parsing.
- No organelle exclusion is inferred from scaffold names. If the source lacks these `region` attributes, check and preprocess the annotation before interpreting it as a nuclear universe.
- Different overlapping gene IDs remain separate loci.
- Identical-coordinate records with the same gene ID are collapsed.
- Repeated IDs with different coordinates are combined into an outer interval only when both records are explicitly marked with `part` and a `trans-splicing` exception, on the same scaffold. Other conflicts stop the analysis.
- The outer interval of a multipart gene can include sequence between its parts; it is an operational locus envelope, not a reconstructed transcript.
- Reading stops at `##FASTA` if a sequence section is embedded in the GFF.

Annotations using `biotype` instead of `gene_biotype`, only transcript/CDS features, or protein identifiers as the sole gene labels need conversion first. Do not silently weaken these filters to obtain more matches.

### Locus matching and audits

Matching uses `ID`, the ID without a leading `gene-`, and the `Name`, `gene`, and `locus_tag` attributes, together with exact scaffold/start/end coordinates. If aliases resolve to candidates, the accepted locus must also match coordinates. If no alias matches, one unique exact-coordinate match is accepted and flagged as `matched_coordinates_only`. Ambiguous or unresolved matches stop the run before randomization and are written to `gene_matching_audit.csv`.

Missing classification labels (including `-`, `NA`, `N/A`, `nan`, and `none`) are recorded in `missing_classifications.csv`; these strings are reserved missing markers, not group names.

## Methods

### Unit of analysis

Each summary row is a species × superfamily × level × group combination. Family, subfamily, clade, and class analyses are separate. A locus may occur at multiple levels, so rows are not necessarily independent. The script does not infer or validate the classification hierarchy.

### Test 1: nearest-neighbor proximity

For two gene intervals `[s_i, e_i]` and `[s_j, e_j]` on the same scaffold:

```text
d(i,j) = max(0, max(s_i,s_j) - min(e_i,e_j) - 1)
d_i    = min d(i,j), over other same-group loci on that scaffold
```

Distances count intervening bases with 1-based inclusive coordinates. Overlap and immediate adjacency both yield zero. All possible partners are examined, including nested intervals. A scaffold singleton has an undefined distance and is excluded from distance statistics, not assigned zero.

One distance is retained per evaluable gene. Reciprocal nearest-neighbor pairs therefore contribute twice. This is intentional and is repeated in the null model; individual gene distances are not treated as independent replicates in a t-test.

The primary statistic is the median of these distances, with **no distance threshold applied to select genes**. The threshold is used only for an additional descriptive nearest-neighbor fraction and for the clustering test.

### Test 2: clustering at a specified threshold

Within each group and scaffold, order genes by start and then end. A gene continues the preceding run when:

```text
current_start - previous_gene_end <= threshold
```

Only runs of two or more genes count as clusters. Other-classification genes do not interrupt a run because they are outside the selected group. Strand is not a clustering criterion. Clusters can span more than the threshold through successive links.

The statistic is `G/N`, where `G` is the number of genes in clusters and `N` is **all** genes of the group, including scaffold singletons. Cluster number is recorded descriptively but does not receive a separate P-value.

**Distance conventions:** clustering uses `start - previous_end` without subtracting one, whereas nearest-neighbor distances use `start - end - 1`, bounded at zero. At the boundary, a nearest-neighbor gap of 35,000 intervening bases corresponds to a clustering difference of 35,001 and is not linked. The clustering algorithm compares the **immediate predecessor**, not the maximum end of all earlier genes. For nested genes, it may split a run that an interval-connectivity method would merge. This behavior is explicitly tested; changing it constitutes a methodological change.

Consequently, the nearest-neighbor threshold fraction and the clustered fraction are not interchangeable: they have different denominators and can differ in interval handling.

### Null model and randomization

For each group, draw uniformly without replacement from the eligible loci **within each scaffold**, preserving the observed number of group members on every scaffold. Original target loci remain eligible. Draws are independent, so the same locus or subset can recur across draws. Both test statistics are calculated on the **same** sampled sets, in a single pass.

Scaffolds with only one group member have no possible within-group pair. Their one-gene contribution is retained in the total `N`, but no random draw is consumed for them: their contribution to both statistics is fixed. This fixed handling is part of the reproducible sampling procedure.

The model conditions on the observed scaffold distribution. It does not test whether genes are unusually concentrated on particular scaffolds, nor match target-specific local density, gene length, annotation confidence, or chromosomal neighborhoods. Real annotated positions are retained in the eligible universe.

### Monte Carlo P-values

For `B` randomizations:

```text
P_proximity  = (1 + count(null_median <= observed_median)) / (B + 1)
P_clustering = (1 + count(null_clustered_fraction >= observed_clustered_fraction)) / (B + 1)
```

Ties are included. Clustering comparisons use integer clustered-gene counts internally, equivalent to comparing fractions because `N` is fixed. The plus-one correction avoids zero Monte Carlo P-values (Phipson & Smyth, 2010). With 9,999 draws, the smallest attainable P-value is 0.0001; it should not be reported as zero or `<0.0001`.

If the observed clustered count is zero, every draw is at least as large and `P_clustering = 1`. This is expected, not a software error. A P-value of one does not prove the null model true.

Groups with no within-scaffold pairs are reported but not tested. Both inferential P-values and distance statistics are blank for these groups; observed cluster counts/fraction are zero. Constant null distributions are flagged but still tested when observed distances exist.

### Multiple comparisons

Benjamini–Hochberg adjustment (Benjamini & Hochberg, 1995) is applied across **all tested combinations in the current invocation**, separately for the proximity and clustering P-values. Values named `q_*_bh` are BH-adjusted P-values, not independently estimated Storey q-values. The reporting threshold used in the study is adjusted P < 0.05.

Filtering a run changes the family of tests and therefore can change adjusted P-values. Compute the full intended testing family before filtering figures or tables; do not adjust again only among significant or displayed rows. Hierarchical groups overlap; do not describe these tests as independent or claim an unconditional FDR guarantee under arbitrary dependence.

### Null summaries and effect sizes

Each distance statistic is computed within each random draw and then summarized across draws by its mean and 2.5th, 50th, and 97.5th percentiles. Thus `null_median_bp_p50` is the **median of simulated medians**, not the median of a pooled list of all sampled distances or the genome-wide intergenic median. Percentiles use NumPy's default linear interpolation.

The interval from P2.5 to P97.5 describes the central 95% of the null statistic distribution; it is not a confidence interval for the observed statistic. The observed-to-null median ratio summarizes proximity effect size. An optional clustering effect size in percentage points is:

```text
100 * (observed_clustered_fraction - null_clustered_fraction_mean)
```

This effect-size expression is not used to compute the Monte Carlo P-value.

## Outputs

All CSVs are UTF-8 with a BOM, comma separators, decimal points, and CSV quoting. Numeric missing values are empty fields. Proportions/fractions are **0–1, not percentages**. The optional `.xlsx` stores numeric cells as numbers and is preferable to opening CSVs by double-clicking in a locale with decimal commas. Do not replace decimal separators or delete punctuation in raw CSVs manually.

| File | Contents |
|---|---|
| `summary.csv` | One row per species/superfamily/level/group, including untested groups. |
| `summary.xlsx` | Optional typed copy of the full summary (`--xlsx`), sheet `Results`; fractions remain 0–1. |
| `simulation_statistics.csv` | One row per tested group per draw, containing both analyses. Suppress with `--no-simulations`. |
| `observed_genes.csv` | Per-group gene coordinates, nearest neighbor, and new cluster membership. |
| `observed_clusters.csv` | Catalogue of observed clusters, members, and spans. |
| `gene_matching_audit.csv` | Input-to-GFF matching evidence and failures. |
| `missing_classifications.csv` | Missing group labels excluded from the selected analysis. |
| `parameters.json` | Input hashes, parameters, versions, GFF audits, timing, and completion flag. |

Only use a completed run: `parameters.json` must contain `"completed": true`. A failure may leave partial output for diagnosis. New output directories are mandatory, preventing overwrites.

### Complete summary column dictionary

The four identity columns form the join key. All distances have units of bp; `*_fraction*` values are dimensionless on the 0–1 scale. Blank is not zero. Every summary column is listed below and in `column_dictionary.csv`.

| Column | Meaning |
|---|---|
| `species` | Species file key; joins to SPECIES_proteome_genomic.gff. |
| `superfamily` | Curated superfamily label. |
| `level` | Classification level: family, subfamily, clade, or class. |
| `group` | Curated group label at that classification level. |
| `n_genes` | N: all observed group loci, including scaffold singletons. |
| `n_scaffolds` | Number of scaffolds containing observed group loci. |
| `n_without_partner` | Number of genes with no other same-group gene on their scaffold. |
| `n_evaluable` | A = N - n_without_partner; denominator of nearest-neighbor fractions. |
| `counts_per_scaffold` | JSON object mapping scaffold IDs to observed gene counts, held fixed in every draw. |
| `threshold_bp` | L in bp (default 35000), applied to both threshold statistics with their distinct gap definitions. |
| `iterations` | B draws actually performed; zero for untested groups. |
| `status` | tested or no_pairs_on_scaffold. Untested inferential fields are blank, not zero. |
| `p_proximity` | Lower-tail Monte Carlo P: (1 + number of null medians <= observed median)/(B + 1). |
| `q_proximity_bh` | BH-adjusted p_proximity across all tested groups in this invocation. |
| `p_clustering` | Upper-tail Monte Carlo P: (1 + number of null clustered counts >= observed count)/(B + 1). |
| `q_clustering_bh` | BH-adjusted p_clustering across all tested groups in this invocation, separately from proximity. |
| `observed_n_clusters` | Observed number of consecutive same-group runs containing at least two genes. |
| `observed_n_clustered_genes` | G: number of observed genes in those runs. |
| `observed_mean_bp` | Observed arithmetic mean of evaluable nearest-neighbor distances (bp). |
| `observed_q1_bp` | Observed 25th percentile of evaluable nearest-neighbor distances (bp). |
| `observed_median_bp` | Observed median of evaluable nearest-neighbor distances (bp). |
| `observed_q3_bp` | Observed 75th percentile of evaluable nearest-neighbor distances (bp). |
| `observed_p90_bp` | Observed 90th percentile of evaluable nearest-neighbor distances (bp). |
| `observed_p95_bp` | Observed 95th percentile of evaluable nearest-neighbor distances (bp). |
| `observed_minimum_bp` | Observed minimum evaluable nearest-neighbor distance (bp). |
| `observed_maximum_bp` | Observed maximum evaluable nearest-neighbor distance (bp), not the maximum pairwise or consecutive gap. |
| `observed_zero_fraction` | Observed fraction of A evaluable genes with zero nearest-neighbor distance (overlap or immediate adjacency). |
| `observed_within_threshold_fraction` | Observed fraction of A evaluable genes with nearest-neighbor distance <= L. |
| `observed_clustered_fraction` | Observed fraction G/N of all group genes belonging to clusters. |
| `null_mean_bp_mean` | Arithmetic mean across B per-draw values of the arithmetic mean of evaluable nearest-neighbor distances (bp). |
| `null_mean_bp_p025` | 2.5th percentile across B per-draw values of the arithmetic mean of evaluable nearest-neighbor distances (bp). |
| `null_mean_bp_p50` | Median across B per-draw values of the arithmetic mean of evaluable nearest-neighbor distances (bp). |
| `null_mean_bp_p975` | 97.5th percentile across B per-draw values of the arithmetic mean of evaluable nearest-neighbor distances (bp). |
| `null_q1_bp_mean` | Arithmetic mean across B per-draw values of the 25th percentile of evaluable nearest-neighbor distances (bp). |
| `null_q1_bp_p025` | 2.5th percentile across B per-draw values of the 25th percentile of evaluable nearest-neighbor distances (bp). |
| `null_q1_bp_p50` | Median across B per-draw values of the 25th percentile of evaluable nearest-neighbor distances (bp). |
| `null_q1_bp_p975` | 97.5th percentile across B per-draw values of the 25th percentile of evaluable nearest-neighbor distances (bp). |
| `null_median_bp_mean` | Arithmetic mean across B per-draw values of the median of evaluable nearest-neighbor distances (bp). |
| `null_median_bp_p025` | 2.5th percentile across B per-draw values of the median of evaluable nearest-neighbor distances (bp). |
| `null_median_bp_p50` | Median across B per-draw values of the median of evaluable nearest-neighbor distances (bp). |
| `null_median_bp_p975` | 97.5th percentile across B per-draw values of the median of evaluable nearest-neighbor distances (bp). |
| `null_q3_bp_mean` | Arithmetic mean across B per-draw values of the 75th percentile of evaluable nearest-neighbor distances (bp). |
| `null_q3_bp_p025` | 2.5th percentile across B per-draw values of the 75th percentile of evaluable nearest-neighbor distances (bp). |
| `null_q3_bp_p50` | Median across B per-draw values of the 75th percentile of evaluable nearest-neighbor distances (bp). |
| `null_q3_bp_p975` | 97.5th percentile across B per-draw values of the 75th percentile of evaluable nearest-neighbor distances (bp). |
| `null_p90_bp_mean` | Arithmetic mean across B per-draw values of the 90th percentile of evaluable nearest-neighbor distances (bp). |
| `null_p90_bp_p025` | 2.5th percentile across B per-draw values of the 90th percentile of evaluable nearest-neighbor distances (bp). |
| `null_p90_bp_p50` | Median across B per-draw values of the 90th percentile of evaluable nearest-neighbor distances (bp). |
| `null_p90_bp_p975` | 97.5th percentile across B per-draw values of the 90th percentile of evaluable nearest-neighbor distances (bp). |
| `null_p95_bp_mean` | Arithmetic mean across B per-draw values of the 95th percentile of evaluable nearest-neighbor distances (bp). |
| `null_p95_bp_p025` | 2.5th percentile across B per-draw values of the 95th percentile of evaluable nearest-neighbor distances (bp). |
| `null_p95_bp_p50` | Median across B per-draw values of the 95th percentile of evaluable nearest-neighbor distances (bp). |
| `null_p95_bp_p975` | 97.5th percentile across B per-draw values of the 95th percentile of evaluable nearest-neighbor distances (bp). |
| `null_minimum_bp_mean` | Arithmetic mean across B per-draw values of the minimum evaluable nearest-neighbor distance (bp). |
| `null_minimum_bp_p025` | 2.5th percentile across B per-draw values of the minimum evaluable nearest-neighbor distance (bp). |
| `null_minimum_bp_p50` | Median across B per-draw values of the minimum evaluable nearest-neighbor distance (bp). |
| `null_minimum_bp_p975` | 97.5th percentile across B per-draw values of the minimum evaluable nearest-neighbor distance (bp). |
| `null_maximum_bp_mean` | Arithmetic mean across B per-draw values of the maximum evaluable nearest-neighbor distance (bp), not the maximum pairwise or consecutive gap. |
| `null_maximum_bp_p025` | 2.5th percentile across B per-draw values of the maximum evaluable nearest-neighbor distance (bp), not the maximum pairwise or consecutive gap. |
| `null_maximum_bp_p50` | Median across B per-draw values of the maximum evaluable nearest-neighbor distance (bp), not the maximum pairwise or consecutive gap. |
| `null_maximum_bp_p975` | 97.5th percentile across B per-draw values of the maximum evaluable nearest-neighbor distance (bp), not the maximum pairwise or consecutive gap. |
| `null_zero_fraction_mean` | Arithmetic mean across B per-draw values of the fraction of A evaluable genes with zero nearest-neighbor distance (overlap or immediate adjacency). |
| `null_zero_fraction_p025` | 2.5th percentile across B per-draw values of the fraction of A evaluable genes with zero nearest-neighbor distance (overlap or immediate adjacency). |
| `null_zero_fraction_p50` | Median across B per-draw values of the fraction of A evaluable genes with zero nearest-neighbor distance (overlap or immediate adjacency). |
| `null_zero_fraction_p975` | 97.5th percentile across B per-draw values of the fraction of A evaluable genes with zero nearest-neighbor distance (overlap or immediate adjacency). |
| `null_within_threshold_fraction_mean` | Arithmetic mean across B per-draw values of the fraction of A evaluable genes with nearest-neighbor distance <= L. |
| `null_within_threshold_fraction_p025` | 2.5th percentile across B per-draw values of the fraction of A evaluable genes with nearest-neighbor distance <= L. |
| `null_within_threshold_fraction_p50` | Median across B per-draw values of the fraction of A evaluable genes with nearest-neighbor distance <= L. |
| `null_within_threshold_fraction_p975` | 97.5th percentile across B per-draw values of the fraction of A evaluable genes with nearest-neighbor distance <= L. |
| `null_clustered_fraction_mean` | Arithmetic mean across B per-draw values of the fraction G/N of all group genes belonging to clusters. |
| `null_clustered_fraction_p025` | 2.5th percentile across B per-draw values of the fraction G/N of all group genes belonging to clusters. |
| `null_clustered_fraction_p50` | Median across B per-draw values of the fraction G/N of all group genes belonging to clusters. |
| `null_clustered_fraction_p975` | 97.5th percentile across B per-draw values of the fraction G/N of all group genes belonging to clusters. |
| `observed_to_null_median_ratio` | Observed nearest-neighbor median / median of null medians; blank when the null denominator is zero. |
| `null_median_variable` | True if the simulated nearest-neighbor medians have nonzero range. |
| `null_clustering_variable` | True if simulated clustered fractions have nonzero range. |

### Simulation table columns

`species`, `superfamily`, `level`, and `group` identify the set. `iteration` is the one-based draw number. `mean_bp`, `q1_bp`, `median_bp`, `q3_bp`, `p90_bp`, `p95_bp`, `minimum_bp`, `maximum_bp`, `zero_fraction`, and `within_threshold_fraction` are the corresponding statistics for that single draw (same definitions as the observed metrics above). `clustered_fraction` is the simulated clustered count divided by total group size. `n_clusters` is the number of simulated clusters and `n_clustered_genes` their total membership.

No per-draw gene IDs or coordinates are saved. Reconstruct them by rerunning with the same inputs, versions, group identity, and seed; the table alone is insufficient to reconstruct membership.

### Observed gene table columns

| Column | Meaning |
|---|---|
| `species`, `superfamily`, `level`, `group` | Group identity. A locus can appear again at another level. |
| `gene_id` | Matched GFF gene ID. |
| `scaffold` | Exact GFF scaffold accession. |
| `start`, `end` | Inclusive locus interval; multipart envelopes may be used. |
| `nearest_gene_id` | Closest same-group locus on the scaffold; ties use the first coordinate-sorted candidate. |
| `nearest_distance_bp` | Intervening bases to that neighbor, bounded at zero. Blank for scaffold singletons. |
| `overlap` | True when the two inclusive intervals overlap; adjacent non-overlapping intervals can have distance zero but overlap=False. Blank if no partner exists. |
| `status` | `evaluable` or `no_partner_on_scaffold`. |
| `cluster_id` | New cluster label within the group; blank for unclustered genes. |
| `previous_group_gene_gap_bp` | Start-minus-end distance to the previous member inside the same cluster. Blank for its first member or for an unclustered gene; may be negative. |

### Observed cluster catalogue columns

| Column | Meaning |
|---|---|
| `species`, `superfamily`, `level`, `group` | Group identity. |
| `cluster_id` | `cluster_001`, etc., numbered within each group across coordinate-sorted scaffolds. Only unique together with all four identity columns. |
| `scaffold` | Scaffold containing the cluster. |
| `n_genes` | Number of cluster members (at least two). |
| `start` | Minimum member start. |
| `end` | Maximum member end. |
| `span_bp` | `end - start + 1`, not the sum of gaps; may exceed the clustering threshold. |
| `maximum_link_gap_bp` | Largest immediate-predecessor start-minus-end link in the run; may be negative for overlapping genes. |
| `gene_ids` | Semicolon-separated ordered GFF IDs. |

**Cluster IDs are unique only within a species × superfamily × level × group combination.** Use all four group identifiers together with `cluster_id` when linking tables. Proximity and clustering P-values apply to whole groups; they are not cluster-specific P-values.

### Matching audit and missing-classification columns

`gene_matching_audit.csv` contains `species`, `superfamily`, `input_gene`, `scaffold`, `start`, `end`, `status`, `gene_id`, `name_candidates`, and `coordinate_candidates`. Candidate fields contain semicolon-separated GFF IDs. Status is `matched_name_and_coordinates`, `matched_coordinates_only`, or `error`. This audit is per input locus/superfamily, not repeated for each hierarchy level.

`missing_classifications.csv` contains `species`, `superfamily`, `level`, `gene`, and `input_location`. The last field identifies the workbook sheet and row, or the CSV row. It reports missing labels within the selected species/superfamily/level scope, even if a group-name filter is also supplied.

### Run metadata

`parameters.json` records:

- `script_version`, `parameters`, `input_file`, `input_sha256`, `script_sha256`;
- `python`, `numpy`, `openpyxl` versions;
- `rng_namespace`, `clustering_rule`, `nearest_rule`, and `bh_scope`;
- `completed` (false until all final outputs have been written);
- `seconds`, `n_groups`, `n_tests` after completion;
- `gffs`, keyed by species. Each entry has `file`, `sha256`, `gene_line_biotypes`, `eligible_genes`, `organelle_genes_excluded`, `scaffolds`, `detected_organelle_scaffolds`, `merged_multipart_ids`, and `duplicate_lines`.

`gene_line_biotypes` counts gene feature lines before eligibility filters/merging, while `eligible_genes` counts final loci. `duplicate_lines` counts same-ID/same-coordinate repetitions; distinct multipart additions are reported separately. All timing and hashes are run-specific.

## Reproducibility

The workflow calculates both statistics in a single pass through the same sampled gene sets. No precomputed analysis results are required.

Each group receives its own random stream, derived from the base seed and its identity. Fixed internal classification-level codes are used in the seed calculation:

```text
family -> familia
subfamily -> subfamilia
clade -> clado
class -> classe
SHA256(json.dumps([seed, species, superfamily, level_seed_code, group]))
```

These internal codes are part of the reproducibility specification; changing them changes the random draws. They do not affect the English output labels.

The first 16 digest bytes, interpreted as a big-endian integer, initialize NumPy `default_rng`. The default base seed is 20261103. Eligible pools are sorted by start, end, and gene ID, and scaffolds are processed in sorted order. Group filtering therefore preserves that group's random stream when inputs, seed, and numerical environment are unchanged; BH-adjusted P-values can still change because the testing family changes.

For reproducible results, retain the input files, their hashes, the software release or commit, command-line parameters, draw count, threshold, and dependency versions. `parameters.json` records these run details except the Git commit, which should be archived separately. Summary and simulation outputs are deterministic under a matching environment; timestamps, paths, timing, and file hashes can differ between runs.

`summary.xlsx` uses sheet `Results` and English machine-readable headers. Fractions remain on the 0–1 scale. Join tables by the four group-identity columns rather than row position. Before comparing analyses based on different curated inputs, check gene memberships, coordinates, classifications, and input hashes.

## Command-line options

| Option | Default / behavior |
|---|---|
| `--genes PATH` | Normalized CSV; mutually exclusive with `--workbook`. |
| `--workbook PATH` | Study workbook input. One input mode is required. |
| `--gff-dir PATH` | Required directory of species GFFs. |
| `--output PATH` | New directory; default `results/run_TIMESTAMP` relative to the working directory. |
| `--species TEXT` | Optional exact species file key filter. |
| `--superfamily TEXT` | Optional exact superfamily filter. |
| `--level TEXT` | Optional `family`, `subfamily`, `clade`, or `class` filter. |
| `--group TEXT` | Optional exact group label; quote labels containing spaces. |
| `--iterations N` | 9999; integer >=1. |
| `--seed N` | 20261103; integer base seed. |
| `--threshold N` | 35000 bp; integer >=0. |
| `--no-simulations` | Save disk space by omitting the per-draw CSV; tests and summary unchanged. |
| `--xlsx` | Also write a typed numeric full summary workbook. |

## Testing

From the repository directory:

```bash
python3 -m unittest discover -s tests -v
```

Tests cover nearest distances against a scalar oracle; threshold boundary and nested-interval behavior; BH adjustment; GFF organelle/pseudogene filtering and multipart conflicts; full CLI execution; singleton denominators; Monte Carlo tails reconstructed from saved draws; filtered-run random-stream reproducibility; numeric Excel cells; and refusal to overwrite an existing output directory.

Runtime is driven mainly by the number of draws and group size on each scaffold. The nearest-neighbor calculation uses an n-by-n distance matrix per scaffold, so it is intended for gene-family-scale groups rather than an all-genome nearest-neighbor matrix. Per-draw CSVs can become large. Use a 99- or 999-draw pilot to estimate runtime, then choose a fixed final draw count. A filtered pilot changes the BH testing family even when raw random streams match.

## Interpretation and limitations

- A significant proximity test indicates a small median relative to the specified random-locus model.
- A significant clustering test indicates excess clustered membership at the chosen threshold relative to that same model.
- A group may be significant for only one statistic: a minority of tight clusters can coexist with distant genes.
- These results support nonrandom organization under the model; they do not establish the optimality of 35 kb or prove tandem duplication, shared function, selection, or a specific evolutionary mechanism.
- Pairwise species differences are not tested by comparing their separate P-values.
- Sampling the entire group before clustering avoids selecting close observed subsets and comparing them to unselected random subsets. Individual observed clusters do not receive independent significance tests.
- Fragmented assemblies, missing annotation, multipart envelopes, unequal annotation quality, and scaffold-conditioned sampling affect interpretation. Distance across scaffold boundaries is unknown.
- Parent-level groups contain child-level genes and are not independent replicates. BH adjustment is applied as documented, but arbitrary dependence is not addressed by an additional correction.
- The two gap conventions are explicitly defined above. For a different interval-connectivity definition, revise the algorithm, test it, version the change, and rerun observed and null analyses consistently.

## Citation and attribution

When using this workflow, cite the software repository and exact release or commit, the associated study when available, and the statistical references below.

The 35-kb default is a study parameter, not an optimized threshold. When reporting an analysis, cite the biological source supporting your threshold choice in addition to the statistical references below.

- Phipson, B., & Smyth, G. K. (2010). Permutation P-values should never be zero: Calculating exact P-values when permutations are randomly drawn. *Statistical Applications in Genetics and Molecular Biology*, 9, Article 39. https://doi.org/10.2202/1544-6115.1585
- Benjamini, Y., & Hochberg, Y. (1995). Controlling the false discovery rate: A practical and powerful approach to multiple testing. *Journal of the Royal Statistical Society: Series B*, 57(1), 289–300. https://doi.org/10.1111/j.2517-6161.1995.tb02031.x

### Suggested concise methods text

> Within-group genomic proximity and clustering were assessed using one-sided Monte Carlo randomization tests. For each species and classification group, 9,999 gene sets were sampled without replacement from eligible protein-coding loci, preserving observed counts per scaffold. The median nearest-neighbor distance and the fraction of genes in consecutive same-group clusters at a 35-kb threshold were compared with their respective null distributions. Monte Carlo P-values used the plus-one correction (Phipson & Smyth, 2010) and were adjusted using Benjamini–Hochberg separately for each statistic across all tested groups (Benjamini & Hochberg, 1995). Genes lacking a same-group partner on their scaffold were excluded from distance statistics but retained in the clustered-fraction denominator. Full coordinate conventions and input eligibility criteria are documented in this repository.
