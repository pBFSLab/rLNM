# Robust Lesion Network Mapping （rLNM） Sensitivity Analysis

MATLAB scripts for observed lesion network mapping (LNM), spatial null lesion
generation, voxel-wise Monte Carlo testing, and Benjamini-Hochberg false
discovery rate (FDR) correction. This directory contains a 43-lesion Amnesia
example on the MNI-8 mm resolution.

## Project Structure

```text
C1_sensitivity_test/
|-- README.md
|-- code_algorithm_latest/
|   |-- S1_true_lesion_r_map_generate.m
|   |-- S2_true_lesion_LNM_generate.m
|   |-- S3_generate_null_model_annealing.m
|   |-- S3_generate_null_model_reposition.m
|   |-- S4_null_model_r_map_generate.m
|   |-- S5_null_model_LNM_generate.m
|   `-- S6_permutation_test_FDR.m
|-- disease_data/
|   |-- Amnesia_mask_voxel_indices_v2.xlsx
|   `-- Amnesia_lesion_mask/
|-- template/
|   |-- rho_z_mean.mat
|   |-- valid_indices_standard_wh_MNI8mm_32_32_32.mat
|   |-- norm_8mm.nii.gz
|   `-- norm.nii.gz
`-- output/
```

All paths are relative to the script location. Keep the scripts in a direct
child directory of the project root. S1-S6 do not read the raw lesion NIfTI 
files or the 1 mm `norm.nii.gz` template.

## Requirements and Inputs

- MATLAB R2019a or later with `readcell` and local functions in scripts.
- Statistics and Machine Learning Toolbox for `prctile` and `kstest2`.
- Bioinformatics Toolbox for `mafdr` in S6.
- FreeSurfer MATLAB utilities, including `MRIread`, for S3.
- Memory for a 9,473-by-9,473 double FC matrix, approximately 685 MiB, plus
  temporary arrays.
- Storage for 10,000 null MAT datasets and the downstream maps.

The input spreadsheet must have a header and one lesion per subsequent row:

| Column | Content |
|---|---|
| 1 | Unique subject/lesion ID suitable for filenames |
| 2 | Number of lesion voxels (`NumInd`) |
| 3 onward | MATLAB 1-based whole-volume linear voxel indices |

Indices must follow the voxel ordering and spatial orientation of
`MRIread('norm_8mm.nii.gz').vol`, which has dimensions 32-by-32-by-32.
Matching dimensions alone does not establish anatomical alignment.

### Current FC Matrix and Spatial Domains

The `rho_z_mean.mat` contains the variable `rho_z`, a
**9,473-by-9,473 double matrix**. The corresponding template `valid_indices`
contains **9,473 unique whole-volume indices**. Matrix rows and columns must
follow this index ordering.


## Execution

```matlab
S1_true_lesion_r_map_generate
S2_true_lesion_LNM_generate
S3_generate_null_model_annealing  % Or S3_generate_null_model_reposition
S4_null_model_r_map_generate
S5_null_model_LNM_generate
S6_permutation_test_FDR
```

### S1: Observed Lesion Connectivity Maps

S1 maps lesion indices to FC matrix rows and averages the selected rows using
`mean(..., 'omitnan')`. 

Output: `output/R0_true_lesion_r_maps_from_meanFC/Amnesia_rmap_meanFC_<ID>.mat`,
containing the 9,473-value vector `lesion_r_map`.

### S2: Observed Signed LNM

For each lesion map, values in the top 25% (at or above the 75th percentile)
receive label +1; values in the bottom 25% (at or below the 25th percentile)
receive label -1. The cutoffs are calculated separately for each map
from percentile ranks.

For the functional-connectivity matrix used in this analysis, selecting the
top and bottom 25% of each map approximately corresponds to an absolute
t-statistic threshold of `|t| > 7`.

Let C+ and C- be the numbers of positive and negative labels at a voxel, and L
the number of loaded lesions. The saved statistic is:

```text
overlap =  C+ / L    if C+ > C-
overlap = -C- / L    if C- > C+
overlap =  0         if C+ = C-
```

The current output is `output/R1_true_lesion_rLNM_map_results/Amnesia_LNM_rmap_25.
mat`, containing `overlap`; update the filename when synchronizing thresholds.

### S3: Spatial Null Lesion Datasets

Both methods preserve lesion voxel count, shape, and orientation through integer
translation. Every synthetic lesion must fit within the brain mask. Different 
lesions may overlap.

**Annealing:** up to 1,000 placement attempts per lesion and 1,000 optimization
steps per candidate. The loss is the RMSE between sorted observed and synthetic
lesion coverage proportions across all placement-mask voxels. The current 
temperature is `T0 = 1e-3`, with cooling coefficient `ALPHA = 0.997`.
The optimizer stops early when `RMSE <= TARGET_RMSE`, currently 0.005, or reaches
the iteration limit. 

A two-sample Kolmogorov-Smirnov test subsequently compares observed and synthetic
coverage frequency distributions. Candidates are retained when p > 0.05.
This constrains the overall distribution, allowing its spatial locations to
vary. It does not require voxel-wise agreement or preserve individual lesion
pairs' overlaps. The KS test is an acceptance check; RMSE is the optimization
objective. The RMSE target is an early-stopping setting, not a hard acceptance
bound. Candidates reaching the step limit still undergo KS testing.

**Repositioning:** up to 10,000 placement attempts per lesion, without frequency
matching or KS filtering. A placement failure stops the script.

Only MAT files are generated:

```text
output/R2_null_model_random_indices_iter10000_sortedfreq/null_<NNNNN>.mat
output/R2_null_model_random_indices_iter10000_reposition/null_<NNNNN>.mat
```

Both contain `null_data` with fields `Subid`, `NumInd`, and `lin_idx`. Annealing
files also contain `null_ratio`, placement-mask `valid_indices`,
`true_sorted_valid`, `final_loss`, `p_ks`, `h_ks`, `final_acceptance_rate`,
`initial_loss`, `optimization_steps`, and `run_config`. The latter records the
generation settings and initial RNG state. Both methods initialize with
`rng('shuffle', 'twister')`; the repositioning script does not save its RNG state.

### S4 and S5: Null Connectivity and LNM Maps

S4 applies S1's averaging procedure to each synthetic lesion and saves
`lesion_r_map` MAT files in an R3 subfolder for each null. S5 applies S2's
percentile labeling and signed overlap calculation to those files.

| Stage | Annealing directory under `output/` | Repositioning directory under `output/` |
|---|---|---|
| S4 | `R3_null10000_generate_r_map_sortedfreq/` | `R3_null10000_generate_r_map_reposition/` |
| S5, current code | `R4_LNM_rmaps_null_10000_results_percent25_sortedfreq/` | `R4_LNM_rmaps_null_10000_results_percent30_reposition/` |

Each S5 file is named `Amnesia_LNM_rmap_null<NNNNN>.mat` and contains `overlap`.
S4/S5 support `START_NULL` and skip existing output files by default. Only resume
with skipping enabled when upstream files, thresholds, and voxel order have
not changed.

### S6: Voxel-wise Null Testing and FDR

For nonnegative observed overlap, equally or more extreme null values satisfy
`null >= observed`. For negative overlap they satisfy `null <= observed`.
Ties count in the extreme tail. With the required complete null set:

```text
p = (number of equally or more extreme null values + 1) / (N_NULL + 1)
```

S6 applies `mafdr(pval, 'BHFDR', true)` across the FC target voxels and retains
nonzero observed values with `qval < 0.005`. This is a direction selected by the
observed sign, not an absolute-value two-sided test.

The current S6 directories are `R5_Permutation_FDR_results_percent25_sortedfreq/`
and `R5_Permutation_FDR_results_percent25_reposition/` under `output/`. For 10,000
nulls, the filename is `Amnesia_true_vs_null_10000_percent25_ratio_pval_FDR.mat`.

| Saved variable | Content |
|---|---|
| `ratio` | Fraction of null values strictly less extreme in the selected direction |
| `pval` | Monte Carlo p-values with the plus-one correction |
| `qval` | Benjamini-Hochberg adjusted p-values |
| `true_overlap_masked` | Observed signed overlap at retained voxels; zero elsewhere |
| `significant_mask_nonzero` | Nonzero observed voxels passing the FDR threshold |
| `Q_THRESHOLD` | Applied FDR threshold |
| `NULL_MODEL_METHOD` | Selected null-model method |

The plus-one correction means `pval` is not exactly `1 - ratio`. The current save
statement omits the unmasked observed map, voxel indices, null count, and file
list. Keep the source observed map and template index file for reference. 

