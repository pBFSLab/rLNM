# rLNM Selective Specificity Testing

This directory implements a **symptom-specificity analysis** for lesion network mapping (LNM). It identifies the connectivity pattern that is unique to each symptom of interest (SOI) after removing connectivity shared across symptoms.

Inputs are two pre-computed sets of files:

- lesion-wise connectivity maps (r-maps) for every lesion of every symptom;
- a binary rLNM sensitivity mask for every symptom (voxels that survive a null model).

The main output is one **symptom-specific map** per SOI, together with all intermediate statistics.

> The layout follows `../C1_sensitivity_test` (rLNM sensitivity analysis): scripts in
> `code_algorithm/`, input data in `disease_data/`, templates in `template/`,
> and all results under lowercase `output/`.

---

## 1. Pipeline

The pipeline has five steps. Each step reads the `.mat` file written by the previous step, so the pipeline can be run as a whole or step by step.

| Step | What it does | Main output |
|---|---|---|
| **S1** `S1_load_inputs` | Loads all lesion r-maps and sensitivity masks, aligned to one grid | `S1_inputs.mat` |
| **S2** `S2_symptom_maps` | Averages within each symptom, restricts to positive connectivity, defines the analysis domain, builds the cross-symptom map | `S2_symptom_maps.mat`, maps/domain files |
| **S3** `S3_lesion_specific` | Builds a leave-one-lesion-out (LOLO) shared-connectivity regressor per lesion and computes the OLS residual | `S3_lesion_specific.mat`, per-lesion maps |
| **S4** `S4_group_comparison` | Each SOI versus all other lesions, voxel-wise Welch t-test and specific effect | `S4_group_comparison.mat`, t / effect maps |
| **S5** `S5_permutation_fdr` | 10,000 symptom-label permutations → empirical P → BH-FDR → retain q < 0.05 | `S5_permutation_fdr.mat`, p / q / final maps |

**Core result:**

```
output/S5_specific_map_<SOI>.nii.gz
```

Significant voxels carry the specific effect of that SOI (standardized target − contrast difference); all other voxels are zero.

---

## 2. Directory layout

```text
C2_specificity_test/
├── README.md                 Chinese documentation
├── README_EN.md              English documentation (this file)
│
├── code_algorithm/    algorithm scripts (same folder name as S1_sensitivity_test)
│   ├── run_all.m             single entry point: adds paths, runs S1–S5
│   ├── S1_load_inputs.m
│   ├── S2_symptom_maps.m
│   ├── S3_lesion_specific.m
│   ├── S4_group_comparison.m
│   ├── S5_permutation_fdr.m
│   └── utils/               methods and statistics
│       ├── specificity_config.m       shared settings
│       ├── specificity_rethreshold.m  offline re-thresholding (no new permutations)
│       ├── welch_t.m / specific_effect_score.m / bh_fdr.m / label_permutation.m
│       ├── io_nifti_read.m / io_nifti_write.m / resample_nearest.m
│       └── write_map.m / write_lesionwise_maps.m
│
├── disease_data/             input data
│   ├── Lesion_rmaps/
│   │   ├── <SOI_1>/*.nii.gz        one 3D r-map per lesion
│   │   ├── <SOI_2>/*.nii.gz
│   │   └── ...
│   └── Lesion_sensitivitymap/
│       ├── <SOI_1>/significant_8mm.nii.gz
│       ├── <SOI_2>/significant_8mm.nii.gz
│       └── ...
│
├── template/                 templates (single copy of the brain mask etc.)
│   ├── norm_8mm.nii.gz        brain mask used for the analysis domain
│   ├── norm.nii.gz
│   └── valid_indices_standard_wh_MNI8mm_32_32_32.mat
│
└── output/                   all results (created on first run)
```

---

## 3. Inputs

| Input | Requirement |
|---|---|
| `disease_data/Lesion_rmaps/<SOI>/*.nii.gz` | Same 3D grid (32×32×32 for the current data); each file is one lesion; the file stem is the lesion ID |
| `disease_data/Lesion_sensitivitymap/<SOI>/significant_8mm.nii.gz` | Binary 0/1; same grid as the r-maps, or resampled by affine |
| `template/norm_8mm.nii.gz` | Optional. Falls back to `../S1_sensitivity_test/template/norm_8mm.nii.gz` |
| Number of symptoms | Detected automatically from the subfolders of `disease_data/Lesion_rmaps/`; adding a symptom requires no code change |

`<SOI>` is the symptom name, e.g. `Face`, `Lan`, `Spatial`. By default a missing sensitivity map is an error (`REQUIRE_SENSITIVITY_MAP = true`).

---

## 4. Running

### Requirements

- **MATLAB** R2020a or newer; NIfTI I/O uses `niftiread` (Image Processing Toolbox) or SPM.

### Full pipeline

```matlab
cd('matlab_specificity/C2_specificity_test/code_algorithm');
run_all
```

`run_all` adds `utils/` to the MATLAB path.
On the current dataset (184 lesions, 5,486 voxels, 10,000 permutations) the run takes about 1–2 minutes.

### Running one step

```matlab
cd('matlab_specificity/C2_specificity_test/code_algorithm');
addpath('utils');

S3_lesion_specific       % e.g. recompute residuals only
S5_permutation_fdr
```

### Parameters

| Parameter | Location | Default |
|---|---|---|
| `REQUIRE_SENSITIVITY_MAP` | `code_algorithm_latest/S1_load_inputs.m` | `true` |
| `BRAIN_MASK_FILE` | `code_algorithm_latest/S2_symptom_maps.m` | `template/norm_8mm.nii.gz` |
| `fdr_level` | `code_algorithm_latest/utils/specificity_config.m` | `0.05` |
| `sided` | `code_algorithm_latest/utils/specificity_config.m` | `'one'` |
| `n_permutations` | `code_algorithm_latest/utils/specificity_config.m` | `10000` |
| `n_workers` | `code_algorithm_latest/utils/specificity_config.m` | `8` |

`sided` is shared by S5 and `specificity_rethreshold`, so it only needs to be set once.

---

## 5. Outputs

```
output/
├── S1_inputs/S1_inputs.mat
├── S2_symptom_maps/S2_symptom_maps.mat
├── S3_lesion_specific/S3_lesion_specific.mat
├── S4_group_comparison/S4_group_comparison.mat
├── S5_permutation_fdr/S5_permutation_fdr.mat
│
├── S1_input_summary.txt
├── S2_symptom_map_<SOI>.nii.gz               raw symptom-level map
├── S2_symptom_map_positive_<SOI>.nii.gz      positive-restricted map
├── S2_domain_mask.nii.gz                     analysis domain
├── S2_sens_intersection.nii.gz               intersection of sensitivity masks
├── S2_cross_symptom_map.nii.gz               cross-symptom map
│
├── S3_lolo/<SOI>/<lesion>.nii.gz             per-lesion LOLO regressor
├── S3_lesion_specific/<SOI>/<lesion>.nii.gz  per-lesion OLS residual
│
├── S4_tstat_<SOI>.nii.gz                     Welch t
├── S4_specific_effect_<SOI>.nii.gz           specific effect
│
├── S5_pvalue_<SOI>.nii.gz                    empirical P
├── S5_qvalue_<SOI>.nii.gz                    BH-FDR q
├── S5_specific_map_<SOI>.nii.gz              ★ core result
├── S5_specific_mask_<SOI>.nii.gz             significance mask
└── S5_summary.txt
```

All volumes share the same 3D grid (32³ for the current data). The `.mat` files keep the complete intermediate state of each step for inspection and recomputation.

### Re-thresholding without re-running permutations

```matlab
addpath('code_algorithm_latest/utils');
specificity_rethreshold
```

Recomputes significance for different FDR levels or specific-effect thresholds. Results are written to `output/S5_rethreshold/`.

---

## 6. Conventions

- **Analysis domain**: all voxels inside the brain mask (5,486 in the current data). The domain is fixed by anatomy, not selected from the data.
- **Significance**: only BH-FDR `q < fdr_level`. The specific effect is reported as an effect size and is not used as a threshold.
- **Regression**: both the lesion r-maps and the LOLO regressors remain signed, so the residual removes shared connectivity of both signs.
- **Sensitivity masks**: used only for the cross-symptom map in S2; they do not enter the regression or the tests.
- **Grid consistency**: S1 checks shape and affine across symptoms; sensitivity maps on a different grid are resampled by nearest neighbour.
- **Voxel order**: column-major (`vol(:)`), consistently across all steps.
- **Paths**: all scripts resolve paths relative to this directory (one level up from `code_algorithm_latest/`), with no hard-coded absolute paths.

---

## 7. Troubleshooting

| Symptom | Fix |
|---|---|
| `Grid mismatch` | Check that all symptom r-maps have the same affine |
| `Missing sensitivity maps` | Add `<SOI>/significant_8mm.nii.gz`, or set `REQUIRE_SENSITIVITY_MAP = false` |
| Empty core result | Check `q_min` in `S5_summary.txt`; relax `fdr_level` or scan thresholds with `specificity_rethreshold` |
| Adding a symptom | Create matching subfolders under `disease_data/Lesion_rmaps/` and `disease_data/Lesion_sensitivitymap/` |
