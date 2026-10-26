# rLNM: Robust Lesion Network Mapping

This repository provides MATLAB code for **robust lesion network mapping (rLNM)**, a null model-based framework designed to identify reliable and symptom-specific brain networks while reducing biases associated with normative functional connectivity.

## Code Structure

The repository contains two main modules:

- **`C1_sensitivity_test/`**: Spatial null model-based sensitivity testing to identify brain regions showing lesion-connectivity overlap beyond chance.
- **`C2_specificity_test/`**: Regression-based selective specificity testing to distinguish symptom-specific networks after accounting for shared connectivity patterns.

## Required Data

The sensitivity analysis requires a normative functional connectivity (FC) matrix (~600 MB).

Due to GitHub's file size limitations, the matrix is provided separately in [Release v1.0](https://github.com/pBFSLab/rLNM/releases/tag/v1.0).

**Instructions:**

1. Download `rho_z_mean.mat` from Release v1.0.
2. Place the file in the `C1_sensitivity_test/template/` directory.
3. Run the corresponding MATLAB analysis scripts.

## Usage

Clone or download this repository:

```bash
git clone https://github.com/pBFSLab/rLNM.git
```

Open MATLAB, navigate to the corresponding analysis directory, and run the scripts following their instructions. Adjust local file paths if necessary.

## License

This repository is distributed under the Apache License 2.0.
