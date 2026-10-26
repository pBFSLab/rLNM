# rLNM: Robust Lesion Network Mapping

This repository provides MATLAB code for **robust lesion network mapping (rLNM)**, a null model-based framework designed to identify reliable and symptom-specific brain networks while reducing biases associated with normative functional connectivity.

## Code Structure

The repository contains two main modules:

- **`C1_sensitivity_test/`**: Spatial null model-based sensitivity testing to identify brain regions showing lesion-connectivity overlap beyond chance.
- **`C2_specificity_test/`**: Regression-based selective specificity testing to distinguish symptom-specific networks after accounting for shared connectivity patterns.

## Required Data

The sensitivity analysis requires a normative functional connectivity (FC) matrix (~600 MB), stored at:

`C1_sensitivity_test/template/rho_z_mean.mat`

This file is managed using **Git Large File Storage (Git LFS)**.

To download the repository with the complete FC matrix, install [Git LFS](https://git-lfs.com/) and run:

```bash
git lfs install
git clone https://github.com/pBFSLab/rLNM.git
```

The FC matrix will be automatically downloaded to the appropriate directory. No separate download is required.

## Usage

Open MATLAB, navigate to the corresponding analysis directory, and run the scripts. Adjust local file paths if necessary.

## License

This project is licensed under the Apache License 2.0 License.
