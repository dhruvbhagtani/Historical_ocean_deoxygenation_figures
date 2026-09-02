# Historical ocean deoxygenation figures

This repository contains Jupyter notebooks and publication-ready figures used to examine historical ocean deoxygenation, uncertainty from observational sampling and mapping, and emergent constraints on global oxygen trends. The analyses compare six gridded observational products and use Earth system model sampling experiments, including GFDL-ESM4, to evaluate how incomplete observations affect reconstructed trends.

## Repository contents

| Notebook | Purpose | Main output |
| --- | --- | --- |
| `Intro_figure_deoxygenation_three_depths.ipynb` | Summarizes oxygen changes, observational coverage, sampled-and-mapped model estimates, and emergent-constraint estimates for the full water column, upper 2000 m, and ocean below 2000 m. | `Figures/Intro_figure_deoxygenation_three_depths.pdf` |
| `Emergent_constraint_examples_and_individual_ec_1965_2021.ipynb` | Shows representative emergent-constraint relationships and reconstruction-specific uncertainty distributions for 1965–2021. | `Figures/emergent_constraint_examples_and_individual_ec_1965_2021.pdf` |
| `GFDL_ESM4_sampling_sensitivity_six_panel.ipynb` | Evaluates GFDL-ESM4 trend errors under uniform, random-global, and enhanced ship-track sampling for objective-interpolation, neural-network, and random-forest reconstructions. | `Figures/GFDL_ESM4_sampling_sensitivity_six_panel.pdf` |
| `Oxygen_sampling_mapping_difference_trends_two_periods.ipynb` | Compares sampling–mapping trend differences across six observational products for 1965–2021 and 1993–2021. | Two PDF variants in `Figures/` |

Rendered PDFs are tracked in [`Figures`](Figures/) so the results can be viewed without rerunning the analyses.

## Requirements

The notebooks use Python 3 and Jupyter, with the following principal packages:

- `numpy`
- `pandas`
- `matplotlib`
- `scipy`
- `xarray`
- `netCDF4`
- `cartopy`

Only the packages needed by a particular notebook have to be installed. For example, the GFDL-ESM4 sampling notebook requires `xarray`, `netCDF4`, and `cartopy`, while the plotting notebooks built from pickle caches do not.

## Data and portability

The notebooks are plotting and analysis products rather than a self-contained data archive. They read precomputed pickle, CSV, and NetCDF inputs from the companion analysis directories `Paper_figures` and `GFDL_ESM4_artificial_subsampling`. These large source data and intermediate caches are not included in this repository.

Several input locations are currently configured as absolute Princeton Research Computing paths near the beginning of each notebook. To run the notebooks elsewhere:

1. Obtain or reproduce the required input caches and model-sampling outputs.
2. Update `SOURCE_DIR`, `ROOT`, `PAPER`, `ARTIFICIAL`, and/or `PYTHON_FUNCTIONS` in the relevant setup cell.
3. Preserve the expected filenames and directory structure, or update the individual path constants accordingly.

The emergent-constraint notebook additionally imports `emergent_constraint_tests` from the external `Python_functions` directory.

## Running the notebooks

From the repository root, start Jupyter:

```bash
jupyter lab
```

Open a notebook and run its cells from top to bottom. Output directories are created by the notebooks where needed. Most final figures are written to `Figures/`; the GFDL-ESM4 workflow also writes intermediate metrics and a cached random-sampling trend map to its configured companion `Paper_figures` directory.

Because some workflows process multi-dimensional model output, memory use and execution time vary substantially by notebook. The checked-in PDFs correspond to completed notebook runs and provide a reference for the expected result.

## License

This project is distributed under the terms in [`LICENSE`](LICENSE).
