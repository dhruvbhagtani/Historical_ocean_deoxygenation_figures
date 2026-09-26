# Paper figures

The repository contains Jupyter notebooks and the twelve PDFs used in the paper, plus a standalone map of the six analysis regions. Finished figures are in `Figures/`. The notebooks below generate them.

| PDF in `Figures/` | Notebook |
| --- | --- |
| `GFDL_ESM4_Ito2024_RF_sampling_sensitivity_full_depth_maps.pdf` | `GFDL_ESM4_sampling_sensitivity_full_depth_maps.ipynb` or `make_gfdl_ito2024_rf_3x3.ipynb` |
| `GFDL_ESM4_Ito2022_sampling_sensitivity_full_depth_maps.pdf` | `GFDL_ESM4_sampling_sensitivity_full_depth_maps.ipynb` or `make_gfdl_ito2022_3x3.ipynb` |
| `GFDL_ESM4_sampling_sensitivity_RF_consolidated.pdf` | `GFDL_ESM4_sampling_sensitivity_six_panel.ipynb` |
| `sampling_measurement_counts_six_panel.pdf` | `sampling_measurement_counts_six_panel.ipynb` |
| `oxygen_sampling_mapping_difference_trends_boxplots_1993_2021_supplement.pdf` | `Oxygen_sampling_mapping_difference_trends_two_periods.ipynb` |
| `oxygen_sampling_mapping_difference_trends_two_periods_boxplots.pdf` | `Oxygen_sampling_mapping_difference_trends_two_periods.ipynb` |
| `GFDL_ESM4_sampling_sensitivity_six_panel.pdf` | `GFDL_ESM4_sampling_sensitivity_six_panel.ipynb` |
| `all_models_multimodel_mean_mapped_minus_full_deoxygenation_difference_maps_1965_2021.pdf` | `plot_olivelli_difference_maps.ipynb` |
| `deoxygenation_EC_six_regions.pdf` | `Deoxygenation_EC_six_regions.ipynb` |
| `emergent_constraint_examples_and_individual_ec_1965_2021.pdf` | `Emergent_constraint_examples_and_individual_ec_1965_2021.ipynb` |
| `Intro_figure_deoxygenation_three_depths.pdf` | `Intro_figure_deoxygenation_three_depths.ipynb` |
| `gfdl_esm4_sampling_and_mapping_pdf_bias_by_period_1965_2021.pdf` | `GFDL_ESM4_sampling_and_mapping_pdf_bias_by_period_1965_2021.ipynb` |
| `six_analysis_regions_map.pdf` | `six_analysis_regions_map.ipynb` |

## Paper reference

Dhruv Bhagtani, Laure Resplandy, Lijing Cheng, Juan Du, Christopher Roach, Arianna Olivelli, and Takamitsu Ito. *Model-derived constraints on observational estimates imply stronger historical ocean deoxygenation*. Submitted to *Global Biogeochemical Cycles*.

## Rebuilding

Open the notebooks from this directory and run their cells in order. `build_gfdl_full_depth_sampling_caches.ipynb` generates the intermediate OI and RF NetCDF files needed by `GFDL_ESM4_sampling_sensitivity_full_depth_maps.ipynb`. `build_multimethod_mean_full_column.ipynb` generates the intermediate data used by `plot_olivelli_difference_maps.ipynb` and `Oxygen_sampling_mapping_difference_trends_two_periods.ipynb`. Those cache files are generated as needed and are not part of the stored figure set.

Figure notebooks display their plots when run, and their saved cell outputs show the figures on GitHub. The two cache-building notebooks produce intermediate data and have no figure output.

Analysis inputs live in the companion `Paper_figures`, `GFDL_ESM4_artificial_subsampling`, and `Cause_model_obs_differences` directories. Paths in the notebooks reflect their current Princeton Research Computing locations. Dependencies include Jupyter, NumPy, pandas, Matplotlib, SciPy, xarray, netCDF4, and Cartopy.

## License

See [LICENSE](LICENSE).
