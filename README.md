# Paper figures

This folder contains the twelve PDFs used in the paper, their plotting code, and code needed to prepare intermediate data. The finished PDFs are in `Figures/`.

| PDF in `Figures/` | Producer |
| --- | --- |
| `GFDL_ESM4_Ito2024_RF_sampling_sensitivity_full_depth_maps.pdf` | `GFDL_ESM4_sampling_sensitivity_full_depth_maps.ipynb` (also `make_gfdl_ito2024_rf_3x3.py`) |
| `GFDL_ESM4_Ito2022_sampling_sensitivity_full_depth_maps.pdf` | `GFDL_ESM4_sampling_sensitivity_full_depth_maps.ipynb` (also `make_gfdl_ito2022_3x3.py`) |
| `GFDL_ESM4_sampling_sensitivity_RF_consolidated.pdf` | `GFDL_ESM4_sampling_sensitivity_six_panel.ipynb` |
| `sampling_measurement_counts_six_panel.pdf` | `sampling_measurement_counts_six_panel.ipynb` |
| `oxygen_sampling_mapping_difference_trends_boxplots_1993_2021_supplement.pdf` | `Oxygen_sampling_mapping_difference_trends_two_periods.ipynb` |
| `oxygen_sampling_mapping_difference_trends_two_periods_boxplots.pdf` | `Oxygen_sampling_mapping_difference_trends_two_periods.ipynb` |
| `GFDL_ESM4_sampling_sensitivity_six_panel.pdf` | `GFDL_ESM4_sampling_sensitivity_six_panel.ipynb` |
| `all_models_multimodel_mean_mapped_minus_full_deoxygenation_difference_maps_1965_2021.pdf` | `plot_olivelli_difference_maps.py` |
| `deoxygenation_EC_six_regions.pdf` | `Deoxygenation_EC_six_regions.ipynb` |
| `emergent_constraint_examples_and_individual_ec_1965_2021.pdf` | `Emergent_constraint_examples_and_individual_ec_1965_2021.ipynb` |
| `Intro_figure_deoxygenation_three_depths.pdf` | `Intro_figure_deoxygenation_three_depths.ipynb` |
| `gfdl_esm4_sampling_and_mapping_pdf_bias_by_period_1965_2021.pdf` | `GFDL_ESM4_sampling_and_mapping_pdf_bias_by_period_1965_2021.ipynb` |

## Rebuilding

Run from this folder. Analysis inputs live in the companion `Paper_figures`, `GFDL_ESM4_artificial_subsampling`, and `Cause_model_obs_differences` directories. Paths in the code reflect their current Princeton Research Computing locations. Typical dependencies are Jupyter, NumPy, pandas, Matplotlib, SciPy, xarray, netCDF4, and Cartopy.

The full-depth sampling notebook reads two generated NetCDF caches. Create them with `build_gfdl_full_depth_sampling_caches.py --method oi` and `--method rf`, or use `build_gfdl_full_depth_sampling_caches.slurm`. The standalone 3×3 plot scripts and their SLURM files offer another route to the same two PDFs.

For the six-method mean map PDF and the two-period boxplot/map PDF, first run `build_multimethod_mean_full_column.py` (or its SLURM file). Then run `plot_olivelli_difference_maps.py` and `Oxygen_sampling_mapping_difference_trends_two_periods.ipynb`. `render_multimethod_mean_boxplot.slurm` runs the boxplot/map cell. Intermediate NetCDF and NPZ files are generated as needed and are omitted from this folder's stored contents.

The other notebooks write their listed PDFs when run from top to bottom. Some notebooks also write analysis caches to companion directories.

## License

See [LICENSE](LICENSE).
