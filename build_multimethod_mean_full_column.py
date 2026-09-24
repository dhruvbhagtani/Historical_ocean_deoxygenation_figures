"""Build direct full-column and Olivelli maps for a six-method spatial mean.

The five-method upper and deep maps already exist. This calculates matching
full-column maps for those methods and all three Olivelli depth regions.
"""

import os
import re
import sys
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib-db9274")

import numpy as np
import xarray as xr


SOURCE = Path(
    "/scratch/gpfs/LRGROUP/db9274/Subsampled_o2_variability/"
    "Cause_model_obs_differences"
)
sys.path.insert(0, str(SOURCE))

import all_models_local_trend_weighted_pdf_median_iqr_heatmaps as global_pdf
import plot_multimodel_mean_difference_maps as difference_maps


OUTPUT = Path(__file__).resolve().parent / "Figures/multimethod_mean_full_column_1965_2021.nc"
CACHE_DIR = Path(__file__).resolve().parent / "Figures/cache_six_method_maps_1965_2021"
RF_CACHE = SOURCE / "Figures/ito2024_rf_mapped_minus_full_deoxygenation_difference_full_water_column_1965_2021.nc"
SEAM_INVENTORY = SOURCE / "Figures/individual_model_seam_diagnostics_1965_2021.csv"
OLIVELLI_NAME = "observation-Olivelli(2026)"
OLIVELLI_LABEL = "Olivelli et al. (2026)"
OLIVELLI_DIR = Path(
    "/scratch/gpfs/GEOCLIM/LRGROUP/db9274/Colab_data/"
    "Arianna_gridded_products/o2_reconstructions"
)
global_pdf.ect.TREND_PRODUCTS[OLIVELLI_NAME] = global_pdf.ect.TrendProduct(
    label="Olivelli(2026)",
    model_pattern=str(OLIVELLI_DIR / "Oxygen_reconstruction_1x1bin_{model}_{ensemble}.nc"),
    observation_path=str(OLIVELLI_DIR / "Oxygen_reconstruction_1x1bin_ORAS5.nc"),
    variable="oxygen",
)
global_pdf.PRODUCT_LABELS[OLIVELLI_NAME] = OLIVELLI_LABEL


def local_trend_in_depth_chunks(path, variable, fill_value, chunk_size=8):
    """Limit memory while computing trends from the monthly Olivelli fields."""
    with xr.open_dataset(path) as ds:
        data = global_pdf.ect.standardize_oxygen_grid(
            ds, variable=variable, fill_value=fill_value
        )
        data = global_pdf.ect._subset_years(data, global_pdf.YEARS)
        chunks = []
        for start in range(0, data.sizes["lev"], chunk_size):
            selected = data.isel(lev=slice(start, start + chunk_size))
            chunks.append(global_pdf.local_trend_per_decade(selected).load())
        return xr.concat(chunks, dim="lev")


def product_mean(product_name, inventory):
    label = global_pdf.PRODUCT_LABELS[product_name]
    depths = ("full", "upper_2000_m", "below_2000_m") if product_name == OLIVELLI_NAME else ("full",)
    fields = {depth: {} for depth in depths}
    for model_index, model in enumerate(global_pdf.ect.expt_name):
        print(f"{label}: {model}", flush=True)
        product = global_pdf.ect.TREND_PRODUCTS[product_name]
        full = global_pdf.load_local_trend(
            global_pdf.ect.fully_sampled_model_path(model_index, 0),
            variable="o2", fill_value=None,
            depth_range=global_pdf.DEPTH_RANGES["Full depth"],
        )
        if product_name == OLIVELLI_NAME:
            mapped = local_trend_in_depth_chunks(
                product.model_path(model_index, 0),
                variable=product.variable, fill_value=product.fill_value,
            )
        else:
            mapped = global_pdf.load_local_trend(
                product.model_path(model_index, 0),
                variable=product.variable, fill_value=product.fill_value,
                depth_range=global_pdf.DEPTH_RANGES["Full depth"],
            )
        full = global_pdf.ecp.align_horizontal_grid(full, mapped).interp_like(
            mapped, kwargs={"fill_value": "extrapolate"}
        )
        common = np.isfinite(full) & np.isfinite(mapped)
        difference = (mapped - full).where(common)
        maps = {"full": difference_maps.depth_weighted_mean(difference)}
        if product_name == OLIVELLI_NAME:
            maps["upper_2000_m"] = difference_maps.depth_weighted_mean(
                difference.sel(lev=slice(0, 2000))
            )
            maps["below_2000_m"] = difference_maps.depth_weighted_mean(
                difference.sel(lev=slice(2000, None))
            )
        seam_lons = sorted(set(
            inventory.get((model, label, "upper_2000_m"), [])
            + inventory.get((model, label, "below_2000_m"), [])
        ))
        if product_name == OLIVELLI_NAME:
            seam_lons = [71.5, 72.5, 73.5, 74.5]
        for depth, field in maps.items():
            fields[depth][model] = difference_maps.interpolate_missing_native_longitudes(
                field, seam_lons
            )

    means = {}
    for depth, by_model in fields.items():
        if difference_maps.MPI_LR_MODEL in by_model and difference_maps.MPI_HR_MODEL in by_model:
            by_model[difference_maps.MPI_LR_MODEL], _ = difference_maps.periodic_2d_fill_to_valid_mask(
                by_model[difference_maps.MPI_LR_MODEL],
                np.isfinite(by_model[difference_maps.MPI_HR_MODEL]),
            )

        stack = xr.concat(
            [by_model[model] for model in global_pdf.ect.expt_name],
            dim=xr.IndexVariable("model", global_pdf.ect.expt_name),
        )
        means[depth] = stack.where(stack.notnull().all("model")).mean(
            "model", skipna=False
        ).astype("float32")
    return xr.Dataset(means)


def main():
    inventory = difference_maps.load_seam_inventory(SEAM_INVENTORY)
    means = []
    labels = []
    olivelli = None
    CACHE_DIR.mkdir(parents=True, exist_ok=True)
    for name in [*global_pdf.PRODUCTS, OLIVELLI_NAME]:
        label = global_pdf.PRODUCT_LABELS[name]
        if name == "observation-Ito(2024; RF)":
            with xr.open_dataset(RF_CACHE) as ds:
                result = xr.Dataset({"full": ds.difference.load()})
        else:
            slug = re.sub(r"[^a-z0-9]+", "_", name.lower()).strip("_")
            cache = CACHE_DIR / f"{slug}.nc"
            if cache.exists():
                print(f"Loading {cache}", flush=True)
                with xr.open_dataset(cache) as ds:
                    result = ds.load()
            else:
                result = product_mean(name, inventory)
                result.to_netcdf(cache, encoding={
                    key: {"zlib": True, "complevel": 4} for key in result.data_vars
                })
                print(f"Wrote {cache}", flush=True)
        means.append(result["full"].reset_coords(drop=True))
        labels.append(label)
        if name == OLIVELLI_NAME:
            olivelli = result

    stack = xr.concat(means, dim=xr.IndexVariable("product", labels))
    common = stack.notnull().all("product")
    output = xr.Dataset({
        "difference": stack.where(common).mean("product", skipna=False).astype("float32"),
        "method_difference": stack.astype("float32"),
        "olivelli_upper": olivelli["upper_2000_m"],
        "olivelli_below": olivelli["below_2000_m"],
    })
    output.attrs.update({
        "title": "Six-method mean full-column mapped-minus-full deoxygenation trend difference",
        "period": "1965-2021",
        "units": "micromol kg-1 decade-1",
        "calculation_order": "levelwise difference, full-column volume weighting, common-model multimodel mean, common-method mean",
    })
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    output.to_netcdf(OUTPUT, encoding={
        name: {"zlib": True, "complevel": 4} for name in output.data_vars
    })
    print(f"Wrote {OUTPUT}", flush=True)


if __name__ == "__main__":
    main()
