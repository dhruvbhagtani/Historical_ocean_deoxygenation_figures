from pathlib import Path
import string
import textwrap
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
from cartopy import crs as ccrs
import cartopy.feature as cfeature
from cartopy.util import add_cyclic_point

import make_gfdl_ito2022_3x3 as common

WORK = Path('/scratch/gpfs/LRGROUP/db9274/GFDL_ESM4_artificial_subsampling')
RF_DIR = WORK / 'workflows/Ito_2024_NN/NETCDF/Artificial_subsampling_RF'
MAP_DIR = WORK / 'trend_comparison_outputs/fully_minus_gap_filled_trend_maps'
ORIGINAL_MAP = MAP_DIR / 'ito2024_rf_original_4x_4xeast_8x.nc'
RANDOM_MAP = Path('/scratch/gpfs/LRGROUP/db9274/Variability_trends_EC/Paper_figures/cache_gfdl_esm4_randobs_trend_error_maps.nc')
PCT30_MAP = MAP_DIR / 'fully_minus_gap_filled_rf_pct30_level_first_volume_weighted_1965_2021.nc'
OUT = Path('Figures/GFDL_ESM4_Ito2024_RF_sampling_sensitivity_full_depth_maps.pdf')

def main():
    z, lat, lon, full = common.trend(common.FULL, 1965, 2021)
    upper_v, lower_v = common.band_volumes()
    with (
        xr.open_dataset(ORIGINAL_MAP) as original,
        xr.open_dataset(RANDOM_MAP) as random,
        xr.open_dataset(PCT30_MAP) as pct30,
    ):
        fields = [
            -common.from_bands(original, 'fully_minus_gap_filled_original', upper_v, lower_v),
            -common.from_bands(random, 'randobs_rf', upper_v, lower_v),
            -common.from_bands(pct30, 'fully_minus_gap_filled_rf30', upper_v, lower_v),
        ]

    cases = [
        'regional_tropical_2000',
        'regional_extratropical_2000', 'regional_southern_ocean_2000',
    ]
    for case in cases:
        path = RF_DIR / f'GFDL_ESM4_{case}_r1i1p1f1_oxygen_gridded_RF.nc'
        zz, yy, xx, mapped = common.trend(path, 1965, 2021)
        if not (np.allclose(z, zz) and np.allclose(lat, yy) and np.allclose(lon, xx)):
            raise ValueError(f'Coordinate mismatch: {path}')
        fields.append(common.volume_mean(mapped, zz) - common.volume_mean(full, z))
        print(f'Completed {case}', flush=True)

    labels = [
        'Historical\ncoverage',
        'Uniform sampling at\nhistorical coverage',
        'Uniform sampling at\n30% coverage',
        '2,000 added profiles:\nTropical Oceans',
        '2,000 added profiles:\nNorth Atlantic +\nNorth Pacific Oceans',
        '2,000 added profiles:\nSouthern Ocean',
    ]
    available = [f for f in fields if f is not None]
    values = np.concatenate([np.abs(f[np.isfinite(f)]) for f in available])
    limit = float(np.nanpercentile(values, 99))
    norm = TwoSlopeNorm(vmin=-limit, vcenter=0, vmax=limit)
    plt.rcParams.update({'font.size': 11, 'pdf.fonttype': 42, 'savefig.dpi': 600, 'savefig.facecolor': 'white'})
    fig, axes = plt.subplots(2, 3, figsize=(10, 6.0), subplot_kw={'projection': ccrs.Robinson(central_longitude=180)}, constrained_layout=True)
    mesh = None
    for i, (ax, field, label) in enumerate(zip(axes.flat, fields, labels)):
        if field is not None:
            cyclic, cyclic_lon = add_cyclic_point(field, coord=lon)
            mesh = ax.pcolormesh(cyclic_lon, lat, cyclic, transform=ccrs.PlateCarree(), cmap='bwr', norm=norm, shading='auto', rasterized=True)
        ax.set_global(); ax.add_feature(cfeature.LAND, facecolor='0.85', edgecolor='none', zorder=2); ax.coastlines(linewidth=.45, color='.15', zorder=3); ax.gridlines(linewidth=.25, color='.55', alpha=.4)
        ax.set_title(f'{string.ascii_lowercase[i]}. {label}', fontsize=12, fontweight='bold', pad=7)
        if field is None:
            ax.text(.5, .5, 'Global RF output\nnot available', transform=ax.transAxes, ha='center', va='center', fontsize=11, color='.25')
    cb = fig.colorbar(mesh, ax=axes, orientation='horizontal', shrink=.62, pad=.08, aspect=35, extend='both'); cb.set_label('Sampling-mapping effect (µmol kg⁻¹ decade⁻¹)', fontsize=11); cb.ax.tick_params(labelsize=11)
    OUT.parent.mkdir(exist_ok=True); fig.savefig(OUT, bbox_inches='tight', pad_inches=.03); print(f'Saved {OUT}', flush=True)

if __name__ == '__main__':
    main()
