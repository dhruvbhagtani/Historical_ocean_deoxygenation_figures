"""Cache full-depth mapped-minus-full trends for five GFDL sampling experiments.

The calculation matches the existing upper/below-2000 m caches: annualize,
calculate a trend at every depth, subtract on common valid levels, then weight
the level differences by model-layer thickness over the full water column.
"""

import argparse
from pathlib import Path

import netCDF4 as nc
import numpy as np
import xarray as xr


WORK = Path('/scratch/gpfs/LRGROUP/db9274/GFDL_ESM4_artificial_subsampling')
FULL = Path('/scratch/gpfs/LRGROUP/db9274/Variability_trends_EC/Paper_figures/NETCDF/Fully_sampled_models_annual_o2/o2_1x1bin_GFDL_ESM4_r1i1p1f1.nc')
OI_DIR = WORK / 'workflows/Ito_2022_OI/NETCDF/processed'
RF_DIR = WORK / 'workflows/Ito_2024_NN/NETCDF/Artificial_subsampling_RF'
OUTPUT_DIR = Path(__file__).resolve().parent / 'Figures'
REGIONS = ('tropical', 'extratropical', 'southern_ocean')
EXPERIMENTS = ('uniform_historical_coverage', 'uniform_30_percent',
               'added_tropical_2000', 'added_north_atlantic_pacific_2000',
               'added_southern_ocean_2000')


def inputs(method):
    if method == 'oi':
        paths = [
            OI_DIR / 'o2_grid_1x1_GFDL_ESM4_randobs.nc',
            OI_DIR / 'o2_grid_1x1_GFDL_ESM4_pct30.nc',
            *(OI_DIR / f'o2_grid_1x1_GFDL_ESM4_regional_{region}_2000.nc'
              for region in REGIONS),
        ]
        return (1967, 2019), paths
    paths = [
        RF_DIR / 'GFDL_ESM4_randobs_r1i1p1f1_oxygen_gridded_RF.nc',
        RF_DIR / 'GFDL_ESM4_pct30_r1i1p1f1_oxygen_gridded_RF.nc',
        *(RF_DIR / f'GFDL_ESM4_regional_{region}_2000_r1i1p1f1_oxygen_gridded_RF.nc'
          for region in REGIONS),
    ]
    return (1965, 2021), paths


def as_nan(values):
    result = np.ma.filled(values, np.nan).astype(np.float64)
    result[result > 1e19] = np.nan
    return result


def level_trend(path, start, end):
    """Read annual fields and require a complete time series at each level."""
    with nc.Dataset(path) as ds:
        time = ds.variables['time']
        years = np.array([date.year for date in nc.num2date(
            time[:], time.units, calendar=getattr(time, 'calendar', 'standard'))])
        zname = 'lev' if 'lev' in ds.variables else 'depth'
        yname = 'y' if 'y' in ds.variables else 'lat'
        xname = 'x' if 'x' in ds.variables else 'lon'
        z, lat, lon = (np.asarray(ds.variables[name][:], dtype=float)
                       for name in (zname, yname, xname))
        lon = np.mod(lon, 360)
        var = ds.variables['o2']
        time_axis = var.dimensions.index('time')
        remaining = [name for name in var.dimensions if name != 'time']
        order = [remaining.index(name) for name in (zname, yname, xname)]
        shape = (len(z), len(lat), len(lon))
        wanted = np.arange(start, end + 1)
        missing = [year for year in wanted if not np.any(years == year)]
        if missing:
            raise ValueError(f'{path}: missing years {missing}')

        # Ito OI stores time last. Longitude slabs avoid slow strided reads.
        if var.dimensions == (xname, yname, zname, 'time'):
            if any(np.sum(years == year) != 1 for year in wanted):
                raise ValueError(f'{path}: expected one OI field per year')
            indices = np.array([np.flatnonzero(years == year)[0] for year in wanted])
            centered = wanted.astype(float) - wanted.mean()
            denominator = np.sum(centered ** 2)
            trend = np.full(shape, np.nan)
            for first in range(0, len(lon), 90):
                last = min(first + 90, len(lon))
                block = as_nan(var[first:last, :, :, indices])
                block = np.transpose(block, (2, 1, 0, 3))
                valid = np.isfinite(block).all(axis=-1)
                slope = np.sum(block * centered[None, None, None, :], axis=-1)
                trend[:, :, first:last] = np.where(valid, slope / denominator * 10, np.nan)
            return z, lat, lon, trend

        n = np.zeros(shape)
        sx = np.zeros(shape)
        sy = np.zeros(shape)
        sxx = np.zeros(shape)
        sxy = np.zeros(shape)
        for year in wanted:
            total = np.zeros(shape)
            count = np.zeros(shape, dtype=np.uint8)
            for index in np.flatnonzero(years == year):
                slicer = [slice(None)] * var.ndim
                slicer[time_axis] = int(index)
                field = np.transpose(as_nan(var[tuple(slicer)]), order)
                valid = np.isfinite(field)
                total += np.where(valid, field, 0)
                count += valid
            annual = np.divide(total, count, out=np.full(shape, np.nan), where=count > 0)
            valid = np.isfinite(annual)
            n += valid
            sx += np.where(valid, year, 0)
            sy += np.where(valid, annual, 0)
            sxx += np.where(valid, year * year, 0)
            sxy += np.where(valid, year * annual, 0)
        denominator = n * sxx - sx * sx
        slope = np.divide(n * sxy - sx * sy, denominator,
                          out=np.full(shape, np.nan),
                          where=(n == len(wanted)) & (denominator != 0))
        return z, lat, lon, slope * 10


def full_depth_difference(mapped, full, z):
    edges = np.empty(len(z) + 1)
    edges[1:-1] = (z[:-1] + z[1:]) / 2
    edges[0] = max(0, z[0] - (z[1] - z[0]) / 2)
    edges[-1] = z[-1] + (z[-1] - z[-2]) / 2
    thickness = np.diff(edges)[:, None, None]
    difference = mapped - full
    valid = np.isfinite(difference)
    numerator = np.sum(np.where(valid, difference * thickness, 0), axis=0)
    denominator = np.sum(np.where(valid, thickness, 0), axis=0)
    return np.divide(numerator, denominator,
                     out=np.full(numerator.shape, np.nan), where=denominator > 0)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--method', choices=('oi', 'rf'), required=True)
    parser.add_argument('--check-inputs', action='store_true')
    args = parser.parse_args()
    (start, end), paths = inputs(args.method)
    for path in [FULL, *paths]:
        if not path.is_file():
            raise FileNotFoundError(path)
    if args.check_inputs:
        print(f'{args.method}: {start}-{end}, {len(paths)} experiments; inputs exist')
        return

    z, lat, lon, full = level_trend(FULL, start, end)
    fields = []
    for name, path in zip(EXPERIMENTS, paths):
        print(f'{args.method}: calculating {name} from {path}', flush=True)
        zz, yy, xx, mapped = level_trend(path, start, end)
        if not (np.allclose(z, zz) and np.allclose(lat, yy)
                and np.allclose(lon, xx)):
            raise ValueError(f'Coordinate mismatch for {path}')
        field = full_depth_difference(mapped, full, z)
        fields.append(field.astype('float32'))
        del mapped

    data = xr.Dataset(
        {'mapped_minus_full': (('experiment', 'lat', 'lon'), np.stack(fields))},
        coords={'experiment': list(EXPERIMENTS), 'lat': lat, 'lon': lon},
        attrs={
            'title': f'GFDL-ESM4 {args.method.upper()} idealized sampling full-depth trend errors',
            'period': f'{start}-{end}',
            'units': 'micromol kg-1 decade-1',
            'difference_sign': 'mapped minus fully sampled',
            'calculation_order': 'annualize; levelwise trend; subtract on common valid levels; full-depth thickness weighting',
            'historical_sampling': 'Use the existing original-sampling map',
            'source_files': '\n'.join(map(str, paths)),
        },
    )
    if args.method == 'oi':
        data.attrs['unit_caveat'] = 'OI files are labeled micromol O2 L-1; no density conversion applied'
    output = OUTPUT_DIR / f'GFDL_ESM4_{args.method}_idealized_sampling_full_depth_1965_2021.nc'
    if args.method == 'oi':
        output = OUTPUT_DIR / 'GFDL_ESM4_oi_idealized_sampling_full_depth_1967_2019.nc'
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    temporary = output.with_suffix('.tmp.nc')
    data.to_netcdf(temporary, encoding={'mapped_minus_full': {'zlib': True, 'complevel': 4}})
    temporary.replace(output)
    print(f'Wrote {output}', flush=True)


if __name__ == '__main__':
    main()
