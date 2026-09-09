"""Create and verify staging NetCDFs; these are NOT CMOR/ESGF products.

One file per cumulative depth avoids treating overlapping stocks as disjoint layers.
Values are preserved; no spatial resampling or time coordinate is introduced.
"""
import argparse
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path

import netCDF4
import numpy as np
import rasterio
from rasterio.windows import Window

ROOT = Path(__file__).resolve().parents[3]
FILL = np.float32(1e20)
SOURCE_MD5 = {
    30: "382205cda96f807f9391645f773dc107",
    100: "54f0a8d4b3fc07be56229842e62cba30",
}


def checksum(path, algorithm="sha256"):
    h = hashlib.new(algorithm)
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(8 * 1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def bounds(centres):
    centres = np.asarray(centres, dtype=np.float64)
    delta = np.diff(centres)
    if len(centres) < 2 or not np.all(np.isfinite(centres)):
        raise ValueError("At least two finite coordinate centres required")
    if not np.all(delta > 0) or not np.allclose(delta, delta[0], rtol=1e-7, atol=1e-10):
        raise ValueError("Coordinates must be increasing and regularly spaced")
    edges = np.r_[centres[0] - delta[0] / 2, (centres[:-1] + centres[1:]) / 2,
                  centres[-1] + delta[-1] / 2]
    return np.column_stack((edges[:-1], edges[1:]))


def clean(values):
    a = np.ma.asarray(values, dtype=np.float32)
    mask = np.ma.getmaskarray(a) | ~np.isfinite(a.data)
    return np.ma.array(a.data, mask=mask)


class Input:
    def __init__(self, path, depth, units):
        self.path, self.depth, self.units = Path(path), depth, units
        self.tiff = self.path.suffix.lower() in (".tif", ".tiff")
        if self.tiff:
            self.ds = rasterio.open(path)
            t = self.ds.transform
            if self.ds.crs != rasterio.crs.CRS.from_epsg(4326) or t.b != 0 or t.d != 0 or t.a <= 0 or t.e >= 0:
                self.ds.close()
                raise ValueError("Expected north-up EPSG:4326 raster with positive longitude spacing")
            if self.ds.count != 1 or self.ds.scales != (1.0,) or self.ds.offsets != (0.0,):
                self.ds.close()
                raise ValueError("Expected one band with no scale/offset")
            self.ny, self.nx = self.ds.height, self.ds.width
            self.lon = t.c + (np.arange(self.nx) + 0.5) * t.a
            self.lat = (t.f + (np.arange(self.ny) + 0.5) * t.e)[::-1]
            self.source_variable = f"band1:0-{depth}cm"
        else:
            self.ds = netCDF4.Dataset(path)
            self.source_variable = f"SOC_{depth}cm_mean"
            v = self.ds[self.source_variable]
            if v.dimensions != ("lat", "lon") or v.units not in ("kg C m-2", "kg m-2"):
                self.ds.close()
                raise ValueError("Expected static (lat, lon) SOC in kg C m-2 or kg m-2")
            self.lat, self.lon = np.array(self.ds['lat'][:]), np.array(self.ds['lon'][:])
            self.ny, self.nx = v.shape
            self.units = v.units
        bounds(self.lat)
        bounds(self.lon)

    def read(self, start, stop):
        if self.tiff:
            a = self.ds.read(1, window=Window(0, self.ny-stop, self.nx, stop-start), masked=True)[::-1]
        else:
            a = self.ds[self.source_variable][start:stop, :]
        return clean(a)

    def close(self):
        self.ds.close()


def prepare(source, output):
    output = Path(output)
    if output.exists():
        raise FileExistsError(f"Refusing to overwrite {output}")
    output.parent.mkdir(parents=True, exist_ok=True)
    temporary = output.with_suffix('.partial.nc')
    sha = checksum(source.path)
    with netCDF4.Dataset(temporary, 'w') as d:
        for dim, size in [('lat', source.ny), ('lon', source.nx), ('bnds', 2)]:
            d.createDimension(dim, size)
        for name, values, units, axis, standard in [
            ('lat', source.lat, 'degrees_north', 'Y', 'latitude'),
            ('lon', source.lon, 'degrees_east', 'X', 'longitude')]:
            v = d.createVariable(name, 'f8', (name,))
            v.setncatts(dict(units=units, standard_name=standard, axis=axis, bounds=name+'_bnds'))
            v[:] = values
            d.createVariable(name+'_bnds', 'f8', (name, 'bnds'))[:] = bounds(values)
        depth = d.createVariable('depth', 'f8', ())
        depth.setncatts(dict(standard_name='depth', long_name='Midpoint of cumulative SOC integration interval',
                           units='m', positive='down', axis='Z', bounds='depth_bnds'))
        depth.assignValue(source.depth / 200)
        d.createVariable('depth_bnds', 'f8', ('bnds',))[:] = [0, source.depth / 100]
        soc = d.createVariable('cSoil', 'f4', ('lat', 'lon'), fill_value=FILL,
                               zlib=True, complevel=2, chunksizes=(min(128, source.ny), min(512, source.nx)))
        soc.setncatts(dict(units='kg m-2', standard_name='soil_mass_content_of_carbon',
                          long_name=f'Soil organic carbon stock integrated from 0 to {source.depth} cm',
                          coordinates='depth', comment='Cumulative stock; do not sum with overlapping depth products.'))
        d.setncatts(dict(
            title=f'SOC 0-{source.depth} cm - staging data for obs4MIPs preparation',
            processing_status='STAGING ONLY: this file is input to CMOR, not an obs4MIPs product',
            source_data_url='https://doi.org/10.5281/zenodo.8057232', source_file=source.path.name,
            source_sha256=sha, source_variable=source.source_variable, source_units=source.units,
            license='Data in this file is licensed under a Creative Commons Attribution 4.0 International License (https://creativecommons.org/licenses/by/4.0/).',
            references='https://doi.org/10.1029/2023JG007702',
            history=f'{datetime.now(timezone.utc).isoformat()}: prepare_soc.py; preserve values and grid; latitude ascending; no resampling or time axis; CV excluded.',
            comment='Full source footprint retained; the cumulative depth intervals are stored separately.'))
        for start in range(0, source.ny, 128):
            stop = min(start+128, source.ny)
            soc[start:stop] = source.read(start, stop)
    temporary.replace(output)
    return validate(source, output)


def validate(source, output):
    """Compare every output value and mask to source, in bounded-memory chunks."""
    count = negative = 0
    minimum, maximum = np.inf, -np.inf
    finite_extent = [np.inf, np.inf, -np.inf, -np.inf]
    with netCDF4.Dataset(output) as d:
        assert 'time' not in d.variables and 'time' not in d.dimensions
        assert not any('CV' in name or 'coefficient' in name for name in d.variables)
        assert d['cSoil'].dimensions == ('lat', 'lon')
        assert d['cSoil'].units == 'kg m-2'
        np.testing.assert_array_equal(d['depth_bnds'][:], [0, source.depth / 100])
        for name, coord in [('lat', source.lat), ('lon', source.lon)]:
            np.testing.assert_array_equal(d[name][:], coord)
            np.testing.assert_array_equal(d[name+'_bnds'][:], bounds(coord))
        for start in range(0, source.ny, 128):
            stop = min(start+128, source.ny)
            expected, actual = source.read(start, stop), clean(d['cSoil'][start:stop])
            np.testing.assert_array_equal(np.ma.getmaskarray(actual), np.ma.getmaskarray(expected))
            np.testing.assert_array_equal(actual.compressed(), expected.compressed())
            vals = actual.compressed()
            if vals.size:
                count += int(vals.size)
                negative += int(np.sum(vals < 0))
                minimum, maximum = min(minimum, float(vals.min())), max(maximum, float(vals.max()))
                rows, cols = np.where(~np.ma.getmaskarray(actual))
                finite_extent = [min(finite_extent[0],float(source.lon[cols.min()])),
                                 min(finite_extent[1],float(source.lat[start+rows.min()])),
                                 max(finite_extent[2],float(source.lon[cols.max()])),
                                 max(finite_extent[3],float(source.lat[start+rows.max()]))]
    if not count:
        raise ValueError('No valid SOC cells')
    return dict(output=str(output), depth_cm=[0,source.depth], shape=[source.ny,source.nx],
                input_sha256=checksum(source.path), output_sha256=checksum(output),
                valid_cells=count, negative_cells=negative, minimum=minimum, maximum=maximum,
                finite_cell_centre_extent_wsen=finite_extent,
                grid_spacing_degrees=[float(np.diff(source.lon)[0]),float(np.diff(source.lat)[0])],
                full_value_and_mask_comparison='PASS', coordinates_and_bounds='PASS',
                publication_ready=False)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    group=p.add_mutually_exclusive_group(required=True)
    group.add_argument('--sample', type=Path, help='Upstream two-depth sample NetCDF')
    group.add_argument('--raw-dir', type=Path, help='Directory with both original GeoTIFFs')
    p.add_argument('--input-units', choices=['kg C m-2','kg m-2'], help='Required declaration for GeoTIFF inputs')
    p.add_argument('--output-dir', type=Path, required=True)
    p.add_argument('--validate-only', action='store_true')
    args=p.parse_args()
    if args.raw_dir and not args.input_units:
        p.error('GeoTIFF input requires --input-units; no unit conversion is inferred')
    results=[]
    for depth in [30,100]:
        path=args.sample or args.raw_dir/f'{depth}cm_SOC_mean.tif'
        if args.raw_dir:
            if checksum(path,'md5') != SOURCE_MD5[depth]:
                raise ValueError(f'Zenodo checksum mismatch: {path}')
        source=Input(path,depth,args.input_units)
        try:
            output=(args.output_dir/f'SOC_0-{depth}cm_staging.nc').resolve()
            results.append(validate(source,output) if args.validate_only else prepare(source,output))
            print(json.dumps(results[-1]), flush=True)
        finally:
            source.close()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    (args.output_dir/'validation.json').write_text(json.dumps(results,indent=2)+'\n')


if __name__=='__main__':
    main()
