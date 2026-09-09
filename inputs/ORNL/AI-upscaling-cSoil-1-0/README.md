# AI-upscaling-cSoil-1-0 preparation

This directory prepares the two soil organic carbon (SOC) stock products from
[Zenodo record 8057232](https://doi.org/10.5281/zenodo.8057232) for obs4MIPs.
The registered source ID is `AI-upscaling-cSoil-1-0` and the preparation is
tracked in issue #593.

## Source data

Only the mean SOC rasters are processed. The coefficient-of-variation rasters
are outside the scope of this contribution.

| File | Cumulative depth | Zenodo MD5 |
| --- | --- | --- |
| `30cm_SOC_mean.tif` | 0--30 cm | `382205cda96f807f9391645f773dc107` |
| `100cm_SOC_mean.tif` | 0--100 cm | `54f0a8d4b3fc07be56229842e62cba30` |

The GeoTIFFs use WGS84 longitude/latitude coordinates at 30 arc-second
resolution. Their complete rectangular grid is retained: 180 W to 55 W and
15 N to 80 N, with 7,800 latitude rows and 15,000 longitude columns. Values,
missing-data masks, and coordinates are not cropped or spatially resampled.

The products are cumulative stocks over overlapping intervals. They are
therefore written as separate files, provisionally distinguished by variant
labels `ORNL-0to30cm` and `ORNL-0to100cm`. They must not be treated as adjacent
layers or summed.

## Reproduce

Python dependencies are `cmor`, `netCDF4`, `numpy`, and `rasterio`. The full
workflow used CMOR 3.15.3.

```bash
python fetch_soc.py raw
python prepare_soc.py --raw-dir raw --input-units 'kg C m-2' --output-dir staging
python runCMOR_soc.py --input-dir staging --output-dir output
python validate_cmor_soc.py --input-dir staging --output-dir output
```

For a small processing test, use the repository sample:

```bash
python prepare_soc.py --sample soc_data_test.nc --output-dir staging-sample
python runCMOR_soc.py --input-dir staging-sample --output-dir output-sample
python validate_cmor_soc.py --input-dir staging-sample --output-dir output-sample
```

`prepare_soc.py` reverses latitude and data rows together to satisfy the
increasing coordinate requirement. It performs no numeric unit conversion:
`kg C m-2` in the source represents the same carbon mass per area encoded as
the CF unit `kg m-2` for `soil_mass_content_of_carbon`.

## Review items

The code uses `grid_label=gn` because it preserves the source grid. The current
source registration says `contiguous_united_states`, while valid cells in the
complete source footprint also cover Alaska, Hawaii, and Puerto Rico. Maintainer
review is requested for the final region value, the `gn`/`gr` choice, and the
depth-specific file identifiers. After those decisions, the full products must
be regenerated from the merged tables before ESGF publication.
