"""Create separate CMOR files for the two cumulative SOC depth products."""

import argparse
import importlib.metadata
import json
from pathlib import Path
import sys

import cmor
import netCDF4
import numpy as np


SCRIPT_DIR = Path(__file__).resolve().parent
REPOSITORY_ROOT = SCRIPT_DIR.parents[2]
TABLES = REPOSITORY_ROOT / "Tables"
sys.path.append(str(REPOSITORY_ROOT / "inputs" / "misc"))
import obs4MIPsLib  # noqa: E402


def processing_code_location():
    git_commit_number = obs4MIPsLib.get_git_revision_hash()
    path_to_code = SCRIPT_DIR.relative_to(REPOSITORY_ROOT).as_posix()
    return (
        f"https://github.com/PCMDI/obs4MIPs-cmor-tables/tree/"
        f"{git_commit_number}/{path_to_code}"
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()

    output_dir = args.output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    base_metadata = json.loads((SCRIPT_DIR / "soc.json").read_text())
    base_metadata.update(
        outpath=str(output_dir),
        has_aux_unc="FALSE",
        processing_code_location=processing_code_location(),
        history=(
            "Original source raster values and masks retained; latitude and "
            "data rows reversed together; no time coordinate or spatial "
            "resampling; converted with CMOR."
        ),
    )

    products = []
    variants = {30: "ORNL-r1", 100: "ORNL-r2"}
    for depth_cm in (30, 100):
        metadata = dict(base_metadata)
        metadata["variant_label"] = variants[depth_cm]
        metadata["variant_info"] = (
            f"Prepared at ORNL from the cumulative 0-{depth_cm} cm SOC mean. "
            "The depth-specific label distinguishes the two overlapping "
            "cumulative products."
        )
        input_json = output_dir / f"cmor_input_{depth_cm}cm.json"
        input_json.write_text(json.dumps(metadata, indent=2) + "\n")

        staging_path = args.input_dir / f"SOC_0-{depth_cm}cm_staging.nc"
        with netCDF4.Dataset(staging_path) as source:
            cmor.setup(
                inpath=str(TABLES),
                netcdf_file_action=cmor.CMOR_REPLACE_4,
                exit_control=cmor.CMOR_EXIT_ON_MAJOR,
                logfile=str(output_dir / f"cmor_{depth_cm}cm.log"),
            )
            cmor.dataset_json(str(input_json))
            cmor.load_table("obs4MIPs_fx.json")
            cmor.set_cur_dataset_attribute(
                "geospatial_coverage",
                (
                    "Original 30 arc-second rectangle: 180 W to 55 W, "
                    "15 N to 80 N; source mask retained"
                ),
            )
            cmor.set_cur_dataset_attribute(
                "source_doi", "10.5281/zenodo.8057232"
            )

            latitude = cmor.axis(
                "latitude",
                coord_vals=source["lat"][:],
                cell_bounds=source["lat_bnds"][:],
                units="degrees_north",
            )
            longitude = cmor.axis(
                "longitude",
                coord_vals=source["lon"][:],
                cell_bounds=source["lon_bnds"][:],
                units="degrees_east",
            )
            depth = cmor.axis(
                "sdepth",
                coord_vals=np.array([depth_cm / 200.0]),
                cell_bounds=np.array([[0.0, depth_cm / 100.0]]),
                units="m",
            )
            variable = cmor.variable(
                "cSoil",
                "kg m-2",
                [depth, latitude, longitude],
                missing_value=1e20,
            )
            cmor.set_deflate(variable, 1, 1, 2)
            values = np.asarray(
                source["cSoil"][:].filled(1e20), dtype="float32"
            )[None, :, :]
            cmor.write(variable, values)
            output = cmor.close(variable, file_name=True)
            cmor.close()

        if isinstance(output, bytes):
            output = output.decode()
        products.append({"depth_cm": [0, depth_cm], "output": output})
        print(output, flush=True)

    report = {
        "cmor_version": importlib.metadata.version("cmor"),
        "products": products,
    }
    (output_dir / "run.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
