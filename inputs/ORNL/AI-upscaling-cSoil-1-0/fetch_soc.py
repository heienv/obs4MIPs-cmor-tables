"""Download and verify the two SOC mean rasters used for this product."""

import argparse
import hashlib
from pathlib import Path
from urllib.request import urlopen


FILES = {
    "30cm_SOC_mean.tif": {
        "url": "https://zenodo.org/api/records/8057232/files/30cm_SOC_mean.tif/content",
        "md5": "382205cda96f807f9391645f773dc107",
    },
    "100cm_SOC_mean.tif": {
        "url": "https://zenodo.org/api/records/8057232/files/100cm_SOC_mean.tif/content",
        "md5": "54f0a8d4b3fc07be56229842e62cba30",
    },
}


def digest(path):
    checksum = hashlib.md5()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(8 * 1024 * 1024), b""):
            checksum.update(block)
    return checksum.hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output_dir", type=Path)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    for name, metadata in FILES.items():
        output = args.output_dir / name
        if output.exists():
            if digest(output) != metadata["md5"]:
                raise ValueError(f"Existing file has an unexpected checksum: {output}")
            print(f"Already verified: {output}")
            continue

        partial = output.with_suffix(output.suffix + ".download")
        checksum = hashlib.md5()
        with urlopen(metadata["url"], timeout=120) as response, partial.open("wb") as stream:
            for block in iter(lambda: response.read(1024 * 1024), b""):
                checksum.update(block)
                stream.write(block)
        if checksum.hexdigest() != metadata["md5"]:
            raise ValueError(f"Downloaded file has an unexpected checksum: {partial}")
        partial.replace(output)
        print(f"Downloaded and verified: {output}")


if __name__ == "__main__":
    main()
