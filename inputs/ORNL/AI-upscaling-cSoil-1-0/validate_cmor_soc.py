"""Validate actual CMOR output, including every cell against staging inputs."""
import argparse
import hashlib
import json
from pathlib import Path
import uuid

import netCDF4
import numpy as np

ROOT=Path(__file__).resolve().parents[3]


def validate(input_dir,output_dir):
    run=json.loads((output_dir/'run.json').read_text())
    cv=json.loads((ROOT/'Tables/obs4MIPs_CV.json').read_text())['CV']
    reports=[]
    tracking=[]
    variants={30:'ORNL-r1',100:'ORNL-r2'}
    for entry in run['products']:
        depth=entry['depth_cm'][1]
        path=Path(entry['output'])
        with netCDF4.Dataset(input_dir/f'SOC_0-{depth}cm_staging.nc') as source, netCDF4.Dataset(path) as target:
            assert target.data_model=='NETCDF4_CLASSIC'
            assert target.cmor_version==run['cmor_version']
            assert target.frequency=='fx' and target.variable_id=='cSoil'
            assert 'time' not in target.variables and 'time' not in target.dimensions
            assert target['cSoil'].dimensions==('depth','lat','lon')
            assert target['cSoil'].units=='kg m-2'
            assert target['cSoil'].standard_name=='soil_mass_content_of_carbon'
            assert target['cSoil'].dtype==np.dtype('float32')
            assert target.has_aux_unc=='FALSE'
            assert target.variant_label==variants[depth]
            assert path.name==f'cSoil_fx_AI-upscaling-cSoil-1-0_{variants[depth]}_gn.nc'
            missing=set(cv['required_global_attributes'])-set(target.ncattrs())
            assert not missing,missing
            tracking.append(str(uuid.UUID(target.tracking_id)))
            np.testing.assert_allclose(target['depth'][:],[depth/200.],atol=1e-12,rtol=0)
            np.testing.assert_allclose(target['depth_bnds'][:],[[0,depth/100.]],atol=1e-12,rtol=0)
            np.testing.assert_allclose(target['lat'][:],source['lat'][:],atol=1e-10,rtol=0)
            np.testing.assert_allclose(target['lat_bnds'][:],source['lat_bnds'][:],atol=1e-10,rtol=0)
            # CMOR normalizes the western-hemisphere longitude coordinates to 0..360.
            np.testing.assert_allclose(np.mod(target['lon'][:],360),np.mod(source['lon'][:],360),atol=1e-10,rtol=0)
            np.testing.assert_allclose(np.mod(target['lon_bnds'][:],360),np.mod(source['lon_bnds'][:],360),atol=1e-10,rtol=0)
            count=0
            for j in range(0,len(source.dimensions['lat']),128):
                expected=source['cSoil'][j:j+128,:]
                actual=target['cSoil'][0,j:j+128,:]
                np.testing.assert_array_equal(np.ma.getmaskarray(actual),np.ma.getmaskarray(expected))
                np.testing.assert_array_equal(actual.compressed(),expected.compressed())
                count+=int(actual.count())
            log=(output_dir/f'cmor_{depth}cm.log').read_text()
            assert 'All files were closed successfully' in log
            assert 'Error:' not in log and 'Warning:' not in log,log
            with path.open('rb') as f:digest=hashlib.file_digest(f,'sha256').hexdigest()
            reports.append(dict(depth_cm=[0,depth],file=str(path),
                shape=list(target['cSoil'].shape),valid_cells=count,sha256=digest,
                all_values_and_masks='PASS: exact equality',coordinates_and_bounds='PASS',
                cmor_log='PASS: no errors or warnings',required_attributes='PASS',
                tracking_id=target.tracking_id,region=target.region,
                cmor_completed=True))
    assert len(set(tracking))==len(tracking)
    (output_dir/'validation.json').write_text(json.dumps(reports,indent=2)+'\n')
    print(json.dumps(reports,indent=2))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--input-dir',type=Path,required=True)
    p.add_argument('--output-dir',type=Path,required=True)
    a=p.parse_args()
    validate(a.input_dir,a.output_dir)
