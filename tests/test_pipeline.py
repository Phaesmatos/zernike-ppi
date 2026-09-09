from pathlib import Path
import numpy as np
from zernike_ppi.pdb import read_pdb
from zernike_ppi.surface import approximate_surface
from zernike_ppi.projection import local_frame, project_patch
from zernike_ppi.pipeline import build_descriptors, save_descriptor_bundle, load_descriptor_bundle
from zernike_ppi.matching import match_descriptors
from zernike_ppi.config import load_config

FIX=Path(__file__).parent/'data/fixture_a.pdb'
def test_pdb_surface_and_frame():
 s=read_pdb(FIX); assert len(s['xyz'])==4
 surf=approximate_surface(s); assert len(surf['xyz'])>0 and surf['normal'].shape[1]==3
 x,y,z=local_frame([0,0,1]); assert np.allclose(np.cross(x,y),z)
 p=project_patch(surf,0,radius=6,grid_size=16); assert p['fields']['shape'].shape==(16,16)
def test_bundle_roundtrip_and_config():
 import shutil
 root=Path('tests/_tmp_bundle')
 try:
  ds,_=build_descriptors(FIX,names=['shape','charge'],sample_every=20,grid_size=16,order=4)
  save_descriptor_bundle(root,ds,{'order':4}); loaded=load_descriptor_bundle(root,['shape','charge'])
  assert len(loaded)==len(ds) and np.allclose(loaded[0].channels['shape'],ds[0].channels['shape'])
 finally:
  shutil.rmtree(root,ignore_errors=True)
 cfg=load_config('examples/configs/shape_charge.yaml'); assert cfg['patches']['radius']==6.0
def test_match_weight_and_disabled_channel():
 ds,_=build_descriptors(FIX,names=['shape','charge'],sample_every=50,grid_size=16,order=4)
 r=match_descriptors(ds,ds,['shape'],top_k=1)[0]; assert 'score_shape' in r and 'score_charge' not in r
