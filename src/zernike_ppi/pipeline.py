import json
from pathlib import Path
import numpy as np
from .pdb import read_pdb, read_surface_csv
from .surface import approximate_surface
from .channels import get_channels
from .projection import project_patch
from .descriptors import describe_patches, PatchDescriptor
from .cache import cache_key, save_descriptors
from .zernike import zernike_descriptor

def load_surface(path):
    p=Path(path)
    if p.suffix.lower()=='.csv': return read_surface_csv(p)
    if p.suffix.lower()=='.npz':
        z=np.load(p,allow_pickle=True); return {k:z[k] for k in z.files}
    return approximate_surface(read_pdb(p))

def surface_from_input(protein=None, surface=None):
    if surface: return load_surface(surface)
    if not protein: raise ValueError('provide --protein or --surface')
    return load_surface(protein)

def add_channel_fields(s, names):
    channels=get_channels(names)
    for n,c in channels.items():
        if n=='shape':
            # Height is represented by the signed local normal coordinate in projection.
            s.setdefault('height', np.zeros(len(s['xyz'])))
        else: s.setdefault('fields',{})[n]=c.values(s)
    return s

def build_descriptors(protein=None, surface=None, names=('shape',), radius=6., sample_every=10, order=20, grid_size=64):
    s=add_channel_fields(surface_from_input(protein,surface), names)
    s['fields']=dict(s.get('fields',{})); s['fields']['shape']=s.get('height',np.zeros(len(s['xyz'])))
    idx=range(0,len(s['xyz']),max(1,int(sample_every))); patches=[]
    for i in idx:
        p=project_patch(s,i,radius,grid_size,1); p['metadata']={'residue_name':str(s.get('residue_name',[''])[i]),'residue_id':str(s.get('residue_id',[''])[i]),'chain':str(s.get('chain',[''])[i])}; p['fields'].update({n:p['fields'].get(n,np.zeros((grid_size,grid_size))) for n in names}); patches.append(p)
    return describe_patches(patches,list(names),order,grid_size), s

def save_surface(path,s):
    p=Path(path); p.parent.mkdir(parents=True,exist_ok=True); np.savez(p,**{k:v for k,v in s.items() if isinstance(v,np.ndarray)})

def save_descriptor_bundle(path, descriptors, params):
    p=Path(path); save_descriptors(p,descriptors,params)
    meta=[]
    for d in descriptors: meta.append({'patch_id':int(d.patch_id),'center':d.center.tolist(),'normal':d.normal.tolist(),'metadata':d.metadata})
    (p/'patches.json').write_text(json.dumps(meta,indent=2)); (p/'parameters.json').write_text(json.dumps(params,indent=2))

def load_descriptor_bundle(path,names):
    p=Path(path); meta=json.loads((p/'patches.json').read_text()); z=np.load(p/'channels.npz'); out=[]
    for i,m in enumerate(meta): out.append(PatchDescriptor(m['patch_id'],np.array(m['center']),np.array(m['normal']),{n:z[f'{i}_{n}'] for n in names},m.get('metadata',{})))
    return out
