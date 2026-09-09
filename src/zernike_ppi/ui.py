"""Thin UI-facing adapters; all scientific work remains in existing modules."""
from dataclasses import dataclass, asdict
from pathlib import Path
from tempfile import TemporaryDirectory
import json
import numpy as np
from .pdb import read_pdb
from .surface import approximate_surface
from .interface_map import run_interface_map, write_interface_map

@dataclass(frozen=True)
class InterfaceMapConfig:
    interface_distance: float = 3.0
    patch_radius: float = 6.0
    sample_every: int = 10
    zernike_order: int = 20
    grid_size: int = 64
    grid_spacing: float = .5
    axis_max_distance: float | None = None

def chains(atoms): return sorted(str(c) for c in set(atoms['chain']) if str(c))
def extract_chain(atoms, chain):
    if chain not in chains(atoms): raise ValueError(f"Chain '{chain}' is not present in this PDB.")
    mask=np.asarray(atoms['chain'])==chain
    return {k:(v[mask] if isinstance(v,np.ndarray) and len(v)==len(mask) else v) for k,v in atoms.items()}
def read_pdb_upload(path): return read_pdb(path)
def analyse_atoms(atoms_a, chain_a, atoms_b, chain_b, config: InterfaceMapConfig, output_dir):
    aa=extract_chain(atoms_a,chain_a); bb=extract_chain(atoms_b,chain_b)
    if not len(aa['xyz']) or not len(bb['xyz']): raise ValueError('The selected chain has no usable atoms.')
    sa=approximate_surface(aa); sb=approximate_surface(bb)
    rows,meta=run_interface_map(sa,sb,config.interface_distance,config.patch_radius,config.sample_every,config.zernike_order,config.grid_size,config.grid_spacing,config.axis_max_distance)
    write_interface_map(output_dir,rows,meta,asdict(config),sa,sb)
    return rows,meta,sa,sb

def result_files(output_dir):
    names=['patch_pairs.csv','projected_complementarity.csv','complementarity_matrix.csv','interface_geometry.json','config.yaml','complementarity_map.png','interface_geometry.png']
    return {n:(Path(output_dir)/n).read_bytes() for n in names}
