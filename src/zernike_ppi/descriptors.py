from dataclasses import dataclass, field
import numpy as np
from .zernike import zernike_descriptor
@dataclass
class PatchDescriptor:
    patch_id: int; center: np.ndarray; normal: np.ndarray
    channels: dict[str,np.ndarray]; metadata: dict = field(default_factory=dict)
def describe_field(field, order=20, grid_size=64): return zernike_descriptor(field,order,grid_size)
def describe_patches(patches, channel_names, order=20, grid_size=64):
    return [PatchDescriptor(p["patch_id"],np.asarray(p["center"]),np.asarray(p["normal"]),{n:describe_field(p["fields"][n],order,grid_size) for n in channel_names},p.get("metadata",{})) for p in patches]
