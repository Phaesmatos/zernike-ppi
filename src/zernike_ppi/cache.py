import hashlib, json
from pathlib import Path
import numpy as np
def cache_key(parameters): return hashlib.sha256(json.dumps(parameters,sort_keys=True,default=str).encode()).hexdigest()[:16]
def save_descriptors(path, descriptors, parameters):
 p=Path(path); p.mkdir(parents=True,exist_ok=True); np.savez(p/"channels.npz",**{f"{i}_{n}":v for i,d in enumerate(descriptors) for n,v in d.channels.items()}); (p/"metadata.json").write_text(json.dumps({"parameters":parameters,"count":len(descriptors)},indent=2))
