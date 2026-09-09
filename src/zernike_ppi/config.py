"""Reproducible YAML experiment configuration helpers."""
from pathlib import Path
import yaml
DEFAULTS={"patches":{"radius":6.0,"sample_every":10},"zernike":{"order":20,"grid_size":64}}
def load_config(path):
    data=yaml.safe_load(Path(path).read_text()) or {}
    out={**DEFAULTS,**data}; out["patches"]={**DEFAULTS["patches"],**data.get("patches",{})}; out["zernike"]={**DEFAULTS["zernike"],**data.get("zernike",{})}; return out
