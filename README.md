# ZernikePPI

ZernikePPI is an experimental continuation of [matmi8/Zernike2D](https://github.com/matmi8/Zernike2D), extending local 2D Zernike representations of protein surfaces from geometric shape complementarity toward configurable multichannel physicochemical descriptors. It does not claim invention of the underlying 2D Zernike method. The continuation is motivated by Milanetti et al. (2021), [doi:10.1016/j.csbj.2020.11.051](https://doi.org/10.1016/j.csbj.2020.11.051).

## Quick start

```bash
python -m venv .venv
pip install -e .[test]
zernike-ppi channels
zernike-ppi demo --channels shape
zernike-ppi demo --channels shape,charge,hydrophobicity
pytest
```

The MVP keeps channel descriptors separate and combines normalized Euclidean distances with user weights. Lower scores indicate better complementarity.

| Channel | Surface field | Expected relation |
|---|---|---|
| shape | local surface height | convex ↔ concave |
| hydrophobicity | Kyte–Doolittle residue scale | similar ↔ similar |
| charge | qualitative residue charge proxy (HIS=+0.25) | positive ↔ negative |
| donor | annotated donor propensity | donor ↔ acceptor |
| acceptor | annotated acceptor propensity | acceptor ↔ donor |

The donor and acceptor channels are intentionally separate. The charge channel is not Poisson–Boltzmann electrostatics; residue annotations and hydrophobicity are approximations.

## Scientific pipeline

PDB or DMS CSV → surface points and normals → shared local patch frame → one registered 2D disk per channel → independent order-20 Zernike invariant vectors → property-specific complementarity → interpretable weighted ranking.

The original-compatible defaults are a 6 Å patch radius, order 20 (121 components), every-tenth-point sampling, orientation-aware shape complementarity corresponding to `verso=+1/-1`, and Euclidean descriptor distance. See [docs/original_method.md](docs/original_method.md). This repository currently exposes a compact synthetic CLI demo; production PDB/DMS matching and benchmark orchestration are the next implementation milestone.

## Limitations and roadmap

Scores are candidate interface complementarity, not binding probabilities or affinity. The MVP ignores flexibility, solvent mediation, side-chain rearrangement, entropy, and full 3D topology. Next, add PDB/DMS CLI workflows, YAML-resolved ablations, ROC/PR benchmark reports, and a validated `zepyros` backend; later learn channel weights and use complex coefficient phases for coarse docking.

## Attribution

Please cite the original Zernike2D repository and Milanetti et al. (2021) when using this continuation. The original student/professor project remains the scientific and software baseline.
