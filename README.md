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

## Graphical interface

Install the optional local interface and start it:

```bash
pip install -e ".[ui]"
streamlit run app.py
```

Upload either one bound complex PDB (then choose two chains) or two PDB files, select the interacting chains, and click **Run analysis**. The page shows the heatmap, 3D diagnostic, patch-pair table, scientific diagnostics, and download controls. It uses the automatic surface generator by default, so DMS is not required.

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

The original-compatible defaults are a 6 Å patch radius, order 20 (121 components), every-tenth-point sampling, orientation-aware shape complementarity corresponding to `verso=+1/-1`, and Euclidean descriptor distance. See [docs/original_method.md](docs/original_method.md).

## What is implemented today

The CLI supports PDB/DMS/NPZ surface input, approximate solvent-accessible surface generation, shared-frame patch projection, descriptor bundle serialization, cache reuse through `--descriptors-a/--descriptors-b`, weighted multichannel matching, CSV/JSON result export, ranking-only benchmark manifests, and patch visualization. The fast `demo` uses the same descriptor and matcher APIs as production workflows.

For a real workflow:

```bash
zernike-ppi match --protein-a A.pdb --protein-b B.pdb --channels shape,charge --top-k 20 --output results/
zernike-ppi visualize --protein-a A.pdb --protein-b B.pdb --matches results/matches.csv --rank 1 --output results/match.png
```

Descriptor bundles produced by `describe` can be reused without recomputation:

```bash
zernike-ppi match --protein-a A.pdb --protein-b B.pdb --descriptors-a descriptors/A --descriptors-b descriptors/B --channels shape,charge --output results/
```

## Limitations and roadmap

Scores are candidate interface complementarity, not binding probabilities or affinity. The MVP ignores flexibility, solvent mediation, side-chain rearrangement, entropy, and full 3D topology. Next, add PDB/DMS CLI workflows, YAML-resolved ablations, ROC/PR benchmark reports, and a validated `zepyros` backend; later learn channel weights and use complex coefficient phases for coarse docking.

## Attribution

Please cite the original Zernike2D repository and Milanetti et al. (2021) when using this continuation. The original student/professor project remains the scientific and software baseline.

## Geometric interface complementarity map

`interface-map` analyzes an already bound protein pair, rather than attempting blind docking. It uses DMS-style surfaces (`x,y,z,nx,ny,nz`) to identify points closer than 3 Å, forms 6 Å local patches, finds the facing B patch along A's outward normal, evaluates shape-channel Zernike complementarity, and projects the patch-pair scores onto a global A→B interface plane. The plane is an analysis coordinate system; it does not assume that the physical interface is flat.

```bash
zernike-ppi interface-map --surface-a chain_A.csv --surface-b chain_B.csv --interface-distance 3 --patch-radius 6 --sample-every 10 --grid-spacing 0.5 --output interface_map/
```

The output includes interface masks, paired patch diagnostics, a complementarity matrix, interface-axis metadata, and 2D/3D figures. Points with similarly oriented opposing normals are flagged as suspicious rather than silently accepted.




## Why Zernike?

A natural question is: **why use 2D Zernike polynomials instead of another basis?**

The main reason is that local protein surface patches are projected onto a **unit disk**, and Zernike polynomials form an orthogonal basis defined exactly on that domain. In addition, under an in-plane rotation by an angle \(\alpha\), the Zernike coefficients transform as

$$
A'_{nm} = A_{nm} e^{-im\alpha},
$$

so their magnitudes remain unchanged:

$$
|A'_{nm}| = |A_{nm}|.
$$

This provides a compact and naturally rotation-invariant representation of local surface patches, which is particularly convenient when many candidate protein interfaces must be compared efficiently.

Other representations are possible — for example Fourier-Bessel functions, pseudo-Zernike moments, Legendre moments, spherical harmonics, or 3D Zernike descriptors — but they come with different domain assumptions, computational costs, and invariance properties. Zernike polynomials are therefore not the only possible choice, but they are a particularly natural one for **local, disk-mapped surface patches**.

### A small historical note

This question has a slightly personal motivation. During the oral discussion of the original Zernike2D student project, I was asked why Zernike polynomials had been chosen over alternative representations. At the time, I had largely accepted that methodological choice as given and could not provide an answer.

This continuation of the project is also an opportunity to revisit that question explicitly: to understand not only **how** the method works, but **why this representation is appropriate**, where its advantages come from, and where alternative bases might eventually perform better.

In that sense, this section closes a small unfinished question from the original project.
