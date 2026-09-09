# Original method and continuity

ZernikePPI is a continuation of `matmi8/Zernike2D`. The upstream protocol reads DMS molecular-surface CSV data, extracts 6 Å local patches, projects oriented surface geometry to a disk, computes order-20 2D Zernike invariants (121 values), compares patches with Euclidean distance, and labels interfaces using a roughly 3 Å surface separation. The upstream scripts sample every tenth point for the expensive binding-propensity calculation and expose `verso=+1/-1` orientation.

This MVP preserves those defaults conceptually while isolating descriptor generation behind `ZernikeBackend`. Multichannel scalar fields, property-specific complementarity, cross donor/acceptor matching, vectorized matching, and descriptor caching are new extensions. The bundled backend is numerically self-contained; `zepyros` can be added as a future drop-in backend after validating its exact API against a pinned version.
