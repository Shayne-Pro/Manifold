# Grassmann Manifold Image-Set Classification

MATLAB research scripts for *Convex Hull Collaborative Representation Learning on Grassmann Manifold for Image Set Classification*.

The implementation forms local image-set models, extracts subspaces, builds Grassmann projection kernels, and performs collaborative-representation classification.

## Files

| Path | Purpose |
| --- | --- |
| `search_extendedyaleB.m` | Example parameter selection and experiment entry |
| `ten_HDC_GLSC.m` | Load the supplied split indices, build train/test sets, and report accuracy statistics |
| `HDC_local_model.m`, `nonlinearity_score.m`, `Geodisic.m` | Local model construction and neighborhood/geodesic calculations |
| `Compute_Subspace.m`, `grassmann_proj.m` | Subspace extraction and projection kernels |
| `GLSC_func.m`, `Classify_ISCRC.m` | Representation and classification |
| `Extended_yaleB.mat`, `random_index_extyaleB.mat` | Historical experiment data and split indices |

## Run in MATLAB

Run from the repository root so relative `.mat` file loads resolve. The code uses `pdist` and `squareform`, which require the Statistics and Machine Learning Toolbox.

```matlab
cd('/path/to/Manifold');
[accuracy, variance] = ten_HDC_GLSC(10, 2, 3);       % one split
[mean_accuracy, variance] = ten_HDC_GLSC(10, 2, 3, 10); % ten supplied splits
```

Arguments are subspace dimension, local-model count, neighbor count, and optional trial count. Trials must not exceed the number of supplied split-index sets. Accuracy is returned as a fraction; variance is calculated over the trials actually executed. A single trial therefore has zero sample variance.

The historical function name includes `ten`, but the old body only executed one split while averaging a ten-element array with unused zero entries. The optional trial count makes the execution and summary consistent. Algorithm parameters and the supplied data are otherwise preserved.

## Reproduction status

MATLAB execution and classification accuracy have not been validated on a new runtime. Comments in several legacy `.m` files use a historical non-UTF-8 encoding; avoid rewriting them as UTF-8 without checking their original encoding. No current accuracy claim is inferred from old comments or the filename.

The bundled data's original terms still apply. This repository currently contains no repository-level license grant or complete publication citation; retain known source attribution when adding either.
