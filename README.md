# Non-Binary Detector Error Models — Stim Patch

This repository contains a patch to [Stim](https://github.com/quantumlib/Stim)
implementing **non-binary detector error models (NB-DEM)** as described in
*Non-Binary Detector Error Models* (Tinelli, Lázaro, García Herrero).

In a non-binary DEM every fault location of the circuit is one random variable over
its native alphabet F_2^r (r = 1, 2, 4 for bit flips, single-qubit and two-qubit Pauli
channels) with its exact prior, instead of being split into independent binary
mechanisms. The patch builds a compact non-binary DEM from any Clifford circuit with
the split-merge-assign-regroup construction of the paper.

## Contents

| file | description |
|---|---|
| `nb_dem.patch` | the Stim patch (git diff, applies with `git apply`) |
| `setup.sh` | clones Stim at the tested commit, applies the patch and builds the CLI |
| `scripts/validate.py` | checks the NB-DEM against Stim's DEM and prints the Tanner graph statistics |
| `circuits/*.stim` | the circuits of the paper at p = 10^-3 (surface and color code d = 5, [[72,12,6]] BB code) |
| `LICENSE` | Apache License 2.0 (same as Stim) |

## Quick start

```bash
./setup.sh                                   # clones Stim into ./Stim, applies the patch, builds
./Stim/build/out/stim analyze_errors_nb --in circuits/surface_d5_r5_p1e-3.stim
python3 scripts/validate.py ./Stim/build/out/stim
```

The patch was tested on Stim commit `ffe0d043598ed5df0851a3e6b4920359a23fafbd`
(the one checked out by `setup.sh`).

## What the patch adds

### Construction (`src/stim/util_top/circuit_to_nb_dem.{h,cc}`)

`circuit_to_nonbinary_dem(circuit, options)` implements the construction of Sec. III of
the paper:

1. **Locations.** The circuit is walked backward with Stim's error analyzer as a
   sensitivity tracker (loops are unrolled). Every application of a noise channel is
   recorded as *one* column with its exact channel prior: single-qubit Pauli channels
   with r = 2 (generators `X@q`, `Z@q`), two-qubit channels with r = 4, bit and
   measurement flips with r = 1, `E`/`ELSE` chains as one column with one generator per
   Pauli component, heralded channels with r = 3 (the herald record is the third
   generator).
2. **Reduce.** Each column is restricted to the span of its support. Columns that flip
   no detector and no observable are dropped.
3. **Split.** Columns whose law factorizes into independent binary parts (solved with the
   character formula `log(1-2q_v) = -2^{1-r} Σ_a (-1)^{<a,v>} log χ_a`, with a numerical
   check) are split into their parts. Columns whose law does not factorize (e.g. a
   single-qubit channel without `Y` component, and in general `PAULI_CHANNEL_2`, `E`/`ELSE`
   chains and heralded channels) are kept whole with their exact prior.
4. **Merge.** Parts with equal symptom (detectors and observables) are merged by XOR of
   their probabilities, exactly as in Stim's binary DEM.
5. **Assign.** Every merged part goes to its origin location of largest rank, ties
   broken by circuit order in favor of the earlier location.
6. **Regroup.** Every location is rebuilt from its assigned parts: its basis spans the
   parts and its prior is the XOR-convolution of the parts. Locations without parts are
   dropped.

For channels that factorize (bit flips, depolarizing channels), the NB-DEM and Stim's
binary DEM have the same distribution of detectors and observables, and the binary
image of the NB-DEM (its merged parts) is Stim's DEM. Columns kept whole are exact also
when the channel does not factorize, whereas Stim's binary DEM approximates them (test
`nb_dem.non_factorizable_channel_kept_whole`).

### Edge labels

Column c is connected to detector d with a label h in F_2^r, and d is flipped by
u_c iff <h, u_c> = 1. When the basis of a column is the native basis of a Pauli channel
(`X@q, Z@q` or `X@a, Z@a, X@b, Z@b`), the label is printed as a Pauli operator in the
Kuo–Lai convention: the detector is flipped iff the error anticommutes with the label.
For instance, a detector with label `Z` is flipped by `X` and `Y` errors. Other bases (e.g.
a two-qubit channel reduced to `X@a, X@b`) print the label as a bit string over the basis,
bit i being the coefficient of generator i. The Python API returns the labels as integer
bit masks.

### Interfaces

- **C++**: `NonbinaryDem circuit_to_nonbinary_dem(const Circuit &, NbDemOptions)`
  with `NbDemMergeMode::{PARTS, COLUMNS, NONE}` (`PARTS` is the construction above,
  `COLUMNS` only merges columns with equal span, `NONE` keeps one column per channel
  application), plus `NonbinaryDem::str()`, `to_json()`, `binary_image()` (the merged
  parts as a `DetectorErrorModel`, in Stim's canonical error order) and
  `num_detector_edges()`.
- **CLI**: `stim analyze_errors_nb [--merge=parts|columns|none] [--format=text|json]
  [--allow_gauge_detectors] [--in=f.stim] [--out=f.txt]`
- **Python**: `circuit.nonbinary_detector_error_model(merge="parts",
  allow_gauge_detectors=False) -> stim.NonbinaryDem` with `.to_dict()`,
  `.to_json()`, `.binary_image()`, `.column(i)`, `.num_columns()`,
  `.num_detector_edges()`, `__str__` and `__len__`.

Supported noise: `X_ERROR`, `Y_ERROR`, `Z_ERROR`, `DEPOLARIZE1`, `DEPOLARIZE2`,
`PAULI_CHANNEL_1`, `PAULI_CHANNEL_2`, `E`/`ELSE_CORRELATED_ERROR` chains, measurement
flips (`M/MX/MY/MR/MRX/MRY/MXX/MYY/MZZ` with a parens argument), `HERALDED_ERASE` and
`HERALDED_PAULI_CHANNEL_1`. Noise arguments of `MPP`/`SPP` raise an error. Columns of
rank up to 16 are supported.

## Building by hand

```bash
git clone https://github.com/quantumlib/Stim.git
cd Stim
git checkout ffe0d043598ed5df0851a3e6b4920359a23fafbd
git apply /path/to/nb_dem.patch
mkdir build && cd build
cmake .. -DCMAKE_BUILD_TYPE=Release
make stim -j4                                  # CLI: build/out/stim
```

- **Tests** need GTest (`find_package(GTest)`, e.g. `-DCMAKE_PREFIX_PATH=/path/to/gtest`):
  `make stim_test -j4`, then run `./build/out/stim_test --gtest_filter='nb_dem.*'`
  from the Stim root directory (some of Stim's own tests read reference files relative to it).
- **Python module** needs pybind11 2.x (Stim does not compile against pybind11 3.x):
  `cmake .. -Dpybind11_DIR=$(python3 -c "import pybind11; print(pybind11.get_cmake_dir())")`,
  then `make stim_python_bindings -j4`, which produces `build/out/stim*.so`.

## Validation

`python3 scripts/validate.py ./Stim/build/out/stim` runs the checks on the three circuits
of the paper (`circuits/`, uniform circuit-level noise p = 10^-3, d rounds for the surface
and color codes, 6 rounds for the BB code):

| | surface d=5 | color d=5 | BB [[72,12,6]] |
|---|---|---|---|
| binary DEM: VN / edges / 4-cycles | 1677 / 5004 / 26046 | 2194 / 8617 / 186974 | 16164 / 91116 / 2628756 |
| NB-DEM: VN / edges / 4-cycles | 284 / 1156 / 2672 | 327 / 1908 / 21089 | 2016 / 17352 / 199332 |
| NB-DEM ranks 1/2/3/4 | 39/97/28/120 | 72/59/0/196 | 36/432/0/1548 |
| binary image == Stim DEM (max \|Δp\|) | ✓ (7e-16) | ✓ (5e-16) | ✓ (3e-16) |
| max \|E[(-1)^⟨f,(s,ℓ)⟩] NB − Stim DEM\|, 200 masks | 4e-14 | 4e-14 | 3e-14 |

The parity check computes the NB-DEM expectations from the printed edge labels, so it
also checks the label convention.

## Tests

- `src/stim/util_top/circuit_to_nb_dem.test.cc`: 18 gtest cases covering merging,
  reduction (invisible directions), measurement flips, `MXX` flips, factorizing and
  non-factorizable `PAULI_CHANNEL`s, `E`/`ELSE` chains, loops, the [[4,2,2]] example of the
  paper, merge modes, output formats, the label convention and the Assign rule. Each
  exactness case compares random parity expectations against Stim's DEM.
- Stim's full test suite passes with the patch (1911/1911).

## Circuits

The surface code circuit is Stim's `surface_code:rotated_memory_z`. The color code
circuit uses the geometry and CNOT schedule of Stim's `color_code:memory_xyz` with one Z-
and one X-ancilla per plaquette. The BB code circuit uses the depth-8 syndrome cycle of
Bravyi et al. and was generated with the circuit builder of Gong, Cammerer and Renes
([SlidingWindowDecoder](https://github.com/gongaa/SlidingWindowDecoder), MIT license).
Noise: two-qubit depolarizing after every CNOT, single-qubit depolarizing on the data
qubits before every round and after every Hadamard, bit flips after resets and before
measurements, all with parameter p.

## Citation

If you use this code, please cite *Non-Binary Detector Error Models* by S. Tinelli,
F. Lázaro and F. García Herrero (reference to be added).

## License

Apache License 2.0, see `LICENSE`. Stim is also distributed under the Apache License 2.0.
