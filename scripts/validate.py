"""Validates the non-binary DEM produced by the patched Stim against Stim's binary DEM.

For every circuit it runs `stim analyze_errors_nb --format json` and checks that
  1. the binary image (the merged parts of all columns) equals Stim's detector error model,
  2. the non-binary DEM has the same law of (detectors, observables) as Stim's DEM, through the
     parity expectations E[(-1)^<f,(s,l)>] on random masks f (computed from the rendered edge
     labels, so this also checks the label convention),
and prints the Tanner graph statistics of both models (variable nodes, edges, 4-cycles, ranks).

Usage: python3 scripts/validate.py path/to/patched/stim [circuit.stim ...]
       (default circuits: circuits/*.stim). Needs numpy and the official `stim` Python package.
"""
import glob
import json
import os
import subprocess
import sys
from collections import Counter

import numpy as np
import stim


def render_to_mask(label, generator_labels, rank):
    """Inverse of the patch's label rendering: Pauli label X^hz Z^hx per qubit, or a bit string."""
    if len(label) == rank and set(label) <= {"0", "1"}:
        return sum(int(b) << i for i, b in enumerate(label))
    mask = 0
    for q, ch in enumerate(label):
        hx = ch in "ZY"  # coefficient of the X generator of qubit q
        hz = ch in "XY"  # coefficient of the Z generator of qubit q
        mask |= (hx << (2 * q)) | (hz << (2 * q + 1))
    return mask


def four_cycles(supports, nd):
    """Number of 4-cycles of a Tanner graph given the detector support of every variable node."""
    pair = Counter()
    for sup in supports:
        sup = sorted(sup)
        for i in range(len(sup)):
            for j in range(i + 1, len(sup)):
                pair[sup[i] * nd + sup[j]] += 1
    return sum(k * (k - 1) // 2 for k in pair.values())


def validate(stim_bin, path, trials=200, seed=0):
    circuit = stim.Circuit.from_file(path)
    out = subprocess.run([stim_bin, "analyze_errors_nb", "--format", "json", "--in", path],
                         capture_output=True, text=True, check=True).stdout
    nb = json.loads(out)
    nd, no = nb["num_detectors"], nb["num_observables"]
    dem = circuit.detector_error_model().flattened()
    errors = []
    for inst in dem:
        if inst.type == "error":
            t = inst.targets_copy()
            errors.append((tuple(sorted(x.val for x in t if x.is_relative_detector_id())),
                           tuple(sorted(x.val for x in t if x.is_logical_observable_id())), inst.args_copy()[0]))

    # 1. binary image == Stim DEM
    parts = [(tuple(p["flipped_detectors"]), tuple(p["flipped_observables"]), p["probability"])
             for col in nb["columns"] for p in col["parts"]]
    key = lambda e: (e[0], e[1])
    a, b = sorted(parts, key=key), sorted(errors, key=key)
    same = len(a) == len(b) and all(key(x) == key(y) for x, y in zip(a, b))
    dp = max(abs(x[2] - y[2]) for x, y in zip(a, b)) if same else float("nan")

    # 2. parity expectations
    rng = np.random.default_rng(seed)
    cols = []
    for col in nb["columns"]:
        r, g = col["rank"], col["generator_labels"]
        det = {d: render_to_mask(lab, g, r) for d, lab in col["detector_edges"]}
        obs = {o: render_to_mask(lab, g, r) for o, lab in col["observable_edges"]}
        cols.append((np.array(col["prior"]), det, obs))
    worst = 0.0
    for _ in range(trials):
        dens = rng.choice([0.02, 0.05, 0.2])
        fd = set(np.flatnonzero(rng.random(nd) < dens).tolist())
        fo = set(np.flatnonzero(rng.random(no) < 0.5).tolist())
        e_dem = 1.0
        for ds, os_, p in errors:
            if (len(fd.intersection(ds)) + len(fo.intersection(os_))) & 1:
                e_dem *= 1 - 2 * p
        e_nb = 1.0
        for prior, det, obs in cols:
            t = 0
            for d, m in det.items():
                if d in fd:
                    t ^= m
            for o, m in obs.items():
                if o in fo:
                    t ^= m
            if t:
                sgn = np.array([1 - 2 * (bin(t & u).count("1") & 1) for u in range(len(prior))])
                e_nb *= float(prior @ sgn)
        worst = max(worst, abs(e_nb - e_dem))

    # Tanner graph statistics
    ranks = Counter(col["rank"] for col in nb["columns"])
    nb_edges = sum(len(col["detector_edges"]) for col in nb["columns"])
    nb_cyc = four_cycles([[d for d, _ in col["detector_edges"]] for col in nb["columns"]], nd)
    bin_edges = sum(len(e[0]) for e in errors)
    bin_cyc = four_cycles([e[0] for e in errors], nd)
    print(f"{os.path.basename(path)}: detectors={nd} observables={no}")
    print(f"  binary DEM     : VN={len(errors)} edges={bin_edges} 4-cycles={bin_cyc}")
    print(f"  non-binary DEM : VN={len(nb['columns'])} edges={nb_edges} 4-cycles={nb_cyc} "
          f"ranks 1/2/3/4={'/'.join(str(ranks.get(r, 0)) for r in (1, 2, 3, 4))}")
    print(f"  binary image == Stim DEM: {same} (max |dp| = {dp:.1e})")
    print(f"  parity expectations, max |NB - Stim DEM| over {trials} masks: {worst:.1e}")
    return same and worst < 1e-10


if __name__ == "__main__":
    if len(sys.argv) < 2:
        sys.exit(__doc__)
    here = os.path.dirname(os.path.abspath(__file__))
    paths = sys.argv[2:] or sorted(glob.glob(os.path.join(here, "..", "circuits", "*.stim")))
    ok = all([validate(sys.argv[1], p) for p in paths])
    print("ALL CHECKS PASSED" if ok else "SOME CHECKS FAILED")
    sys.exit(0 if ok else 1)
