"""Compare scripts/popstats.py to tskit truth on msprime loci with varying n.

Ten independent 1 kb loci, each simulated with a different sample size n_w.
They are concatenated into one 10 kb alignment with max(n) rows; rows beyond
n_w are 'N' in that locus, mimicking indel-driven missing data. Within each
window n is constant, so tskit's per-locus statistics are the exact truth for
the script's window output (including Tajima's D).

Run with pytest, or directly:  python tests/test_popstats_msprime.py
"""
import subprocess
import sys
import tempfile
from pathlib import Path

import msprime
import numpy as np
from scipy.stats import hypergeom

SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "popstats.py"
L = 1000
NS = [8, 10, 12, 15, 18, 20, 24, 28, 32, 40]
X = 8  # = min(NS): SFS downsampling target and --min-seqs


def simulate(seed=42):
    loci = []
    for w, n in enumerate(NS):
        ts = msprime.sim_ancestry(
            samples=n, ploidy=1, sequence_length=L, population_size=1e4,
            recombination_rate=1e-8, random_seed=seed + w)
        # binary model -> strictly biallelic sites, so the folded SFS is well defined
        ts = msprime.sim_mutations(ts, rate=5e-7, model=msprime.BinaryMutationModel(),
                                   random_seed=seed + 100 + w)
        loci.append(ts)
    return loci


def to_fasta(loci, path):
    nmax = max(ts.num_samples for ts in loci)
    mat = np.full((nmax, L * len(loci)), "N", dtype="<U1")
    for w, ts in enumerate(loci):
        mat[: ts.num_samples, w * L:(w + 1) * L] = "A"
        for var in ts.variants():
            col = w * L + int(var.site.position)
            mat[: ts.num_samples, col] = np.where(var.genotypes == 0, "A", "T")
    with open(path, "w") as fh:
        for i in range(nmax):
            fh.write(f">s{i}\n{''.join(mat[i])}\n")


def truth(ts):
    n = ts.num_samples
    a1 = (1.0 / np.arange(1, n)).sum()
    S = ts.segregating_sites(mode="site", span_normalise=False)
    t = {
        "n_sites": L,
        "S": S,
        "pi": ts.diversity(mode="site"),  # per-site (span normalised)
        "theta_w": S / (L * a1),
        "tajima_d": ts.Tajimas_D(mode="site"),
    }
    afs = ts.allele_frequency_spectrum(mode="site", polarised=True,
                                       span_normalise=False)
    proj = np.zeros(X + 1)
    for k in range(1, n):  # derived count k of n -> X chromosomes, w/o replacement
        proj += afs[k] * hypergeom(n, k, X).pmf(np.arange(X + 1))
    folded = [proj[0] + proj[X]] + [
        proj[j] + (proj[X - j] if X - j != j else 0.0) for j in range(1, X // 2 + 1)]
    t["sfs"] = np.array(folded)
    return t


def read_tsv(path):
    rows = [l.split("\t") for l in Path(path).read_text().splitlines()]
    return [dict(zip(rows[0], r)) for r in rows[1:]]


def run_popstats(fa, out, *extra):
    subprocess.run([sys.executable, str(SCRIPT), str(fa), "--out", str(out),
                    "--window", str(L), "--overlap", "0", *extra],
                   check=True, stderr=subprocess.DEVNULL)
    return read_tsv(f"{out}.stats.tsv"), read_tsv(f"{out}.sfs.tsv")


def test_windows_match_tskit(tmp_path):
    loci = simulate()
    fa = tmp_path / "sim.fa"
    to_fasta(loci, fa)
    stats, sfs = run_popstats(fa, tmp_path / "o", "--min-seqs", str(X))
    assert len(stats) == len(NS)
    assert sum(int(r["S"]) for r in stats) > 100  # simulation isn't trivially empty
    for ts, st, sf in zip(loci, stats, sfs):
        t = truth(ts)
        assert int(st["n_sites"]) == L
        assert int(st["S"]) == int(t["S"])
        for key in ("pi", "theta_w", "tajima_d"):
            assert np.isclose(float(st[key]), t[key], rtol=1e-5, atol=1e-9), (ts.num_samples, key)
        got = np.array([float(sf["proj_invariant"])] +
                       [float(sf[f"k{j}"]) for j in range(1, X // 2 + 1)])
        assert np.allclose(got, t["sfs"], atol=1e-5), (ts.num_samples, got, t["sfs"])


def test_min_seqs_filters_windows(tmp_path):
    loci = simulate()
    fa = tmp_path / "sim.fa"
    to_fasta(loci, fa)
    stats, _ = run_popstats(fa, tmp_path / "o", "--min-seqs", "20", "--sfs-n", "8")
    for ts, st in zip(loci, stats):
        assert int(st["n_sites"]) == (L if ts.num_samples >= 20 else 0)


if __name__ == "__main__":
    for f in (test_windows_match_tskit, test_min_seqs_filters_windows):
        f(Path(tempfile.mkdtemp()))
        print(f.__name__, "ok")
