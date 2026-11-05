#!/usr/bin/env python3
"""Pi, Tajima's D and the (folded) SFS from an aligned FASTA with missing data.

Input: aligned FASTA (all sequences same length). A/C/G/T (any case) are
called bases; everything else (N, -, ?, IUPAC codes, ...) is missing.

Missing-data handling
---------------------
* --min-seqs X : only sites with >= X called sequences are used (for every
  statistic; the retained-site count is the denominator for per-site values).
* pi: per site, the mean over the n_i *called* pairs, 1 - sum_a c_a(c_a-1)/(n_i(n_i-1)),
  summed over retained sites (monomorphic retained sites contribute 0).
  pi_per_site = sum / n_retained_sites.
* Watterson: S / sum_i a1(n_i), where a1(n)=sum_{k<n} 1/k, summed over all
  retained sites (so the per-site expectation uses each site's own n_i).
* Tajima's D: (pi_sum - S/A) / sqrt(e1*S + e2*S*(S-1)) using A, and the
  variance coefficients, averaged over retained sites (n -> mean n_i).
  This is an approximation when n_i varies; it is exact when n_i is constant.
* SFS: folded (no outgroup), biallelic sites only, each site projected from
  its n_i called sequences down to --sfs-n (<= --min-seqs) chromosomes, so
  every retained site contributes total weight 1 (expected counts). Projection
  is hypergeometric (sampling X of the n_i called chromosomes without
  replacement). Output bins are minor-allele counts 1..X//2; mass projected to 0 (site looks
  monomorphic after downsampling) is reported in the `proj_invariant` column.

Windows: --window W --overlap Y gives windows of W alignment columns that
advance by W-Y. Windows are in alignment coordinates (1-based in output).
Without --window the whole alignment is one window.

Outputs: <out>.stats.tsv and <out>.sfs.tsv
"""
import argparse
import sys
from math import comb

import numpy as np

CODES = np.full(256, 255, dtype=np.uint8)
for i, b in enumerate("ACGT"):
    CODES[ord(b)] = i
    CODES[ord(b.lower())] = i


def read_fasta(path):
    names, seqs, cur = [], [], None
    opener = open
    if path.endswith(".gz"):
        import gzip
        opener = lambda p: gzip.open(p, "rt")
    with opener(path) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if cur is not None:
                    seqs.append("".join(cur))
                names.append(line[1:].split()[0])
                cur = []
            else:
                if cur is None:
                    sys.exit("error: FASTA does not start with '>'")
                cur.append(line)
    if cur is not None:
        seqs.append("".join(cur))
    if not seqs:
        sys.exit("error: no sequences read")
    lens = {len(s) for s in seqs}
    if len(lens) != 1:
        sys.exit(f"error: sequences are not aligned (lengths {sorted(lens)[:5]}...)")
    arr = np.frombuffer("".join(seqs).encode("ascii"), dtype=np.uint8)
    arr = arr.reshape(len(seqs), -1)
    return names, CODES[arr]


def site_counts(codes):
    """(4, L) allele counts of called bases."""
    return np.stack([(codes == a).sum(axis=0) for a in range(4)]).astype(np.int64)


def a_coeffs(n):
    """a1, a2 for integer array n (>=2)."""
    maxn = int(n.max()) if n.size else 2
    k = np.arange(1, max(maxn, 2) + 1, dtype=float)
    a1 = np.concatenate([[0.0], np.cumsum(1.0 / k)])  # a1[m] = sum_{k<=m}
    a2 = np.concatenate([[0.0], np.cumsum(1.0 / k**2)])
    return a1[n - 1], a2[n - 1]  # a(n) = sum_{k<n}


def tajima_d(pi_sum, S, n):
    """n: array of per-site called counts at retained sites."""
    if S == 0 or n.size == 0:
        return float("nan")
    a1s, a2s = a_coeffs(n)
    a1, a2, nm = a1s.mean(), a2s.mean(), n.mean()
    b1 = (nm + 1) / (3 * (nm - 1))
    b2 = 2 * (nm**2 + nm + 3) / (9 * nm * (nm - 1))
    c1 = b1 - 1 / a1
    c2 = b2 - (nm + 2) / (a1 * nm) + a2 / a1**2
    e1 = c1 / a1
    e2 = c2 / (a1**2 + a2)
    var = e1 * S + e2 * S * (S - 1)
    if var <= 0:
        return float("nan")
    return (pi_sum - S / a1) / np.sqrt(var)


def proj_table(nmax, X):
    """P[n][k] -> vector over j=0..X of P(j | k of n), for 0<=k<=n, n<=nmax."""
    tab = {}
    for n in range(X, nmax + 1):
        rows = np.zeros((n + 1, X + 1))
        for k in range(n + 1):
            den = comb(n, X)
            for j in range(max(0, X - (n - k)), min(k, X) + 1):
                rows[k, j] = comb(k, j) * comb(n - k, X - j) / den
        tab[n] = rows
    return tab


def window_stats(cnt, lo, hi, min_seqs, X, tab):
    c = cnt[:, lo:hi]
    n = c.sum(axis=0)
    keep = n >= min_seqs
    c, n = c[:, keep], n[keep]
    nsites = int(keep.sum())
    out = {"n_sites": nsites}
    if nsites == 0:
        out.update(pi_sum=np.nan, pi=np.nan, S=0, theta_w=np.nan, tajima_d=np.nan)
        return out, np.full(X // 2 + 1, np.nan)
    pi_site = 1.0 - (c * (c - 1)).sum(axis=0) / (n * (n - 1))
    pi_sum = float(pi_site.sum())
    nalle = (c > 0).sum(axis=0)
    S = int((nalle > 1).sum())
    a1, _ = a_coeffs(n)
    out.update(
        pi_sum=pi_sum,
        pi=pi_sum / nsites,
        S=S,
        theta_w=S / a1.sum(),
        tajima_d=tajima_d(pi_sum, S, n),
        mean_n=float(n.mean()),
    )
    # SFS: biallelic sites, folded, projected to X
    bi = nalle == 2
    sfs = np.zeros(X + 1)
    sfs_inv = float((nalle == 1).sum())
    cb, nb = c[:, bi], n[bi]
    ks = np.sort(cb, axis=0)[-1]  # count of the more common allele
    kmin = nb - ks
    for nn in np.unique(nb):
        sel = nb == nn
        ks_ = np.bincount(kmin[sel], minlength=nn + 1)
        sfs += ks_ @ tab[int(nn)]
    folded = np.zeros(X // 2 + 1)
    for j in range(X // 2 + 1):
        folded[j] = sfs[j] + (sfs[X - j] if X - j != j else 0.0)
    out["n_biallelic"] = int(bi.sum())
    return out, folded


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("fasta", help="aligned FASTA (optionally .gz)")
    ap.add_argument("--min-seqs", type=int, required=True,
                    help="X: keep only sites with >= X called sequences")
    ap.add_argument("--sfs-n", type=int, default=None,
                    help="downsample SFS to this many chromosomes (default: --min-seqs)")
    ap.add_argument("--window", type=int, default=None, help="W: window size (alignment columns)")
    ap.add_argument("--overlap", type=int, default=0, help="Y: overlap between windows (< W)")
    ap.add_argument("--out", default="popstats", help="output prefix")
    a = ap.parse_args()

    X = a.sfs_n if a.sfs_n is not None else a.min_seqs
    if a.min_seqs < 2:
        sys.exit("error: --min-seqs must be >= 2")
    if X > a.min_seqs or X < 2:
        sys.exit("error: --sfs-n must satisfy 2 <= sfs-n <= --min-seqs")
    if a.window is not None and not (0 <= a.overlap < a.window):
        sys.exit("error: need 0 <= --overlap < --window")

    names, codes = read_fasta(a.fasta)
    nseq, L = codes.shape
    if nseq < a.min_seqs:
        sys.exit(f"error: only {nseq} sequences but --min-seqs {a.min_seqs}")
    cnt = site_counts(codes)
    tab = proj_table(nseq, X)

    if a.window is None:
        wins = [(0, L)]
    else:
        step = a.window - a.overlap
        wins = [(s, min(s + a.window, L)) for s in range(0, L, step)]
        # drop trailing windows fully contained in the previous one
        wins = [w for i, w in enumerate(wins) if i == 0 or w[1] > wins[i - 1][1]]

    cols = ["start", "end", "n_sites", "n_biallelic", "mean_n", "S",
            "pi_sum", "pi", "theta_w", "tajima_d"]
    with open(f"{a.out}.stats.tsv", "w") as fs, open(f"{a.out}.sfs.tsv", "w") as ff:
        fs.write("\t".join(cols) + "\n")
        ff.write("\t".join(["start", "end", "proj_invariant"] +
                           [f"k{j}" for j in range(1, X // 2 + 1)]) + "\n")
        for lo, hi in wins:
            st, folded = window_stats(cnt, lo, hi, a.min_seqs, X, tab)
            st.setdefault("n_biallelic", 0)
            st.setdefault("mean_n", float("nan"))
            row = [lo + 1, hi] + [st[c] for c in cols[2:]]
            fs.write("\t".join(f"{v:.6g}" if isinstance(v, float) else str(v)
                               for v in row) + "\n")
            ff.write("\t".join([str(lo + 1), str(hi)] +
                               [f"{v:.6g}" for v in folded]) + "\n")
    print(f"wrote {a.out}.stats.tsv and {a.out}.sfs.tsv "
          f"({len(wins)} window(s), {nseq} seqs, {L} columns)", file=sys.stderr)


if __name__ == "__main__":
    main()
