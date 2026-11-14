#!/usr/bin/env python3
"""Compare ARGprep's SNP calls with ``wgatools call`` on pairwise MAFs.

Reproduces the numbers in WGATOOLS_COMPARISON.md. For each sample MAF it runs
``wgatools call --snp`` (one single-sample VCF per sample), runs ARGprep's
``scripts/maf_to_sites.py`` once per missingness setting on all samples
together, and compares the per-sample non-reference SNP calls by
(position, REF, ALT).

Two ARGprep settings are compared, both with indel-adjacent masking off:

* ``strict``     -- ``--max-missing-fraction 0.0`` (the example/default)
* ``permissive`` -- ``--max-missing-count <n_samples>`` (no missingness filter)

Usage (bundled example; run on a compute node):

    python scripts/compare_wgatools.py \\
        --wgatools /path/to/wgatools \\
        --out-dir /tmp/argprep-wgatools

Writes ``comparison.json`` and ``comparison.md`` to ``--out-dir`` and prints
the Markdown summary. Exits non-zero if any shared call disagrees on REF/ALT.
"""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
from collections import Counter
from datetime import date
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]


def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--wgatools", default="wgatools", help="wgatools executable (default: on PATH)")
    ap.add_argument("--maf-dir", default=str(REPO_ROOT / "example_data" / "example.maf"))
    ap.add_argument("--reference-fasta", default=str(REPO_ROOT / "example_data" / "example.reference.fa"))
    ap.add_argument("--contig", default="reference", help="Reference contig name (default: reference)")
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--svlen", type=int, default=50, help="wgatools --svlen (default: 50, its own default)")
    return ap.parse_args()


def run(cmd: list[str]) -> str:
    result = subprocess.run(cmd, check=True, capture_output=True, text=True)
    return result.stdout


def tool_versions(wgatools: str) -> dict[str, str]:
    wga = run([wgatools, "--version"]).strip()
    try:
        argprep = run(["git", "-C", str(REPO_ROOT), "describe", "--tags", "--always", "--dirty"]).strip()
        commit = run(["git", "-C", str(REPO_ROOT), "rev-parse", "--short", "HEAD"]).strip()
    except (subprocess.CalledProcessError, FileNotFoundError):
        argprep, commit = "unknown", "unknown"
    return {"wgatools": wga, "argprep_describe": argprep, "argprep_commit": commit, "python": sys.version.split()[0]}


def vcf_records(path: Path):
    with open(path, encoding="utf-8") as handle:
        header: list[str] = []
        for line in handle:
            if line.startswith("##"):
                continue
            fields = line.rstrip("\n").split("\t")
            if line.startswith("#"):
                header = fields
                continue
            yield header, fields


def wgatools_calls(vcf: Path) -> tuple[dict[int, tuple[str, str]], Counter]:
    """SNPs keyed by 1-based POS -> (REF, ALT); plus counts of every record type."""
    snps: dict[int, tuple[str, str]] = {}
    kinds: Counter = Counter()
    for _, f in vcf_records(vcf):
        pos, ref, alt, info = int(f[1]), f[3].upper(), f[4].upper(), f[7]
        svtype = next((kv.split("=", 1)[1] for kv in info.split(";") if kv.startswith("SVTYPE=")), None)
        if svtype:
            kinds[svtype] += 1
        elif len(ref) == 1 and len(alt) == 1 and not alt.startswith("<"):
            kinds["SNP"] += 1
            if pos in snps:
                kinds["duplicate_snp_pos"] += 1
            snps[pos] = (ref, alt)
        else:
            kinds["other"] += 1
    return snps, kinds


def argprep_calls(vcf: Path, samples: list[str]) -> tuple[dict[str, dict[int, tuple[str, str]]], int]:
    """Per-sample non-reference SNP calls keyed by POS -> (REF, ALT allele carried)."""
    calls: dict[str, dict[int, tuple[str, str]]] = {s: {} for s in samples}
    n_records = 0
    for header, f in vcf_records(vcf):
        n_records += 1
        pos, ref = int(f[1]), f[3].upper()
        alts = [a.upper() for a in f[4].split(",")]
        for sample in samples:
            gt = f[header.index(sample)].split(":")[0]
            if gt in (".", "0"):
                continue
            calls[sample][pos] = (ref, alts[int(gt) - 1])
    return calls, n_records


def run_argprep(args, samples: list[str], out_dir: Path, name: str, extra: list[str]) -> Path:
    prefix = out_dir / f"argprep_{name}" / args.contig
    prefix.parent.mkdir(parents=True, exist_ok=True)
    cmd = [
        sys.executable, str(REPO_ROOT / "scripts" / "maf_to_sites.py"),
        "--maf-dir", args.maf_dir,
        "--reference-fasta", args.reference_fasta,
        "--contig", args.contig,
        "--out-prefix", str(prefix),
        "--samples", *samples,
        "--keep-indel-adjacent-snps",
        *extra,
    ]
    run(cmd)
    return Path(str(prefix) + ".vcf")


def compare(argprep: dict[int, tuple[str, str]], wga: dict[int, tuple[str, str]]) -> dict[str, int]:
    shared_pos = argprep.keys() & wga.keys()
    mismatched = sum(1 for p in shared_pos if argprep[p] != wga[p])
    return {
        "argprep": len(argprep),
        "wgatools": len(wga),
        "shared": len(shared_pos) - mismatched,
        "mismatched_alleles": mismatched,
        "argprep_only": len(argprep.keys() - wga.keys()),
        "wgatools_only": len(wga.keys() - argprep.keys()),
    }


def main() -> int:
    args = parse_args()
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    maf_dir = Path(args.maf_dir)
    mafs = sorted(maf_dir.glob("*.maf")) + sorted(maf_dir.glob("*.maf.gz"))
    samples = [m.name.removesuffix(".gz").removesuffix(".maf") for m in mafs]
    if not samples:
        raise SystemExit(f"No .maf/.maf.gz files in {maf_dir}")

    versions = tool_versions(args.wgatools)

    wga_dir = out_dir / "wgatools"
    wga_dir.mkdir(exist_ok=True)
    wga: dict[str, dict[int, tuple[str, str]]] = {}
    wga_kinds: Counter = Counter()
    for sample, maf in zip(samples, mafs):
        vcf = wga_dir / f"{sample}.vcf"
        run([args.wgatools, "call", "--snp", "--svlen", str(args.svlen), "-n", sample, "-r", "-o", str(vcf), str(maf)])
        wga[sample], kinds = wgatools_calls(vcf)
        wga_kinds.update(kinds)

    settings = {
        "strict": ["--max-missing-fraction", "0.0"],
        "permissive": ["--max-missing-count", str(len(samples))],
    }
    results: dict[str, dict] = {}
    for name, extra in settings.items():
        vcf = run_argprep(args, samples, out_dir, name, extra)
        calls, n_records = argprep_calls(vcf, samples)
        per_sample = {s: compare(calls[s], wga[s]) for s in samples}
        totals = Counter()
        for row in per_sample.values():
            totals.update(row)
        results[name] = {"argprep_options": extra, "variant_records": n_records,
                         "per_sample": per_sample, "total": dict(totals)}

    report = {"date": date.today().isoformat(), "versions": versions, "samples": samples,
              "wgatools_svlen": args.svlen, "wgatools_record_counts": dict(wga_kinds), "settings": results}
    (out_dir / "comparison.json").write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")

    lines = [
        f"Date: {report['date']}",
        f"wgatools: {versions['wgatools']} (call --snp --svlen {args.svlen})",
        f"ARGprep: {versions['argprep_describe']} (commit {versions['argprep_commit']}), Python {versions['python']}",
        f"wgatools record counts (all samples): {dict(wga_kinds)}",
        "",
    ]
    for name, res in results.items():
        t = res["total"]
        lines += [
            f"## {name}: {' '.join(res['argprep_options'])}, indel-adjacent masking off",
            f"ARGprep multi-sample variant records: {res['variant_records']:,}",
            "",
            "| Sample | ARGprep non-ref SNPs | wgatools SNPs | Shared | REF/ALT mismatch | ARGprep-only | wgatools-only |",
            "|---|---:|---:|---:|---:|---:|---:|",
        ]
        for s, r in res["per_sample"].items():
            lines.append(f"| {s} | {r['argprep']:,} | {r['wgatools']:,} | {r['shared']:,} | "
                         f"{r['mismatched_alleles']:,} | {r['argprep_only']:,} | {r['wgatools_only']:,} |")
        lines.append(f"| **Total** | **{t['argprep']:,}** | **{t['wgatools']:,}** | **{t['shared']:,}** | "
                     f"**{t['mismatched_alleles']:,}** | **{t['argprep_only']:,}** | **{t['wgatools_only']:,}** |")
        lines.append("")
    text = "\n".join(lines)
    (out_dir / "comparison.md").write_text(text, encoding="utf-8")
    print(text)

    mismatches = sum(res["total"]["mismatched_alleles"] for res in results.values())
    return 1 if mismatches else 0


if __name__ == "__main__":
    raise SystemExit(main())
