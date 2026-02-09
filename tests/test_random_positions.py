import os
import random
from bisect import bisect_right
from pathlib import Path

import pytest


def _require_integration():
    if os.getenv("RUN_INTEGRATION") != "1":
        pytest.skip("Set RUN_INTEGRATION=1 to run integration tests.")


def _read_fai_length(fai_path: Path, contig: str) -> int:
    with fai_path.open("r", encoding="utf-8") as fh:
        for line in fh:
            if not line.strip():
                continue
            parts = line.rstrip().split("\t")
            if parts[0] == contig:
                return int(parts[1])
    raise ValueError(f"Contig {contig} not found in {fai_path}")


def _read_vcf_positions(path: Path) -> set[int]:
    positions = set()
    with path.open("r", encoding="utf-8") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.rstrip().split("\t")
            if len(parts) < 2:
                continue
            try:
                positions.add(int(parts[1]))
            except ValueError:
                continue
    return positions


def _read_bed_intervals(path: Path) -> list[tuple[int, int]]:
    intervals = []
    with path.open("r", encoding="utf-8") as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip().split("\t")
            if len(parts) < 3:
                continue
            try:
                start = int(parts[1])
                end = int(parts[2])
            except ValueError:
                continue
            if end > start:
                intervals.append((start, end))
    intervals.sort()
    return intervals


def _interval_contains(intervals: list[tuple[int, int]], pos1: int) -> bool:
    # BED intervals are 0-based half-open; convert 1-based pos to 0-based.
    pos0 = pos1 - 1
    idx = bisect_right(intervals, (pos0, float("inf"))) - 1
    if idx < 0:
        return False
    start, end = intervals[idx]
    return start <= pos0 < end


def test_random_position_membership():
    _require_integration()

    contig = "1"
    base = Path("results") / "split" / f"combined.{contig}"
    inv_path = base.with_suffix(base.suffix + ".inv")
    clean_path = base.with_suffix(base.suffix + ".clean")
    filtered_bed = base.with_suffix(base.suffix + ".filtered.bed")
    missing_bed = base.with_suffix(base.suffix + ".missing.bed")
    fai_path = Path("results") / "refs" / "reference_gvcf.fa.fai"

    for path in (inv_path, clean_path, filtered_bed, missing_bed, fai_path):
        if not path.exists():
            pytest.skip(f"Missing required file: {path}")

    contig_len = _read_fai_length(fai_path, contig)
    inv_positions = _read_vcf_positions(inv_path)
    clean_positions = _read_vcf_positions(clean_path)
    filtered_intervals = _read_bed_intervals(filtered_bed)
    missing_intervals = _read_bed_intervals(missing_bed)

    rng = random.Random(1337)
    positions = [rng.randint(1, contig_len) for _ in range(10000)]
    positions.sort()

    out_path = Path("tests") / "position_membership.tsv"
    with out_path.open("w", encoding="utf-8") as fh:
        fh.write("pos\tmembership\n")
        for pos in positions:
            membership = []
            if pos in inv_positions:
                membership.append("inv")
            if pos in clean_positions:
                membership.append("clean")
            if _interval_contains(filtered_intervals, pos):
                membership.append("filtered")
            if _interval_contains(missing_intervals, pos):
                membership.append("missing")
            if not membership:
                membership.append("none")
            fh.write(f"{pos}\t{','.join(membership)}\n")
            assert len(membership) == 1, f"Position {pos} has memberships: {membership}"

    assert out_path.exists()
    lines = out_path.read_text(encoding="utf-8").strip().splitlines()
    assert len(lines) == 10001
