"""Shared MAF-projection helpers used by maf_to_sites.py and maf_to_fasta.py.

Nucleotide codes, quality masks, sample/MAF path resolution, reference FASTA
access and reference/sample alignment-row selection. Public API: import from
here rather than from either CLI module.
"""

from __future__ import annotations

import bisect
import sys
from pathlib import Path

try:
    from scripts.common import MafRecord, merge_intervals, normalize_contig, open_text
except ModuleNotFoundError:
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
    from scripts.common import MafRecord, merge_intervals, normalize_contig, open_text


NUC_TO_CODE = {
    "A": 1,
    "C": 2,
    "G": 3,
    "T": 4,
    "-": 5,
    "?": 7,
}
CODE_TO_BASE = {
    1: "A",
    2: "C",
    3: "G",
    4: "T",
}
VALID_BASES = {"A", "C", "G", "T"}
MISSING_CODES = {0, 5, 7}


class QualityMask:
    """Low-quality intervals per source sequence, in that sequence's own
    forward-strand coordinates.

    A per-sample BED file lists ``chrom start end score`` rows (score is a
    0-1 assembly-quality value). Any base whose score is *below* the configured
    threshold is considered low quality; aligned bases at those positions are
    treated as missing during parsing. Intervals with a passing score do not
    need to be listed.
    """

    def __init__(self, low_quality: dict[str, list[tuple[int, int]]]):
        # Merged, sorted low-quality intervals per source sequence, plus a
        # parallel list of interval starts for bisecting.
        self._intervals = low_quality
        self._starts = {
            src: [start for start, _ in intervals]
            for src, intervals in low_quality.items()
        }

    def is_low(self, src: str, pos: int) -> bool:
        intervals = self._intervals.get(src)
        if not intervals:
            return False
        starts = self._starts[src]
        i = bisect.bisect_right(starts, pos) - 1
        if i < 0:
            return False
        start, end = intervals[i]
        return start <= pos < end


def load_quality_mask(bed_path: Path, quality_min: float) -> QualityMask:
    """Build a :class:`QualityMask` from a per-sample quality BED file.

    Rows are ``chrom start end score`` (whitespace-separated); ``track``/
    ``browser``/comment lines are ignored. Malformed data rows raise an error.
    Only intervals whose score is below ``quality_min`` are retained.
    """
    raw: dict[str, list[tuple[int, int]]] = {}
    with open_text(bed_path, "rt", errors="ignore") as handle:
        for line_number, line in enumerate(handle, start=1):
            stripped = line.strip()
            if not stripped or stripped.startswith(("#", "track", "browser")):
                continue
            parts = stripped.split()
            if len(parts) < 4:
                raise ValueError(
                    f"Malformed quality BED row in {bed_path} at line {line_number}: "
                    "expected at least 4 fields"
                )
            try:
                start = int(parts[1])
                end = int(parts[2])
                score = float(parts[3])
            except ValueError as exc:
                raise ValueError(
                    f"Malformed quality BED row in {bed_path} at line {line_number}: "
                    "start and end must be integers and score must be numeric"
                ) from exc
            if score < quality_min:
                raw.setdefault(parts[0], []).append((start, end))
    merged = {src: merge_intervals(intervals) for src, intervals in raw.items()}
    return QualityMask(merged)


def quality_bed_for_sample(quality_bed_dir: Path, sample: str) -> Path | None:
    plain = quality_bed_dir / f"{sample}.bed"
    gz = quality_bed_dir / f"{sample}.bed.gz"
    if plain.exists():
        return plain
    if gz.exists():
        return gz
    return None


def discover_samples(maf_dir: Path) -> list[str]:
    samples = set()
    for path in maf_dir.iterdir():
        if path.name.endswith(".maf.gz"):
            samples.add(path.name[: -len(".maf.gz")])
        elif path.name.endswith(".maf"):
            samples.add(path.name[: -len(".maf")])
    return sorted(samples)


def maf_path_for_sample(maf_dir: Path, sample: str) -> Path:
    plain = maf_dir / f"{sample}.maf"
    gz = maf_dir / f"{sample}.maf.gz"
    if plain.exists():
        return plain
    if gz.exists():
        return gz
    raise FileNotFoundError(f"Missing MAF for sample '{sample}' under {maf_dir}")


def parse_maf_path_map(values: list[str] | None) -> dict[str, Path]:
    paths: dict[str, Path] = {}
    for value in values or []:
        if "=" not in value:
            raise ValueError(f"Invalid --maf-paths entry {value!r}; expected SAMPLE=PATH")
        sample, path = value.split("=", 1)
        if not sample:
            raise ValueError(f"Invalid --maf-paths entry {value!r}; sample is empty")
        paths[sample] = Path(path)
    return paths


def maf_path_for_sample_with_map(
    maf_dir: Path,
    sample: str,
    maf_paths: dict[str, Path],
) -> Path:
    if sample in maf_paths:
        return maf_paths[sample]
    return maf_path_for_sample(maf_dir, sample)


def _read_fai_entry(fai_path: Path, contig: str) -> tuple[int, int, int, int] | None:
    """Return (length, offset, linebases, linewidth) for a contig from a .fai file."""
    with fai_path.open("r", encoding="utf-8") as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 5 and parts[0] == contig:
                return int(parts[1]), int(parts[2]), int(parts[3]), int(parts[4])
    return None


def read_contig_sequence(reference_fasta: Path, contig: str) -> str:
    fai_path = Path(str(reference_fasta) + ".fai")
    if fai_path.exists() and not str(reference_fasta).endswith(".gz"):
        entry = _read_fai_entry(fai_path, contig)
        if entry is not None:
            length, offset, linebases, linewidth = entry
            if length == 0:
                return ""
            full_lines, remainder = divmod(length, linebases)
            total_bytes = full_lines * linewidth + remainder
            with open(reference_fasta, "rb") as fh:
                fh.seek(offset)
                raw = fh.read(total_bytes)
            seq = raw.replace(b"\r", b"").replace(b"\n", b"")
            if len(seq) != length:
                raise ValueError(
                    f"FAI seek returned {len(seq)} bases for '{contig}', expected {length}"
                )
            return seq.decode("ascii").upper()
    # Fallback: linear scan (handles .gz or missing .fai)
    seq_parts: list[str] = []
    current: str | None = None
    with open_text(reference_fasta, "rt", errors="ignore") as handle:
        for raw in handle:
            if raw.startswith(">"):
                name = raw[1:].strip().split()[0]
                current = name
                continue
            if current == contig:
                seq_parts.append(raw.strip())
    if not seq_parts:
        raise ValueError(f"Contig '{contig}' not found in {reference_fasta}")
    return "".join(seq_parts).upper()


def read_contig_length(reference_fasta: Path, contig: str) -> int:
    """Length of a reference contig, from the .fai index when one exists."""
    fai_path = Path(str(reference_fasta) + ".fai")
    if fai_path.exists() and not str(reference_fasta).endswith(".gz"):
        entry = _read_fai_entry(fai_path, contig)
        if entry is not None:
            return entry[0]
    return len(read_contig_sequence(reference_fasta, contig))


def read_contig_region(reference_fasta: Path, contig: str, start: int, end: int) -> str:
    """Reference bases over the 0-based half-open interval ``[start, end)``.

    With a .fai index this seeks straight to the region, so pulling a kilobase
    window does not cost a whole-contig read; without one it falls back to
    :func:`read_contig_sequence` and slices.
    """
    if start < 0 or end < start:
        raise ValueError(f"Invalid region [{start}, {end}) for '{contig}'")
    if start == end:
        return ""
    fai_path = Path(str(reference_fasta) + ".fai")
    if fai_path.exists() and not str(reference_fasta).endswith(".gz"):
        entry = _read_fai_entry(fai_path, contig)
        if entry is not None:
            length, offset, linebases, linewidth = entry
            if end > length:
                raise ValueError(
                    f"Region [{start}, {end}) exceeds length {length} of '{contig}'"
                )
            first = offset + (start // linebases) * linewidth + (start % linebases)
            last = offset + ((end - 1) // linebases) * linewidth + ((end - 1) % linebases)
            with open(reference_fasta, "rb") as fh:
                fh.seek(first)
                raw = fh.read(last - first + 1)
            seq = raw.replace(b"\r", b"").replace(b"\n", b"")
            if len(seq) != end - start:
                raise ValueError(
                    f"FAI seek returned {len(seq)} bases for '{contig}':[{start}, {end}), "
                    f"expected {end - start}"
                )
            return seq.decode("ascii").upper()
    seq = read_contig_sequence(reference_fasta, contig)
    if end > len(seq):
        raise ValueError(
            f"Region [{start}, {end}) exceeds length {len(seq)} of '{contig}'"
        )
    return seq[start:end]


def choose_sample_record(block: list[MafRecord], contig: str) -> tuple[MafRecord, MafRecord] | None:
    if not block:
        return None
    ref = block[0]
    if normalize_contig(ref.src) != normalize_contig(contig):
        return None
    if ref.strand != "+":
        # Callers walk reference columns from ref.start as forward-strand
        # positions; a minus-strand reference row would be placed silently at
        # the wrong coordinates. AnchorWave never writes one, so reject it.
        raise ValueError(
            f"MAF reference row for '{ref.src}' at start {ref.start} is on strand "
            f"'{ref.strand}'; only '+' reference rows are supported"
        )
    if len(block) >= 2:
        # Some pairwise MAFs use the same contig name for both the reference
        # and query rows, so the second alignment row is still the sample.
        return ref, block[1]
    return None


def assign_code(calls: bytearray, idx: int, code: int) -> None:
    existing = calls[idx]
    if existing == 0:
        calls[idx] = code
        return
    if existing == code:
        return
    calls[idx] = NUC_TO_CODE["?"]
