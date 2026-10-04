# Plan: separate indel gaps from other missing data in filters

**TODO item:** "Separate indel gaps from other missing data in filters"
**Scope:** Medium, ~1–2 hours. Mechanically straightforward; touches CLI, filter
loop, output schemas, Snakefile config, tests, docs.
**Verified against code:** 2026-07-05 (line numbers below reflect current
`scripts/maf_to_sites.py` / `scripts/summary_report.py`).

## Goal

Let users threshold indel-caused missingness (`-`, code `5`) independently of
other missing data (unaligned `0`, ambiguous/quality-masked `?` `7`). Add
`--max-indel-count` / `--max-indel-fraction` and break out indel-missing vs.
other-missing in the site summary TSV and HTML report.

## Design decisions

1. **Stacking, not replacing.** The new indel thresholds apply *in addition to*
   the existing `max_missing_*` total-missingness thresholds. Both axes must
   pass for a site to be retained. This is the more expressive default and falls
   out naturally from adding a second check.
2. **Partition of existing `MISSING_CODES = {0, 5, 7}`:**
   - `indel_missing` = code `5` (the `-` gap)
   - `other_missing` = codes `0` (unaligned) + `7` (ambiguous / quality-masked)
   - `total_missing` = all three (unchanged; existing `max_missing_*` still
     governs this)
3. **Known limitation to document, not fix here:** code `7` is overloaded —
   it covers both genuinely ambiguous `?` bases *and* v1.6 quality-mask output.
   So `other_missing` lumps unaligned + ambiguous + quality-masked. Breaking
   quality-masked out separately would need a distinct code and is out of scope.
4. **Do not physically change the alphabet.** No new codes; purely new
   accounting + a second threshold gate.

## Implementation steps

### 1. CLI args — `scripts/maf_to_sites.py:149-150`
Add next to `--max-missing-count` / `--max-missing-fraction`:
```python
ap.add_argument("--max-indel-count", type=int, default=None)
ap.add_argument("--max-indel-fraction", type=float, default=None)
```

### 2. Threshold helper — `scripts/maf_to_sites.py:416-425`
Generalize `missing_threshold` into a reusable helper (or add a parallel one).
Recommended: rename the logic to a generic
`count_threshold(sample_count, max_count, max_fraction, *, label)` that raises
with the right flag name in its message, then call it twice:
```python
allowed_missing = count_threshold(n, args.max_missing_count,
                                  args.max_missing_fraction, label="missing")
allowed_indel   = count_threshold(n, args.max_indel_count,
                                  args.max_indel_fraction, label="indel")
```
Keep the "no limit" sentinel behavior (returns `sample_count`, i.e. never
filtered) so the feature is disabled by default. Wire both at the call site
`maf_to_sites.py:584-587`.

### 3. Per-site accounting — `scripts/maf_to_sites.py:680-682`
Currently:
```python
missing_mask = (block == 0) | (block == 5) | (block == 7)
missing_counts = missing_mask.sum(axis=0).tolist()
```
Add a parallel indel count (and derive other-missing for the summary):
```python
indel_counts = (block == 5).sum(axis=0).tolist()
# other_missing = missing_counts[j] - indel_counts[j]
```

### 4. Filter gate — `scripts/maf_to_sites.py:720-727`
Insert an indel-axis check *before* the existing total-missingness check so an
indel-driven mask is attributed correctly:
```python
if indel > allowed_indel:
    mask_intervals.add(idx); masked_total += 1
    counts["masked_indel_missing"] += 1
    continue
if missing > allowed_missing:
    ...  # existing no_alignment / missingness split unchanged
```
Add a `counts["masked_indel_missing"]` counter. (Precedent: the existing
`has_unaligned` split of `masked_no_alignment` vs `masked_missingness` at
`maf_to_sites.py:704-707` shows the pattern.)

### 5. INFO field (optional but cheap) — `maf_to_sites.py:735`
Consider adding an `IS=<indel_missing>` tag alongside the existing
`MS=<missing>` so per-site indel missingness is visible in the VCF. Low effort;
include if it doesn't complicate the all-sites/variant record formatting.

### 6. Site summary TSV — `maf_to_sites.py:794-806`
Add `allowed_indel` next to `allowed_missing`, and emit the new
`masked_indel_missing` count (plus, if useful, aggregate indel-missing vs.
other-missing totals) in the metrics block.

### 7. HTML report — `scripts/summary_report.py:753`
Add a metric row mirroring the existing `("allowed_missing", "Allowed missing")`
tuple, e.g. `("allowed_indel", "Allowed indel-missing")`, and surface the
`masked_indel_missing` mask category wherever masked categories are tabulated
(`masked_counts` handling around `summary_report.py:486-535`).

### 8. Snakefile + config plumbing
- Add `max_indel_count` / `max_indel_fraction` config keys (default null =
  disabled) to `config.yaml` / example config.
- Thread them into the `maf_to_sites` rule's shell command (parallel to how
  `max_missing_*` are passed).

### 9. Tests
- Fixture MAF where a site has enough `-` gaps to exceed `max_indel_count` but
  would pass `max_missing_*` → expect it masked with `masked_indel_missing`.
- Fixture where non-indel missing (unaligned/`?`) exceeds `max_missing_*` but
  indels are under `max_indel_count` → existing missingness masking still fires.
- Default (both unset) → no behavior change vs. current output (regression).
- `count_threshold` unit test for the fraction/count/none cases + validation
  error messages.

### 10. Docs
- README: document the two new options and the stacking semantics.
- `changelog.md`: new entry (next version bump).

## Open question to confirm with user before coding
- Include the optional `IS=` INFO field (step 5), or keep the VCF INFO
  untouched and expose indel-missing only via the summary outputs?
