# ARGprep implementation notes

Detail not needed for typical use, but useful when auditing edge cases or reasoning about output. See [README.md](README.md) for the user-facing pipeline documentation.

## Indel masking semantics

Deletion characters (`-`) in any sample contribute to that sample's missing-sample count at the affected reference position, so indel-overlapped sites are dropped whenever they cross the missingness threshold. To reproduce the old `mask_indels: true` behavior of dropping every site overlapped by *any* deletion, set `max_missing_count: 0`.

`mask_indel_adjacent_snps` (default `false`) masks SNPs whose reference position immediately flanks a deletion (in any sample) or an insertion (a non-`-` sample base aligned to a reference gap) when enabled. Invariant adjacent positions are not masked — only sites that would otherwise be called as variants. The adjacency flags OR across samples, so a single sample's indel is enough to mask a flanking SNP.

## Multi-block alignment conflicts

A single sample's MAF may contain multiple alignment blocks that cover the same reference position (for example overlapping or supplementary alignments). When that happens:

- Agreement on a non-missing base — the base is kept.
- Disagreement on non-missing bases — the position is set to `?` and counts as missing.
- A gap call from one block does not overwrite a base call from another.

Implemented in `_assign_code` in [scripts/maf_to_sites.py](scripts/maf_to_sites.py).

## Site classification in `site_summary.tsv`

The masked-site categories (`masked_ref_non_acgt`, `masked_indel_adjacent`, `masked_multiallelic`, `masked_no_alignment`, `masked_missingness`) are mutually exclusive: each masked position contributes to exactly one. When a site fails multiple filters, the category recorded reflects the cause checked first by the pipeline; see [scripts/maf_to_sites.py](scripts/maf_to_sites.py) for the order.

`masked_no_alignment` specifically marks sites where at least one sample had no alignment block covering that position. `masked_missingness` is used when the missingness threshold was exceeded but every sample had an alignment block (the missing calls came from gaps, `N`s, or ambiguity codes rather than absent alignments).

## Genotype encoding and ploidy

ARGprep is built for haploid samples. Each `GT` field in both VCFs is a single integer (`0` for the reference allele, `1`/`2`/... for ALTs in `ALT` order, `.` for missing). With `add_ref: true`, the synthetic `REF` sample is emitted as `0` at every retained site. Because each MAF is one haplotype, a diploid individual is represented by two samples (one per phased haplotype assembly, aligned separately); ARGprep does not phase or split diploid assemblies itself.
