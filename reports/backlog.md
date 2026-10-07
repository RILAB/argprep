# ARGprep backlog

Open work is tracked in [GitHub issues](https://github.com/RILAB/argprep/issues), not here. On 2026-10-06 this file's open items moved to issues #34–#44. Decisions not to act ("deferred by decision") became closed issues labelled `wontfix` (#45–#47), so a search turns them up before anyone re-reports them.

The dated reports in [code_review/](code_review/) are historical records of what each review found. Don't edit them; file issues for any follow-ups.

## Closed since the 2026-07 review

Recorded so these are not reopened from the old lists:

- Zero-length contig handling: v1.7 added an empty-output path; v1.9 replaced it with an explicit rejection.
- Leading-insertion indel-adjacent flagging: fixed in v1.8.
- `intervals_from_positions`, `summarize_site_and_mask_coverage`: deleted.
- "Unaligned" genome-overview segment always 0: the segment is gone.
- gzip per-contig MAF chunks (half of 2026-07 efficiency #6): shipped in v1.8.
- Test gaps for `parse_maf_path_map` error branches, split `.maf.gz` input, and `read_sample_missing_bp` longest-suffix matching: all now covered.
- Config propagation: every config key that changes which sites are kept has an end-to-end workflow test of both states (closed in v1.9).
- Hardcoded personal paths in the ad-hoc analysis scripts: generalized to argparse; `chr1_variable_plot.py` deleted in v1.9.
- "Real alignment FASTA from the MAFs": shipped as `scripts/maf_to_fasta.py` (v1.10).
