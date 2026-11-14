## Repo Rules

- Do not ask whether to create a git commit or push after making changes.
- When the user says "commit and push", immediately run the necessary git commands (including `git add`, `git status`, `git commit`, `git push`) without any further confirmation or questions.
- Never ask whether to run `git commit` or `git push` after the user has already said "commit and push".
- For Snakemake runs, use a longer command timeout to avoid premature termination.
- Never commit or push SLURM output files (e.g., `slurm-*.out`) or Snakemake log files (e.g., `*.snakemake.log`, `.snakemake/log/*`).
- `direct_maf` is the default branch on GitHub (not `main`; a stale `main` also exists). Push, open PRs and tag releases against `direct_maf`; do not suggest merging into `main`.
