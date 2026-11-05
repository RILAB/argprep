import subprocess, sys
from pathlib import Path

SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "popstats.py"


def run(tmp_path, fasta, *args):
    f = tmp_path / "a.fa"
    f.write_text(fasta)
    out = tmp_path / "o"
    subprocess.run([sys.executable, str(SCRIPT), str(f), "--out", str(out), *args], check=True)
    rows = [l.split("\t") for l in (tmp_path / "o.stats.tsv").read_text().splitlines()]
    return dict(zip(rows[0], rows[1])), (tmp_path / "o.sfs.tsv").read_text().splitlines()


FA = ">a\nAAA\n>b\nAAA\n>c\nTAA\n>d\nTAN\n"


def test_missing_data_pi(tmp_path):
    # min-seqs 4: only sites 1,2 kept. site1 AATT -> 2/3, site2 0
    st, _ = run(tmp_path, FA, "--min-seqs", "4")
    assert st["n_sites"] == "2" and abs(float(st["pi_sum"]) - 2 / 3) < 1e-5
    # min-seqs 3: site 3 (n=3, all A) kept too
    st, _ = run(tmp_path, FA, "--min-seqs", "3")
    assert st["n_sites"] == "3" and abs(float(st["pi"]) - (2 / 3) / 3) < 1e-5


def test_sfs_mass_and_windows(tmp_path):
    st, sfs = run(tmp_path, FA, "--min-seqs", "3")
    vals = [float(x) for x in sfs[1].split("\t")[2:]]
    assert abs(sum(vals) - 1.0) < 1e-6  # one biallelic site
    run(tmp_path, FA, "--min-seqs", "3", "--window", "2", "--overlap", "1")
    assert len((tmp_path / "o.stats.tsv").read_text().splitlines()) == 3  # 2 windows + header (trailing contained window dropped)
