#!/usr/bin/env python3
"""
Cross-trait genetic correlation with LDSC
- Project: Genetic determinants of Lewy body pathology
- Version: Python/3.11 (standard library only); LDSC v1.0.1
- Last Updated: 16-SEP-2026

Overview
- Runs `ldsc.py --rg` for one phenotype (LBD) against several traits in a single call,
  so all comparisons share the same LD reference and merged SNP set, then parses the
  log into a CSV with QC checks (usable SNPs per trait, failed traits, gcov_int, |rg| > 1).
- In the CSV, h2_obs / h2_int are each p2 trait's; phenotype 1's h2 is printed separately.

Usage
    python 05_LDSC.py --ldsc ldsc/ldsc.py --ld-ref eur_w_ld_chr/ \\
        --p1 LBD.sumstats.gz --p2 PD.sumstats.gz AD_Kunkle.sumstats.gz AD_Bellenguez.sumstats.gz \\
        --out results/LBD_rg
    python 05_LDSC.py --parse results/LBD_rg.log
  Inputs are munged .sumstats.gz files (LDSC munge_sumstats.py) and per-chromosome
  LD scores (e.g. eur_w_ld_chr/). See --help for all options.

"""
import argparse
import csv
import glob
import math
import os
import re
import subprocess
import sys

COLUMNS = ["Comparison", "rg", "SE", "Z", "P", "h2_obs", "h2_obs_se",
           "h2_int", "h2_int_se", "gcov_int", "gcov_int_se"]
LDSC_BANNER = "* LD Score Regression (LDSC)"


# ── running LDSC ─────────────────────────────────────────────────────────────
def ld_prefix(path):
    """LDSC appends the chromosome number to --ref-ld-chr/--w-ld-chr, so a
    directory must end with a separator. Returns the absolute prefix."""
    path = os.path.abspath(path) + (os.sep if os.path.isdir(path) else "")
    if not glob.glob(path + "1.l2.ldscore*"):
        sys.exit(f"ERROR: no LD scores found at {path}1.l2.ldscore*")
    return path


def run_ldsc(args):
    for f in [args.ldsc, args.p1, *args.p2]:
        if not os.path.isfile(f):
            sys.exit(f"ERROR: file not found: {f}")
        if "," in os.path.abspath(f):
            sys.exit(f"ERROR: LDSC --rg splits on commas, so paths cannot contain one: {f}")
    ref = ld_prefix(args.ld_ref)
    w = ld_prefix(args.w_ld) if args.w_ld else ref
    sumstats = [os.path.abspath(f) for f in [args.p1, *args.p2]]
    out = os.path.abspath(args.out)
    os.makedirs(os.path.dirname(out), exist_ok=True)
    if os.path.exists(out + ".log"):
        print(f"NOTE: {out}.log exists; the parser reads only the last run in it.",
              file=sys.stderr)

    cmd = [args.python, os.path.abspath(args.ldsc),
           "--rg", ",".join(sumstats),
           "--ref-ld-chr", ref, "--w-ld-chr", w,
           "--out", out]
    print("+ " + " ".join(cmd), file=sys.stderr, flush=True)
    try:
        r = subprocess.run(cmd)
    except OSError as e:
        sys.exit(f"ERROR: could not run {args.python}: {e}")
    if r.returncode:
        sys.exit(f"ERROR: LDSC exited with status {r.returncode}; see {out}.log")
    return out + ".log"


# ── parsing the log ──────────────────────────────────────────────────────────
def last_run(text):
    """Text of the last LDSC run in a log, and the number of runs it contains."""
    n = text.count(LDSC_BANNER)
    return (text[text.rfind(LDSC_BANNER):] if n else text), n


def to_number(s):
    """float for numeric fields; None for LDSC's NA and for nan/inf."""
    try:
        v = float(s)
    except (TypeError, ValueError):
        return None
    return v if math.isfinite(v) else None


READ_RE = re.compile(r"Reading summary statistics from (.+?) \.\.\.\s*$", re.M)


def parse_traits(run):
    """Phenotype-1 path, and one dict per p2 trait in run order:
    {path, n_valid}. Taken from the per-trait log blocks rather than the summary
    table, so paths containing spaces are kept intact."""
    head, *blocks = re.split(r"Computing rg for phenotype \d+/\d+", run)
    p1 = READ_RE.search(head)
    traits = []
    for block in blocks:
        f = READ_RE.search(block)
        v = re.search(r"(\d+) SNPs with valid alleles", block)
        traits.append({"path": f.group(1) if f else None,
                       "n_valid": int(v.group(1)) if v else None})
    return (p1.group(1) if p1 else None), traits


def parse_summary(run, n_traits):
    """Header and numeric fields of the 'Summary of Genetic Correlation Results'
    table. The p1/p2 path columns are skipped: values are read from the right,
    so spaces in paths cannot shift the columns."""
    lines = run.splitlines()
    try:
        start = next(i for i, ln in enumerate(lines)
                     if ln.strip().startswith("Summary of Genetic Correlation Results"))
        header = lines[start + 1].split()
    except (StopIteration, IndexError):
        sys.exit("ERROR: no 'Summary of Genetic Correlation Results' table. "
                 "Did LDSC finish? Check the end of the log.")
    if header[:2] != ["p1", "p2"]:
        sys.exit(f"ERROR: unexpected summary header: {' '.join(header)}")
    value_cols = header[2:]
    rows = []
    for ln in lines[start + 2:]:
        if not ln.strip():
            break
        rows.append(dict(zip(value_cols, ln.split()[-len(value_cols):])))
    if n_traits and len(rows) != n_traits:
        sys.exit(f"ERROR: summary table has {len(rows)} rows but the log shows "
                 f"{n_traits} traits")
    return value_cols, rows


def parse_p1_h2(run):
    """Phenotype 1's h2 (and intercept, if estimated) from its own log block, or None."""
    m = re.search(r"Heritability of phenotype 1\n(.*?)(?=\nHeritability of phenotype|"
                  r"\nGenetic Covariance|\Z)", run, flags=re.S)
    if not m:
        return None
    block = m.group(1)
    h2 = re.search(r"Total (Observed|Liability) scale h2:\s*(\S+)\s*\((\S+)\)", block)
    icpt = re.search(r"Intercept:\s*(\S+)\s*\((\S+)\)", block)
    if not h2:
        return None
    return {"scale": h2.group(1).lower(), "h2": to_number(h2.group(2)),
            "h2_se": to_number(h2.group(3)),
            "int": to_number(icpt.group(1)) if icpt else None,
            "int_se": to_number(icpt.group(2)) if icpt else None}


def default_label(path):
    return re.sub(r"\.sumstats(\.gz)?$", "", os.path.basename(path))


def fmt(v, spec=".4f"):
    return "NA" if v is None else format(v, spec)


def parse_log(logfile, csv_out, labels=None, min_snp_frac=0.9):
    with open(logfile) as fh:
        run, n_runs = last_run(fh.read())
    if n_runs > 1:
        print(f"NOTE: {logfile} holds {n_runs} LDSC runs; parsing the last one.",
              file=sys.stderr)

    p1_path, traits = parse_traits(run)
    value_cols, rows = parse_summary(run, len(traits))
    if labels and len(labels) != len(rows):
        sys.exit(f"ERROR: {len(labels)} labels for {len(rows)} traits in the log")

    # LDSC reports liability-scale h2 instead when run with --samp-prev/--pop-prev
    h2_col = "h2_liab" if "h2_liab" in value_cols else "h2_obs"
    columns = [c.replace("h2_obs", h2_col) for c in COLUMNS]
    keys = ["rg", "se", "z", "p", h2_col, h2_col + "_se",
            "h2_int", "h2_int_se", "gcov_int", "gcov_int_se"]
    missing = [k for k in keys if k not in value_cols]
    if missing:
        sys.exit(f"ERROR: summary table lacks columns {missing}: {' '.join(value_cols)}")

    table = []
    for i, r in enumerate(rows):
        path = traits[i]["path"] if traits else None
        out = {"Comparison": labels[i] if labels else default_label(path or f"trait{i + 2}"),
               "n_valid": traits[i]["n_valid"] if traits else None}
        for col, key in zip(columns[1:], keys):
            out[col] = to_number(r[key])
        table.append(out)

    with open(csv_out, "w", newline="") as fh:
        w = csv.writer(fh, lineterminator="\n")
        w.writerow(columns)
        for r in table:
            w.writerow(["" if r[c] is None else r[c] for c in columns])

    # ── report ───────────────────────────────────────────────────────────────
    h2 = parse_p1_h2(run)
    if h2:
        z = h2["h2"] / h2["h2_se"] if h2["h2"] is not None and h2["h2_se"] else None
        icpt = (f"intercept = {fmt(h2['int'])} (SE {fmt(h2['int_se'])})"
                if h2["int"] is not None else "intercept constrained or not reported")
        print(f"\n[phenotype 1] {p1_path or ''}\n"
              f"  {h2['scale']}-scale h2 = {fmt(h2['h2'])} (SE {fmt(h2['h2_se'])}), "
              f"z = {fmt(z, '.2f')}; {icpt}\n"
              f"  {h2_col} / h2_int in the table are each p2 trait's, not phenotype 1's.",
              file=sys.stderr)

    def cell(v):
        return v if isinstance(v, str) else "NA" if v is None else format(v, ".4g")

    print()
    print("  ".join(f"{c:>13s}" for c in columns))
    for r in table:
        print("  ".join(f"{cell(r[c]):>13s}" for c in columns))

    m = re.search(r"After merging with regression SNP LD, (\d+) SNPs remain", run)
    n_p1 = int(m.group(1)) if m else None
    print("\nQC" + (f" (phenotype-1 SNPs after merging with LD scores: {n_p1:,})"
                    if n_p1 else ""), file=sys.stderr)
    for r in table:
        flags = []
        n = r["n_valid"]
        frac = n / n_p1 if n and n_p1 else None
        if r["rg"] is None:
            flags.append("LDSC FAILED for this trait; see the traceback in the log")
        if frac is not None and frac < min_snp_frac:
            flags.append(f"only {frac:.0%} of SNPs usable; check rsID/build matching at munging")
        if r["gcov_int"] is not None and r["gcov_int_se"] and \
                abs(r["gcov_int"] / r["gcov_int_se"]) > 2:
            flags.append("gcov_int away from 0: possible sample overlap")
        if r["rg"] is not None and abs(r["rg"]) > 1:
            flags.append("|rg| > 1: interpret with its SE")
        n_txt = f"{n:,} SNPs" if n else "SNPs unknown"
        frac_txt = f" ({frac:.1%})" if frac is not None else ""
        print(f"  {r['Comparison']:<20s} {n_txt}{frac_txt}"
              + ("".join(f"\n      ! {f}" for f in flags)), file=sys.stderr)

    print(f"\nWrote {csv_out}", file=sys.stderr)
    return sum(r["rg"] is None for r in table)


# ── CLI ──────────────────────────────────────────────────────────────────────
def main():
    ap = argparse.ArgumentParser(
        description="Run LDSC --rg for one phenotype against several traits and parse "
                    "the log into a CSV. See the module docstring for details.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    run = ap.add_argument_group("run LDSC")
    run.add_argument("--ldsc", help="path to ldsc.py")
    run.add_argument("--ld-ref", help="LD score prefix or directory, e.g. eur_w_ld_chr/")
    run.add_argument("--w-ld", help="regression-weight LD score prefix (default: --ld-ref)")
    run.add_argument("--p1", help="phenotype 1 .sumstats.gz")
    run.add_argument("--p2", nargs="+", help="traits to correlate with phenotype 1")
    run.add_argument("--out", help="output prefix for <out>.log and <out>.csv")
    run.add_argument("--python", default=sys.executable,
                     help="interpreter used to run ldsc.py")

    ap.add_argument("--parse", metavar="LOG", help="parse an existing LDSC --rg log only")
    ap.add_argument("--csv", help="CSV path (default: <out>.csv, or the log path with .csv)")
    ap.add_argument("--labels", nargs="+",
                    help="Comparison labels, one per p2 trait in order "
                         "(default: sumstats file names)")
    ap.add_argument("--min-snp-frac", type=float, default=0.9,
                    help="warn when a trait keeps fewer than this fraction of phenotype-1 SNPs")
    args = ap.parse_args()

    if args.parse:
        extra = [o for o in ("ldsc", "ld_ref", "w_ld", "p1", "p2", "out") if getattr(args, o)]
        if extra:
            ap.error("--parse only reads a log; drop " +
                     ", ".join("--" + o.replace("_", "-") for o in extra))
        logfile = args.parse
        csv_out = args.csv or re.sub(r"\.log$", "", logfile) + ".csv"
    else:
        missing = [o for o in ("ldsc", "ld_ref", "p1", "p2", "out") if not getattr(args, o)]
        if missing:
            ap.error("to run LDSC, provide " + ", ".join("--" + m.replace("_", "-")
                                                        for m in missing)
                     + " (or use --parse LOG)")
        if args.labels and len(args.labels) != len(args.p2):
            ap.error(f"--labels has {len(args.labels)} names for {len(args.p2)} --p2 traits")
        logfile = run_ldsc(args)
        csv_out = args.csv or args.out + ".csv"

    n_failed = parse_log(logfile, csv_out, args.labels, args.min_snp_frac)
    if n_failed:
        sys.exit(f"ERROR: LDSC failed for {n_failed} trait(s); their CSV rows are empty")


if __name__ == "__main__":
    main()
