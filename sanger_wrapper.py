#!/usr/bin/env python3
"""
Sanger wrapper: pairing of forward/reverse files, batch processing,
summary files (FASTA, CSV, JSON, HTML report) and ZIP packaging.

Recognised naming (case-insensitive), e.g.
  sample_forward.ab1 / sample_reverse.ab1      sample_F.ab1 / sample_R.ab1
  sample_fwd.ab1     / sample_rev.ab1          sample_27F.ab1 / sample_1492R.ab1
  sample_M13F(-21).ab1 / sample_M13R.ab1       sample_forward_read1.ab1

CLI:
  python sanger_wrapper.py -i all_ab1_files/ -o results [--reference ref.fasta] [--dry-run]
"""

import argparse
import csv
import glob
import html
import json
import os
import re
import shutil
import sys

import sanger_workflow as wf

_PAT = re.compile(
    r"^(?P<s>.+?)[_\-. ](?P<d>forward|fwd|fw|reverse|rev|rv|[A-Za-z]?\d*[FR])(?:[^A-Za-z].*)?\.ab1$",
    re.IGNORECASE,
)


def _direction(token: str) -> str:
    t = token.lower()
    if t in ("forward", "fwd", "fw"):
        return "F"
    if t in ("reverse", "rev", "rv"):
        return "R"
    return "F" if t[-1] == "f" else "R"


def clean_name(name: str) -> str:
    return re.sub(r"[^\w.\-]+", "_", name).strip("_") or "sample"


def sample_name_from_file(filename: str) -> str:
    """'plasmid3_F.ab1' -> 'plasmid3' (falls back to the file stem)."""
    m = _PAT.match(os.path.basename(filename))
    return clean_name(m.group("s") if m else os.path.splitext(os.path.basename(filename))[0])


def find_ab1(paths):
    out = []
    for p in paths:
        if os.path.isdir(p):
            out += sorted(glob.glob(os.path.join(p, "*.ab1")) + glob.glob(os.path.join(p, "*.AB1")))
        elif os.path.isfile(p) and p.lower().endswith(".ab1"):
            out.append(p)
    return [os.path.abspath(p) for p in out]


def pair_files(paths):
    """
    Group AB1 files into samples.
    Returns (jobs, unmatched) where jobs = [(sample, forward|None, reverse|None)].
    A sample with only one direction is still returned (single-read mode).
    """
    fwd, rev, unmatched = {}, {}, []
    for p in paths:
        m = _PAT.match(os.path.basename(p))
        if not m:
            unmatched.append(p)
            continue
        sample = clean_name(m.group("s"))
        table = fwd if _direction(m.group("d")) == "F" else rev
        if sample in table:
            unmatched.append(p)           # duplicate name
        else:
            table[sample] = p
    jobs = [(s, fwd.get(s), rev.get(s)) for s in sorted(set(fwd) | set(rev))]
    return jobs, unmatched


def run_batch(jobs, output_dir, trim_quality=0.05, reference=None, progress=None, min_base_q=10):
    """
    Run every job. Returns (stats_list, errors) with errors = [(sample, message)].
    One failing sample never stops the others. Writes summary files.
    """
    os.makedirs(output_dir, exist_ok=True)
    stats_list, errors = [], []
    for i, (sample, fwd, rev) in enumerate(jobs, 1):
        if progress:
            progress(i - 1, len(jobs), sample)
        try:
            stats_list.append(wf.run_sample(sample, fwd, rev, os.path.join(output_dir, sample),
                                            trim_quality, reference, min_base_q))
        except Exception as e:  # noqa: BLE001 - report and continue
            errors.append((sample, str(e)))
    if progress:
        progress(len(jobs), len(jobs), "")
    if stats_list:
        write_summary(output_dir, stats_list)
    return stats_list, errors


# ---------------------------------------------------------------------------
# Summary files
# ---------------------------------------------------------------------------
_CSV_COLS = ["sample_name", "mode", "consensus_length", "gc_content", "ambiguity_pct", "n_pct",
             "q20_pct", "mean_consensus_quality", "overlap_length", "overlap_identity", "mismatches",
             "forward_avg_quality", "reverse_avg_quality", "quality_rating"]


def n_variants(s):
    return len(s["reference"]["variants"]) if s.get("reference") else None


def write_summary(output_dir, stats_list):
    with open(os.path.join(output_dir, "all_consensus.fasta"), "w") as f:
        for s in stats_list:
            f.write(f">{s['sample_name']}\n{s['consensus_seq']}\n")
    with open(os.path.join(output_dir, "summary.csv"), "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(_CSV_COLS + ["variants_vs_reference"])
        for s in stats_list:
            w.writerow([s[c] for c in _CSV_COLS] + [n_variants(s) if n_variants(s) is not None else ""])
    slim = [{k: v for k, v in s.items() if k not in ("consensus_quals", "consensus_src", "consensus_fmap", "consensus_rmap")} for s in stats_list]
    with open(os.path.join(output_dir, "all_sample_stats.json"), "w") as f:
        json.dump(slim, f, indent=1)
    write_html_report(output_dir, stats_list)


def write_html_report(output_dir, stats_list):
    e = html.escape
    n = len(stats_list)
    avg = lambda k: sum(s[k] for s in stats_list) / n
    has_ref = any(s.get("reference") for s in stats_list)
    rows = ""
    for s in stats_list:
        nv = n_variants(s)
        rows += (
            f"<tr><td><b>{e(s['sample_name'])}</b></td><td>{e(s['mode'])}</td>"
            f"<td>{s['consensus_length']}</td><td>{s['mean_consensus_quality']}</td><td>{s['q20_pct']}</td>"
            f"<td>{s['ambiguity_pct']}</td><td>{s['mismatches']}</td>"
            + (f"<td>{'' if nv is None else nv}</td>" if has_ref else "")
            + f"<td><span class='r {e(s['quality_rating'])}'>{e(s['quality_rating'])}</span></td></tr>"
        )
    details = ""
    for s in stats_list:
        seq = s["consensus_seq"]
        wrapped = "\n".join(seq[i:i + 80] for i in range(0, len(seq), 80))
        warn = "".join(f"<p class='w'>⚠ {e(w)}</p>" for w in s["warnings"])
        vr = ""
        if s.get("reference"):
            r = s["reference"]
            vr = (f"<p>vs <b>{e(r['reference_name'])}</b>: identity {r['identity']}%, "
                  f"{len(r['variants'])} variant(s)</p><ul>"
                  + "".join(f"<li>{e(v['type'])} {e(v['change'])} (Q{v['quality']})</li>" for v in r["variants"])
                  + "</ul>")
        details += (f"<div class='d'><h3>{e(s['sample_name'])} <span class='r {e(s['quality_rating'])}'>"
                    f"{e(s['quality_rating'])}</span></h3>{warn}{vr}<pre>&gt;{e(s['sample_name'])}\n{e(wrapped)}</pre></div>")
    page = f"""<!DOCTYPE html><html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1"><title>Sanger Quality Report</title>
<style>body{{font-family:system-ui,sans-serif;background:#f5f7fa;color:#222;margin:0;padding:24px}}
.c{{max-width:1100px;margin:auto}}h1{{margin:0 0 4px}}.s{{color:#666;margin-bottom:20px}}
.cards{{display:grid;grid-template-columns:repeat(auto-fit,minmax(160px,1fr));gap:12px;margin-bottom:24px}}
.card{{background:#fff;border-radius:10px;padding:16px;text-align:center;box-shadow:0 1px 6px #0001}}
.card b{{display:block;font-size:26px}}.card span{{color:#888;font-size:12px}}
table{{width:100%;border-collapse:collapse;background:#fff;border-radius:10px;overflow:hidden;box-shadow:0 1px 6px #0001}}
th{{background:#1a1a2e;color:#fff;text-align:left;padding:10px;font-size:12px}}td{{padding:8px 10px;border-bottom:1px solid #eee;font-size:14px}}
.r{{padding:2px 9px;border-radius:10px;font-size:12px;font-weight:600}}.Excellent{{background:#d4edda;color:#155724}}
.Good{{background:#d1ecf1;color:#0c5460}}.Fair{{background:#fff3cd;color:#856404}}.Poor{{background:#f8d7da;color:#721c24}}
.d{{background:#fff;border-radius:10px;padding:16px;margin-top:16px;box-shadow:0 1px 6px #0001}}
pre{{background:#f0f2f5;padding:10px;border-radius:6px;overflow-x:auto;font-size:13px}}.w{{color:#b45309}}
.t{{overflow-x:auto}}</style></head><body><div class="c">
<h1>🧬 Sanger Sequencing – Quality Report</h1><div class="s">{n} sample(s)</div>
<div class="cards"><div class="card"><b>{n}</b><span>Samples</span></div>
<div class="card"><b>{avg('consensus_length'):.0f}</b><span>Avg consensus (bp)</span></div>
<div class="card"><b>{avg('mean_consensus_quality'):.1f}</b><span>Avg consensus Phred</span></div>
<div class="card"><b>{avg('q20_pct'):.1f}%</b><span>Bases ≥ Q20</span></div>
<div class="card"><b>{avg('ambiguity_pct'):.1f}%</b><span>Avg ambiguity</span></div></div>
<div class="t"><table><thead><tr><th>Sample</th><th>Mode</th><th>Length</th><th>Mean Q</th><th>Q20 %</th>
<th>Ambig %</th><th>F/R disagreements</th>{'<th>Variants</th>' if has_ref else ''}<th>Rating</th></tr></thead>
<tbody>{rows}</tbody></table></div>
<h2>Consensus sequences</h2>{details}</div></body></html>"""
    with open(os.path.join(output_dir, "quality_report.html"), "w", encoding="utf-8") as f:
        f.write(page)


def make_zip(output_dir) -> str:
    """Zip the whole output directory next to it; returns the zip path."""
    return shutil.make_archive(output_dir.rstrip("/\\") + "_results", "zip", output_dir)


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------
def main():
    p = argparse.ArgumentParser(description="Sanger batch wrapper")
    p.add_argument("-i", "--input", nargs="+", required=True, help="AB1 files and/or directories")
    p.add_argument("-o", "--output", default="output")
    p.add_argument("-q", "--trim-quality", type=float, default=0.05)
    p.add_argument("--reference", help="Reference FASTA (optional)")
    p.add_argument("--dry-run", action="store_true")
    a = p.parse_args()

    files = find_ab1(a.input)
    jobs, unmatched = pair_files(files)
    print(f"{len(files)} AB1 files -> {len(jobs)} samples")
    for s, f, r in jobs:
        print(f"  {s}: F={os.path.basename(f) if f else '-'}  R={os.path.basename(r) if r else '-'}")
    for u in unmatched:
        print(f"  ? not recognised: {os.path.basename(u)}")
    if a.dry_run or not jobs:
        return
    ref = wf.parse_reference(open(a.reference).read()) if a.reference else None
    stats, errors = run_batch(jobs, a.output, a.trim_quality, ref,
                              lambda i, n, s: print(f"  [{i}/{n}] {s}") if s else None)
    for s, m in errors:
        print(f"  ERROR {s}: {m}")
    print(f"Done: {len(stats)} ok, {len(errors)} failed -> {a.output}/ (open quality_report.html)")
    sys.exit(1 if errors and not stats else 0)


if __name__ == "__main__":
    main()
