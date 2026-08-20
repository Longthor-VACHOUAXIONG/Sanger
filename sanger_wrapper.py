#!/usr/bin/env python3
"""
Sanger Workflow Wrapper
=======================
Convenient wrapper to run the Sanger sequencing workflow on custom AB1 files.

Features:
  - Auto-detects forward/reverse pairs from naming conventions
  - Supports individual files, directories, or mixed input
  - Batch mode for processing multiple samples at once
  - Flexible output naming

Naming conventions recognized for auto-pairing:
  *_forward.ab1 / *_reverse.ab1
  *_Forward.ab1 / *_Reverse.ab1
  *_fwd.ab1     / *_rev.ab1
  *_F.ab1       / *_R.ab1
  *_forward*.ab1 / *_reverse*.ab1  (e.g. sample_forward_read1.ab1)

Usage examples:
  # Auto-detect from directories
  python3 sanger_wrapper.py -d F/ R/ -o my_output

  # Pass individual paired files
  python3 sanger_wrapper.py -f sample1_forward.ab1 -r sample1_reverse.ab1

  # Pass multiple pairs
  python3 sanger_wrapper.py -p sample1_forward.ab1 sample1_reverse.ab1 -p sample2_forward.ab1 sample2_reverse.ab1

  # Auto-detect from flat directory with mixed files
  python3 sanger_wrapper.py --auto-dir all_ab1_files/ -o results

  # Dry run (show what would be processed)
  python3 sanger_wrapper.py --auto-dir all_ab1_files/ --dry-run
"""

import os
import sys
import re
import glob
import argparse
import shutil
import tempfile
from pathlib import Path
from collections import defaultdict


# ---------------------------------------------------------------------------
# Auto-detection patterns
# ---------------------------------------------------------------------------
FORWARD_PATTERNS = [
    r"(.+)_forward(?:_.*)?\.ab1$",
    r"(.+)_Forward(?:_.*)?\.ab1$",
    r"(.+)_fwd(?:_.*)?\.ab1$",
    r"(.+)_F\.ab1$",
    r"(.+)_f\.ab1$",
    r"(.+)_(\d+)_forward\.ab1$",      # e.g. sample_01_forward.ab1
    r"(.+)_(\d+)_Forward\.ab1$",
    r"(.+)_(\d+)_F\.ab1$",
]

REVERSE_PATTERNS = [
    r"(.+)_reverse(?:_.*)?\.ab1$",
    r"(.+)_Reverse(?:_.*)?\.ab1$",
    r"(.+)_rev(?:_.*)?\.ab1$",
    r"(.+)_R\.ab1$",
    r"(.+)_r\.ab1$",
    r"(.+)_(\d+)_reverse\.ab1$",
    r"(.+)_(\d+)_Reverse\.ab1$",
    r"(.+)_(\d+)_R\.ab1$",
]


def find_ab1_files(paths):
    """Collect all .ab1 files from given paths (files or directories)."""
    files = []
    for p in paths:
        p = str(p)
        if os.path.isfile(p) and p.endswith(".ab1"):
            files.append(os.path.abspath(p))
        elif os.path.isdir(p):
            for f in sorted(glob.glob(os.path.join(p, "*.ab1"))):
                files.append(os.path.abspath(f))
    return files


def extract_sample_name(filepath, patterns):
    """Try to extract a sample base name from a filepath using patterns."""
    basename = os.path.basename(filepath)
    for pattern in patterns:
        m = re.match(pattern, basename, re.IGNORECASE)
        if m:
            return m.group(1)
    return None


def auto_pair_files(ab1_files):
    """
    Auto-detect forward/reverse pairs from a flat list of AB1 files.
    Returns (forward_files, reverse_files, unmatched_files).
    """
    forward = {}
    reverse = {}
    unmatched = []

    for f in ab1_files:
        basename = os.path.basename(f)

        # Try forward patterns
        paired = False
        for pattern in FORWARD_PATTERNS:
            m = re.match(pattern, basename, re.IGNORECASE)
            if m:
                sample = m.group(1)
                forward[sample] = f
                paired = True
                break

        if paired:
            continue

        # Try reverse patterns
        for pattern in REVERSE_PATTERNS:
            m = re.match(pattern, basename, re.IGNORECASE)
            if m:
                sample = m.group(1)
                reverse[sample] = f
                paired = True
                break

        if not paired:
            unmatched.append(f)

    return forward, reverse, unmatched


def validate_pairs(forward, reverse):
    """Check that forward and reverse sets match and report issues."""
    fwd_samples = set(forward.keys())
    rev_samples = set(reverse.keys())

    paired = fwd_samples & rev_samples
    fwd_only = fwd_samples - rev_samples
    rev_only = rev_samples - fwd_samples

    return paired, fwd_only, rev_only


def print_detection_report(forward, reverse, unmatched):
    """Print a detailed report of auto-detected pairs."""
    paired, fwd_only, rev_only = validate_pairs(forward, reverse)

    print("\n  Auto-detection report:")
    print(f"  {'─' * 60}")

    if paired:
        print(f"  ✓ Matched pairs ({len(paired)}):")
        for sample in sorted(paired):
            fwd_name = os.path.basename(forward[sample])
            rev_name = os.path.basename(reverse[sample])
            print(f"    {sample}")
            print(f"      Forward: {fwd_name}")
            print(f"      Reverse: {rev_name}")

    if fwd_only:
        print(f"\n  ⚠ Forward-only (no reverse match) ({len(fwd_only)}):")
        for sample in sorted(fwd_only):
            print(f"    {sample} → {os.path.basename(forward[sample])}")

    if rev_only:
        print(f"\n  ⚠ Reverse-only (no forward match) ({len(rev_only)}):")
        for sample in sorted(rev_only):
            print(f"    {sample} → {os.path.basename(reverse[sample])}")

    if unmatched:
        print(f"\n  ✗ Unmatched files ({len(unmatched)}):")
        for f in unmatched:
            print(f"    {os.path.basename(f)}")

    print(f"  {'─' * 60}")
    return paired


def prepare_temp_dirs(forward_files, reverse_files, output_dir):
    """
    Create temporary F/ and R/ directories with symlinks/copies
    so the workflow can find them.
    """
    tmp_forward = os.path.join(output_dir, "_tmp_forward")
    tmp_reverse = os.path.join(output_dir, "_tmp_reverse")

    os.makedirs(tmp_forward, exist_ok=True)
    os.makedirs(tmp_reverse, exist_ok=True)

    for f in forward_files:
        dest = os.path.join(tmp_forward, os.path.basename(f))
        if not os.path.exists(dest):
            shutil.copy2(f, dest)

    for f in reverse_files:
        dest = os.path.join(tmp_reverse, os.path.basename(f))
        if not os.path.exists(dest):
            shutil.copy2(f, dest)

    return tmp_forward, tmp_reverse


def cleanup_temp_dirs(*dirs):
    """Remove temporary directories."""
    for d in dirs:
        if os.path.exists(d):
            shutil.rmtree(d, ignore_errors=True)


def run_workflow(forward_dir, reverse_dir, output_dir, trim_quality=0.05):
    """Run the sanger_workflow.py script with given parameters."""
    script_dir = os.path.dirname(os.path.abspath(__file__))
    workflow_script = os.path.join(script_dir, "sanger_workflow.py")

    cmd = [
        sys.executable, workflow_script,
        "--forward", forward_dir,
        "--reverse", reverse_dir,
        "--output", output_dir,
        "--trim-quality", str(trim_quality),
    ]

    print(f"\n  Running: {' '.join(cmd)}\n")
    result = os.spawnvp(os.P_WAIT, sys.executable, [sys.executable] + cmd[1:])
    return result


def collect_sample_stats(output_dir):
    """Collect sample_stats.json from all subdirectories."""
    import json
    stats_list = []
    for root, dirs, files in os.walk(output_dir):
        if "sample_stats.json" in files:
            with open(os.path.join(root, "sample_stats.json")) as f:
                stats = json.load(f)
                stats["_path"] = root
                stats_list.append(stats)
    return sorted(stats_list, key=lambda s: s.get("sample_name", ""))


def combine_consensus_files(output_dir):
    """
    Combine all per-sample consensus.fasta files into a single all_consensus.fasta.
    Equivalent to: cat */consensus.fasta > all_consensus.fasta
    """
    import glob
    consensus_files = sorted(glob.glob(os.path.join(output_dir, "*", "consensus.fasta")))

    if not consensus_files:
        print("  No per-sample consensus.fasta files found, skipping.")
        return None

    combined_path = os.path.join(output_dir, "all_consensus.fasta")
    total_seqs = 0

    with open(combined_path, "w") as out:
        for cf in consensus_files:
            with open(cf) as f:
                content = f.read().strip()
                if content:
                    out.write(content + "\n")
                    total_seqs += 1

    print(f"  Combined {total_seqs} consensus sequences → {combined_path}")
    return combined_path


def generate_html_report(output_dir, stats_list):
    """Generate an HTML quality report from collected stats."""
    import json

    if not stats_list:
        print("  No sample_stats.json files found, skipping report.")
        return None

    n = len(stats_list)

    # Compute summary stats
    avg_ambiguity = sum(s["ambiguity_pct"] for s in stats_list) / n
    avg_gc = sum(s["gc_content"] for s in stats_list) / n
    avg_coverage = sum(s["coverage_pct"] for s in stats_list) / n
    avg_consensus_len = sum(s["consensus_length"] for s in stats_list) / n

    # Quality distribution
    rating_counts = {}
    for s in stats_list:
        r = s["quality_rating"]
        rating_counts[r] = rating_counts.get(r, 0) + 1

    # Bar chart data
    bar_data = json.dumps([
        {"name": s["sample_name"], "ambiguity": s["ambiguity_pct"], "gc": s["gc_content"], "coverage": s["coverage_pct"]}
        for s in stats_list
    ])

    # Build HTML
    html = f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">
<title>Sanger Workflow - Quality Report</title>
<style>
  * {{ margin: 0; padding: 0; box-sizing: border-box; }}
  body {{ font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, sans-serif; background: #f5f7fa; color: #333; padding: 20px; }}
  .container {{ max-width: 1200px; margin: 0 auto; }}
  h1 {{ color: #1a1a2e; margin-bottom: 8px; font-size: 28px; }}
  h2 {{ color: #16213e; margin: 30px 0 15px; font-size: 20px; border-bottom: 2px solid #e0e0e0; padding-bottom: 8px; }}
  .subtitle {{ color: #666; margin-bottom: 30px; font-size: 14px; }}
  .cards {{ display: grid; grid-template-columns: repeat(auto-fit, minmax(200px, 1fr)); gap: 16px; margin-bottom: 30px; }}
  .card {{ background: white; border-radius: 10px; padding: 20px; box-shadow: 0 2px 8px rgba(0,0,0,0.08); text-align: center; }}
  .card .value {{ font-size: 32px; font-weight: 700; color: #1a1a2e; }}
  .card .label {{ font-size: 13px; color: #888; margin-top: 4px; }}
  .card.excellent .value {{ color: #27ae60; }}
  .card.good .value {{ color: #2ecc71; }}
  .card.fair .value {{ color: #f39c12; }}
  .card.poor .value {{ color: #e74c3c; }}
  table {{ width: 100%; border-collapse: collapse; background: white; border-radius: 10px; overflow: hidden; box-shadow: 0 2px 8px rgba(0,0,0,0.08); margin-bottom: 30px; }}
  th {{ background: #1a1a2e; color: white; padding: 12px 16px; text-align: left; font-size: 13px; text-transform: uppercase; letter-spacing: 0.5px; }}
  td {{ padding: 10px 16px; border-bottom: 1px solid #f0f0f0; font-size: 14px; }}
  tr:hover {{ background: #f8f9ff; }}
  .rating {{ padding: 4px 10px; border-radius: 12px; font-size: 12px; font-weight: 600; display: inline-block; }}
  .rating.Excellent {{ background: #d4edda; color: #155724; }}
  .rating.Good {{ background: #d1ecf1; color: #0c5460; }}
  .rating.Fair {{ background: #fff3cd; color: #856404; }}
  .rating.Poor {{ background: #f8d7da; color: #721c24; }}
  .bar-chart {{ background: white; border-radius: 10px; padding: 20px; box-shadow: 0 2px 8px rgba(0,0,0,0.08); margin-bottom: 30px; }}
  .bar-group {{ margin-bottom: 12px; }}
  .bar-label {{ font-size: 13px; color: #555; margin-bottom: 4px; }}
  .bar-track {{ height: 24px; background: #eee; border-radius: 4px; position: relative; overflow: hidden; }}
  .bar-fill {{ height: 100%; border-radius: 4px; transition: width 0.5s ease; display: flex; align-items: center; padding-left: 8px; font-size: 11px; color: white; font-weight: 600; min-width: 30px; }}
  .bar-fill.ambiguity {{ background: linear-gradient(90deg, #e74c3c, #c0392b); }}
  .bar-fill.gc {{ background: linear-gradient(90deg, #3498db, #2980b9); }}
  .bar-fill.coverage {{ background: linear-gradient(90deg, #27ae60, #229954); }}
  .legend {{ display: flex; gap: 20px; margin-bottom: 15px; font-size: 13px; }}
  .legend-item {{ display: flex; align-items: center; gap: 6px; }}
  .legend-dot {{ width: 12px; height: 12px; border-radius: 3px; }}
  .detail {{ background: white; border-radius: 10px; padding: 20px; box-shadow: 0 2px 8px rgba(0,0,0,0.08); margin-bottom: 20px; }}
  .detail h3 {{ margin-bottom: 12px; font-size: 16px; color: #1a1a2e; }}
  .detail-grid {{ display: grid; grid-template-columns: repeat(auto-fit, minmax(180px, 1fr)); gap: 10px; }}
  .detail-item {{ padding: 8px 12px; background: #f8f9fa; border-radius: 6px; }}
  .detail-item .dl {{ font-size: 12px; color: #888; }}
  .detail-item .dv {{ font-size: 15px; font-weight: 600; color: #333; }}
  .footer {{ text-align: center; color: #aaa; font-size: 12px; margin-top: 40px; padding: 20px; }}
</style>
</head>
<body>
<div class="container">
  <h1>🧬 Sanger Sequencing - Quality Report</h1>
  <p class="subtitle">Generated by Sanger Workflow | {n} sample{'s' if n != 1 else ''} processed</p>

  <div class="cards">
    <div class="card"><div class="value">{n}</div><div class="label">Total Samples</div></div>
    <div class="card"><div class="value">{avg_consensus_len:.0f}</div><div class="label">Avg Consensus (bp)</div></div>
    <div class="card"><div class="value">{avg_gc:.1f}%</div><div class="label">Avg GC Content</div></div>
    <div class="card"><div class="value">{avg_ambiguity:.1f}%</div><div class="label">Avg Ambiguity</div></div>
    <div class="card"><div class="value">{avg_coverage:.1f}%</div><div class="label">Avg Coverage</div></div>
  </div>

  <h2>📊 Quality Distribution</h2>
  <div class="cards">
"""

    for rating in ["Excellent", "Good", "Fair", "Poor"]:
        count = rating_counts.get(rating, 0)
        cls = rating.lower()
        html += f'    <div class="card {cls}"><div class="value">{count}</div><div class="label">{rating}</div></div>\n'

    html += """  </div>

  <h2>📋 Sample Comparison</h2>
  <table>
    <thead>
      <tr>
        <th>Sample</th>
        <th>Fwd (bp)</th>
        <th>Rev (bp)</th>
        <th>Trimmed F/R</th>
        <th>Alignment</th>
        <th>Consensus</th>
        <th>GC %</th>
        <th>Ambiguity %</th>
        <th>Coverage %</th>
        <th>Fwd QV</th>
        <th>Rev QV</th>
        <th>Rating</th>
      </tr>
    </thead>
    <tbody>
"""

    for s in stats_list:
        rating = s["quality_rating"]
        html += f"""      <tr>
        <td><strong>{s['sample_name']}</strong></td>
        <td>{s['forward_raw_length']}</td>
        <td>{s['reverse_raw_length']}</td>
        <td>{s['forward_trimmed_length']} / {s['reverse_trimmed_length']}</td>
        <td>{s['alignment_length']}</td>
        <td>{s['consensus_length']}</td>
        <td>{s['gc_content']:.1f}</td>
        <td>{s['ambiguity_pct']:.1f}</td>
        <td>{s['coverage_pct']:.1f}</td>
        <td>{s['forward_avg_quality']:.1f}</td>
        <td>{s['reverse_avg_quality']:.1f}</td>
        <td><span class="rating {rating}">{rating}</span></td>
      </tr>
"""

    html += """    </tbody>
  </table>

  <h2>📈 Metrics Comparison</h2>
  <div class="bar-chart">
    <div class="legend">
      <div class="legend-item"><div class="legend-dot" style="background:#e74c3c"></div> Ambiguity %</div>
      <div class="legend-item"><div class="legend-dot" style="background:#3498db"></div> GC Content %</div>
      <div class="legend-item"><div class="legend-dot" style="background:#27ae60"></div> Coverage %</div>
    </div>
"""

    for s in stats_list:
        name = s["sample_name"]
        amb = s["ambiguity_pct"]
        gc = s["gc_content"]
        cov = s["coverage_pct"]
        html += f"""    <div class="bar-group">
      <div class="bar-label"><strong>{name}</strong></div>
      <div class="bar-track"><div class="bar-fill ambiguity" style="width:{max(amb, 3):.1f}%">{amb:.1f}%</div></div>
      <div class="bar-track"><div class="bar-fill gc" style="width:{max(gc, 3):.1f}%">{gc:.1f}%</div></div>
      <div class="bar-track"><div class="bar-fill coverage" style="width:{max(cov, 3):.1f}%">{cov:.1f}%</div></div>
    </div>
"""

    html += """  </div>

  <h2>🔍 Per-Sample Details</h2>
"""

    for s in stats_list:
        html += f"""  <div class="detail">
    <h3>{s['sample_name']} <span class="rating {s['quality_rating']}">{s['quality_rating']}</span></h3>
    <div class="detail-grid">
      <div class="detail-item"><div class="dl">Forward File</div><div class="dv">{s['forward_file']}</div></div>
      <div class="detail-item"><div class="dl">Reverse File</div><div class="dv">{s['reverse_file']}</div></div>
      <div class="detail-item"><div class="dl">Raw Length (F/R)</div><div class="dv">{s['forward_raw_length']} / {s['reverse_raw_length']} bp</div></div>
      <div class="detail-item"><div class="dl">Trimmed Length (F/R)</div><div class="dv">{s['forward_trimmed_length']} / {s['reverse_trimmed_length']} bp</div></div>
      <div class="detail-item"><div class="dl">Alignment Length</div><div class="dv">{s['alignment_length']} bp</div></div>
      <div class="detail-item"><div class="dl">Consensus Length</div><div class="dv">{s['consensus_length']} bp</div></div>
      <div class="detail-item"><div class="dl">Valid Bases</div><div class="dv">{s['valid_bases']} ({s['coverage_pct']:.1f}%)</div></div>
      <div class="detail-item"><div class="dl">Ambiguous Bases</div><div class="dv">{s['ambiguous_bases']} ({s['ambiguity_pct']:.1f}%)</div></div>
      <div class="detail-item"><div class="dl">N Bases</div><div class="dv">{s['n_bases']} ({s['n_pct']:.1f}%)</div></div>
      <div class="detail-item"><div class="dl">GC Content</div><div class="dv">{s['gc_content']:.1f}%</div></div>
      <div class="detail-item"><div class="dl">Forward Avg Quality</div><div class="dv">Phred {s['forward_avg_quality']:.1f}</div></div>
      <div class="detail-item"><div class="dl">Reverse Avg Quality</div><div class="dv">Phred {s['reverse_avg_quality']:.1f}</div></div>
      <div class="detail-item"><div class="dl">Sequences Aligned</div><div class="dv">{s['num_sequences_aligned']}</div></div>
    </div>
  </div>
"""

    html += f"""
  <div class="footer">
    Report generated by Sanger Workflow
  </div>
</div>
</body>
</html>
"""

    report_path = os.path.join(output_dir, "quality_report.html")
    with open(report_path, "w") as f:
        f.write(html)

    # Also save raw JSON
    json_path = os.path.join(output_dir, "all_sample_stats.json")
    with open(json_path, "w") as f:
        json.dump(stats_list, f, indent=2, default=str)

    print(f"\n  📊 Quality report: {report_path}")
    print(f"  📄 Raw stats JSON: {json_path}")
    return report_path


# ===========================================================================
# MAIN
# ===========================================================================
def main():
    parser = argparse.ArgumentParser(
        description="Sanger Workflow Wrapper - Run Sanger workflow on custom AB1 files",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Auto-detect pairs from directories (per-sample output)
  %(prog)s -d forward_files/ reverse_files/ -o results
  #   Creates: results/sample1/, results/sample2/, ...

  # Auto-detect from a mixed directory (per-sample output)
  %(prog)s --auto-dir all_ab1_files/ -o results
  #   Creates: results/sample1/, results/sample2/, ...

  # Single pair (flat output)
  %(prog)s -f sample_forward.ab1 -r sample_reverse.ab1 -o my_output

  # Multiple pairs (per-sample output)
  %(prog)s -p fwd1.ab1 rev1.ab1 -p fwd2.ab1 rev2.ab1 -o results
  #   Creates: results/fwd1/, results/fwd2/, ...

  # Force flat output for multi-pair (all in one dir)
  %(prog)s --auto-dir all_files/ --flat -o results

  # Dry run
  %(prog)s --auto-dir all_files/ --dry-run
        """,
    )

    # Input modes (mutually exclusive)
    input_group = parser.add_argument_group("Input modes (choose one)")

    input_group.add_argument(
        "-d", "--dirs", nargs=2, metavar=("FORWARD_DIR", "REVERSE_DIR"),
        help="Two directories: first for forward AB1 files, second for reverse"
    )
    input_group.add_argument(
        "--auto-dir",
        help="Single directory with mixed forward/reverse AB1 files (auto-detect pairs)"
    )
    input_group.add_argument(
        "-f", "--forward",
        help="Single forward AB1 file (use with -r)"
    )
    input_group.add_argument(
        "-r", "--reverse",
        help="Single reverse AB1 file (use with -f)"
    )
    input_group.add_argument(
        "-p", "--pair", nargs=2, action="append", metavar=("FORWARD", "REVERSE"),
        help="Explicit forward/reverse pair (can be used multiple times)"
    )

    # Options
    parser.add_argument("-o", "--output", default="output", help="Output directory (default: output/)")
    parser.add_argument("-q", "--trim-quality", type=float, default=0.05, help="Trim quality threshold (default: 0.05)")
    parser.add_argument("--dry-run", action="store_true", help="Show what would be processed without running")
    parser.add_argument("--no-cleanup", action="store_true", help="Keep temporary files")
    parser.add_argument("--flat", action="store_true", help="Put all results in one directory instead of per-sample subdirs")
    parser.add_argument("--batch", action="store_true", help="Process ALL files together. Produces 1 consensus from all reads.")

    args = parser.parse_args()

    # Check that at least one input mode is specified
    if not any([args.dirs, args.auto_dir, args.forward, args.pair]):
        parser.error("Specify input with -d, --auto-dir, -f/-r, or -p")

    print("=" * 70)
    print("SANGER WORKFLOW WRAPPER")
    print("=" * 70)

    # ---- BATCH MODE: Process all files together (matches Galaxy) ----
    if args.batch:
        print(f"\nMode: BATCH (process all files together)")

        all_fwd = []
        all_rev = []

        if args.dirs:
            fwd_dir, rev_dir = args.dirs
            all_fwd = find_ab1_files([fwd_dir])
            all_rev = find_ab1_files([rev_dir])
        elif args.auto_dir:
            all_files = find_ab1_files([args.auto_dir])
            forward, reverse, unmatched = auto_pair_files(all_files)
            for sample in sorted(forward.keys()):
                all_fwd.append(forward[sample])
            for sample in sorted(reverse.keys()):
                all_rev.append(reverse[sample])
        elif args.pair:
            for f, r in args.pair:
                all_fwd.append(f)
                all_rev.append(r)
        elif args.forward and args.reverse:
            all_fwd = [args.forward]
            all_rev = [args.reverse]
        else:
            print("  ERROR: --batch requires input files. Use -d, --auto-dir, -f/-r, or -p")
            sys.exit(1)

        if not all_fwd or not all_rev:
            print("  ERROR: No forward/reverse files found")
            sys.exit(1)

        print(f"  Found: {len(all_fwd)} forward, {len(all_rev)} reverse")
        print(f"  Processing ALL files together → 1 consensus")

        if args.dry_run:
            print("\n  Dry run - would process:")
            for f in all_fwd:
                print(f"    F: {os.path.basename(f)}")
            for f in all_rev:
                print(f"    R: {os.path.basename(f)}")
            return

        # Create temp dirs with ALL files
        tmp_fwd, tmp_rev = prepare_temp_dirs(all_fwd, all_rev, args.output)
        try:
            run_workflow(tmp_fwd, tmp_rev, args.output, args.trim_quality)
        finally:
            if not args.no_cleanup:
                cleanup_temp_dirs(tmp_fwd, tmp_rev)

        # Generate report
        print("\n" + "=" * 70)
        print("GENERATING QUALITY REPORT")
        print("=" * 70)
        stats = collect_sample_stats(args.output)
        if stats:
            generate_html_report(args.output, stats)
        print("\n" + "=" * 70)
        print("DONE!")
        print("=" * 70)
        return

    # ---- Mode 1: Two directories ----
    if args.dirs:
        fwd_dir, rev_dir = args.dirs
        print(f"\nMode: Directory pair")
        print(f"  Forward: {fwd_dir}/")
        print(f"  Reverse: {rev_dir}/")

        fwd_files = find_ab1_files([fwd_dir])
        rev_files = find_ab1_files([rev_dir])

        if not fwd_files:
            print(f"\n  ERROR: No .ab1 files found in {fwd_dir}/")
            sys.exit(1)
        if not rev_files:
            print(f"\n  ERROR: No .ab1 files found in {rev_dir}/")
            sys.exit(1)

        print(f"  Found: {len(fwd_files)} forward, {len(rev_files)} reverse")

        # Pair files by matching names
        fwd_map = {os.path.basename(f).replace("_forward", "").replace("_Forward", "").replace("_fwd", "").replace("_F", ""): f for f in fwd_files}
        rev_map = {os.path.basename(f).replace("_reverse", "").replace("_Reverse", "").replace("_rev", "").replace("_R", ""): f for f in rev_files}
        common = sorted(set(fwd_map.keys()) & set(rev_map.keys()))

        if not common:
            print("  WARNING: No matching pairs found by name, running all forward + all reverse together")
            if args.dry_run:
                print("\n  Dry run - would process all files together")
                return
            run_workflow(fwd_dir, rev_dir, args.output, args.trim_quality)
        else:
            print(f"  Matched {len(common)} pairs")
            if args.dry_run:
                print("\n  Dry run - would process:")
                for s in common:
                    print(f"    {s}: {os.path.basename(fwd_map[s])} + {os.path.basename(rev_map[s])}")
                return
            for sample in common:
                sample_dir = os.path.join(args.output, sample)
                tmp_fwd, tmp_rev = prepare_temp_dirs([fwd_map[sample]], [rev_map[sample]], sample_dir)
                try:
                    run_workflow(tmp_fwd, tmp_rev, sample_dir, args.trim_quality)
                finally:
                    if not args.no_cleanup:
                        cleanup_temp_dirs(tmp_fwd, tmp_rev)

    # ---- Mode 2: Auto-detect from single directory ----
    elif args.auto_dir:
        print(f"\nMode: Auto-detect from {args.auto_dir}/")

        all_files = find_ab1_files([args.auto_dir])
        if not all_files:
            print(f"\n  ERROR: No .ab1 files found in {args.auto_dir}/")
            sys.exit(1)

        print(f"  Found {len(all_files)} AB1 files")

        forward, reverse, unmatched = auto_pair_files(all_files)
        paired_samples = print_detection_report(forward, reverse, unmatched)

        if not paired_samples:
            print("\n  ERROR: No valid forward/reverse pairs found!")
            print("  Tip: Use naming convention like sample_forward.ab1 / sample_reverse.ab1")
            sys.exit(1)

        if args.dry_run:
            print(f"\n  Dry run - would process {len(paired_samples)} pairs into {args.output}/")
            for sample in sorted(paired_samples):
                print(f"    {args.output}/{sample}/")
            return

        for sample in sorted(paired_samples):
            sample_dir = os.path.join(args.output, sample)
            tmp_fwd, tmp_rev = prepare_temp_dirs([forward[sample]], [reverse[sample]], sample_dir)
            try:
                run_workflow(tmp_fwd, tmp_rev, sample_dir, args.trim_quality)
            finally:
                if not args.no_cleanup:
                    cleanup_temp_dirs(tmp_fwd, tmp_rev)

    # ---- Mode 3: Single pair ----
    elif args.forward and args.reverse:
        print(f"\nMode: Single pair")
        print(f"  Forward: {args.forward}")
        print(f"  Reverse: {args.reverse}")

        if not os.path.isfile(args.forward):
            print(f"\n  ERROR: Forward file not found: {args.forward}")
            sys.exit(1)
        if not os.path.isfile(args.reverse):
            print(f"\n  ERROR: Reverse file not found: {args.reverse}")
            sys.exit(1)

        if args.dry_run:
            print(f"\n  Dry run - would process into {args.output}/")
            return

        tmp_fwd, tmp_rev = prepare_temp_dirs([args.forward], [args.reverse], args.output)

        try:
            run_workflow(tmp_fwd, tmp_rev, args.output, args.trim_quality)
        finally:
            if not args.no_cleanup:
                cleanup_temp_dirs(tmp_fwd, tmp_rev)

    # ---- Mode 4: Multiple explicit pairs ----
    elif args.pair:
        print(f"\nMode: {len(args.pair)} explicit pair(s)")

        for i, (f, r) in enumerate(args.pair, 1):
            if not os.path.isfile(f):
                print(f"\n  ERROR: Forward file {i} not found: {f}")
                sys.exit(1)
            if not os.path.isfile(r):
                print(f"\n  ERROR: Reverse file {i} not found: {r}")
                sys.exit(1)
            print(f"  Pair {i}: {os.path.basename(f)} + {os.path.basename(r)}")

        if args.dry_run:
            print(f"\n  Dry run - would process {len(args.pair)} pairs into {args.output}/")
            for i, (f, r) in enumerate(args.pair, 1):
                name = Path(f).stem.replace("_forward", "").replace("_Forward", "").replace("_fwd", "").replace("_F", "")
                print(f"    {args.output}/{name}/")
            return

        for i, (f, r) in enumerate(args.pair, 1):
            # Derive sample name from forward file
            name = Path(f).stem.replace("_forward", "").replace("_Forward", "").replace("_fwd", "").replace("_F", "")
            if not name:
                name = f"pair_{i}"
            sample_dir = os.path.join(args.output, name)
            tmp_fwd, tmp_rev = prepare_temp_dirs([f], [r], sample_dir)
            try:
                run_workflow(tmp_fwd, tmp_rev, sample_dir, args.trim_quality)
            finally:
                if not args.no_cleanup:
                    cleanup_temp_dirs(tmp_fwd, tmp_rev)

    # Generate quality report from collected stats
    print("\n" + "=" * 70)
    print("GENERATING QUALITY REPORT")
    print("=" * 70)

    stats = collect_sample_stats(args.output)
    if stats:
        generate_html_report(args.output, stats)
    else:
        print("  No sample_stats.json found - report skipped.")

    # Combine all per-sample consensus files
    print("\n  Combining all consensus sequences...")
    combine_consensus_files(args.output)

    print("\n" + "=" * 70)
    print("DONE!")
    print("=" * 70)


if __name__ == "__main__":
    main()
