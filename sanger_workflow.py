#!/usr/bin/env python3
"""
Sanger Sequencing Workflow (v2)

For each sample (forward and/or reverse AB1):
  1. Read bases + Phred qualities with Biopython
  2. Quality trimming (Mott algorithm, same as `seqtk trimfq -q`)
  3. Reverse-complement the reverse read
  4. Overlap-align forward and reverse reads (pure Python, no MAFFT needed)
  5. Quality-aware consensus:
       - agree            -> base, qualities combined
       - disagree         -> higher quality wins if clearly better,
                             otherwise an IUPAC ambiguity code (R, Y, ...)
       - single-read ends -> that read's base
  6. Optional: compare consensus to a reference and list mutations
  7. Statistics + quality rating

Only dependency: biopython.
"""

import argparse
import csv
import json
import math
import os
import re
import shutil
import sys
import time
from collections import Counter

from Bio import SeqIO
from Bio.Align import PairwiseAligner

# ---------------------------------------------------------------------------
# IUPAC helpers
# ---------------------------------------------------------------------------
IUPAC_EXPANSION = {
    "A": "A", "C": "C", "G": "G", "T": "T",
    "R": "AG", "Y": "CT", "W": "AT", "S": "CG", "K": "GT", "M": "AC",
    "B": "CGT", "D": "AGT", "H": "ACT", "V": "ACG", "N": "ACGT",
}
IUPAC_CODE = {frozenset(v): k for k, v in IUPAC_EXPANSION.items()}
AMBIGUOUS = set("RYSWKMBDHVN")

_COMP = str.maketrans("ACGTRYSWKMBDHVNacgtryswkmbdhvn",
                      "TGCAYRSWMKVHDBNtgcayrswmkvhdbn")


def revcomp(seq: str) -> str:
    return seq.translate(_COMP)[::-1]


def _expand(base: str) -> str:
    return IUPAC_EXPANSION.get(base, "ACGT")


def _compatible(a: str, b: str) -> bool:
    return bool(set(_expand(a)) & set(_expand(b)))


# ---------------------------------------------------------------------------
# Reading + trimming
# ---------------------------------------------------------------------------
def read_ab1(path: str) -> dict:
    """Read an AB1 file -> {seq, quals, quality_missing}."""
    rec = SeqIO.read(path, "abi")
    seq = str(rec.seq).upper()
    quals = list(rec.letter_annotations.get("phred_quality") or [])
    missing = len(quals) != len(seq)
    if missing:
        quals = [20] * len(seq)
    return {"seq": seq, "quals": quals, "quality_missing": missing}


def mott_trim(quals, q: float = 0.05):
    """
    Modified Mott trimming (identical to seqtk trimfq).
    Returns (start, end) of the best-scoring segment; (0, 0) if none.
    """
    best, best_start, best_end = 0.0, 0, 0
    run, start = 0.0, 0
    for i, Q in enumerate(quals):
        run += q - min(1.0, 10 ** (-Q / 10))
        if run < 0:
            run, start = 0.0, i + 1
        elif run > best:
            best, best_start, best_end = run, start, i + 1
    return best_start, best_end


# ---------------------------------------------------------------------------
# Pairwise alignment helpers
# ---------------------------------------------------------------------------
def _make_aligner(mode):
    a = PairwiseAligner()
    a.mode = mode
    a.match_score = 2
    a.mismatch_score = -3
    a.open_gap_score = -5
    a.extend_gap_score = -2
    a.wildcard = "N"
    if mode == "global":  # free end gaps -> overlap alignment
        a.end_gap_score = 0.0
    return a


_OVERLAP = _make_aligner("global")
_LOCAL = _make_aligner("local")


def _columns(aln):
    """Alignment -> list of (target_idx|None, query_idx|None) per column."""
    c = aln.coordinates
    cols = []
    for k in range(c.shape[1] - 1):
        a0, a1 = int(c[0][k]), int(c[0][k + 1])
        b0, b1 = int(c[1][k]), int(c[1][k + 1])
        da, db = a1 - a0, b1 - b0
        if da == db:
            cols += [(a0 + i, b0 + i) for i in range(da)]
        elif db == 0:
            cols += [(a0 + i, None) for i in range(da)]
        elif da == 0:
            cols += [(None, b0 + i) for i in range(db)]
    return cols


# ---------------------------------------------------------------------------
# Consensus from two reads
# ---------------------------------------------------------------------------
def merge_reads(fseq, fq, fidx, rseq, rq, ridx, min_q=10, margin=10):
    """
    Merge forward read and (already reverse-complemented) reverse read.
    fidx / ridx map each base back to its 0-based position in the original AB1.
    Returns None if the reads cannot be aligned at all, otherwise a dict.
    """
    aln = _OVERLAP.align(fseq, rseq)[0]
    cols = _columns(aln)
    paired = [k for k, (a, b) in enumerate(cols) if a is not None and b is not None]
    if not paired:
        return {"ok": False, "overlap_len": 0, "identity": 0.0}

    first, last = paired[0], paired[-1]
    matches = sum(1 for k in paired if _compatible(fseq[cols[k][0]], rseq[cols[k][1]]))
    identity = matches / len(paired)
    overlap_len = last - first + 1
    if len(paired) < 20 or identity < 0.8:
        return {"ok": False, "overlap_len": overlap_len, "identity": identity}

    cb, cq, src = [], [], []
    fmap, rmap = [], []                    # original AB1 index per consensus base (-1 = none)
    frow, rrow, crow = [], [], []
    disagreements, mism, indels = [], 0, 0

    for k, (a, b) in enumerate(cols):
        fb = fseq[a] if a is not None else "-"
        rb = rseq[b] if b is not None else "-"
        frow.append(fb)
        rrow.append(rb)

        if a is not None and b is not None:
            q1, q2 = fq[a], rq[b]
            inter = set(_expand(fb)) & set(_expand(rb))
            if fb == rb:
                base, q, s = fb, min(q1 + q2, 60), "both"
            elif fb == "N":
                base, q, s = rb, q2, "R"
            elif rb == "N":
                base, q, s = fb, q1, "F"
            elif len(inter) == 1:
                base, q, s = next(iter(inter)), min(q1 + q2, 60), "both"
            else:
                mism += 1
                if q1 >= q2 + margin:
                    base, q, s, res = fb, q1, "F", "forward"
                elif q2 >= q1 + margin:
                    base, q, s, res = rb, q2, "R", "reverse"
                elif min(q1, q2) < min_q:
                    base, q, s, res = (fb, q1, "F", "forward") if q1 >= q2 else (rb, q2, "R", "reverse")
                else:
                    union = frozenset(set(_expand(fb)) | set(_expand(rb)))
                    base, q, s, res = IUPAC_CODE.get(union, "N"), min(q1, q2), "conflict", "IUPAC"
                disagreements.append({
                    "consensus_pos": len(cb) + 1, "call": base, "resolution": res,
                    "forward_pos": fidx[a] + 1, "forward_base": fb, "forward_q": q1,
                    "reverse_pos": ridx[b] + 1, "reverse_base": rb, "reverse_q": q2,
                })
        elif k < first or k > last:           # single-read overhang
            if a is not None:
                base, q, s = fb, fq[a], "F"
            else:
                base, q, s = rb, rq[b], "R"
        else:                                  # indel inside the overlap
            indels += 1
            if a is not None:
                base, q, s = fb, fq[a], "F"
            else:
                base, q, s = rb, rq[b], "R"
            if q < 20:                         # unreliable insertion -> drop
                crow.append("-")
                continue

        if q < min_q:
            base = "N"
        cb.append(base)
        cq.append(q)
        src.append(s)
        crow.append(base)
        fmap.append(fidx[a] if a is not None else -1)
        rmap.append(ridx[b] if b is not None else -1)

    return {
        "ok": True, "seq": "".join(cb), "quals": cq, "src": src,
        "rows": ("".join(frow), "".join(rrow), "".join(crow)), "fmap": fmap, "rmap": rmap,
        "overlap_len": overlap_len, "identity": identity,
        "disagreements": disagreements, "mismatches": mism, "indels": indels,
        "alignment_length": len(cols),
    }


# ---------------------------------------------------------------------------
# Reference comparison (mutation calling)
# ---------------------------------------------------------------------------
def parse_reference(text):
    """Accept FASTA or a raw sequence. Returns (name, seq) or None if empty."""
    text = (text or "").strip()
    if not text:
        return None
    name, lines = "reference", []
    for line in text.splitlines():
        line = line.strip()
        if line.startswith(">"):
            if lines:
                break
            if line[1:].strip():
                name = line[1:].split()[0]
        elif line:
            lines.append(line)
    seq = re.sub(r"[^A-Za-z]", "", "".join(lines)).upper().replace("U", "T")
    if len(seq) < 20 or sum(c in "ACGTN" for c in seq) / len(seq) < 0.9:
        raise ValueError("needs at least 20 nucleotides (A, C, G, T)")
    return name, seq


def compare_to_reference(cons: str, quals, ref_name: str, ref: str) -> dict:
    """Local-align consensus to the reference (both strands) and list variants."""
    n = len(cons)
    best = None
    for strand, s in (("+", cons), ("-", revcomp(cons))):
        aln = _LOCAL.align(ref, s)[0]
        if best is None or aln.score > best[0].score:
            best = (aln, strand)
    aln, strand = best
    cols = _columns(aln)
    r0 = int(aln.coordinates[0][0])   # coordinates are already absolute indices
    qry = cons if strand == "+" else revcomp(cons)

    def cons_pos(j):       # 0-based index in the query -> 1-based consensus position
        return j + 1 if strand == "+" else n - j

    def qual(j):
        return quals[j if strand == "+" else n - 1 - j]

    variants, matches, aligned = [], 0, 0
    i = 0
    while i < len(cols):
        a, b = cols[i]
        if a is not None and b is not None:
            aligned += 1
            rb, qb = ref[a], qry[b]
            if rb == qb or (qb in AMBIGUOUS and rb in _expand(qb)):
                matches += 1 if rb == qb else 0
            elif qb == "N":
                pass
            else:
                kind = "AMBIGUOUS" if qb in AMBIGUOUS else "SNP"
                variants.append({"type": kind, "ref_pos": a + 1, "ref": rb, "alt": qb,
                                 "consensus_pos": cons_pos(b), "quality": qual(b),
                                 "change": f"{rb}{a + 1}{qb}"})
            i += 1
        else:
            j = i
            while j < len(cols) and ((a is None) == (cols[j][0] is None)) and ((b is None) == (cols[j][1] is None)):
                j += 1
            run = cols[i:j]
            if a is None:   # insertion in sample
                bases = "".join(qry[x[1]] for x in run)
                p = next((c[0] for c in reversed(cols[:i]) if c[0] is not None), r0 - 1)
                variants.append({"type": "INS", "ref_pos": p + 1, "ref": "-", "alt": bases,
                                 "consensus_pos": cons_pos(run[0][1]),
                                 "quality": min(qual(x[1]) for x in run),
                                 "change": f"ins {bases} after {p + 1}"})
            else:           # deletion in sample
                bases = "".join(ref[x[0]] for x in run)
                p = run[0][0]
                nb = next((c[1] for c in cols[i:] if c[1] is not None), None)
                variants.append({"type": "DEL", "ref_pos": p + 1, "ref": bases, "alt": "-",
                                 "consensus_pos": cons_pos(nb) if nb is not None else 0,
                                 "quality": qual(nb) if nb is not None else 0,
                                 "change": f"del {bases} at {p + 1}"})
            i = j

    ref_cols = [c[0] for c in cols if c[0] is not None]
    ref_start, ref_end = (min(ref_cols) + 1, max(ref_cols) + 1) if ref_cols else (0, 0)
    return {
        "reference_name": ref_name, "reference_length": len(ref), "strand": strand,
        "identity": round(100 * matches / aligned, 2) if aligned else 0.0,
        "ref_start": ref_start, "ref_end": ref_end,
        "ref_coverage_pct": round(100 * (ref_end - ref_start + 1) / len(ref), 1) if ref_cols else 0.0,
        "variants": variants,
    }


# ---------------------------------------------------------------------------
# Rating + outputs
# ---------------------------------------------------------------------------
def rate(amb_pct, n_pct, q20_pct):
    if amb_pct < 1 and n_pct < 1 and q20_pct >= 90:
        return "Excellent"
    if amb_pct < 3 and n_pct < 3 and q20_pct >= 75:
        return "Good"
    if amb_pct < 5 and n_pct < 5 and q20_pct >= 50:
        return "Fair"
    return "Poor"


def _write_fastq(path, name, seq, quals):
    with open(path, "w") as f:
        f.write(f"@{name}\n{seq}\n+\n{''.join(chr(min(q, 93) + 33) for q in quals)}\n")


def _mean(x):
    return sum(x) / len(x) if x else 0.0


# ---------------------------------------------------------------------------
# Main per-sample entry point
# ---------------------------------------------------------------------------
def run_sample(sample, forward, reverse, out_dir, trim_quality=0.05,
               reference=None, min_base_q=10):
    """
    Process one sample. `forward` / `reverse` are AB1 paths (either may be None).
    `reference` is (name, seq) or None. Returns the stats dict (also saved as
    sample_stats.json). Raises ValueError if no usable read is available.
    """
    t0 = time.perf_counter()
    if not forward and not reverse:
        raise ValueError("No AB1 file given")
    os.makedirs(os.path.join(out_dir, "ab1"), exist_ok=True)
    warnings, reads = [], {}

    for key, path in (("forward", forward), ("reverse", reverse)):
        if not path:
            continue
        shutil.copy2(path, os.path.join(out_dir, "ab1", os.path.basename(path)))
        raw = read_ab1(path)
        if raw["quality_missing"]:
            warnings.append(f"{key}: no quality scores in file, assumed Q20")
        s, e = mott_trim(raw["quals"], trim_quality)
        if e - s < 20:
            warnings.append(f"{key}: read is too short/poor after trimming ({e - s} bp) - ignored")
            reads[key] = {"file": os.path.basename(path), "raw": raw, "trim": (s, e), "usable": False}
            continue
        seq, quals = raw["seq"][s:e], raw["quals"][s:e]
        idx = list(range(s, e))
        if key == "reverse":
            seq, quals, idx = revcomp(seq), quals[::-1], idx[::-1]
        reads[key] = {"file": os.path.basename(path), "raw": raw, "trim": (s, e), "usable": True,
                      "seq": seq, "quals": quals, "idx": idx}
        _write_fastq(os.path.join(out_dir, f"{key}_trimmed.fastq"), f"{sample}_{key}",
                     raw["seq"][s:e], raw["quals"][s:e])

    usable = [k for k in ("forward", "reverse") if reads.get(k, {}).get("usable")]
    if not usable:
        raise ValueError("; ".join(warnings) or "No usable reads")

    merged, mode, disagreements = None, "", []
    overlap_len, identity, mism, indels, aln_len = 0, 0.0, 0, 0, 0
    if len(usable) == 2:
        f, r = reads["forward"], reads["reverse"]
        merged = merge_reads(f["seq"], f["quals"], f["idx"], r["seq"], r["quals"], r["idx"], min_base_q)
        overlap_len, identity = merged["overlap_len"], merged["identity"]
        if merged["ok"]:
            mode = "paired"
            cons, cq, src = merged["seq"], merged["quals"], merged["src"]
            fmap, rmap = merged["fmap"], merged["rmap"]
            disagreements, mism, indels = merged["disagreements"], merged["mismatches"], merged["indels"]
            aln_len = merged["alignment_length"]
            with open(os.path.join(out_dir, "alignment.fasta"), "w") as fh:
                for name, row in zip((f"{sample}_forward", f"{sample}_reverse_rc", f"{sample}_consensus"), merged["rows"]):
                    fh.write(f">{name}\n{row}\n")
        else:
            warnings.append(
                f"Forward and reverse reads do not overlap well (overlap {overlap_len} bp, "
                f"identity {identity * 100:.0f}%). Check that they belong to the same sample/amplicon. "
                "Using the better single read instead.")
            usable = [max(usable, key=lambda k: len(reads[k]["seq"]) * _mean(reads[k]["quals"]))]
            mode = "fallback_single"
    if len(usable) == 1:
        only = reads[usable[0]]
        cons = only["seq"]
        cq = list(only["quals"])
        cons = "".join(b if q >= min_base_q else "N" for b, q in zip(cons, cq))
        src = ["F" if usable[0] == "forward" else "R"] * len(cons)
        fmap = list(only["idx"]) if usable[0] == "forward" else [-1] * len(cons)
        rmap = list(only["idx"]) if usable[0] == "reverse" else [-1] * len(cons)
        mode = mode or f"{usable[0]}_only"
        aln_len = len(cons)

    # ---- statistics ----
    L = len(cons)
    counts = Counter(cons)
    valid = sum(counts[b] for b in "ACGT")
    amb = sum(counts[b] for b in "RYSWKMBDHV")
    n_b = counts["N"]
    amb_pct = amb / L * 100 if L else 0.0
    n_pct = n_b / L * 100 if L else 0.0
    q20_pct = sum(q >= 20 for q in cq) / L * 100 if L else 0.0
    gc = (counts["G"] + counts["C"]) / valid * 100 if valid else 0.0

    def info(key):
        r = reads.get(key)
        if not r:
            return {"file": "", "raw": 0, "trim": [0, 0], "trimmed": 0, "avg": 0.0}
        return {"file": r["file"], "raw": len(r["raw"]["seq"]), "trim": list(r["trim"]),
                "trimmed": r["trim"][1] - r["trim"][0], "avg": round(_mean(r["raw"]["quals"]), 1)}

    fi, ri = info("forward"), info("reverse")
    stats = {
        "sample_name": sample, "mode": mode,
        "forward_file": fi["file"], "reverse_file": ri["file"],
        "forward_raw_length": fi["raw"], "reverse_raw_length": ri["raw"],
        "forward_trimmed_length": fi["trimmed"], "reverse_trimmed_length": ri["trimmed"],
        "forward_trim": fi["trim"], "reverse_trim": ri["trim"],
        "forward_avg_quality": fi["avg"], "reverse_avg_quality": ri["avg"],
        "alignment_length": aln_len, "consensus_length": L,
        "valid_bases": valid, "ambiguous_bases": amb, "n_bases": n_b,
        "gc_content": round(gc, 2), "ambiguity_pct": round(amb_pct, 2), "n_pct": round(n_pct, 2),
        "coverage_pct": round(valid / L * 100, 2) if L else 0.0,
        "q20_pct": round(q20_pct, 1), "mean_consensus_quality": round(_mean(cq), 1),
        "overlap_length": overlap_len, "overlap_identity": round(identity * 100, 1),
        "mismatches": mism, "indels": indels, "disagreements": disagreements,
        "num_sequences_aligned": len(usable),
        "quality_rating": rate(amb_pct, n_pct, q20_pct),
        "warnings": warnings, "reference": None, "min_base_q": min_base_q,
    }

    # ---- files ----
    with open(os.path.join(out_dir, "consensus.fasta"), "w") as fh:
        fh.write(f">{sample}\n{cons}\n")
    with open(os.path.join(out_dir, "consensus.tsv"), "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["position", "base", "quality", "source"])
        for i, (b, q, s) in enumerate(zip(cons, cq, src), 1):
            w.writerow([i, b, q, s])
    if disagreements:
        with open(os.path.join(out_dir, "disagreements.tsv"), "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(disagreements[0].keys()), delimiter="\t")
            w.writeheader()
            w.writerows(disagreements)

    if reference:
        cmp = compare_to_reference(cons, cq, *reference)
        stats["reference"] = cmp
        with open(os.path.join(out_dir, "variants.tsv"), "w", newline="") as fh:
            w = csv.writer(fh, delimiter="\t")
            w.writerow(["type", "ref_pos", "ref", "alt", "consensus_pos", "quality", "change"])
            for v in cmp["variants"]:
                w.writerow([v[k] for k in ("type", "ref_pos", "ref", "alt", "consensus_pos", "quality", "change")])

    stats["consensus_quals"] = cq
    stats["consensus_src"] = "".join(src)
    stats["consensus_fmap"], stats["consensus_rmap"] = fmap, rmap
    stats["consensus_seq"] = cons
    stats["elapsed_s"] = round(time.perf_counter() - t0, 2)
    with open(os.path.join(out_dir, "sample_stats.json"), "w") as fh:
        json.dump(stats, fh, indent=1)
    return stats


# ---------------------------------------------------------------------------
def main():
    p = argparse.ArgumentParser(description="Sanger workflow for one sample")
    p.add_argument("-f", "--forward")
    p.add_argument("-r", "--reverse")
    p.add_argument("-o", "--output", default="output")
    p.add_argument("-q", "--trim-quality", type=float, default=0.05)
    p.add_argument("--sample-name", default=None)
    p.add_argument("--reference", help="Reference FASTA file")
    a = p.parse_args()
    ref = None
    if a.reference:
        ref = parse_reference(open(a.reference).read())
    name = a.sample_name or os.path.basename(os.path.abspath(a.output))
    try:
        s = run_sample(name, a.forward, a.reverse, a.output, a.trim_quality, ref)
    except ValueError as e:
        print(f"ERROR: {e}")
        sys.exit(1)
    print(f"{name}: {s['consensus_length']} bp, rating {s['quality_rating']}, "
          f"{s['mismatches']} F/R disagreements -> {a.output}/")


if __name__ == "__main__":
    main()
