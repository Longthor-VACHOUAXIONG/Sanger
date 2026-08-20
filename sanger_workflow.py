#!/usr/bin/env python3
"""
Sanger Sequencing Workflow

Pipeline:
  1. AB1 to FASTQ conversion (forward & reverse collections)
  2. Quality trimming (seqtk trimfq equivalent, q=0.05 → Phred 13)
  3. Sort collections alphabetically
  4. FASTQ Groomer (ensure fastqsanger encoding)
  5. Reverse complement (reverse reads)
  6. Sort reversed collection alphabetically
  7. Merge paired-end reads (interleaved FASTQ)
  8. FASTQ Groomer (post-merge encoding fix)
  9. FASTQ to Tabular conversion
  10. Tabular to FASTA conversion
  11. Align sequences (MAFFT)
  12. Consensus from aligned FASTA (IUPAC ambiguity codes)
  13. Regex header formatting
  14. Final MAFFT alignment
"""

import os
import sys
import re
import math
import json
import glob
import subprocess
import argparse
from collections import OrderedDict, Counter

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment


# ---------------------------------------------------------------------------
# IUPAC ambiguity dictionary (used for consensus)
# ---------------------------------------------------------------------------
IUPAC = {
    frozenset(["A"]): "A",
    frozenset(["T"]): "T", frozenset(["U"]): "T",
    frozenset(["C"]): "C",
    frozenset(["G"]): "G",
    frozenset(["A", "T"]): "W",
    frozenset(["C", "G"]): "S",
    frozenset(["A", "G"]): "R",
    frozenset(["C", "T"]): "Y",
    frozenset(["G", "T"]): "K",
    frozenset(["A", "C"]): "M",
    frozenset(["A", "C", "G"]): "V",
    frozenset(["A", "C", "T"]): "H",
    frozenset(["A", "G", "T"]): "D",
    frozenset(["C", "G", "T"]): "B",
    frozenset(["A", "C", "G", "T"]): "N",
}

IUPAC_EXPANSION = {
    "A": ["A"],
    "T": ["T"],
    "C": ["C"],
    "G": ["G"],
    "R": ["A", "G"],
    "Y": ["C", "T"],
    "W": ["A", "T"],
    "S": ["G", "C"],
    "K": ["G", "T"],
    "M": ["A", "C"],
    "B": ["C", "G", "T"],
    "D": ["A", "G", "T"],
    "H": ["A", "C", "T"],
    "V": ["A", "C", "G"],
    "N": ["A", "C", "G", "T"],
}


# ===========================================================================
# STEP 1 – AB1 to FASTQ
# ===========================================================================
def read_ab1_quality(ab1_path: str) -> list:
    """
    Read quality scores directly from AB1 file using the PHRD tag.
    This matches the ab1_fastq_converter which reads PHRD directly.
    """
    try:
        import struct
        with open(ab1_path, "rb") as f:
            data = f.read()

        # AB1 file format: look for PHRD tag (Phred quality scores)
        # The PHRD tag contains quality scores as signed 32-bit integers
        phrd_pos = data.find(b"PHRD")
        if phrd_pos == -1:
            # Try alternative tag names
            for tag in [b"P1Acq", b"P9Acq", b"QUAL"]:
                phrd_pos = data.find(tag)
                if phrd_pos != -1:
                    break

        if phrd_pos == -1:
            return None

        # Skip tag name (4 bytes) and find data
        # PHRD tag structure: tag_name(4) + data_size(4) + data...
        # Data starts after the tag and size field
        data_start = phrd_pos + 4  # Skip tag name

        # Read data size (32-bit little-endian integer)
        if data_start + 4 > len(data):
            return None
        data_size = struct.unpack("<I", data[data_start:data_start + 4])[0]

        # Quality scores are stored as 32-bit signed integers
        qual_start = data_start + 4
        num_quals = data_size // 4  # Each quality is 4 bytes

        if qual_start + num_quals * 4 > len(data):
            return None

        quals = []
        for i in range(num_quals):
            offset = qual_start + i * 4
            qual = struct.unpack("<i", data[offset:offset + 4])[0]
            # Clamp to valid Phred range (0-93)
            qual = max(0, min(93, qual))
            quals.append(qual)

        return quals if quals else None

    except Exception:
        return None


def ab1_to_fastq(ab1_path: str) -> SeqRecord:
    """
    Convert a single .ab1 chromatogram file to a SeqRecord with quality scores.
    Matches the 'ab1 to FASTQ converter'.
    Reads PHRD quality scores directly from the AB1 file.
    """
    record = SeqIO.read(ab1_path, "abi")

    # Try to read quality scores directly from PHRD tag
    phrd_quals = read_ab1_quality(ab1_path)

    if phrd_quals and len(phrd_quals) == len(record.seq):
        # Use PHRD quality scores directly
        record.letter_annotations["phred_quality"] = phrd_quals
    elif not record.letter_annotations.get("phred_quality"):
        # Fallback: use BioPython's extracted quality or default
        # BioPython's quality extraction might differ
        record.letter_annotations["phred_quality"] = [20] * len(record.seq)

    return record


def records_to_fastq_string(records: list, tag: str = "") -> str:
    """Write a list of SeqRecords to FASTQ string format."""
    lines = []
    for rec in records:
        quals = rec.letter_annotations.get("phred_quality", [20] * len(rec.seq))
        qual_str = "".join(chr(q + 33) for q in quals)
        lines.append(f"@{rec.id}")
        lines.append(str(rec.seq))
        lines.append("+")
        lines.append(qual_str)
    return "\n".join(lines) + "\n"


def parse_fastq_string(fastq_str: str) -> list:
    """Parse a FASTQ string into a list of SeqRecord objects."""
    records = []
    lines = fastq_str.strip().split("\n")
    i = 0
    while i < len(lines) - 3:
        if lines[i].startswith("@"):
            header = lines[i][1:].strip()
            seq = lines[i + 1].strip()
            # lines[i+2] should be "+"
            qual_str = lines[i + 3].strip()
            quals = [ord(c) - 33 for c in qual_str]
            rec = SeqRecord(Seq(seq), id=header, description="")
            rec.letter_annotations["phred_quality"] = quals
            records.append(rec)
            i += 4
        else:
            i += 1
    return records


# ===========================================================================
# STEP 2 – Quality Trimming (seqtk trimfq)
# ===========================================================================
def trimfq_quality(records: list, q: float = 0.05) -> list:
    """
    Quality-based trimming using seqtk trimfq.
    Falls back to a Python implementation if seqtk is not available.
    """
    import shutil
    import tempfile
    import io as _io

    # Try to use real seqtk if available
    if shutil.which("seqtk"):
        # Write records to temp FASTQ
        with tempfile.NamedTemporaryFile(mode="w", suffix=".fastq", delete=False) as f:
            SeqIO.write(records, f, "fastq")
            input_file = f.name

        output_file = input_file + ".trimmed"

        try:
            result = subprocess.run(
                ["seqtk", "trimfq", "-q", str(q), input_file],
                capture_output=True, text=True, timeout=60
            )

            if result.returncode != 0:
                print(f"  WARNING: seqtk failed: {result.stderr.strip()}")
                raise RuntimeError(result.stderr)

            # Parse seqtk output
            trimmed = list(SeqIO.parse(_io.StringIO(result.stdout), "fastq"))
            if not trimmed:
                raise RuntimeError("seqtk produced no output")

            return trimmed

        except (FileNotFoundError, subprocess.TimeoutExpired) as e:
            print(f"  WARNING: seqtk not available ({e}), falling back to Python implementation")
        finally:
            if os.path.exists(input_file):
                os.unlink(input_file)
            if os.path.exists(output_file):
                os.unlink(output_file)

    # Fallback: Python implementation
    trimmed = []
    for rec in records:
        quals = list(rec.letter_annotations.get("phred_quality", [20] * len(rec.seq)))
        seq = str(rec.seq)
        n = len(quals)

        if n == 0:
            trimmed.append(rec)
            continue

        qual_threshold = -10 * math.log10(q) if q > 0 else 20
        window_size = max(1, n // 5)

        left_trim = 0
        for i in range(n - window_size + 1):
            window_mean = sum(quals[i:i + window_size]) / window_size
            if window_mean >= qual_threshold:
                left_trim = i
                break

        right_trim = n
        for i in range(n - 1, window_size - 1, -1):
            window_mean = sum(quals[i - window_size + 1:i + 1]) / window_size
            if window_mean >= qual_threshold:
                right_trim = i + 1
                break

        for i in range(left_trim, right_trim):
            if quals[i] >= qual_threshold:
                left_trim = i
                break

        for i in range(right_trim - 1, left_trim - 1, -1):
            if quals[i] >= qual_threshold:
                right_trim = i + 1
                break

        min_len = min(50, n)
        if right_trim - left_trim < min_len and n >= min_len:
            left_trim = 0
            right_trim = min_len

        new_seq = Seq(seq[left_trim:right_trim])
        new_quals = quals[left_trim:right_trim]
        new_rec = SeqRecord(new_seq, id=rec.id, description="")
        new_rec.letter_annotations["phred_quality"] = new_quals
        trimmed.append(new_rec)

    return trimmed


# ===========================================================================
# STEP 3 – Sort Collection (alphabetical)
# ===========================================================================
def sort_records_alpha(records: list) -> list:
    """Sort SeqRecords alphabetically by ID."""
    return sorted(records, key=lambda r: r.id)


# ===========================================================================
# STEP 4 – FASTQ Groomer (ensure fastqsanger encoding)
# ===========================================================================
def fastq_groomer(records: list) -> list:
    """
    Ensures quality values are in Sanger (Phred+33) encoding
    and fixes IDs to match the sequence.
    """
    groomed = []
    for rec in records:
        quals = list(rec.letter_annotations.get("phred_quality", [20] * len(rec.seq)))
        # Clamp to 0-93 range for Sanger encoding
        quals = [max(0, min(93, q)) for q in quals]
        new_rec = SeqRecord(Seq(str(rec.seq)), id=rec.id, description="")
        new_rec.letter_annotations["phred_quality"] = quals
        groomed.append(new_rec)
    return groomed


# ===========================================================================
# STEP 5 – Reverse Complement
# ===========================================================================
def reverse_complement_records(records: list) -> list:
    """Reverse complement each SeqRecord – equivalent to FASTX Reverse-Complement."""
    rc_records = []
    for rec in records:
        rc_seq = rec.seq.reverse_complement()
        quals = list(rec.letter_annotations.get("phred_quality", []))
        quals.reverse()  # Reverse quality scores too
        new_rec = SeqRecord(rc_seq, id=rec.id, description="")
        if quals:
            new_rec.letter_annotations["phred_quality"] = quals
        rc_records.append(new_rec)
    return rc_records


# ===========================================================================
# STEP 6 – Merge Paired-End (interleaved FASTQ)
# ===========================================================================
def mergepe_interleaved(forward_records: list, reverse_records: list) -> list:
    """
    Merge forward and reverse into interleaved FASTQ.
    Equivalent to seqtk mergepe.
    Forward and reverse collections must be sorted the same way.
    """
    merged = []
    for fwd, rev in zip(forward_records, reverse_records):
        merged.append(fwd)
        merged.append(rev)
    return merged


# ===========================================================================
# STEP 7 – FASTQ to Tabular
# ===========================================================================
def fastq_to_tabular(records: list) -> str:
    """
    Convert FASTQ to tabular format.
    Columns: sequence identifier, sequence
    Equivalent to fastq_to_tabular.
    """
    lines = []
    for rec in records:
        # Description columns = 1 (from workflow: descr_columns="1")
        lines.append(f"{rec.id}\t{str(rec.seq)}")
    return "\n".join(lines) + "\n"


# ===========================================================================
# STEP 8 – Tabular to FASTA
# ===========================================================================
def tabular_to_fasta(tabular_str: str) -> list:
    """
    Convert tabular to FASTA.
    seq_col=2 (sequence), title_col=1 (identifier)
    Equivalent to tabular_to_fasta.
    """
    records = []
    for line in tabular_str.strip().split("\n"):
        if not line.strip():
            continue
        parts = line.split("\t")
        if len(parts) >= 2:
            header = parts[0].strip()
            seq = parts[1].strip()
            rec = SeqRecord(Seq(seq), id=header, description="")
            records.append(rec)
    return records


# ===========================================================================
# STEP 9 – MAFFT Alignment
# ===========================================================================
def mafft_align(records: list, method: str = "auto") -> tuple:
    """
    Run MAFFT alignment on a list of SeqRecords.
    Returns (aligned_records, log_string).
    """
    # Write sequences to temp file
    import tempfile
    with tempfile.NamedTemporaryFile(mode="w", suffix=".fasta", delete=False) as f:
        SeqIO.write(records, f, "fasta")
        input_file = f.name

    with tempfile.NamedTemporaryFile(mode="w", suffix=".fasta", delete=False) as f:
        output_file = f.name

    try:
        # Run MAFFT
        result = subprocess.run(
            ["mafft", "--auto", input_file],
            capture_output=True,
            text=True,
            timeout=300
        )
        aligned_str = result.stdout
        log_str = result.stderr

        # Parse aligned sequences
        import io
        aligned_records = list(SeqIO.parse(io.StringIO(aligned_str), "fasta"))

        return aligned_records, log_str
    except subprocess.TimeoutExpired:
        raise RuntimeError("MAFFT alignment timed out after 300 seconds")
    except FileNotFoundError:
        raise RuntimeError("MAFFT not found. Please install MAFFT: conda install -c bioconda mafft")
    finally:
        os.unlink(input_file)
        if os.path.exists(output_file):
            os.unlink(output_file)


# ===========================================================================
# STEP 10 – Consensus from Aligned FASTA (chr_ambiguity method)
# ===========================================================================
def aligned_to_consensus(aligned_records: list, aligned_with_quals: list = None, sample_name: str = "consensus") -> SeqRecord:
    """
    Build consensus from aligned FASTA.
    method=chr_ambiguity, gaps=false, seqtype=DNA

    Algorithm:
    1. For each alignment column, collect non-gap bases + their quality scores
    2. Expand IUPAC ambiguity codes to constituent bases (e.g. Y→C/T)
    3. Weight each base call by its Phred quality score
    4. The base with the highest total quality weight wins
    5. Skip columns where all sequences have gaps
    """
    if not aligned_records:
        raise ValueError("No sequences to build consensus from")

    max_len = max(len(r.seq) for r in aligned_records)
    n_seqs = len(aligned_records)
    use_quality = aligned_with_quals is not None and len(aligned_with_quals) == n_seqs

    raw_consensus = []
    raw_coverage = []

    for i in range(max_len):
        base_weights = {}
        nongap_count = 0

        for idx, rec in enumerate(aligned_records):
            if i >= len(rec.seq):
                continue
            base = str(rec.seq[i]).upper()
            if base in ("-", "."):
                continue

            nongap_count += 1

            if use_quality and idx < len(aligned_with_quals):
                quals = aligned_with_quals[idx]
                weight = max(quals[i], 1) if i < len(quals) else 10
            else:
                weight = 10

            # Expand IUPAC ambiguity codes before voting
            if base in IUPAC_EXPANSION:
                for eb in IUPAC_EXPANSION[base]:
                    base_weights[eb] = base_weights.get(eb, 0) + weight
            else:
                base_weights[base] = base_weights.get(base, 0) + weight

        raw_coverage.append(nongap_count)

        if not base_weights:
            raw_consensus.append("N")
            continue

        best_base = max(base_weights, key=base_weights.get)
        raw_consensus.append(best_base)

    consensus_seq = Seq("".join(raw_consensus))
    return SeqRecord(consensus_seq, id=sample_name, description="")


# ===========================================================================
# STEP 11 – Regex Find And Replace (format FASTA headers)
# ===========================================================================
def regex_find_replace_fasta(fasta_str: str) -> str:
    """
    Regex find and replace on FASTA file.
    Pattern: ([A-Z-])>  →  \1\n>
    This ensures there's a newline before the '>' header when the previous
    line ends with a base character, effectively formatting multi-line FASTA
    into one-sequence-per-line format.
    """
    result = re.sub(r"([A-Z\-])>", r"\1\n>", fasta_str)
    return result


# ===========================================================================
# MAIN WORKFLOW
# ===========================================================================
def main():
    parser = argparse.ArgumentParser(
        description="Sanger Sequencing Workflow"
    )
    parser.add_argument(
        "--forward", "-f",
        default="F",
        help="Directory containing forward AB1 files (default: F/)"
    )
    parser.add_argument(
        "--reverse", "-r",
        default="R",
        help="Directory containing reverse AB1 files (default: R/)"
    )
    parser.add_argument(
        "--output", "-o",
        default="output",
        help="Output directory (default: output/)"
    )
    parser.add_argument(
        "--trim-quality", "-q",
        type=float,
        default=0.05,
        help="Quality trimming threshold for seqtk trimfq (default: 0.05, ~Phred 13)"
    )
    args = parser.parse_args()

    # Create output directory
    os.makedirs(args.output, exist_ok=True)
    work = args.output  # shorthand

    print("=" * 70)
    print("SANGER SEQUENCING WORKFLOW")
    print("Python Sanger Sequencing Data Analysis")
    print("=" * 70)

    # ---------------------------------------------------------------
    # STEP 1 & 2: AB1 → FASTQ conversion (Forward & Reverse)
    # ---------------------------------------------------------------
    print("\n[Step 1] Reading AB1 files and converting to FASTQ...")

    forward_ab1_files = sorted(glob.glob(os.path.join(args.forward, "*.ab1")))
    reverse_ab1_files = sorted(glob.glob(os.path.join(args.reverse, "*.ab1")))

    if not forward_ab1_files:
        print(f"ERROR: No .ab1 files found in {args.forward}/")
        sys.exit(1)
    if not reverse_ab1_files:
        print(f"ERROR: No .ab1 files found in {args.reverse}/")
        sys.exit(1)

    print(f"  Forward files: {len(forward_ab1_files)}")
    print(f"  Reverse files: {len(reverse_ab1_files)}")

    forward_records_raw = [ab1_to_fastq(f) for f in forward_ab1_files]
    reverse_records_raw = [ab1_to_fastq(f) for f in reverse_ab1_files]

    # Write raw FASTQ
    with open(os.path.join(work, "forward_raw.fastq"), "w") as f:
        f.write(records_to_fastq_string(forward_records_raw))
    with open(os.path.join(work, "reverse_raw.fastq"), "w") as f:
        f.write(records_to_fastq_string(reverse_records_raw))

    print(f"  → forward_raw.fastq ({len(forward_records_raw)} sequences)")
    print(f"  → reverse_raw.fastq ({len(reverse_records_raw)} sequences)")

    # ---------------------------------------------------------------
    # STEP 3: Quality Trimming (seqtk trimfq -q 0.05)
    # ---------------------------------------------------------------
    print(f"\n[Step 2] Quality trimming (q={args.trim_quality}, Phred≥{int(-10*math.log10(args.trim_quality))})...")

    forward_trimmed = trimfq_quality(forward_records_raw, q=args.trim_quality)
    reverse_trimmed = trimfq_quality(reverse_records_raw, q=args.trim_quality)

    print(f"  Forward: {len(forward_records_raw[0].seq)} → {len(forward_trimmed[0].seq)} bases (avg)")
    print(f"  Reverse: {len(reverse_records_raw[0].seq)} → {len(reverse_trimmed[0].seq)} bases (avg)")

    with open(os.path.join(work, "forward_trimmed.fastq"), "w") as f:
        f.write(records_to_fastq_string(forward_trimmed))
    with open(os.path.join(work, "reverse_trimmed.fastq"), "w") as f:
        f.write(records_to_fastq_string(reverse_trimmed))

    # ---------------------------------------------------------------
    # STEP 4: Sort collections alphabetically
    # ---------------------------------------------------------------
    print("\n[Step 3] Sorting collections alphabetically...")

    forward_sorted = sort_records_alpha(forward_trimmed)
    reverse_sorted = sort_records_alpha(reverse_trimmed)

    print(f"  Forward order: {[r.id for r in forward_sorted]}")
    print(f"  Reverse order: {[r.id for r in reverse_sorted]}")

    # ---------------------------------------------------------------
    # STEP 5: FASTQ Groomer (reverse reads → fastqsanger)
    # ---------------------------------------------------------------
    print("\n[Step 4] FASTQ Groomer (ensuring Sanger encoding)...")

    reverse_groomed = fastq_groomer(reverse_sorted)

    with open(os.path.join(work, "reverse_groomed.fastq"), "w") as f:
        f.write(records_to_fastq_string(reverse_groomed))

    # ---------------------------------------------------------------
    # STEP 6: Reverse Complement (reverse reads)
    # ---------------------------------------------------------------
    print("\n[Step 5] Reverse Complementing reverse reads...")

    reverse_rc = reverse_complement_records(reverse_groomed)

    with open(os.path.join(work, "reverse_complement.fastq"), "w") as f:
        f.write(records_to_fastq_string(reverse_rc))

    print(f"  RC first seq: {str(reverse_rc[0].seq[:50])}...")

    # ---------------------------------------------------------------
    # STEP 7: Sort reversed collection
    # ---------------------------------------------------------------
    print("\n[Step 6] Sorting reversed collection alphabetically...")

    reverse_rc_sorted = sort_records_alpha(reverse_rc)

    # ---------------------------------------------------------------
    # STEP 8: Merge Paired-End
    # ---------------------------------------------------------------
    print("\n[Step 7] Merging paired-end reads (interleaved)...")

    merged = mergepe_interleaved(forward_sorted, reverse_rc_sorted)

    with open(os.path.join(work, "merged_interleaved.fastq"), "w") as f:
        f.write(records_to_fastq_string(merged))

    print(f"  Merged {len(merged)} records ({len(merged)//2} pairs)")

    # ---------------------------------------------------------------
    # STEP 9: FASTQ Groomer (post-merge)
    # ---------------------------------------------------------------
    print("\n[Step 8] FASTQ Groomer (post-merge encoding fix)...")

    merged_groomed = fastq_groomer(merged)

    with open(os.path.join(work, "merged_groomed.fastq"), "w") as f:
        f.write(records_to_fastq_string(merged_groomed))

    # ---------------------------------------------------------------
    # STEP 10: FASTQ → Tabular
    # ---------------------------------------------------------------
    print("\n[Step 9] FASTQ to Tabular conversion...")

    tabular = fastq_to_tabular(merged_groomed)

    with open(os.path.join(work, "merged.tabular"), "w") as f:
        f.write(tabular)

    tab_lines = tabular.strip().split("\n")
    print(f"  {len(tab_lines)} tabular records created")

    # ---------------------------------------------------------------
    # STEP 11: Tabular → FASTA
    # ---------------------------------------------------------------
    print("\n[Step 10] Tabular to FASTA conversion...")

    fasta_records = tabular_to_fasta(tabular)

    with open(os.path.join(work, "merged.fasta"), "w") as f:
        SeqIO.write(fasta_records, f, "fasta")

    print(f"  {len(fasta_records)} FASTA sequences created")

    # ---------------------------------------------------------------
    # STEP 12: MAFFT Alignment (first pass)
    # ---------------------------------------------------------------
    print("\n[Step 11] MAFFT alignment (first pass)...")

    aligned_records_1, log1 = mafft_align(fasta_records)

    with open(os.path.join(work, "aligned.fasta"), "w") as f:
        SeqIO.write(aligned_records_1, f, "fasta")
    with open(os.path.join(work, "alignment_log_1.txt"), "w") as f:
        f.write(log1)

    print(f"  Aligned {len(aligned_records_1)} sequences")
    print(f"  Alignment length: {len(aligned_records_1[0].seq)} bp")

    # ---------------------------------------------------------------
    # STEP 13: Consensus from Aligned FASTA
    # ---------------------------------------------------------------
    print("\n[Step 12] Building consensus sequence (quality-weighted chr_ambiguity)...")

    # Pass quality scores from the original records for weighting
    # The aligned records should have quality annotations from the input
    aligned_quals = []
    for rec in aligned_records_1:
        quals = rec.letter_annotations.get("phred_quality", [])
        if quals:
            aligned_quals.append(quals)
        else:
            # Use default quality if not available
            aligned_quals.append([20] * len(rec.seq))

    sample_name = os.path.basename(os.path.abspath(work))
    consensus = aligned_to_consensus(aligned_records_1, aligned_with_quals=aligned_quals, sample_name=sample_name)

    with open(os.path.join(work, "consensus.fasta"), "w") as f:
        SeqIO.write([consensus], f, "fasta")

    print(f"  Consensus length: {len(consensus.seq)} bp")
    print(f"  First 80 bp: {str(consensus.seq[:80])}")

    # ---------------------------------------------------------------
    # STEP 14: Merge.files (merge consensus with aligned sequences)
    # ---------------------------------------------------------------
    print("\n[Step 13] Merge.files (combining consensus + aligned sequences)...")

    # Merge consensus with alignment
    all_sequences = aligned_records_1 + [consensus]

    with open(os.path.join(work, "all_sequences.fasta"), "w") as f:
        SeqIO.write(all_sequences, f, "fasta")

    print(f"  Total sequences: {len(all_sequences)} ({len(aligned_records_1)} aligned + 1 consensus)")

    # ---------------------------------------------------------------
    # STEP 15: Regex Find And Replace (format headers)
    # ---------------------------------------------------------------
    print("\n[Step 14] Regex Find And Replace (formatting FASTA headers)...")

    # Read the FASTA file
    fasta_str = ""
    for rec in all_sequences:
        fasta_str += f">{rec.id}\n{str(rec.seq)}\n"

    formatted_fasta = regex_find_replace_fasta(fasta_str)

    with open(os.path.join(work, "formatted_sequences.fasta"), "w") as f:
        f.write(formatted_fasta)

    # ---------------------------------------------------------------
    # STEP 16: Final MAFFT Alignment
    # ---------------------------------------------------------------
    print("\n[Step 15] MAFFT alignment (final pass on formatted sequences)...")

    # Parse the formatted FASTA back
    import io
    formatted_records = list(SeqIO.parse(io.StringIO(formatted_fasta), "fasta"))

    aligned_records_final, log2 = mafft_align(formatted_records)

    with open(os.path.join(work, "final_aligned.fasta"), "w") as f:
        SeqIO.write(aligned_records_final, f, "fasta")
    with open(os.path.join(work, "alignment_log_2.txt"), "w") as f:
        f.write(log2)

    print(f"  Final alignment: {len(aligned_records_final)} sequences, {len(aligned_records_final[0].seq)} bp")

    # ---------------------------------------------------------------
    # QUALITY STATISTICS
    # ---------------------------------------------------------------
    print("\n[Step 16] Computing quality statistics...")

    consensus_seq = str(consensus.seq).upper()
    consensus_len = len(consensus_seq)

    base_counts = Counter(consensus_seq)
    valid_bases = sum(base_counts[b] for b in "ATCG")
    ambiguous_bases = sum(base_counts[b] for b in "RYSWKMBDHVN")
    n_bases = base_counts.get("N", 0)
    gap_bases = base_counts.get("-", 0) + base_counts.get(".", 0)

    gc_count = base_counts.get("G", 0) + base_counts.get("C", 0)
    gc_content = (gc_count / valid_bases * 100) if valid_bases > 0 else 0.0
    ambiguity_pct = (ambiguous_bases / consensus_len * 100) if consensus_len > 0 else 0.0
    n_pct = (n_bases / consensus_len * 100) if consensus_len > 0 else 0.0
    gap_pct = (gap_bases / consensus_len * 100) if consensus_len > 0 else 0.0
    coverage_pct = (valid_bases / consensus_len * 100) if consensus_len > 0 else 0.0

    # Per-read quality stats from raw reads
    all_fwd_quals = []
    for rec in forward_records_raw:
        all_fwd_quals.extend(rec.letter_annotations.get("phred_quality", []))
    all_rev_quals = []
    for rec in reverse_records_raw:
        all_rev_quals.extend(rec.letter_annotations.get("phred_quality", []))

    fwd_avg_qual = sum(all_fwd_quals) / len(all_fwd_quals) if all_fwd_quals else 0
    rev_avg_qual = sum(all_rev_quals) / len(all_rev_quals) if all_rev_quals else 0
    fwd_min_qual = min(all_fwd_quals) if all_fwd_quals else 0
    rev_min_qual = min(all_rev_quals) if all_rev_quals else 0

    # Quality rating
    if ambiguity_pct < 1 and n_pct < 1:
        quality_rating = "Excellent"
    elif ambiguity_pct < 3 and n_pct < 3:
        quality_rating = "Good"
    elif ambiguity_pct < 5 and n_pct < 5:
        quality_rating = "Fair"
    else:
        quality_rating = "Poor"

    stats = {
        "sample_name": os.path.basename(work),
        "forward_file": os.path.basename(forward_ab1_files[0]) if forward_ab1_files else "",
        "reverse_file": os.path.basename(reverse_ab1_files[0]) if reverse_ab1_files else "",
        "forward_raw_length": len(forward_records_raw[0].seq) if forward_records_raw else 0,
        "reverse_raw_length": len(reverse_records_raw[0].seq) if reverse_records_raw else 0,
        "forward_trimmed_length": len(forward_trimmed[0].seq) if forward_trimmed else 0,
        "reverse_trimmed_length": len(reverse_trimmed[0].seq) if reverse_trimmed else 0,
        "alignment_length": len(aligned_records_1[0].seq) if aligned_records_1 else 0,
        "consensus_length": consensus_len,
        "valid_bases": valid_bases,
        "ambiguous_bases": ambiguous_bases,
        "n_bases": n_bases,
        "gap_bases": gap_bases,
        "gc_content": round(gc_content, 2),
        "ambiguity_pct": round(ambiguity_pct, 2),
        "n_pct": round(n_pct, 2),
        "gap_pct": round(gap_pct, 2),
        "coverage_pct": round(coverage_pct, 2),
        "forward_avg_quality": round(fwd_avg_qual, 1),
        "reverse_avg_quality": round(rev_avg_qual, 1),
        "forward_min_quality": fwd_min_qual,
        "reverse_min_quality": rev_min_qual,
        "quality_rating": quality_rating,
        "num_sequences_aligned": len(aligned_records_1),
    }

    stats_path = os.path.join(work, "sample_stats.json")
    with open(stats_path, "w") as f:
        json.dump(stats, f, indent=2)

    print(f"  Quality rating: {quality_rating}")
    print(f"  Consensus: {consensus_len} bp, GC: {gc_content:.1f}%, Ambiguity: {ambiguity_pct:.1f}%")
    print(f"  → sample_stats.json")

    # ---------------------------------------------------------------
    # SUMMARY
    # ---------------------------------------------------------------
    print("\n" + "=" * 70)
    print("WORKFLOW COMPLETE!")
    print("=" * 70)
    print(f"\nOutput directory: {os.path.abspath(work)}/")
    print("\nFiles produced:")
    for fname in sorted(os.listdir(work)):
        fpath = os.path.join(work, fname)
        size = os.path.getsize(fpath)
        print(f"  {fname:40s} ({size:,} bytes)")

    print("\n" + "=" * 70)
    print("WORKFLOW STEP SUMMARY:")
    print("=" * 70)
    print(f"  {'Step':12s} | {'Tool':30s} | {'Output'}")
    print(f"  {'-'*12}-+-{'-'*30}-+-{'-'*30}")
    steps = [
        ("Steps 1-2", "AB1 to FASTQ", "forward_raw.fastq, reverse_raw.fastq"),
        ("Step 3", "Quality trimming", "forward_trimmed.fastq, reverse_trimmed.fastq"),
        ("Step 4", "Sort collection", "sorted in memory"),
        ("Step 5", "FASTQ Groomer", "reverse_groomed.fastq"),
        ("Step 6", "Reverse-Complement", "reverse_complement.fastq"),
        ("Step 7", "Sort collection", "sorted in memory"),
        ("Step 8", "Merge paired-end", "merged_interleaved.fastq"),
        ("Step 9", "FASTQ Groomer", "merged_groomed.fastq"),
        ("Step 10", "FASTQ to Tabular", "merged.tabular"),
        ("Step 11", "Tabular-to-FASTA", "merged.fasta"),
        ("Step 12", "MAFFT alignment", "aligned.fasta"),
        ("Step 13", "Consensus", "consensus.fasta"),
        ("Steps 14-15", "Merge files", "all_sequences.fasta"),
        ("Step 16", "Regex format", "formatted_sequences.fasta"),
        ("Step 17", "MAFFT alignment", "final_aligned.fasta"),
    ]
    for step, tool, output in steps:
        print(f"  {step:12s} | {tool:30s} | {output}")

    print("\nDone!")


if __name__ == "__main__":
    main()
