# 🧬 Sanger AB1 → Consensus

**Free, no install.** Upload forward/reverse `.ab1` files and get a quality-weighted
consensus, an interactive chromatogram, forward/reverse disagreement review, optional
mutation report against a reference, and an HTML report + ZIP of all results.
Available in English and ລາວ.

## Features
- Quality trimming (Mott algorithm) and overlap alignment of forward + reverse reads
- Consensus uses Phred quality; unresolved disagreements become IUPAC codes (R, Y, …) – nothing is hidden
- Interactive chromatogram with trimmed regions, quality bars, and "jump to disagreement"
- Optional reference sequence → SNPs / insertions / deletions with quality
- Batch mode: auto-pairs `_F/_R`, `_forward/_reverse`, `_27F/_1492R`, `M13F/M13R` …
- Downloads: ZIP, HTML report, FASTA, CSV

## Run locally
```
pip install -r requirements.txt
streamlit run app.py
```
CLI: `python sanger_wrapper.py -i ab1_folder/ -o results --reference ref.fasta`

Add a real pair as `examples/example_forward.ab1` + `examples/example_reverse.ab1` to enable the demo button.
