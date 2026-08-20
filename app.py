#!/usr/bin/env python3
"""
Sanger Sequencing Data analysis workflow
Streamlit interface for uploading AB1 files, running the pipeline, and viewing results.
"""

import os
import sys
import json
import shutil
import tempfile
import subprocess
import base64
from pathlib import Path
from datetime import datetime

import streamlit as st

# ---------------------------------------------------------------------------
# Config
# ---------------------------------------------------------------------------
WORKFLOW_SCRIPT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "sanger_workflow.py")
WRAPPER_SCRIPT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "sanger_wrapper.py")
TEMP_UPLOAD_DIR = os.path.join(tempfile.gettempdir(), "sanger_uploads")
OUTPUT_BASE = os.path.join(tempfile.gettempdir(), "sanger_output")

# Ensure temp dirs exist
os.makedirs(TEMP_UPLOAD_DIR, exist_ok=True)
os.makedirs(OUTPUT_BASE, exist_ok=True)


# ---------------------------------------------------------------------------
# Page config
# ---------------------------------------------------------------------------
st.set_page_config(
    page_title="Sanger Workflow",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded",
)

# Custom CSS
st.markdown("""
<style>
    .stMetric .metric-value { font-size: 1.5rem; }
    .block-container { padding-top: 1rem; }
    .quality-excellent { color: #27ae60; font-weight: bold; }
    .quality-good { color: #2ecc71; font-weight: bold; }
    .quality-fair { color: #f39c12; font-weight: bold; }
    .quality-poor { color: #e74c3c; font-weight: bold; }
    div[data-testid="stExpander"] { border: 1px solid #e0e0e0; border-radius: 8px; }
</style>
""", unsafe_allow_html=True)


# ---------------------------------------------------------------------------
# Helper functions
# ---------------------------------------------------------------------------

def save_uploaded_files(files, subfolder):
    """Save uploaded files to a temporary directory."""
    upload_dir = os.path.join(TEMP_UPLOAD_DIR, subfolder)
    os.makedirs(upload_dir, exist_ok=True)
    saved = []
    for f in files:
        path = os.path.join(upload_dir, f.name)
        with open(path, "wb") as out:
            out.write(f.getbuffer())
        saved.append(path)
    return upload_dir, saved


def stream_process(cmd, status_container, progress_bar=None, total_steps=16, timeout=3600):
    """Run a subprocess and stream output line by line into a Streamlit status container."""
    import re
    import time
    proc = subprocess.Popen(
        cmd,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        encoding="utf-8",
        bufsize=1,
    )
    log_lines = []
    last_update = 0
    current_step = 0
    start_time = time.time()

    for line in iter(proc.stdout.readline, ""):
        line = line.rstrip()
        if not line:
            continue
        log_lines.append(line)
        # Write line to status container
        status_container.write(f"`{line}`")

        # Detect step number from [Step X] pattern
        m = re.search(r"\[Step (\d+)\]", line)
        if m:
            step_num = int(m.group(1))
            if step_num != current_step:
                current_step = step_num
                pct = min(current_step / total_steps, 1.0)
                if progress_bar:
                    progress_bar.progress(pct, text=f"Step {current_step}/{total_steps}")

        # Update label periodically
        now = time.time()
        if now - last_update >= 1.0:
            elapsed = int(now - start_time)
            status_container.update(label=f"Running... ({elapsed}s elapsed, step {current_step}/{total_steps})", state="running")
            last_update = now

    # Final progress update
    if progress_bar:
        progress_bar.progress(1.0, text="Complete!")

    proc.wait(timeout=timeout)
    elapsed = int(time.time() - start_time)
    return proc.returncode, "\n".join(log_lines), elapsed


def run_workflow(forward_dir, reverse_dir, output_dir, trim_quality=0.05, sample_name=None, status_container=None, progress_bar=None):
    """Run sanger_workflow.py as a subprocess with streaming output."""
    cmd = [
        sys.executable, WORKFLOW_SCRIPT,
        "--forward", forward_dir,
        "--reverse", reverse_dir,
        "--output", output_dir,
        "--trim-quality", str(trim_quality),
    ]
    if sample_name:
        cmd.extend(["--sample-name", sample_name])
    if status_container:
        return stream_process(cmd, status_container, progress_bar=progress_bar, total_steps=16, timeout=600)
    result = subprocess.run(cmd, capture_output=True, text=True, timeout=600)
    return result.returncode, result.stdout + result.stderr, 0


def run_wrapper_batch(upload_dir, output_dir, trim_quality=0.05, status_container=None, progress_bar=None):
    """Run sanger_wrapper.py in auto-dir mode with streaming output."""
    cmd = [
        sys.executable, WRAPPER_SCRIPT,
        "--auto-dir", upload_dir,
        "-o", output_dir,
        "-q", str(trim_quality),
    ]
    if status_container:
        return stream_process(cmd, status_container, progress_bar=progress_bar, total_steps=16, timeout=3600)
    result = subprocess.run(cmd, capture_output=True, text=True, timeout=3600)
    return result.returncode, result.stdout + result.stderr, 0


def load_sample_stats(output_dir):
    """Load sample_stats.json from output directory."""
    stats_path = os.path.join(output_dir, "sample_stats.json")
    if os.path.exists(stats_path):
        with open(stats_path) as f:
            return json.load(f)
    return None


def load_all_sample_stats(output_dir):
    """Load all_sample_stats.json from output directory."""
    stats_path = os.path.join(output_dir, "all_sample_stats.json")
    if os.path.exists(stats_path):
        with open(stats_path) as f:
            return json.load(f)
    return None


def read_fasta(filepath):
    """Read a FASTA file and return header + sequence."""
    if not os.path.exists(filepath):
        return None, None
    with open(filepath) as f:
        lines = f.readlines()
    header = ""
    seq = ""
    for line in lines:
        if line.startswith(">"):
            header = line.strip()[1:]
        else:
            seq += line.strip()
    return header, seq


def get_download_link(filepath, filename=None):
    """Generate a download link for a file."""
    if not os.path.exists(filepath):
        return None
    with open(filepath, "rb") as f:
        data = f.read()
    b64 = base64.b64encode(data).decode()
    if filename is None:
        filename = os.path.basename(filepath)
    return f'<a href="data:file/octet-stream;base64,{b64}" download="{filename}">⬇️ Download {filename}</a>'


def cleanup_session(session_id):
    """Clean up temp files for a session."""
    upload_dir = os.path.join(TEMP_UPLOAD_DIR, session_id)
    output_dir = os.path.join(OUTPUT_BASE, session_id)
    for d in [upload_dir, output_dir]:
        if os.path.exists(d):
            shutil.rmtree(d, ignore_errors=True)


def init_session_state():
    """Initialize session state variables."""
    defaults = {
        "run_returncode": None,
        "run_output": "",
        "run_elapsed": 0,
        "run_output_dir": "",
        "run_session_id": "",
        "run_mode": "Single pair",
    }
    for key, value in defaults.items():
        if key not in st.session_state:
            st.session_state[key] = value


# ---------------------------------------------------------------------------
# Sidebar
# ---------------------------------------------------------------------------

with st.sidebar:
    st.title("🧬 Sanger Sequencing")
    st.markdown("---")
    st.markdown("**Pipeline:** Sanger Sequencing Data Analysis")
    st.markdown("**Steps:**")
    st.markdown("""
    1. AB1 → FASTQ conversion
    2. Quality trimming (seqtk)
    3. Reverse complement
    4. Paired-end merge
    5. MAFFT alignment
    6. Consensus generation
    7. Quality report
    """)
    st.markdown("---")
    st.markdown("### Settings")
    trim_quality = st.slider(
        "Trim quality threshold",
        min_value=0.01,
        max_value=0.20,
        value=0.05,
        step=0.01,
        help="Lower = more aggressive trimming. Default 0.05 ≈ Phred 13."
    )
    st.markdown("---")
    st.markdown("### Links")
    st.markdown("[📖 Usage Guide](./USAGE.md)")


# ---------------------------------------------------------------------------
# Main content
# ---------------------------------------------------------------------------

st.title("🧬 Sanger Sequencing Data Analysis")
st.markdown("Upload your AB1 files and run the Sanger pipeline.")

# ---- Mode selection ----
mode = st.radio(
    "Select mode:",
    ["Single pair", "Multiple samples (auto-detect pairs)"],
    horizontal=True,
)

# ---- File upload ----
st.markdown("### 📁 Upload AB1 Files")

if mode == "Single pair":
    col1, col2 = st.columns(2)
    with col1:
        forward_files = st.file_uploader(
            "Forward reads (.ab1)",
            type=["ab1"],
            accept_multiple_files=True,
            key="fwd",
        )
    with col2:
        reverse_files = st.file_uploader(
            "Reverse reads (.ab1)",
            type=["ab1"],
            accept_multiple_files=True,
            key="rev",
        )
else:
    all_files = st.file_uploader(
        "Upload ALL AB1 files (forward + reverse, mixed)",
        type=["ab1"],
        accept_multiple_files=True,
        key="mixed",
        help="Pairs are auto-detected by filename (e.g., *_forward.ab1 / *_reverse.ab1)"
    )

# ---- Run button ----
st.markdown("---")

run_disabled = False
if mode == "Single pair":
    if not forward_files or not reverse_files:
        run_disabled = True
        st.info("⬆️ Upload forward and reverse AB1 files to start.")
else:
    if not all_files:
        run_disabled = True
        st.info("⬆️ Upload AB1 files to start.")

if st.button("🚀 Run Workflow", type="primary", disabled=run_disabled):
    session_id = f"run_{datetime.now().strftime('%Y%m%d_%H%M%S')}"
    output_dir = os.path.join(OUTPUT_BASE, session_id)

    with st.spinner("Setting up..."):
        cleanup_session(session_id)

    returncode = None
    output = ""
    elapsed = 0

    if mode == "Single pair":
        # Save uploaded files
        fwd_dir, fwd_paths = save_uploaded_files(forward_files, f"{session_id}/forward")
        rev_dir, rev_paths = save_uploaded_files(reverse_files, f"{session_id}/reverse")

        st.success(f"📁 Saved {len(fwd_paths)} forward + {len(rev_paths)} reverse files")

        sample_name = None
        if forward_files:
            stem = Path(forward_files[0].name).stem
            for suffix in ["_forward", "_Forward", "_fwd", "_F"]:
                if stem.endswith(suffix):
                    stem = stem[: -len(suffix)]
            sample_name = stem

        # Run workflow with streaming output
        with st.status("🚀 Running workflow...", expanded=True) as status:
            st.write("Setting up pipeline...")
            progress_bar = st.progress(0, text="Starting...")
            returncode, output, elapsed = run_workflow(
                fwd_dir, rev_dir, output_dir, trim_quality,
                sample_name=sample_name,
                status_container=status, progress_bar=progress_bar,
            )
            if returncode == 0:
                status.update(label=f"✅ Workflow completed in {elapsed}s", state="complete")
            else:
                status.update(label=f"❌ Workflow failed after {elapsed}s", state="error")

    else:
        # Save all files to one directory
        upload_dir, saved_paths = save_uploaded_files(all_files, session_id)
        st.success(f"📁 Saved {len(saved_paths)} files")

        # Run wrapper with streaming output
        with st.status("🚀 Running batch processing...", expanded=True) as status:
            st.write("Auto-detecting pairs...")
            progress_bar = st.progress(0, text="Starting...")
            returncode, output, elapsed = run_wrapper_batch(
                upload_dir, output_dir, trim_quality,
                status_container=status, progress_bar=progress_bar,
            )
            if returncode == 0:
                status.update(label=f"✅ Batch completed in {elapsed}s", state="complete")
            else:
                status.update(label=f"❌ Batch failed after {elapsed}s", state="error")

    # Persist results in session state so they survive reruns
    st.session_state.run_returncode = returncode
    st.session_state.run_output = output
    st.session_state.run_elapsed = elapsed
    st.session_state.run_output_dir = output_dir
    st.session_state.run_session_id = session_id
    st.session_state.run_mode = mode

init_session_state()
returncode = st.session_state.run_returncode
output = st.session_state.run_output
elapsed = st.session_state.run_elapsed
output_dir = st.session_state.run_output_dir
session_id = st.session_state.run_session_id
run_mode = st.session_state.run_mode

# ---- Show results ----
if returncode is not None:
    if returncode == 0:
        st.success(f"✅ Workflow completed in {elapsed}s!")
    else:
        st.error("❌ Workflow failed!")

    # Show output log
    with st.expander("📋 Full workflow output", expanded=False):
        st.code(output, language="text")

    # ---- Display results ----
    st.markdown("---")
    st.markdown("## 📊 Results")

    if run_mode == "Single pair":
        # Single sample results
        stats = load_sample_stats(output_dir)
        if stats:
            # Quality metrics cards
            col1, col2, col3, col4, col5 = st.columns(5)
            with col1:
                st.metric("Consensus", f"{stats['consensus_length']} bp")
            with col2:
                st.metric("GC Content", f"{stats['gc_content']:.1f}%")
            with col3:
                st.metric("Ambiguity", f"{stats['ambiguity_pct']:.1f}%")
            with col4:
                st.metric("Coverage", f"{stats['coverage_pct']:.1f}%")
            with col5:
                rating = stats["quality_rating"]
                st.metric("Rating", rating)

            # Detailed stats
            with st.expander("Detailed Statistics"):
                col1, col2 = st.columns(2)
                with col1:
                    st.markdown(f"**Forward file:** {stats['forward_file']}")
                    st.markdown(f"**Reverse file:** {stats['reverse_file']}")
                    st.markdown(f"**Raw length (F/R):** {stats['forward_raw_length']} / {stats['reverse_raw_length']} bp")
                    st.markdown(f"**Trimmed length (F/R):** {stats['forward_trimmed_length']} / {stats['reverse_trimmed_length']} bp")
                    st.markdown(f"**Alignment length:** {stats['alignment_length']} bp")
                with col2:
                    st.markdown(f"**Valid bases:** {stats['valid_bases']} ({stats['coverage_pct']:.1f}%)")
                    st.markdown(f"**Ambiguous bases:** {stats['ambiguous_bases']} ({stats['ambiguity_pct']:.1f}%)")
                    st.markdown(f"**N bases:** {stats['n_bases']} ({stats['n_pct']:.1f}%)")
                    st.markdown(f"**Forward avg quality:** Phred {stats['forward_avg_quality']:.1f}")
                    st.markdown(f"**Reverse avg quality:** Phred {stats['reverse_avg_quality']:.1f}")

        # Consensus sequence
        consensus_path = os.path.join(output_dir, "consensus.fasta")
        header, seq = read_fasta(consensus_path)
        if seq:
            with st.expander("🧬 Consensus Sequence", expanded=True):
                st.code(f">{header}\n{seq}", language="text")
                st.markdown(get_download_link(consensus_path, "consensus.fasta"), unsafe_allow_html=True)

        # Download all results
        st.markdown("### ⬇️ Download Results")
        col1, col2, col3 = st.columns(3)
        with col1:
            st.markdown(get_download_link(os.path.join(output_dir, "final_aligned.fasta"), "final_aligned.fasta"), unsafe_allow_html=True)
        with col2:
            st.markdown(get_download_link(os.path.join(output_dir, "consensus.fasta"), "consensus.fasta"), unsafe_allow_html=True)
        with col3:
            st.markdown(get_download_link(os.path.join(output_dir, "sample_stats.json"), "sample_stats.json"), unsafe_allow_html=True)

    else:
        # Multi-sample results
        all_stats = load_all_sample_stats(output_dir)
        if all_stats:
            # Summary cards
            n = len(all_stats)
            avg_amb = sum(s["ambiguity_pct"] for s in all_stats) / n
            avg_gc = sum(s["gc_content"] for s in all_stats) / n
            avg_cov = sum(s["coverage_pct"] for s in all_stats) / n
            avg_len = sum(s["consensus_length"] for s in all_stats) / n

            col1, col2, col3, col4, col5 = st.columns(5)
            with col1:
                st.metric("Samples", n)
            with col2:
                st.metric("Avg Consensus", f"{avg_len:.0f} bp")
            with col3:
                st.metric("Avg GC", f"{avg_gc:.1f}%")
            with col4:
                st.metric("Avg Ambiguity", f"{avg_amb:.1f}%")
            with col5:
                st.metric("Avg Coverage", f"{avg_cov:.1f}%")

            # Quality distribution
            rating_counts = {}
            for s in all_stats:
                r = s["quality_rating"]
                rating_counts[r] = rating_counts.get(r, 0) + 1

            st.markdown("### Quality Distribution")
            col1, col2, col3, col4 = st.columns(4)
            with col1:
                st.metric("Excellent", rating_counts.get("Excellent", 0))
            with col2:
                st.metric("Good", rating_counts.get("Good", 0))
            with col3:
                st.metric("Fair", rating_counts.get("Fair", 0))
            with col4:
                st.metric("Poor", rating_counts.get("Poor", 0))

            # Sample table
            st.markdown("### Sample Results")
            st.dataframe(
                [
                    {
                        "Sample": s["sample_name"],
                        "Consensus (bp)": s["consensus_length"],
                        "GC %": s["gc_content"],
                        "Ambiguity %": s["ambiguity_pct"],
                        "Coverage %": s["coverage_pct"],
                        "Fwd QV": s["forward_avg_quality"],
                        "Rev QV": s["reverse_avg_quality"],
                        "Rating": s["quality_rating"],
                    }
                    for s in all_stats
                ],
                use_container_width=True,
                hide_index=True,
            )

        # Download combined consensus
        combined_path = os.path.join(output_dir, "all_consensus.fasta")
        if os.path.exists(combined_path):
            st.markdown("### ⬇️ Download Results")
            col1, col2, col3 = st.columns(3)
            with col1:
                st.markdown(get_download_link(combined_path, "all_consensus.fasta"), unsafe_allow_html=True)
            with col2:
                st.markdown(get_download_link(os.path.join(output_dir, "quality_report.html"), "quality_report.html"), unsafe_allow_html=True)
            with col3:
                st.markdown(get_download_link(os.path.join(output_dir, "all_sample_stats.json"), "all_sample_stats.json"), unsafe_allow_html=True)

        # Show individual consensus sequences
        st.markdown("### 🧬 Individual Consensus Sequences")
        if all_stats:
            for s in all_stats:
                sample_name = s["sample_name"]
                consensus_path = os.path.join(output_dir, sample_name, "consensus.fasta")
                header, seq = read_fasta(consensus_path)
                if seq:
                    with st.expander(f"{sample_name} — {len(seq)} bp — {s['quality_rating']}"):
                        st.code(f">{header}\n{seq}", language="text")

# Cleanup button
st.markdown("---")
if returncode is not None and st.button("🗑️ Clean up temp files"):
    cleanup_session(session_id)
    st.success("Cleaned up!")

# ---- Footer ----
st.markdown("---")
st.markdown(
    "Sanger Sequencing Data Analysis",
    help="Upload AB1 files, run the pipeline, and download results."
)
