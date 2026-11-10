#!/usr/bin/env python3
"""
Sanger Sequencing Data analysis workflow
Streamlit interface for uploading AB1 files, running the pipeline, and viewing results.
Supports English and Lao (Noto Sans Lao).
"""

import os
import json
import shutil
import tempfile
import base64
import traceback
from pathlib import Path
from datetime import datetime

import streamlit as st

from sanger_workflow import run_workflow as workflow_run
from sanger_wrapper import run_batch_workflow
from i18n import LANGUAGES, set_lang, t, inject_font
from donate import render_donate

# Global exception handler to show full traceback in Streamlit
def _excepthook(exc_type, exc_value, exc_tb):
    st.error(f"Unhandled exception: {exc_type.__name__}: {exc_value}")
    st.code("".join(traceback.format_exception(exc_type, exc_value, exc_tb)), language="text")

import sys
sys.excepthook = _excepthook

# ---------------------------------------------------------------------------
# Config
# ---------------------------------------------------------------------------
OUTPUT_BASE = os.path.join(tempfile.gettempdir(), "sanger_output")
os.makedirs(OUTPUT_BASE, exist_ok=True)

# Internal (language-independent) mode keys
MODE_SINGLE = "single"
MODE_MULTI = "multi"


# ---------------------------------------------------------------------------
# Page config
# ---------------------------------------------------------------------------
st.set_page_config(
    page_title="Sanger Workflow",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded",
)

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

def run_workflow(forward_dir, reverse_dir, output_dir, trim_quality=0.05, sample_name=None, status_container=None, progress_bar=None):
    """Run sanger_workflow.py directly (no subprocess)."""
    result = workflow_run(forward_dir, reverse_dir, output_dir, trim_quality=trim_quality, sample_name=sample_name)
    if isinstance(result, tuple) and len(result) == 2:
        returncode, output = result
    elif isinstance(result, int):
        returncode, output = result, ""
    else:
        returncode, output = 0, str(result)
    return returncode, output, 0


def run_wrapper_batch(upload_dir, output_dir, trim_quality=0.05, status_container=None, progress_bar=None):
    """Run sanger_wrapper.py batch mode directly (no subprocess)."""
    try:
        result = run_batch_workflow(upload_dir, output_dir, trim_quality=trim_quality)
    except Exception as e:
        return 1, f"Batch wrapper error: {e}\n" + traceback.format_exc(), 0
    if isinstance(result, tuple) and len(result) == 2:
        returncode, output = result
    elif isinstance(result, int):
        returncode, output = result, ""
    else:
        returncode, output = 0, str(result)
    return returncode, output, 0


def load_json(path):
    if os.path.exists(path):
        with open(path, encoding="utf-8") as f:
            return json.load(f)
    return None


def load_sample_stats(output_dir):
    return load_json(os.path.join(output_dir, "sample_stats.json"))


def load_all_sample_stats(output_dir):
    return load_json(os.path.join(output_dir, "all_sample_stats.json"))


def read_fasta(filepath):
    """Read a FASTA file and return header + sequence."""
    if not os.path.exists(filepath):
        return None, None
    with open(filepath) as f:
        lines = f.readlines()
    header, seq = "", ""
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
    return f'<a href="data:file/octet-stream;base64,{b64}" download="{filename}">⬇️ {filename}</a>'


def show_download(filepath, filename):
    link = get_download_link(filepath, filename)
    if link:
        st.markdown(link, unsafe_allow_html=True)


def cleanup_session(session_id):
    output_dir = os.path.join(OUTPUT_BASE, session_id)
    if os.path.exists(output_dir):
        shutil.rmtree(output_dir, ignore_errors=True)


def init_session_state():
    defaults = {
        "run_returncode": None,
        "run_output": "",
        "run_elapsed": 0,
        "run_output_dir": "",
        "run_session_id": "",
        "run_mode": MODE_SINGLE,
    }
    for key, value in defaults.items():
        if key not in st.session_state:
            st.session_state[key] = value


# ---------------------------------------------------------------------------
# Sidebar
# ---------------------------------------------------------------------------

with st.sidebar:
    lang = st.radio(
        "Language / ພາສາ",
        options=list(LANGUAGES.keys()),
        format_func=lambda k: LANGUAGES[k],
        horizontal=True,
        key="lang",
    )
    set_lang(lang)
    inject_font()  # Noto Sans Lao

    st.title(t("title"))
    st.markdown("---")
    st.markdown(t("sidebar_pipeline"))
    st.markdown(t("steps_title"))
    st.markdown(t("steps"))
    st.markdown("---")
    st.markdown(f"### {t('settings')}")
    trim_quality = st.slider(
        t("trim_label"),
        min_value=0.01,
        max_value=0.20,
        value=0.05,
        step=0.01,
        help=t("trim_help"),
    )
    st.markdown("---")
    st.markdown(f"### {t('links')}")
    st.markdown(f"[{t('usage')}](./USAGE.md)")
    st.markdown("---")
    st.markdown(f"<small>{t('cleanup_hint')}</small>", unsafe_allow_html=True)


# ---------------------------------------------------------------------------
# Main content
# ---------------------------------------------------------------------------

st.title(t("title"))
st.markdown(t("intro"))

render_donate(expanded=True)  # donate box at the top, open by default

mode = st.radio(
    t("mode_label"),
    [MODE_SINGLE, MODE_MULTI],
    format_func=lambda m: t("mode_single") if m == MODE_SINGLE else t("mode_multi"),
    horizontal=True,
)

st.markdown(t("upload_header"))

if mode == MODE_SINGLE:
    col1, col2 = st.columns(2)
    with col1:
        forward_files = st.file_uploader(
            t("fwd_label"), type=["ab1"], accept_multiple_files=True, key="fwd",
        )
    with col2:
        reverse_files = st.file_uploader(
            t("rev_label"), type=["ab1"], accept_multiple_files=True, key="rev",
        )
else:
    all_files = st.file_uploader(
        t("mixed_label"), type=["ab1"], accept_multiple_files=True,
        key="mixed", help=t("mixed_help"),
    )

st.markdown("---")

run_disabled = False
if mode == MODE_SINGLE:
    if not forward_files or not reverse_files:
        run_disabled = True
        st.info(t("info_single"))
else:
    if not all_files:
        run_disabled = True
        st.info(t("info_multi"))

if st.button(t("run"), type="primary", disabled=run_disabled):
    session_id = f"run_{datetime.now().strftime('%Y%m%d_%H%M%S')}"
    output_dir = os.path.join(OUTPUT_BASE, session_id)
    os.makedirs(output_dir, exist_ok=True)

    with tempfile.TemporaryDirectory(prefix="sanger_uploads_") as upload_tmp:
        returncode = None
        output = ""
        elapsed = 0

        if mode == MODE_SINGLE:
            fwd_dir = os.path.join(upload_tmp, "forward")
            rev_dir = os.path.join(upload_tmp, "reverse")
            os.makedirs(fwd_dir, exist_ok=True)
            os.makedirs(rev_dir, exist_ok=True)

            for f in forward_files:
                with open(os.path.join(fwd_dir, f.name), "wb") as out:
                    out.write(f.getbuffer())
            for f in reverse_files:
                with open(os.path.join(rev_dir, f.name), "wb") as out:
                    out.write(f.getbuffer())

            st.success(t("saved_pair", nf=len(forward_files), nr=len(reverse_files)))

            sample_name = None
            if forward_files:
                stem = Path(forward_files[0].name).stem
                for suffix in ["_forward", "_Forward", "_fwd", "_F"]:
                    if stem.endswith(suffix):
                        stem = stem[: -len(suffix)]
                sample_name = stem

            with st.status(t("running"), expanded=True) as status:
                st.write(t("setup"))
                progress_bar = st.progress(0, text=t("starting"))
                returncode, output, elapsed = run_workflow(
                    fwd_dir, rev_dir, output_dir, trim_quality,
                    sample_name=sample_name,
                    status_container=status, progress_bar=progress_bar,
                )
                if returncode == 0:
                    status.update(label=t("done", s=elapsed), state="complete")
                else:
                    status.update(label=t("failed", s=elapsed), state="error")

        else:
            upload_dir = os.path.join(upload_tmp, session_id)
            os.makedirs(upload_dir, exist_ok=True)

            for f in all_files:
                with open(os.path.join(upload_dir, f.name), "wb") as out:
                    out.write(f.getbuffer())

            st.success(t("saved_n", n=len(all_files)))

            try:
                with st.status(t("batch_running"), expanded=True) as status:
                    st.write(t("detect_pairs"))
                    progress_bar = st.progress(0, text=t("starting"))
                    returncode, output, elapsed = run_wrapper_batch(
                        upload_dir, output_dir, trim_quality,
                        status_container=status, progress_bar=progress_bar,
                    )
                    if returncode == 0:
                        status.update(label=t("batch_done", s=elapsed), state="complete")
                    else:
                        status.update(label=t("batch_failed", s=elapsed), state="error")
            except Exception as e:
                returncode = 1
                output = f"Error: {e}\n" + traceback.format_exc()
                elapsed = 0

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
        st.success(t("final_ok", s=elapsed))
    else:
        st.error(t("final_fail"))

    with st.expander(t("full_output"), expanded=False):
        st.code(output, language="text")

    st.markdown("---")
    st.markdown(t("results"))

    if run_mode == MODE_SINGLE:
        stats = load_sample_stats(output_dir)
        if stats:
            col1, col2, col3, col4, col5 = st.columns(5)
            with col1:
                st.metric(t("m_consensus"), f"{stats['consensus_length']} bp")
            with col2:
                st.metric(t("m_gc"), f"{stats['gc_content']:.1f}%")
            with col3:
                st.metric(t("m_amb"), f"{stats['ambiguity_pct']:.1f}%")
            with col4:
                st.metric(t("m_cov"), f"{stats['coverage_pct']:.1f}%")
            with col5:
                st.metric(t("m_rating"), t(stats["quality_rating"]))

            with st.expander(t("detailed")):
                col1, col2 = st.columns(2)
                with col1:
                    st.markdown(f"**{t('fwd_file')}:** {stats['forward_file']}")
                    st.markdown(f"**{t('rev_file')}:** {stats['reverse_file']}")
                    st.markdown(f"**{t('raw_len')}:** {stats['forward_raw_length']} / {stats['reverse_raw_length']} bp")
                    st.markdown(f"**{t('trim_len')}:** {stats['forward_trimmed_length']} / {stats['reverse_trimmed_length']} bp")
                    st.markdown(f"**{t('aln_len')}:** {stats['alignment_length']} bp")
                with col2:
                    st.markdown(f"**{t('valid_bases')}:** {stats['valid_bases']} ({stats['coverage_pct']:.1f}%)")
                    st.markdown(f"**{t('amb_bases')}:** {stats['ambiguous_bases']} ({stats['ambiguity_pct']:.1f}%)")
                    st.markdown(f"**{t('n_bases')}:** {stats['n_bases']} ({stats['n_pct']:.1f}%)")
                    st.markdown(f"**{t('fwd_q')}:** Phred {stats['forward_avg_quality']:.1f}")
                    st.markdown(f"**{t('rev_q')}:** Phred {stats['reverse_avg_quality']:.1f}")

        consensus_path = os.path.join(output_dir, "consensus.fasta")
        header, seq = read_fasta(consensus_path)
        if seq:
            with st.expander(t("consensus_seq"), expanded=True):
                st.code(f">{header}\n{seq}", language="text")
                show_download(consensus_path, "consensus.fasta")

        st.markdown(t("download"))
        col1, col2, col3 = st.columns(3)
        with col1:
            show_download(os.path.join(output_dir, "final_aligned.fasta"), "final_aligned.fasta")
        with col2:
            show_download(consensus_path, "consensus.fasta")
        with col3:
            show_download(os.path.join(output_dir, "sample_stats.json"), "sample_stats.json")

    else:
        all_stats = load_all_sample_stats(output_dir)
        if all_stats:
            n = len(all_stats)
            avg_amb = sum(s["ambiguity_pct"] for s in all_stats) / n
            avg_gc = sum(s["gc_content"] for s in all_stats) / n
            avg_cov = sum(s["coverage_pct"] for s in all_stats) / n
            avg_len = sum(s["consensus_length"] for s in all_stats) / n

            col1, col2, col3, col4, col5 = st.columns(5)
            with col1:
                st.metric(t("samples"), n)
            with col2:
                st.metric(t("avg_consensus"), f"{avg_len:.0f} bp")
            with col3:
                st.metric(t("avg_gc"), f"{avg_gc:.1f}%")
            with col4:
                st.metric(t("avg_amb"), f"{avg_amb:.1f}%")
            with col5:
                st.metric(t("avg_cov"), f"{avg_cov:.1f}%")

            rating_counts = {}
            for s in all_stats:
                r = s["quality_rating"]
                rating_counts[r] = rating_counts.get(r, 0) + 1

            st.markdown(t("quality_dist"))
            cols = st.columns(4)
            for c, name in zip(cols, ["Excellent", "Good", "Fair", "Poor"]):
                with c:
                    st.metric(t(name), rating_counts.get(name, 0))

            st.markdown(t("sample_results"))
            st.dataframe(
                [
                    {
                        t("col_sample"): s["sample_name"],
                        t("col_consensus"): s["consensus_length"],
                        t("col_gc"): s["gc_content"],
                        t("col_amb"): s["ambiguity_pct"],
                        t("col_cov"): s["coverage_pct"],
                        t("col_fwdqv"): s["forward_avg_quality"],
                        t("col_revqv"): s["reverse_avg_quality"],
                        t("col_rating"): t(s["quality_rating"]),
                    }
                    for s in all_stats
                ],
                use_container_width=True,
                hide_index=True,
            )

        combined_path = os.path.join(output_dir, "all_consensus.fasta")
        if os.path.exists(combined_path):
            st.markdown(t("download"))
            col1, col2, col3 = st.columns(3)
            with col1:
                show_download(combined_path, "all_consensus.fasta")
            with col2:
                show_download(os.path.join(output_dir, "quality_report.html"), "quality_report.html")
            with col3:
                show_download(os.path.join(output_dir, "all_sample_stats.json"), "all_sample_stats.json")

        st.markdown(t("individual"))
        if all_stats:
            for s in all_stats:
                sample_name = s["sample_name"]
                consensus_path = os.path.join(output_dir, sample_name, "consensus.fasta")
                header, seq = read_fasta(consensus_path)
                if seq:
                    with st.expander(f"{sample_name} — {len(seq)} bp — {t(s['quality_rating'])}"):
                        st.code(f">{header}\n{seq}", language="text")

st.markdown("---")
if returncode is not None and st.button(t("cleanup")):
    cleanup_session(session_id)
    st.success(t("cleaned"))

st.markdown("---")
st.markdown(t("footer"), help=t("footer_help"))
