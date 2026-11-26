#!/usr/bin/env python3
"""
Sanger Sequencing Data Analysis - Streamlit app (English + Lao)
Upload AB1 files -> consensus, chromatogram viewer, F/R disagreement review,
optional mutation report vs a reference, ZIP / HTML report download.
"""

import os
import shutil
import tempfile
import time
import traceback
from datetime import datetime
from pathlib import Path

import streamlit as st

import chromatogram as chrom
import sanger_workflow as wf
import sanger_wrapper as wr
import i18n_extra  # noqa: F401  (registers the extra translations)
from donate import render_donate
from i18n import LANGUAGES, inject_font, set_lang, t

OUTPUT_BASE = os.path.join(tempfile.gettempdir(), "sanger_output")
os.makedirs(OUTPUT_BASE, exist_ok=True)
EXAMPLE_DIR = Path(__file__).parent / "examples"
FEEDBACK_URL = ""  # e.g. "https://forms.gle/..." or "mailto:you@example.com"

MODE_SINGLE, MODE_MULTI = "single", "multi"

st.set_page_config(page_title="Sanger Workflow", page_icon="🧬", layout="wide",
                   initial_sidebar_state="expanded")
st.markdown("<style>.block-container{padding-top:1.2rem}</style>", unsafe_allow_html=True)

for k, v in {"stats": None, "errors": [], "unmatched": [], "out_dir": "", "elapsed": 0.0}.items():
    st.session_state.setdefault(k, v)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
def consensus_html(seq, quals, width=60):
    rows = []
    for start in range(0, len(seq), width):
        cells = []
        for j, (b, q) in enumerate(zip(seq[start:start + width], quals[start:start + width])):
            if b in "NRYSWKMBDHV":
                cell = f'<span style="background:#ffd54f;color:#000;border-radius:2px">{b}</span>'
            elif q < 20:
                cell = f'<span style="color:#9e9e9e">{b}</span>'
            else:
                cell = b
            cells.append(cell + ("&thinsp;" if (j + 1) % 10 == 0 else ""))
        rows.append(f'<div><span style="color:#888;display:inline-block;width:3.6em">{start + 1}</span>'
                    f'{"".join(cells)}</div>')
    return ('<div style="font-family:\'Source Code Pro\',monospace;font-size:15px;'
            f'line-height:1.7;letter-spacing:1px;overflow-x:auto">{"".join(rows)}</div>')


@st.cache_data(show_spinner=False)
def get_trace(path):
    return chrom.load_trace(path)


def file_bytes(path):
    with open(path, "rb") as f:
        return f.read()


def run_jobs(jobs, unmatched, trim_quality, reference):
    """Run all jobs, store results in session_state."""
    old = st.session_state.get("out_dir")
    if old:
        shutil.rmtree(old, ignore_errors=True)
    out_dir = os.path.join(OUTPUT_BASE, f"run_{datetime.now().strftime('%Y%m%d_%H%M%S_%f')}")
    prog = st.progress(0.0)

    def cb(i, n, name):
        if name:
            prog.progress(i / n, text=t("progress", name=name, i=i + 1, n=n))

    t0 = time.perf_counter()
    try:
        stats, errors = wr.run_batch(jobs, out_dir, trim_quality, reference, cb)
    except Exception as e:  # noqa: BLE001
        stats, errors = [], [("run", f"{e}\n{traceback.format_exc()}")]
    prog.empty()
    st.session_state.update(stats=stats, errors=errors, unmatched=unmatched, out_dir=out_dir,
                            elapsed=round(time.perf_counter() - t0, 1))


# ---------------------------------------------------------------------------
# Sidebar
# ---------------------------------------------------------------------------
with st.sidebar:
    lang = st.radio("Language / ພາສາ", list(LANGUAGES), format_func=lambda k: LANGUAGES[k],
                    horizontal=True, key="lang")
    set_lang(lang)
    inject_font()

    st.title(t("title"))
    st.markdown("---")
    st.markdown(f"### {t('settings')}")
    trim_quality = st.slider(t("trim_label"), 0.01, 0.20, 0.05, 0.01, help=t("trim_help"))

    st.markdown("---")
    st.markdown(f"### {t('ref_title')}")
    st.caption(t("ref_help"))
    ref_text = st.text_area(t("ref_text"), height=110, label_visibility="collapsed",
                            placeholder=">my_gene\nATGGCC...")
    ref_file = st.file_uploader(t("ref_upload"), type=["fasta", "fa", "fna", "txt"], key="ref_file")
    if ref_file is not None:
        ref_text = ref_file.getvalue().decode("utf-8", "ignore")
    reference = None
    try:
        reference = wf.parse_reference(ref_text)
        if reference:
            st.success(t("ref_ok", name=reference[0], n=len(reference[1])))
    except ValueError as e:
        st.error(t("ref_invalid", err=e))

    st.markdown("---")
    st.markdown(f"[{t('usage')}](./USAGE.md)")
    if FEEDBACK_URL:
        st.link_button(t("feedback"), FEEDBACK_URL)
    st.markdown(f"<small>{t('cleanup_hint')}</small>", unsafe_allow_html=True)


# ---------------------------------------------------------------------------
# Header + input
# ---------------------------------------------------------------------------
st.title(t("title"))
st.subheader(t("hero"))
st.caption(t("privacy"))

example_jobs = []
if EXAMPLE_DIR.exists():
    example_jobs, _ = wr.pair_files(wr.find_ab1([str(EXAMPLE_DIR)]))
if example_jobs and st.session_state.stats is None:
    c1, c2 = st.columns([1, 3])
    with c1:
        try_example = st.button(t("try_example"), type="secondary")
    with c2:
        st.caption(t("example_hint"))
    if try_example:
        run_jobs(example_jobs, [], trim_quality, reference)
        st.rerun()

mode = st.radio(t("mode_label"), [MODE_SINGLE, MODE_MULTI], horizontal=True,
                format_func=lambda m: t("mode_single") if m == MODE_SINGLE else t("mode_multi"))
st.markdown(t("upload_header"))

fwd_up = rev_up = all_up = None
if mode == MODE_SINGLE:
    st.caption(t("single_hint"))
    c1, c2 = st.columns(2)
    fwd_up = c1.file_uploader(t("fwd_one"), type=["ab1"], key="fwd1")
    rev_up = c2.file_uploader(t("rev_one"), type=["ab1"], key="rev1")
    disabled = not (fwd_up or rev_up)
    if disabled:
        st.info(t("info_single"))
else:
    all_up = st.file_uploader(t("mixed_label"), type=["ab1"], accept_multiple_files=True,
                              key="mixed", help=t("mixed_help"))
    disabled = not all_up
    if disabled:
        st.info(t("info_multi"))

if st.button(t("run"), type="primary", disabled=disabled):
    with tempfile.TemporaryDirectory(prefix="sanger_uploads_") as tmp:
        def save(upload):
            path = os.path.join(tmp, os.path.basename(upload.name))
            with open(path, "wb") as fh:
                fh.write(upload.getbuffer())
            return path

        if mode == MODE_SINGLE:
            f = save(fwd_up) if fwd_up else None
            r = save(rev_up) if rev_up else None
            name = wr.sample_name_from_file((fwd_up or rev_up).name)
            jobs, unmatched = [(name, f, r)], []
        else:
            jobs, unmatched = wr.pair_files([save(u) for u in all_up])
        if jobs:
            run_jobs(jobs, [os.path.basename(u) for u in unmatched], trim_quality, reference)
        else:
            st.session_state.update(stats=[], errors=[("-", t("no_files_ok"))], unmatched=[],
                                    out_dir="", elapsed=0.0)


# ---------------------------------------------------------------------------
# Results
# ---------------------------------------------------------------------------
def tab_consensus(s, out_dir):
    st.markdown(consensus_html(s["consensus_seq"], s["consensus_quals"]), unsafe_allow_html=True)
    st.caption(f"{t('legend')}: ◼ {t('legend_hq')} · "
               f"<span style='color:#9e9e9e'>◼ {t('legend_lq')}</span> · "
               f"<span style='background:#ffd54f;color:#000'>&nbsp;◼ {t('legend_amb')}&nbsp;</span>",
               unsafe_allow_html=True)
    st.code(f">{s['sample_name']}\n{s['consensus_seq']}", language="text")
    st.download_button(t("dl_fasta"), f">{s['sample_name']}\n{s['consensus_seq']}\n",
                       file_name=f"{s['sample_name']}_consensus.fasta", key=f"dlc_{s['sample_name']}")


def tab_chromatogram(s, out_dir):
    name = s["sample_name"]
    ab1 = os.path.join(out_dir, name, "ab1")
    reads = {}
    if s["forward_file"]:
        reads["forward"] = (os.path.join(ab1, s["forward_file"]), s["forward_trim"], "forward_pos")
    if s["reverse_file"]:
        reads["reverse"] = (os.path.join(ab1, s["reverse_file"]), s["reverse_trim"], "reverse_pos")
    if not reads:
        return
    choice = st.radio(t("chrom_read"), list(reads), horizontal=True, key=f"cr_{name}",
                      format_func=lambda k: t("chrom_fwd") if k == "forward" else t("chrom_rev"))
    path, trim, pos_key = reads[choice]
    try:
        trace = get_trace(path)
    except Exception as e:  # noqa: BLE001
        st.error(t("chrom_error", err=e))
        return
    n = len(trace["bases"])
    dis = s["disagreements"]
    c1, c2, c3 = st.columns([2, 1, 3])
    with c1:
        sel = st.selectbox(t("chrom_goto"), [None] + list(range(len(dis))), key=f"cg_{name}_{choice}",
                           format_func=lambda i: t("chrom_none") if i is None else
                           f"#{dis[i]['consensus_pos']}  F:{dis[i]['forward_base']}(Q{dis[i]['forward_q']}) "
                           f"R:{dis[i]['reverse_base']}(Q{dis[i]['reverse_q']}) → {dis[i]['call']}")
    with c2:
        window = st.select_slider(t("chrom_window"), [30, 60, 100, 150], value=60, key=f"cw_{name}_{choice}")
    with c3:
        default = min(n, max(1, trim[0] + 30))
        center = st.slider(t("chrom_center"), 1, max(2, n), default, key=f"cc_{name}_{choice}")
    highlight = None
    if sel is not None:
        highlight = center = dis[sel][pos_key]
    fig = chrom.make_figure(trace, center, window, trim=tuple(trim), highlight=highlight)
    st.plotly_chart(fig, key=f"fig_{name}_{choice}_{center}_{window}")
    st.caption(t("chrom_hint"))


def tab_differences(s):
    st.markdown(f"#### {t('disagree_title')}")
    dis = s["disagreements"]
    if not dis:
        st.success(t("disagree_none"))
    else:
        st.caption(t("disagree_hint"))
        st.dataframe([{
            t("col_cpos"): d["consensus_pos"], t("col_call"): d["call"], t("col_how"): d["resolution"],
            t("col_fpos"): d["forward_pos"], t("col_fbase"): d["forward_base"], t("col_fq"): d["forward_q"],
            t("col_rpos"): d["reverse_pos"], t("col_rbase"): d["reverse_base"], t("col_rq"): d["reverse_q"],
        } for d in dis], hide_index=True)

    st.markdown(f"#### {t('var_title')}")
    ref = s.get("reference")
    if not ref:
        st.info(t("var_none_ref"))
        return
    st.caption(t("var_summary", name=ref["reference_name"], ident=ref["identity"],
                 cov=ref["ref_coverage_pct"], strand=ref["strand"]))
    if not ref["variants"]:
        st.success(t("var_zero"))
    else:
        st.dataframe([{
            t("col_type"): v["type"], t("col_refpos"): v["ref_pos"], t("col_ref"): v["ref"],
            t("col_alt"): v["alt"], t("col_cpos"): v["consensus_pos"], t("col_q"): v["quality"],
        } for v in ref["variants"]], hide_index=True)


def tab_details(s):
    c1, c2 = st.columns(2)
    with c1:
        st.markdown(f"**{t('fwd_file')}:** {s['forward_file'] or '—'}")
        st.markdown(f"**{t('rev_file')}:** {s['reverse_file'] or '—'}")
        st.markdown(f"**{t('raw_len')}:** {s['forward_raw_length']} / {s['reverse_raw_length']} bp")
        st.markdown(f"**{t('trim_len')}:** {s['forward_trimmed_length']} / {s['reverse_trimmed_length']} bp")
        st.markdown(f"**{t('aln_len')}:** {s['alignment_length']} bp")
    with c2:
        st.markdown(f"**{t('valid_bases')}:** {s['valid_bases']} ({s['coverage_pct']:.1f}%)")
        st.markdown(f"**{t('amb_bases')}:** {s['ambiguous_bases']} ({s['ambiguity_pct']:.1f}%)")
        st.markdown(f"**{t('n_bases')}:** {s['n_bases']} ({s['n_pct']:.1f}%)")
        st.markdown(f"**{t('fwd_q')}:** Phred {s['forward_avg_quality']:.1f}")
        st.markdown(f"**{t('rev_q')}:** Phred {s['reverse_avg_quality']:.1f}")
        st.markdown(f"**{t('details_meanq')}:** Phred {s['mean_consensus_quality']:.1f}")
        st.markdown(f"**{t('details_time')}:** {s['elapsed_s']} s")


def show_sample(s, out_dir):
    c = st.columns(6)
    c[0].metric(t("m_consensus"), f"{s['consensus_length']} bp")
    c[1].metric(t("m_gc"), f"{s['gc_content']:.1f}%")
    c[2].metric(t("m_q20"), f"{s['q20_pct']:.0f}%")
    c[3].metric(t("m_amb"), f"{s['ambiguity_pct']:.1f}%")
    c[4].metric(t("m_disagree"), len(s["disagreements"]))
    c[5].metric(t("m_rating"), t(s["quality_rating"]))
    st.caption(f"{t('col_mode')}: {t('mode_' + s['mode'])}"
               + (f" · {t('m_overlap')}: {s['overlap_length']} bp ({s['overlap_identity']:.0f}%)"
                  if s["mode"] == "paired" else ""))
    for w in s["warnings"]:
        st.warning(w)

    tabs = st.tabs([t("tab_consensus"), t("tab_chrom"), t("tab_variants"), t("tab_details")])
    with tabs[0]:
        tab_consensus(s, out_dir)
    with tabs[1]:
        tab_chromatogram(s, out_dir)
    with tabs[2]:
        tab_differences(s)
    with tabs[3]:
        tab_details(s)


stats = st.session_state.stats
if stats is not None:
    out_dir = st.session_state.out_dir
    st.markdown("---")
    st.markdown(t("results"))

    for sample, msg in st.session_state.errors:
        st.error(f"{t('errors_title')}: **{sample}** — {msg}")
    if st.session_state.unmatched:
        st.warning(t("unmatched", n=len(st.session_state.unmatched),
                     names=", ".join(st.session_state.unmatched[:8])))

    if stats:
        st.success(t("run_done", n=len(stats), s=st.session_state.elapsed))

        # downloads
        zip_path = wr.make_zip(out_dir)
        d = st.columns(4)
        d[0].download_button(t("dl_zip"), file_bytes(zip_path), file_name="sanger_results.zip", type="primary")
        d[1].download_button(t("dl_report"), file_bytes(os.path.join(out_dir, "quality_report.html")),
                             file_name="quality_report.html", mime="text/html")
        d[2].download_button(t("dl_fasta"), file_bytes(os.path.join(out_dir, "all_consensus.fasta")),
                             file_name="all_consensus.fasta")
        d[3].download_button(t("dl_csv"), file_bytes(os.path.join(out_dir, "summary.csv")),
                             file_name="summary.csv", mime="text/csv")

        if len(stats) > 1:
            st.markdown(t("sample_results"))
            has_ref = any(s.get("reference") for s in stats)
            st.dataframe([{
                t("col_sample"): s["sample_name"], t("col_mode"): t("mode_" + s["mode"]),
                t("col_consensus"): s["consensus_length"], t("col_gc"): s["gc_content"],
                t("m_q20"): s["q20_pct"], t("col_amb"): s["ambiguity_pct"],
                t("m_disagree"): len(s["disagreements"]),
                **({t("col_nvar"): wr.n_variants(s)} if has_ref else {}),
                t("col_rating"): t(s["quality_rating"]),
            } for s in stats], hide_index=True)
            names = [s["sample_name"] for s in stats]
            chosen = st.selectbox(t("select_sample"), names, key="sample_pick")
            current = stats[names.index(chosen)]
        else:
            current = stats[0]
        show_sample(current, out_dir)

    if st.button(t("cleanup")):
        shutil.rmtree(out_dir, ignore_errors=True)
        if out_dir:
            Path(out_dir.rstrip("/\\") + "_results.zip").unlink(missing_ok=True)
        st.session_state.update(stats=None, errors=[], unmatched=[], out_dir="")
        st.success(t("cleaned"))
        st.rerun()

st.markdown("---")
with st.expander(t("how_title")):
    st.markdown(t("how_text"))
    st.markdown(t("steps_title"))
    st.markdown(t("steps"))
render_donate(expanded=False)
st.markdown(t("footer"), help=t("footer_help"))
