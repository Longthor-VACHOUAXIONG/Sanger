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


def run_jobs(jobs, unmatched, trim_quality, reference, min_base_q=10):
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
        stats, errors = wr.run_batch(jobs, out_dir, trim_quality, reference, cb, min_base_q)
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
    min_base_q = st.slider(t("minq_label"), 0, 30, 10, 1, help=t("minq_help"))

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
        run_jobs(example_jobs, [], trim_quality, reference, min_base_q)
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
            run_jobs(jobs, [os.path.basename(u) for u in unmatched], trim_quality, reference, min_base_q)
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


@st.cache_data(show_spinner=False)
def get_trace_rc(path):
    return chrom.reverse_complement_trace(chrom.load_trace(path))


def _set_pos(name, pos, L):
    st.session_state[f"cpos_{name}"] = max(1, min(L, pos))


def _jump(name):
    i = st.session_state.get(f"cg_{name}")
    if i is not None:
        st.session_state[f"cpos_{name}"] = st.session_state[f"dis_{name}"][i]


def _next_dis(name, direction):
    cur, dis = st.session_state[f"cpos_{name}"], sorted(st.session_state[f"dis_{name}"])
    cand = [p for p in dis if p > cur] if direction > 0 else [p for p in dis if p < cur][::-1]
    if cand:
        st.session_state[f"cpos_{name}"] = cand[0]


def tab_chromatogram(s, out_dir):
    name, seq, quals, L = s["sample_name"], s["consensus_seq"], s["consensus_quals"], s["consensus_length"]
    dis = s["disagreements"]
    dis_pos = [d["consensus_pos"] for d in dis]
    pos_key, ctr_key = f"cpos_{name}", f"ctr_{name}"
    st.session_state[f"dis_{name}"] = dis_pos
    st.session_state.setdefault(pos_key, dis_pos[0] if dis_pos else 1)
    st.session_state.setdefault(ctr_key, 0)
    st.session_state[pos_key] = max(1, min(L, st.session_state[pos_key]))

    # 1) clickable consensus map
    st.caption(t("chrom_click_hint"))
    ev = st.plotly_chart(
        chrom.make_consensus_map(seq, quals, set(dis_pos), st.session_state[pos_key]),
        key=f"map_{name}_{st.session_state[ctr_key]}", on_select="rerun",
        selection_mode="points", config={"displayModeBar": False})
    try:
        pts = ev.selection.points
    except Exception:  # noqa: BLE001
        pts = []
    if pts and pts[0].get("point_index") is not None:
        st.session_state[pos_key] = int(pts[0]["point_index"]) + 1
        st.session_state[ctr_key] += 1          # fresh chart -> same base can be clicked again
        st.rerun()
    st.caption(f"🟩 {t('legend_hq')} · ⬜ {t('legend_lq')} · 🟨 {t('legend_amb')} · 🟧 {t('chrom_disagree_key')}")

    # 2) navigation
    c = st.columns([0.5, 0.5, 0.5, 0.5, 1.6, 2.2, 1.6])
    c[0].button("◀", key=f"pv_{name}", on_click=_set_pos, args=(name, st.session_state[pos_key] - 1, L))
    c[1].button("▶", key=f"nx_{name}", on_click=_set_pos, args=(name, st.session_state[pos_key] + 1, L))
    c[2].button("⏮", key=f"pd_{name}", on_click=_next_dis, args=(name, -1), disabled=not dis_pos,
                help=t("chrom_prev_dis"))
    c[3].button("⏭", key=f"nd_{name}", on_click=_next_dis, args=(name, 1), disabled=not dis_pos,
                help=t("chrom_next_dis"))
    c[4].number_input(t("chrom_cpos"), 1, max(1, L), key=pos_key)
    c[5].selectbox(t("chrom_goto"), [None] + list(range(len(dis))), key=f"cg_{name}", on_change=_jump,
                   args=(name,),
                   format_func=lambda i: t("chrom_none") if i is None else
                   f"#{dis[i]['consensus_pos']}  F:{dis[i]['forward_base']}(Q{dis[i]['forward_q']}) "
                   f"R:{dis[i]['reverse_base']}(Q{dis[i]['reverse_q']}) → {dis[i]['call']}")
    window = c[6].select_slider(t("chrom_window"), [20, 40, 60, 100], value=40, key=f"cw_{name}")

    # 3) both reads around the selected consensus base
    cpos = st.session_state[pos_key]
    st.markdown(f"**{t('chrom_cpos')} {cpos}:** `{seq[cpos - 1]}` · Q{quals[cpos - 1]} · "
                f"{t('chrom_source')}: {s['consensus_src'][cpos - 1]}")
    ab1 = os.path.join(out_dir, name, "ab1")
    fm, rm = s["consensus_fmap"][cpos - 1], s["consensus_rmap"][cpos - 1]
    panels = []
    try:
        if s["forward_file"]:
            if fm >= 0:
                tr = get_trace(os.path.join(ab1, s["forward_file"]))
                c1 = fm + 1
                panels.append({"label": f"{t('chrom_fwd')} — {t('chrom_base')} {c1}: "
                                        f"{tr['bases'][c1 - 1]} (Q{tr['quals'][c1 - 1]})",
                               "trace": tr, "center": c1, "trim": tuple(s["forward_trim"]),
                               "called": (tr["bases"][c1 - 1], tr["quals"][c1 - 1])})
            else:
                st.info(t("chrom_nocover", read=t("chrom_fwd")))
        if s["reverse_file"]:
            if rm >= 0:
                tr = get_trace_rc(os.path.join(ab1, s["reverse_file"]))
                n = len(tr["bases"])
                c1 = n - rm                      # 1-based position in the reverse-complemented trace
                e0, e1 = s["reverse_trim"]
                panels.append({"label": f"{t('chrom_rev')} — {t('chrom_base')} {c1}: "
                                        f"{tr['bases'][c1 - 1]} (Q{tr['quals'][c1 - 1]})",
                               "trace": tr, "center": c1, "trim": (n - e1, n - e0),
                               "called": (tr["bases"][c1 - 1], tr["quals"][c1 - 1])})
            else:
                st.info(t("chrom_nocover", read=t("chrom_rev")))
    except Exception as e:  # noqa: BLE001
        st.error(t("chrom_error", err=e))
        return
    if seq[cpos - 1] == "N" and any(p["called"][0] != "N" for p in panels):
        calls = ", ".join(f"{p['called'][0]} (Q{p['called'][1]})" for p in panels)
        st.info(t("chrom_masked", calls=calls, q=s.get("min_base_q", 10)))
    if panels:
        st.plotly_chart(chrom.make_stacked_figure(panels, window), key=f"fig_{name}_{cpos}_{window}")
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

    # A radio (not st.tabs) so the chosen view survives reruns, e.g. after clicking a base
    views = {"consensus": (t("tab_consensus"), tab_consensus), "chrom": (t("tab_chrom"), tab_chromatogram),
             "diff": (t("tab_variants"), lambda s, d: tab_differences(s)),
             "details": (t("tab_details"), lambda s, d: tab_details(s))}
    view = st.radio("view", list(views), format_func=lambda k: views[k][0], horizontal=True,
                    key=f"view_{s['sample_name']}", label_visibility="collapsed")
    views[view][1](s, out_dir)


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
