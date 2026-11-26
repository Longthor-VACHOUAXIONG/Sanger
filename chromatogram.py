"""
Chromatogram viewer for AB1 files (interactive Plotly).

- make_stacked_figure(): forward and reverse stacked, each centred on the same
  consensus base, with quality bars, trimmed regions greyed and the base highlighted.
- make_consensus_map(): clickable map of the whole consensus (colour = quality/ambiguity).
"""

import math

from Bio import SeqIO
import plotly.graph_objects as go
from plotly.subplots import make_subplots

COLORS = {"A": "#2ca02c", "C": "#1f77b4", "G": "#7f7f7f", "T": "#d62728"}
OTHER = "#9467bd"
_COMP = dict(zip("ACGTRYSWKMBDHVN", "TGCAYRSWMKVHDBN"))   # full IUPAC complement


def _text(v):
    return v.decode("ascii", "ignore") if isinstance(v, (bytes, bytearray)) else str(v)


def load_trace(path: str) -> dict:
    """Return {channels:{base:[...]}, ploc:[...], bases:str, quals:[...]}."""
    rec = SeqIO.read(path, "abi")
    raw = rec.annotations["abif_raw"]
    order = _text(raw.get("FWO_1", "GATC"))[:4]
    channels = {}
    for i, base in enumerate(order):
        data = raw.get(f"DATA{9 + i}")
        if data is not None:
            channels[base] = [int(x) for x in data]
    ploc = [int(x) for x in (raw.get("PLOC2") or raw.get("PLOC1") or [])]
    bases = str(rec.seq).upper()
    quals = list(rec.letter_annotations.get("phred_quality") or [])
    n = min(len(bases), len(ploc))
    if len(quals) != len(bases):
        quals = [0] * len(bases)
    return {"channels": channels, "ploc": ploc[:n], "bases": bases[:n], "quals": quals[:n]}


def reverse_complement_trace(trace: dict) -> dict:
    """Mirror a trace so it reads in the reverse-complement orientation."""
    length = max(len(d) for d in trace["channels"].values())
    return {
        "channels": {_COMP.get(b, "N"): d[::-1] for b, d in trace["channels"].items()},
        "ploc": [length - 1 - p for p in reversed(trace["ploc"])],
        "bases": "".join(_COMP.get(b, "N") for b in reversed(trace["bases"])),
        "quals": list(reversed(trace["quals"])),
    }


def _window(center, window, n):
    window = max(10, min(window, n))
    start = max(0, center - 1 - window // 2)
    end = min(n, start + window)
    return max(0, end - window), end


def _draw(fig, p, window, rt, rq, show_legend):
    trace = p["trace"]
    ploc, bases, quals = trace["ploc"], trace["bases"], trace["quals"]
    n = len(ploc)
    if n == 0:
        return
    start, end = _window(p["center"], window, n)

    def edge(i):                      # scan position of the boundary before base i
        if i <= 0:
            return ploc[0] - 6
        if i >= n:
            return ploc[-1] + 6
        return (ploc[i - 1] + ploc[i]) / 2

    x0, x1 = max(0, int(edge(start))), int(edge(end))
    ymax = 1
    for base, data in trace["channels"].items():
        seg = data[x0:x1 + 1]
        ymax = max(ymax, max(seg) if seg else 1)
        fig.add_trace(go.Scatter(x=list(range(x0, x0 + len(seg))), y=seg, mode="lines", name=base,
                                 legendgroup=base, showlegend=show_legend,
                                 line=dict(color=COLORS.get(base, OTHER), width=1.6), hoverinfo="skip"),
                      row=rt, col=1)
    idx = list(range(start, end))
    fig.add_trace(go.Scatter(
        x=[ploc[i] for i in idx], y=[ymax * 1.08] * len(idx), mode="text", text=[bases[i] for i in idx],
        textfont=dict(size=15, color=[COLORS.get(bases[i], OTHER) for i in idx]),
        customdata=[[i + 1, quals[i]] for i in idx], showlegend=False,
        hovertemplate="Base %{customdata[0]}: %{text}<br>Phred Q%{customdata[1]}<extra></extra>"),
        row=rt, col=1)
    bar_col = ["#2ca02c" if quals[i] >= 30 else "#ff7f0e" if quals[i] >= 20 else "#d62728" for i in idx]
    fig.add_trace(go.Bar(x=[ploc[i] for i in idx], y=[quals[i] for i in idx], marker_color=bar_col,
                         showlegend=False, hovertemplate="Q%{y}<extra></extra>"), row=rq, col=1)

    for row in (rt, rq):
        if p.get("trim"):
            s_, e_ = p["trim"]
            if s_ > start:
                fig.add_vrect(x0=x0, x1=edge(min(s_, end)), fillcolor="grey", opacity=0.25, line_width=0,
                              row=row, col=1)
            if e_ < end:
                fig.add_vrect(x0=edge(max(e_, start)), x1=x1, fillcolor="grey", opacity=0.25, line_width=0,
                              row=row, col=1)
        h = p["center"]
        if start < h <= end:
            fig.add_vrect(x0=edge(h - 1), x1=edge(h), fillcolor="gold", opacity=0.4, line_width=0,
                          row=row, col=1)

    ticks = [i for i in idx if (i + 1) % 10 == 0]
    fig.update_xaxes(range=[x0, x1], showgrid=False, showticklabels=False, row=rt, col=1)
    fig.update_xaxes(range=[x0, x1], showgrid=False, tickvals=[ploc[i] for i in ticks],
                     ticktext=[str(i + 1) for i in ticks], row=rq, col=1)
    fig.update_yaxes(range=[0, ymax * 1.2], showticklabels=False, showgrid=False, row=rt, col=1)
    fig.update_yaxes(range=[0, 60], title_text="Q", row=rq, col=1)


def make_stacked_figure(panels, window=40):
    """panels: list of {label, trace, center (1-based, in that trace), trim (start,end)|None}."""
    k = max(1, len(panels))
    fig = make_subplots(rows=2 * k, cols=1, shared_xaxes=False, row_heights=[0.8, 0.2] * k,
                        vertical_spacing=0.05,
                        subplot_titles=[x for p in panels for x in (p["label"], "")] or None)
    for i, p in enumerate(panels):
        _draw(fig, p, window, 2 * i + 1, 2 * i + 2, show_legend=(i == 0))
    fig.update_layout(height=380 * k + 30, margin=dict(l=10, r=10, t=40, b=10), bargap=0.15,
                      legend=dict(orientation="h", y=1.0 + 0.08 / k, x=0), hovermode="closest")
    fig.update_annotations(font_size=13, xanchor="left", x=0)
    return fig


def make_figure(trace, center=1, window=60, trim=None, highlight=None):
    """Single-read convenience wrapper."""
    return make_stacked_figure([{"label": "", "trace": trace, "center": highlight or center, "trim": trim}], window)


def make_consensus_map(seq, quals, flagged=(), selected=None, per_row=50):
    """Clickable consensus map. Point index + 1 == consensus position."""
    n = len(seq)
    rows = max(1, math.ceil(n / per_row))
    xs = [i % per_row for i in range(n)]
    ys = [-(i // per_row) for i in range(n)]
    fills = []
    for i, b in enumerate(seq):
        if b in "RYSWKMBDHVN":
            fills.append("#ffd54f")
        elif (i + 1) in flagged:
            fills.append("#ffab91")
        elif quals[i] < 20:
            fills.append("#d9d9d9")
        else:
            fills.append("#c8e6c9")
    fig = go.Figure(go.Scatter(
        x=xs, y=ys, mode="markers+text", text=list(seq),
        textfont=dict(size=11, color="#111", family="monospace"),
        marker=dict(symbol="square", size=15, color=fills),
        customdata=[[i + 1, quals[i]] for i in range(n)],
        hovertemplate="Position %{customdata[0]}: %{text}<br>Q%{customdata[1]}<extra></extra>"))
    if selected:
        i = selected - 1
        fig.add_shape(type="rect", x0=xs[i] - 0.5, x1=xs[i] + 0.5, y0=ys[i] - 0.45, y1=ys[i] + 0.45,
                      line=dict(color="#1565c0", width=3))
    fig.update_layout(
        height=min(60 + rows * 26, 640), margin=dict(l=50, r=10, t=5, b=5), dragmode=False,
        xaxis=dict(visible=False, range=[-1, per_row]), yaxis=dict(visible=False, range=[-rows + 0.4, 0.8]),
        annotations=[dict(x=-1.2, y=-r, text=str(r * per_row + 1), showarrow=False, xanchor="right",
                          font=dict(size=10, color="#888")) for r in range(rows)])
    return fig
