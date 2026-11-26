"""
Chromatogram viewer for AB1 files (interactive Plotly figure).
Shows the 4 trace channels, base calls, a quality bar track, the regions that
were trimmed away, and an optional highlighted base.
"""

from Bio import SeqIO
import plotly.graph_objects as go
from plotly.subplots import make_subplots

COLORS = {"A": "#2ca02c", "C": "#1f77b4", "G": "#7f7f7f", "T": "#d62728", "N": "#9467bd"}


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


def make_figure(trace: dict, center: int = 1, window: int = 60, trim=None, highlight=None):
    """
    center / highlight are 1-based base positions in the original read.
    trim = (start, end) 0-based half-open range that was KEPT; the rest is greyed.
    """
    ploc, bases, quals = trace["ploc"], trace["bases"], trace["quals"]
    n = len(ploc)
    if n == 0:
        return go.Figure()
    window = max(10, min(window, n))
    start = max(0, center - 1 - window // 2)
    end = min(n, start + window)
    start = max(0, end - window)

    def edge(i):  # scan position of the boundary before base i
        if i <= 0:
            return ploc[0] - 6
        if i >= n:
            return ploc[-1] + 6
        return (ploc[i - 1] + ploc[i]) / 2

    x0, x1 = int(edge(start)), int(edge(end))
    fig = make_subplots(rows=2, cols=1, shared_xaxes=True, row_heights=[0.8, 0.2],
                        vertical_spacing=0.03)

    ymax = 1
    for base, data in trace["channels"].items():
        seg = data[max(0, x0):x1 + 1]
        ymax = max(ymax, max(seg) if seg else 1)
        fig.add_trace(go.Scatter(x=list(range(max(0, x0), max(0, x0) + len(seg))), y=seg,
                                 mode="lines", name=base, line=dict(color=COLORS.get(base, "#999"), width=1.6),
                                 hoverinfo="skip"), row=1, col=1)

    idx = list(range(start, end))
    fig.add_trace(go.Scatter(
        x=[ploc[i] for i in idx], y=[ymax * 1.08] * len(idx), mode="text",
        text=[bases[i] for i in idx], textfont=dict(size=15, color=[COLORS.get(bases[i], "#999") for i in idx]),
        customdata=[[i + 1, quals[i]] for i in idx], showlegend=False,
        hovertemplate="Base %{customdata[0]}: %{text}<br>Phred Q%{customdata[1]}<extra></extra>"),
        row=1, col=1)

    bar_col = ["#2ca02c" if quals[i] >= 30 else "#ff7f0e" if quals[i] >= 20 else "#d62728" for i in idx]
    fig.add_trace(go.Bar(x=[ploc[i] for i in idx], y=[quals[i] for i in idx], marker_color=bar_col,
                         showlegend=False, hovertemplate="Q%{y}<extra></extra>"), row=2, col=1)

    if trim:
        s, e = trim
        if s > start:
            fig.add_vrect(x0=x0, x1=edge(min(s, end)), fillcolor="grey", opacity=0.25, line_width=0)
        if e < end:
            fig.add_vrect(x0=edge(max(e, start)), x1=x1, fillcolor="grey", opacity=0.25, line_width=0)
    if highlight and start < highlight <= end:
        fig.add_vrect(x0=edge(highlight - 1), x1=edge(highlight), fillcolor="gold", opacity=0.35, line_width=0)

    ticks = [i for i in idx if (i + 1) % 10 == 0]
    fig.update_xaxes(tickvals=[ploc[i] for i in ticks], ticktext=[str(i + 1) for i in ticks],
                     range=[x0, x1], showgrid=False, row=2, col=1)
    fig.update_xaxes(range=[x0, x1], showgrid=False, row=1, col=1)
    fig.update_yaxes(range=[0, ymax * 1.2], showticklabels=False, showgrid=False, row=1, col=1)
    fig.update_yaxes(range=[0, 60], title_text="Q", row=2, col=1)
    fig.update_layout(height=430, margin=dict(l=10, r=10, t=10, b=10), bargap=0.15,
                      legend=dict(orientation="h", y=1.08, x=0), hovermode="closest")
    return fig
