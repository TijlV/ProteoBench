"""
Plots of the immunopeptidomics identification modules.
"""

from __future__ import annotations

import json
import math
from typing import Dict, List, Optional

import pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots

from proteobench.plotting.plot_generator_base import PlotGeneratorBase

# ProteoBench software colours (plot_generator_entrapment.py); rescored workflows are provisional.
SOFTWARE_COLORS = {
    "Sage": "#671100",
    "Sage + MS2Rescore": "#44aa9a",
    "FragPipe": "#342288",
    "DIA-NN": "#999934",
    "Spectronaut": "#ff7f0e",
    "Custom": "#000000",
}
SOFTWARE_MARKERS = {"Sage": "triangle-down", "Sage + MS2Rescore": "pentagon", "FragPipe": "x", "DIA-NN": "star"}

# Ground-truth group colours.
# Subject to change if dataset is updated.
GROUP_COLORS = {
    "True": "#1baf7a",
    "Truncated": "#4a3aa7",
    "Co-present group": "#008300",
    "Background": "#eda100",
    "False": "#e87ba4",
}

# Precursor charge categories of the charge distribution.
CHARGE_CATEGORIES = ["1", "2", "3", "≥4"]

# Colours of the two results in the comparison plots.
COMPARE_COLORS = {"A": "#1f77b4", "B": "#d62728"}

# Validity of the FDP interval at the threshold, in the colours of the entrapment module.
CATEGORY_COLORS = {
    "valid": "#2ecc71",
    "inconclusive": "#f39c12",
    "invalid": "#e74c3c",
    "not estimable": "#8c8c8c",
}


class ImmunoPlotGenerator(PlotGeneratorBase):
    """
    Plot generator for the immunopeptidomics identification modules.
    """

    def plot_main_metric(
        self,
        result_df: pd.DataFrame,
        hide_annot: bool = False,
        fdp_bound: str = "upper",
        y_metric: str = "n_sample_peptides",
        highlight_color: str = "#d30067",
        label: str = "None",
        qvalue_threshold: Optional[float] = None,
        colorblind_mode: bool = False,
        mapping: Optional[Dict[str, int]] = None,
        **kwargs,
    ) -> go.Figure:
        """
        Plot the pool peptides identified against the estimated FDP, one trace per software tool.

        Workflows without an estimable upper bound (spectral library with r below the dataset's
        min_r) are drawn at their lower bound with an open marker.

        With ``qvalue_threshold`` the metrics are those at that q-value threshold
        (``results["by_qvalue"]``); workflows whose file does not reach the threshold are left out.

        Parameters
        ----------
        result_df : pd.DataFrame
            Datapoints, one row per submission.
        hide_annot : bool, optional
            Hide the line and label of the declared FDR.
        fdp_bound : str, optional
            "upper" or "lower" bound of the FDP on the x axis.
        y_metric : str, optional
            Datapoint column on the y axis; "n_sample_peptides" (pool peptides identified, mean over runs).
        highlight_color : str, optional
            Fill colour of highlighted (``Highlight`` column) datapoints.
        label : str, optional
            Column whose values are written next to the points; "None" for no labels.
        qvalue_threshold : float, optional
            Q-value threshold to show; the declared FDR (the headline metrics) when not given.
        colorblind_mode : bool, optional
            Give each software tool its own marker shape instead of circles for all.
        mapping : Dict[str, int], optional
            Marker size per ``old_new`` value; ``{"old": 10, "new": 20}`` as in the other modules.
            Highlighted datapoints are drawn at twice their size.
        **kwargs : dict
            Ignored; accepted for interface compatibility.

        Returns
        -------
        go.Figure
            The scatter plot.
        """
        if qvalue_threshold is not None:
            result_df = metrics_at_qvalue(result_df, qvalue_threshold)
        x_col = f"fdp_{fdp_bound}"
        y_labels = {"n_sample_peptides": "Pool peptides identified per run"}
        mapping = mapping or {"old": 10, "new": 20}
        fig = go.Figure()
        labels = result_df.apply(workflow_label, axis=1) if len(result_df) else pd.Series(dtype=str)
        for software, d in result_df.groupby(labels, sort=True):
            x = d[x_col].astype(float)
            missing = x.isna()
            x = x.where(~missing, d["fdp_lower"].astype(float))
            symbol = SOFTWARE_MARKERS.get(software, "circle") if colorblind_mode else "circle"
            color = SOFTWARE_COLORS.get(software, "#000000")
            highlight = d.get("Highlight", pd.Series(False, index=d.index)).fillna(False).astype(bool)
            sizes = [mapping.get(str(v), 10) for v in d.get("old_new", pd.Series("old", index=d.index))]
            show_labels = label != "None" and label in d.columns
            fig.add_trace(
                go.Scatter(
                    x=x * 100,
                    y=d[y_metric].astype(float),
                    mode="markers+text" if show_labels else "markers",
                    name=software,
                    marker=dict(
                        size=[s * 2 if h else s for s, h in zip(sizes, highlight)],
                        color=[highlight_color if h else color for h in highlight],
                        symbol=[symbol + "-open" if m else symbol for m in missing],
                        # Open markers are drawn by their outline alone, so they keep a visible one.
                        line=dict(width=[2 if m else 0.5 for m in missing], color="rgba(0,0,0,0.25)"),
                    ),
                    hovertext=[_hover(row) for _, row in d.iterrows()],
                    hoverinfo="text",
                    text=d[label].astype(str) if show_labels else None,
                    textposition="top center",
                )
            )
        if not hide_annot and qvalue_threshold is not None:
            fdr = qvalue_threshold * 100
            fig.add_vline(x=fdr, line_dash="dash", line_color="#8c8c8c", annotation_text=f"q-value threshold {fdr:g}%")
        elif not hide_annot and "declared_fdr" in result_df and len(result_df):
            fdr = float(result_df["declared_fdr"].iloc[0]) * 100
            fig.add_vline(x=fdr, line_dash="dash", line_color="#8c8c8c", annotation_text=f"declared FDR {fdr:g}%")
        fig.update_layout(
            width=None,
            height=700,
            autosize=True,
            xaxis=dict(
                title=f"Estimated FDP, {fdp_bound} bound (%)",
                rangemode="tozero",
                gridcolor="lightgray",
                gridwidth=1,
                linecolor="black",
            ),
            yaxis=dict(
                title=y_labels.get(y_metric, y_metric),
                gridcolor="lightgray",
                gridwidth=1,
                linecolor="black",
            ),
            legend_title="Software",
            margin=dict(l=80, r=20, t=50, b=80),
            clickmode="event+select",
        )
        return fig

    def plot_fdp_interval(
        self,
        result_df: pd.DataFrame,
        qvalue_threshold: Optional[float] = None,
        sort_ascending: bool = False,
        **kwargs,
    ) -> go.Figure:
        """
        Plot a forest / interval plot of the estimated FDP bounds per workflow.

        Laid out as the forest plot of the entrapment module. Each row is one workflow. A thick
        bar with end caps spans the lower to the upper bound of the FDP; open circles mark both
        ends. The bar colour is the validity category at the threshold. A workflow without an
        estimable upper bound shows its lower bound as a grey open marker. Without a threshold,
        a diamond marks the declared FDR of each workflow (a star for the uploaded result). The
        bars on the right give the pool peptides identified, which also sort the rows.

        Parameters
        ----------
        result_df : pd.DataFrame
            Datapoints, one row per submission.
        qvalue_threshold : float, optional
            Q-value threshold to show; the declared FDR (the headline metrics) when not given.
            Workflows whose file does not reach the threshold are left out.
        sort_ascending : bool, optional
            Sort the rows by pool peptides identified, ascending (True) or descending (False).
        **kwargs : dict
            Ignored; accepted for interface compatibility.

        Returns
        -------
        go.Figure
            The forest plot.
        """
        df = metrics_at_qvalue(result_df, qvalue_threshold) if qvalue_threshold is not None else result_df.copy()
        for col in ("fdp_lower", "fdp_upper", "n_sample_peptides", "declared_fdr"):
            df[col] = pd.to_numeric(df[col], errors="coerce") if col in df else math.nan
        df = df.sort_values("n_sample_peptides", ascending=sort_ascending, na_position="last").reset_index(drop=True)
        n = len(df)
        y_pos = list(range(n))
        is_new = (df["old_new"] == "new").tolist() if "old_new" in df else [False] * n
        hover = [_hover(row) for _, row in df.iterrows()]
        cap_half = 0.15

        fig = go.Figure()
        for category, color in CATEGORY_COLORS.items():
            d = df[df["category"].astype(str) == category]
            if d.empty:
                continue
            bar_x, bar_y, cap_x, cap_y, ep_x, ep_y, ep_hover = [], [], [], [], [], [], []
            open_x, open_y, open_hover = [], [], []
            for i, row in d.iterrows():
                lo, hi = row["fdp_lower"] * 100, row["fdp_upper"] * 100
                if pd.isna(lo):
                    continue
                if pd.isna(hi):
                    open_x.append(lo)
                    open_y.append(i)
                    open_hover.append(hover[i])
                    continue
                bar_x += [lo, hi, None]
                bar_y += [i, i, None]
                for x in (lo, hi):
                    cap_x += [x, x, None]
                    cap_y += [i - cap_half, i + cap_half, None]
                ep_x += [lo, hi]
                ep_y += [i, i]
                ep_hover += [hover[i], hover[i]]
            name = category.capitalize()
            if bar_x:
                fig.add_trace(
                    go.Scatter(
                        x=bar_x,
                        y=bar_y,
                        mode="lines",
                        name=name,
                        line=dict(color=color, width=10),
                        hoverinfo="skip",
                        legendgroup=category,
                    )
                )
                fig.add_trace(
                    go.Scatter(
                        x=cap_x,
                        y=cap_y,
                        mode="lines",
                        showlegend=False,
                        line=dict(color=color, width=2),
                        hoverinfo="skip",
                        legendgroup=category,
                    )
                )
                fig.add_trace(
                    go.Scatter(
                        x=ep_x,
                        y=ep_y,
                        mode="markers",
                        showlegend=False,
                        marker=dict(symbol="circle", color="white", size=10, line=dict(color=color, width=2.5)),
                        hovertext=ep_hover,
                        hoverinfo="text",
                        legendgroup=category,
                    )
                )
            if open_x:
                fig.add_trace(
                    go.Scatter(
                        x=open_x,
                        y=open_y,
                        mode="markers",
                        name=f"{name} (lower bound only)" if bar_x else name,
                        marker=dict(symbol="circle-open", size=12, color=color, line=dict(width=2.5)),
                        hovertext=open_hover,
                        hoverinfo="text",
                        legendgroup=category,
                    )
                )

        if qvalue_threshold is None:
            fdr = df["declared_fdr"].fillna(0.01) * 100
            fig.add_trace(
                go.Scatter(
                    x=fdr,
                    y=y_pos,
                    mode="markers",
                    name="Declared FDR threshold",
                    marker=dict(
                        symbol=["star" if new else "diamond" for new in is_new],
                        color="rgba(20,20,20,0.9)",
                        size=[14 if new else 11 for new in is_new],
                        line=dict(width=[2.5 if new else 1.5 for new in is_new], color="black"),
                    ),
                    hoverinfo="skip",
                )
            )

        fig.add_trace(
            go.Bar(
                x=df["n_sample_peptides"],
                y=y_pos,
                orientation="h",
                name="Pool peptides identified",
                marker=dict(color=[CATEGORY_COLORS.get(str(c), "#cccccc") for c in df["category"]], opacity=0.35),
                xaxis="x2",
                hovertemplate="%{x:,.1f} pool peptides<extra></extra>",
                showlegend=False,
            )
        )

        fig.update_layout(
            xaxis=dict(
                title="Estimated FDP, lower to upper bound (%)",
                gridcolor="lightgray",
                linecolor="black",
                rangemode="tozero",
                domain=[0, 0.55],
            ),
            xaxis2=dict(
                title="Pool peptides identified per run",
                side="top",
                domain=[0.60, 1.0],
                gridcolor="lightgray",
                linecolor="black",
                tickformat=",d",
            ),
            yaxis=dict(
                tickmode="array",
                tickvals=y_pos,
                ticktext=[f"<b>{i}</b>" if new else str(i) for i, new in zip(df["id"], is_new)],
                autorange="reversed",
                gridcolor="lightgray",
                linecolor="black",
            ),
            template="plotly_white",
            height=max(400, n * 30 + 140),
            margin=dict(l=180, r=20, t=80, b=80),
        )
        if qvalue_threshold is not None:
            fig.add_vline(
                x=qvalue_threshold * 100,
                line_dash="dash",
                line_color="gray",
                annotation_text=f"Q ≤ {qvalue_threshold:g}",
                annotation_position="top right",
            )
        return fig

    def generate_in_depth_plots(self, performance_data: pd.DataFrame, **kwargs) -> Dict[str, go.Figure]:
        """
        Generate the plots of a single result from its intermediate table.

        Parameters
        ----------
        performance_data : pd.DataFrame
            Intermediate table (one row per run and peptide, with ``run``, ``peptide``, ``group``, ``n_psms``).
        **kwargs : dict
            ``pool_lengths`` (pd.Series, optional): lengths of the designed pool peptides, drawn as grey bars
            in the length distribution. Other keys are ignored (e.g. parse_settings).

        Returns
        -------
        Dict[str, go.Figure]
            ``composition`` (peptides per ground-truth group and run), ``length_charge`` (length and charge
            distributions of the identified pool peptides) and ``false_peptides`` (false peptides by number
            of runs and PSMs; recurring ones suggest carryover).
        """
        counts = performance_data.groupby(["run", "group"]).size().unstack(fill_value=0)
        composition = go.Figure()
        for group, color in GROUP_COLORS.items():
            if group in counts:
                composition.add_trace(
                    go.Bar(x=counts.index.astype(str), y=counts[group], name=group, marker_color=color)
                )
        composition.update_layout(
            barmode="stack",
            xaxis=dict(title="Run", linecolor="black"),
            yaxis=dict(title="Unique peptides", gridcolor="lightgray", gridwidth=1, linecolor="black"),
        )

        false = performance_data[performance_data["group"] == "False"]
        support = (
            false.groupby(["peptide", "class"])
            .agg(runs=("run", "nunique"), psms=("n_psms", "sum"))
            .reset_index()
            .sort_values(["runs", "psms"], ascending=False)
            .head(30)
        )
        false_fig = go.Figure(
            go.Bar(
                x=support["peptide"],
                y=support["psms"],
                marker_color=GROUP_COLORS["False"],
                text=support["runs"].map(lambda n: f"{n} run(s)"),
                hovertext=support["class"],
            )
        )
        false_fig.update_layout(
            xaxis=dict(title="False peptide (most recurrent first)", linecolor="black"),
            yaxis=dict(title="PSMs over all runs", gridcolor="lightgray", gridwidth=1, linecolor="black"),
        )
        series = [
            {"name": "Pool peptides identified", "data": pool_only(performance_data), "color": COMPARE_COLORS["A"]}
        ]
        length_charge = length_charge_figure(series, pool_lengths=kwargs.get("pool_lengths"))
        return {"composition": composition, "length_charge": length_charge, "false_peptides": false_fig}

    def plot_length_charge_comparison(
        self, results: Dict[str, pd.DataFrame], pool_lengths: Optional[pd.Series] = None
    ) -> go.Figure:
        """
        Overlay the length and charge distributions of the identified pool peptides of several results.

        Parameters
        ----------
        results : Dict[str, pd.DataFrame]
            Label -> intermediate table, in legend order.
        pool_lengths : pd.Series, optional
            Lengths of the designed pool peptides, drawn as grey bars.

        Returns
        -------
        go.Figure
            Two panels: length of the pool peptides and charge of their precursors.
        """
        colors = list(COMPARE_COLORS.values())
        series = [
            {"name": label, "data": pool_only(data), "color": colors[i % len(colors)]}
            for i, (label, data) in enumerate(results.items())
        ]
        return length_charge_figure(series, pool_lengths=pool_lengths)

    def get_in_depth_plot_layout(self) -> list:
        """
        Return the layout of the in-depth plots.

        Returns
        -------
        list
            Plot configurations.
        """
        return [
            {"plots": ["composition"], "columns": 1, "title": "Ground-truth composition per run"},
            {"plots": ["length_charge"], "columns": 1, "title": "Length and charge distributions"},
            {"plots": ["false_peptides"], "columns": 1, "title": "False identifications"},
        ]

    def get_in_depth_plot_descriptions(self) -> Dict[str, str]:
        """
        Return the descriptions of the in-depth plots.

        Returns
        -------
        Dict[str, str]
            Plot name -> description.
        """
        return {
            "composition": "Unique accepted peptides per run, by ground-truth group.",
            "length_charge": LENGTH_CHARGE_DESCRIPTION,
            "false_peptides": (
                "False peptides, most recurrent first. A false peptide that is found in many runs with many PSMs is "
                "more consistent with carryover than with a random match."
            ),
        }

    def get_metrics_help_markdown(self) -> Optional[str]:
        """
        Explain how the metrics of the main plot are calculated.

        Returns
        -------
        Optional[str]
            Markdown text.
        """
        return (
            "**Pool peptides identified** is the number of unique designed pool peptides accepted in a run "
            "(mean over runs).\n\n"
            "**FDP** uses the combined entrapment estimator (Wen et al., Nat. Methods 2025): lower bound "
            "N_false / N_all, upper bound N_false × (1 + 1/r) / N_all, with r the entrapment-to-sample ratio "
            "of the search space (from the library for spectral-library searches). Truncated sample peptides, "
            "co-present library groups and background (host, construct, contaminants) are not false."
        )


LENGTH_CHARGE_DESCRIPTION = (
    "Number of identified pool peptides per length (left) and of their precursors (peptide and charge) per "
    "charge state (right), per run and averaged over the runs. The grey bars give the number of designed pool "
    "peptides of each length, the most that can be identified."
)


def pool_only(data: pd.DataFrame) -> pd.DataFrame:
    """
    Rows of the identified designed pool peptides (class ``sample``).

    Parameters
    ----------
    data : pd.DataFrame
        Intermediate table with ``class``.

    Returns
    -------
    pd.DataFrame
        The pool peptide rows.
    """
    return data[data["class"] == "sample"]


def length_distribution(data: pd.DataFrame) -> pd.Series:
    """
    Number of unique peptides per length, per run and averaged over the runs.

    Parameters
    ----------
    data : pd.DataFrame
        Intermediate table with ``run`` and ``peptide``.

    Returns
    -------
    pd.Series
        Length -> mean count per run; empty for an empty table.
    """
    unique = data.drop_duplicates(["run", "peptide"])
    if unique.empty:
        return pd.Series(dtype=float)
    return pd.crosstab(unique["run"], unique["peptide"].str.len()).mean()


def charge_distribution(data: pd.DataFrame) -> pd.Series:
    """
    Number of unique precursors per charge state (``CHARGE_CATEGORIES``), per run and averaged over the runs.

    Parameters
    ----------
    data : pd.DataFrame
        Intermediate table with ``run``, ``peptide`` and ``charges`` ("2" or "2;3").

    Returns
    -------
    pd.Series
        Charge category -> mean count per run; empty when no charge is reported.
    """
    if "charges" not in data.columns:
        return pd.Series(dtype=float)
    charges = data["charges"].map(lambda c: "" if pd.isna(c) else str(c).removesuffix(".0"))
    precursors = data.assign(charge=charges.str.split(";")).explode("charge")
    precursors = precursors[precursors["charge"].str.strip() != ""].drop_duplicates(["run", "peptide", "charge"])
    if precursors.empty:
        return pd.Series(dtype=float)
    category = precursors["charge"].astype(float).astype(int).map(lambda z: "≥4" if z >= 4 else str(z))
    counts = pd.crosstab(precursors["run"], category)
    return counts.reindex(columns=CHARGE_CATEGORIES, fill_value=0).mean()


def length_charge_figure(series: List[dict], pool_lengths: Optional[pd.Series] = None) -> go.Figure:
    """
    Length and charge distributions of one or more sets of accepted peptides, side by side.

    Parameters
    ----------
    series : List[dict]
        One entry per set of bars, with ``name``, ``data`` (intermediate table) and ``color``.
    pool_lengths : pd.Series, optional
        Lengths of the designed pool peptides, drawn as grey bars.

    Returns
    -------
    go.Figure
        The figure.
    """
    lengths = [length_distribution(s["data"]) for s in series]
    observed = {int(k) for d in lengths for k in d.index}
    if pool_lengths is not None and len(pool_lengths):
        observed |= {int(k) for k in pool_lengths}
    axis = [str(k) for k in range(min(observed), max(observed) + 1)] if observed else []

    fig = make_subplots(
        rows=1,
        cols=2,
        horizontal_spacing=0.1,
        subplot_titles=["Length of identified pool peptides", "Charge of their precursors"],
    )
    if pool_lengths is not None and len(pool_lengths):
        pool = pool_lengths.astype(int).value_counts()
        fig.add_trace(
            go.Bar(
                x=axis,
                y=[pool.get(int(k), 0) for k in axis],
                name=f"Designed pool ({len(pool_lengths):,} peptides)",
                marker_color="lightgray",
                hovertemplate="Designed pool, %{x} residues: %{y:,.0f}<extra></extra>",
            ),
            row=1,
            col=1,
        )
    any_charge = False
    for s, length in zip(series, lengths):
        n = s["data"].drop_duplicates(["run", "peptide"]).groupby("run").size().mean()
        label = f"{s['name']} ({n:,.0f} per run)" if not math.isnan(n) else s["name"]
        fig.add_trace(
            go.Bar(
                x=axis,
                y=[length.get(int(k), 0) for k in axis],
                name=label,
                marker_color=s["color"],
                legendgroup=s["name"],
                hovertemplate=f"{s['name']}<br>%{{x}} residues: %{{y:,.1f}} per run<extra></extra>",
            ),
            row=1,
            col=1,
        )
        charge = charge_distribution(s["data"])
        if charge.empty:
            continue
        any_charge = True
        fig.add_trace(
            go.Bar(
                x=CHARGE_CATEGORIES,
                y=charge.to_numpy(),
                name=label,
                marker_color=s["color"],
                legendgroup=s["name"],
                showlegend=False,
                hovertemplate=f"{s['name']}<br>charge %{{x}}: %{{y:,.1f}} per run<extra></extra>",
            ),
            row=1,
            col=2,
        )
    if not any_charge:
        fig.add_annotation(
            text="No precursor charge in this result",
            xref="x2 domain",
            yref="y2 domain",
            x=0.5,
            y=0.5,
            showarrow=False,
            font=dict(color="#8c8c8c"),
        )
    fig.update_xaxes(title_text="Peptide length (residues)", type="category", linecolor="black", row=1, col=1)
    fig.update_xaxes(
        title_text="Precursor charge",
        type="category",
        categoryorder="array",
        categoryarray=CHARGE_CATEGORIES,
        linecolor="black",
        row=1,
        col=2,
    )
    fig.update_yaxes(gridcolor="lightgray", linecolor="black", rangemode="tozero")
    fig.update_yaxes(title_text="Pool peptides per run", row=1, col=1)
    fig.update_yaxes(title_text="Precursors per run", row=1, col=2)
    fig.update_layout(
        template="plotly_white",
        height=500,
        margin=dict(l=60, r=20, t=50, b=130),
        legend=dict(orientation="h", yanchor="top", y=-0.2, x=0),
        barmode="group",
        bargap=0.25,
        bargroupgap=0.05,
    )
    return fig


def _pct(value) -> str:
    """
    Format a fraction as a percentage.

    Parameters
    ----------
    value : float or None
        The fraction.

    Returns
    -------
    str
        The percentage, or "not estimable" for a missing value.
    """
    return (
        "not estimable" if value is None or (isinstance(value, float) and math.isnan(value)) else f"{float(value):.2%}"
    )


def _count(value) -> str:
    """
    Format a peptide count.

    Parameters
    ----------
    value : float or None
        The count, possibly a mean over runs.

    Returns
    -------
    str
        The count with thousands separators, one decimal when it is not whole; "n/a" when missing.
    """
    try:
        value = float(value)
    except (TypeError, ValueError):
        return "n/a"
    if math.isnan(value):
        return "n/a"
    return f"{value:,.0f}" if value.is_integer() else f"{value:,.1f}"


def metrics_at_qvalue(result_df: pd.DataFrame, threshold: float) -> pd.DataFrame:
    """
    Replace the headline metrics of every datapoint with those at a q-value threshold.

    A datapoint without per-threshold metrics (no ``results["by_qvalue"]``) is shown at its declared FDR only.

    Parameters
    ----------
    result_df : pd.DataFrame
        Datapoints, one row per submission.
    threshold : float
        The q-value threshold.

    Returns
    -------
    pd.DataFrame
        The datapoints that have metrics at this threshold, with ``n_sample_peptides``, ``recall``,
        ``recall_all``, ``fdp_lower``, ``fdp_upper`` and ``category`` taken at it.
    """
    rows = []
    for _, row in result_df.iterrows():
        results = row.get("results")
        if isinstance(results, str):
            results = json.loads(results)
        at = (results or {}).get("by_qvalue", {}).get(f"{threshold:g}")
        if at is None and math.isclose(float(row.get("declared_fdr", math.nan)), threshold):
            rows.append(row)  # no per-threshold metrics: its headline metrics are those at the declared FDR
        elif at is not None:
            row = row.copy()
            row["n_sample_peptides"] = at.get("n_sample", math.nan)
            for key in ("recall", "recall_all", "fdp_lower", "fdp_upper", "category"):
                row[key] = at.get(key, math.nan)
            rows.append(row)
    return pd.DataFrame(rows, columns=result_df.columns).reset_index(drop=True)


def _filled(value) -> bool:
    """
    Whether a metadata value is set (not None, NaN or blank).

    Parameters
    ----------
    value : Any
        The value.

    Returns
    -------
    bool
        True if the value is set.
    """
    return value is not None and not (isinstance(value, float) and math.isnan(value)) and str(value).strip() != ""


def workflow_label(row: pd.Series) -> str:
    """
    Name of the workflow of one datapoint: the software, followed by the rescoring tool if any.

    Parameters
    ----------
    row : pd.Series
        The datapoint.

    Returns
    -------
    str
        For example "Sage" or "Sage + MS2Rescore".
    """
    rescoring = row.get("rescoring")
    if isinstance(rescoring, str) and rescoring.strip():
        return f"{row.get('software_name')} + {rescoring.strip()}"
    return str(row.get("software_name"))


def _hover(row: pd.Series) -> str:
    """
    Build the hover text of one datapoint.

    Parameters
    ----------
    row : pd.Series
        The datapoint.

    Returns
    -------
    str
        HTML hover text.
    """
    pct = _pct
    rescoring = ""
    if isinstance(row.get("rescoring"), str) and row.get("rescoring").strip():
        details = ", ".join(str(row.get(k)) for k in ("rescoring_version", "rescoring_model") if _filled(row.get(k)))
        rescoring = f"<br>Rescoring: {row.get('rescoring')}" + (f" ({details})" if details else "")

    return (
        f"ProteoBench ID: {row.get('id')}<br>Software: {row.get('software_name')} {row.get('software_version', '')}"
        f"{rescoring}"
        f"<br>Pool peptides: {_count(row.get('n_sample_peptides'))}<br>FDP lower: {pct(row.get('fdp_lower'))}"
        f"<br>FDP upper: {pct(row.get('fdp_upper'))}<br>r: {float(row.get('r', math.nan)):.3g} ({row.get('r_source')})"
        f"<br>Runs: {row.get('n_runs')}<br>Category: {row.get('category')}"
    )
