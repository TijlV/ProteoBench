"""
Scores for the immunopeptidomics identification modules.

Accepted peptides are classified against the ground truth of an :class:`ImmunoDataset`; recall and
the false discovery proportion (FDP) follow from the class counts per run. The FDP uses the
combined entrapment estimator of Wen et al., Nat. Methods 2025:

    FDP_lower = N_false / N_all
    FDP_upper = N_false * (1 + 1 / r) / N_all

with r the entrapment-to-sample ratio of the search space.
"""

from __future__ import annotations

import math
from typing import Any, Dict, Optional

import pandas as pd

from proteobench.modules.immuno.dataset import FALSE_CLASSES, SAMPLE, ImmunoDataset
from proteobench.score.entrapmentscores import EntrapmentScores
from proteobench.score.score_base import ScoreBase

#: Q-value thresholds at which every datapoint is also scored, for the threshold selector of the main plot
#: (the same set as the entrapment module).
QVALUE_THRESHOLDS = (0.001, 0.01, 0.05, 0.1, 1.0)
#: A file covers a threshold when its largest q-value reaches this fraction of it: a tool that exported at q <= 0.05
#: has a largest q-value just below 0.05, and an unfiltered target-decoy export rarely reaches q = 1.
QVALUE_COVERAGE_TOLERANCE = 0.9
#: Summary metrics stored per q-value threshold.
THRESHOLD_METRICS = ("n_all", "n_false", "n_sample", "recall", "recall_all", "fdp_lower", "fdp_upper", "category")


class ImmunoScores(ScoreBase):
    """
    Classify accepted peptides against a ground truth and compute recall and FDP.

    Parameters
    ----------
    dataset : ImmunoDataset
        Ground-truth definition of the benchmark sample.
    declared_fdr : float
        The FDR threshold the workflows are filtered at.
    """

    def __init__(self, dataset: ImmunoDataset, declared_fdr: float = 0.01):
        """
        Initialize the scorer.

        Parameters
        ----------
        dataset : ImmunoDataset
            Ground-truth definition of the benchmark sample.
        declared_fdr : float
            The FDR threshold the workflows are filtered at.
        """
        self.dataset = dataset
        self.declared_fdr = declared_fdr

    def generate_intermediate(self, filtered_df: pd.DataFrame, replicate_to_raw: Optional[dict] = None) -> pd.DataFrame:
        """
        Classify the unique peptides of every run.

        Parameters
        ----------
        filtered_df : pd.DataFrame
            Accepted identifications with the columns ``Raw file`` (run label), ``Sequence`` and ``Q-Value``.
        replicate_to_raw : dict, optional
            Unused; accepted for interface compatibility with ScoreBase.

        Returns
        -------
        pd.DataFrame
            One row per run and peptide with the columns ``run``, ``peptide``, ``class``, ``group``,
            ``n_psms``, ``best_qvalue`` and ``charges`` (the precursor charges of the peptide in that run,
            sorted and separated by ";"; empty when the tool reports no charge), sorted by run and peptide.
        """
        df = filtered_df.assign(peptide=filtered_df["Sequence"].map(self.dataset.normalise))
        out = (
            df.groupby(["Raw file", "peptide"], sort=True)
            .agg(n_psms=("Q-Value", "size"), best_qvalue=("Q-Value", "min"))
            .reset_index()
            .rename(columns={"Raw file": "run"})
        )
        classes = {p: self.dataset.classify(p) for p in out["peptide"].unique()}
        group_of = {c: g for g, members in self.dataset.display_groups.items() for c in members}
        out["class"] = out["peptide"].map(classes)
        out["group"] = out["class"].map(group_of)
        out["charges"] = ""
        if "Charge" in df.columns:
            charge = pd.to_numeric(df["Charge"], errors="coerce")
            with_charge = df.assign(Charge=charge)[charge.notna()]
            if not with_charge.empty:
                charges = with_charge.groupby(["Raw file", "peptide"])["Charge"].agg(
                    lambda z: ";".join(str(int(c)) for c in sorted(set(z)))
                )
                keys = pd.MultiIndex.from_frame(out[["run", "peptide"]])
                out["charges"] = charges.reindex(keys).fillna("").to_numpy()
        return out[["run", "peptide", "class", "group", "n_psms", "best_qvalue", "charges"]]

    def _fdp(self, n_false: int, n_all: int, r: float) -> tuple[float, float]:
        """
        Compute the lower and upper bound of the FDP.

        Parameters
        ----------
        n_false : int
            Number of false peptides.
        n_all : int
            Number of accepted peptides.
        r : float
            Entrapment-to-sample ratio.

        Returns
        -------
        tuple[float, float]
            Lower and upper bound; the upper bound is NaN when r is below the dataset's min_r.
        """
        if n_all == 0:
            return math.nan, math.nan
        lower = n_false / n_all
        if not r or math.isnan(r) or r < self.dataset.min_r:
            return lower, math.nan
        return lower, n_false * (1 + 1 / r) / n_all

    def _category(self, lower: float, upper: float, threshold: Optional[float] = None) -> str:
        """
        Categorise the FDP bounds at a q-value threshold, as in the entrapment module.

        Parameters
        ----------
        lower : float
            Lower bound of the FDP.
        upper : float
            Upper bound of the FDP, NaN when not estimable.
        threshold : float, optional
            The threshold the bounds are compared with; the declared FDR when not given.

        Returns
        -------
        str
            "valid", "inconclusive", "invalid", "not estimable" or "no identifications".
        """
        threshold = self.declared_fdr if threshold is None else threshold
        if math.isnan(lower):
            return "no identifications"
        if math.isnan(upper):
            return "invalid" if lower > threshold else "not estimable"
        return EntrapmentScores.categorise_metric(lower, upper, threshold)

    def metrics_by_qvalue(
        self, targets: pd.DataFrame, r: float, max_qvalue: Optional[float] = None
    ) -> Dict[str, Dict[str, Any]]:
        """
        Compute the summary metrics at every q-value threshold the submitted file covers.

        A threshold is kept when the file reaches it (``max_qvalue``, within ``QVALUE_COVERAGE_TOLERANCE``),
        so the file holds every identification accepted at that threshold. The declared FDR is always kept: it is the benchmark rule,
        and a tool that already filtered at it exports q-values just below it.
        Thresholds without any accepted identification are left out.

        Parameters
        ----------
        targets : pd.DataFrame
            Rank-1 target identifications with the columns ``Raw file``, ``Sequence`` and ``Q-Value``.
        r : float
            Entrapment-to-sample ratio of the search space.
        max_qvalue : float, optional
            Q-value up to which the file is complete; the largest target q-value when not given.

        Returns
        -------
        Dict[str, Dict[str, Any]]
            Threshold (as a string, e.g. "0.01") -> the metrics in ``THRESHOLD_METRICS``, with the
            category taken at that threshold.
        """
        max_q = float(targets["Q-Value"].max()) if max_qvalue is None else max_qvalue
        covered = {t for t in QVALUE_THRESHOLDS if max_q >= QVALUE_COVERAGE_TOLERANCE * t}
        thresholds = sorted(covered | {self.declared_fdr})
        intermediate = self.generate_intermediate(targets[targets["Q-Value"] <= thresholds[-1]])
        out = {}
        for t in thresholds:
            # A peptide is accepted at t when its best q-value in that run is at most t.
            accepted = intermediate[intermediate["best_qvalue"] <= t]
            if accepted.empty:
                continue  # e.g. 0.001 for target-decoy q-values, which rarely go that low
            summary = self.calculate_metrics(accepted, r, threshold=t)["summary"]
            out[f"{t:g}"] = {k: summary.get(k, math.nan) for k in THRESHOLD_METRICS}
        return out

    def calculate_metrics(
        self, intermediate: pd.DataFrame, r: float, threshold: Optional[float] = None
    ) -> Dict[str, Any]:
        """
        Compute the per-run metrics and their mean over runs.

        Parameters
        ----------
        intermediate : pd.DataFrame
            Output of :meth:`generate_intermediate`.
        r : float
            Entrapment-to-sample ratio of the search space (FASTA or spectral library).
        threshold : float, optional
            Q-value threshold the category is taken at; the declared FDR when not given.

        Returns
        -------
        Dict[str, Any]
            ``summary`` (mean over runs, plus reproducibility over runs) and ``per_run`` (run label ->
            class counts, recall and FDP bounds).
        """
        ds = self.dataset
        in_scope, n_sample = ds.sample_peptides_in_scope, len(ds.sample_peptides)
        per_run = {}
        for run, d in intermediate.groupby("run", sort=True):
            counts = d["class"].value_counts().reindex(ds.classes, fill_value=0)
            sample_found = set(d.loc[d["class"] == SAMPLE, "peptide"])
            n_all, n_false = int(counts.sum()), int(counts[FALSE_CLASSES].sum())
            lower, upper = self._fdp(n_false, n_all, r)
            per_run[str(run)] = {
                **{f"n_{c}": int(counts[c]) for c in ds.classes},
                **{f"n_group_{g}": int(counts[m].sum()) for g, m in ds.display_groups.items()},
                "n_all": n_all,
                "n_false": n_false,
                "recall": len(sample_found & in_scope) / len(in_scope) if in_scope else math.nan,
                "recall_all": len(sample_found) / n_sample if n_sample else math.nan,
                "fdp_lower": lower,
                "fdp_upper": upper,
            }

        table = pd.DataFrame(per_run).T.astype(float)
        summary = table.mean().to_dict() if len(table) else {}
        summary["fdp_upper"] = float(table["fdp_upper"].mean(skipna=False)) if len(table) else math.nan
        sets = [set(d.loc[d["class"] == SAMPLE, "peptide"]) for _, d in intermediate.groupby("run")]
        union = set.union(*sets) if sets else set()
        summary.update(
            {
                "n_runs": len(per_run),
                "n_sample_min": float(table[f"n_{SAMPLE}"].min()) if len(table) else math.nan,
                "n_sample_max": float(table[f"n_{SAMPLE}"].max()) if len(table) else math.nan,
                "n_sample_all_runs": len(set.intersection(*sets)) if sets else 0,
                "jaccard_sample": len(set.intersection(*sets)) / len(union) if union else math.nan,
                "r": r,
                "declared_fdr": self.declared_fdr,
                "n_sample_designed": n_sample,
                "n_sample_in_scope": len(in_scope),
            }
        )
        summary["category"] = self._category(summary.get("fdp_lower", math.nan), summary["fdp_upper"], threshold)
        return {"summary": summary, "per_run": per_run}
