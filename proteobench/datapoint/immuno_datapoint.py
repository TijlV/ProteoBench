"""
Datapoint of the immunopeptidomics identification modules.
"""

from __future__ import annotations

import dataclasses
import hashlib
import logging
import math
from dataclasses import dataclass
from datetime import datetime
from typing import Any, Dict

import pandas as pd

import proteobench
from proteobench.datapoint.datapoint_base import DatapointBase
from proteobench.score.immunoscores import ImmunoScores

# Input formats that rescore the PSMs of a search engine, with the search engine whose output the
# module's parse settings expect. Their datapoints are named after the search engine, with the tool
# in ``rescoring`` (see ``name_rescored_workflow``).
RESCORING_TOOLS = {"MS2Rescore": "Sage"}


def name_rescored_workflow(fields: Dict[str, Any], input_format: str = "") -> Dict[str, Any]:
    """
    Name a rescored workflow after its search engine.

    Applies when ``input_format`` or ``fields["rescoring"]`` is a key of ``RESCORING_TOOLS``; other
    fields are returned unchanged.

    Parameters
    ----------
    fields : Dict[str, Any]
        Datapoint fields.
    input_format : str, optional
        The input format of the submission.

    Returns
    -------
    Dict[str, Any]
        A copy of ``fields`` with ``rescoring`` set to the tool, ``search_engine`` defaulting to the
        search engine in ``RESCORING_TOOLS``, ``software_name`` and ``software_version`` equal to the
        search engine and its version, and ``rescoring_version`` defaulting to the given
        ``software_version`` (the version of the tool that wrote the file).
    """
    tool = input_format if input_format in RESCORING_TOOLS else fields.get("rescoring")
    if tool not in RESCORING_TOOLS:
        return dict(fields)
    out = dict(fields)
    given_version = str(out.get("software_version") or "").strip()
    if given_version and given_version != str(out.get("search_engine_version") or "").strip():
        out["rescoring_version"] = out.get("rescoring_version") or given_version
    out["rescoring"] = tool
    out["search_engine"] = out.get("search_engine") or RESCORING_TOOLS[tool]
    out["software_name"] = out["search_engine"]
    out["software_version"] = out.get("search_engine_version") or ""
    return out


@dataclass
class ImmunoDatapoint(DatapointBase):
    """
    Results of one immunopeptidomics benchmark run.

    The headline metrics are the mean over runs; ``results`` holds the per-run values and the
    full summary. ``dataset`` and ``dataset_version`` identify the ground truth, so results scored
    against different samples are never compared.

    Attributes
    ----------
    id : str
        Unique identifier of the run.
    software_name : str
        Name of the software that identified the peptides: the input format, or for the output
        of a rescoring tool (``RESCORING_TOOLS``) the search engine it rescored.
    software_version : str
        Version of the software.
    search_engine : str
        Search engine used.
    search_engine_version : str
        Version of the search engine.
    rescoring : str
        Rescoring tool, if any.
    rescoring_version : str
        Version of the rescoring tool.
    rescoring_model : str
        Model used by the rescoring tool (e.g. a timsTOF fragmentation model).
    ident_fdr_psm, ident_fdr_peptide : float
        FDR thresholds set in the workflow.
    precursor_mass_tolerance, fragment_mass_tolerance : str
        Mass tolerances with unit.
    enzyme : str
        Enzyme setting (non-specific for HLA class I).
    min_peptide_length, max_peptide_length : int
        Peptide length range.
    max_mods : int
        Maximum number of variable modifications.
    enable_match_between_runs : bool
        Whether match-between-runs was enabled.
    is_temporary : bool
        Whether the datapoint is not yet submitted.
    intermediate_hash : str
        SHA1 of the intermediate table.
    dataset, dataset_version : str
        Ground-truth definition the result was scored against.
    r : float
        Entrapment-to-sample ratio used.
    r_source : str
        "fasta" or "library".
    declared_fdr : float
        FDR threshold of the benchmark rule.
    max_qvalue_in_input : float
        Q-value up to which the submitted file is complete, i.e. the threshold the tool applied (for
        DIA-NN the stricter of the run-level and global cut-off).
    n_below_min_length : int
        Unique peptides accepted at the declared FDR that are shorter than the dataset's
        ``min_length``; they are left out of the scoring.
    n_runs : int
        Number of runs scored.
    n_sample_peptides : float
        Mean number of identified sample peptides per run.
    recall, recall_all : float
        Mean recall over runs, within the recall scope and over all designed peptides.
    fdp_lower, fdp_upper : float
        Mean bounds of the false discovery proportion.
    category : str
        "valid", "inconclusive", "invalid" or "not estimable" at the declared FDR.
    results : dict
        ``summary`` and ``per_run`` metrics, and ``by_qvalue``: the summary metrics at each q-value
        threshold the submitted file covers (``ImmunoScores.metrics_by_qvalue``).
    comments : str
        Comments for plotting.
    proteobench_version : str
        ProteoBench version used.
    """

    id: str = None
    software_name: str = None
    software_version: str = ""
    search_engine: str = ""
    search_engine_version: str = ""
    rescoring: str = ""
    rescoring_version: str = ""
    rescoring_model: str = ""
    ident_fdr_psm: float = math.nan
    ident_fdr_peptide: float = math.nan
    precursor_mass_tolerance: str = ""
    fragment_mass_tolerance: str = ""
    enzyme: str = ""
    min_peptide_length: int = 0
    max_peptide_length: int = 0
    max_mods: int = 0
    enable_match_between_runs: bool = False
    is_temporary: bool = True
    intermediate_hash: str = ""
    dataset: str = ""
    dataset_version: str = ""
    r: float = math.nan
    r_source: str = ""
    declared_fdr: float = math.nan
    max_qvalue_in_input: float = math.nan
    n_below_min_length: int = 0
    n_runs: int = 0
    n_sample_peptides: float = math.nan
    recall: float = math.nan
    recall_all: float = math.nan
    fdp_lower: float = math.nan
    fdp_upper: float = math.nan
    category: str = ""
    results: dict = None
    comments: str = ""
    proteobench_version: str = ""

    def generate_id(self) -> None:
        """
        Generate a unique ID from the software name, its version and a timestamp with microseconds.
        """
        time_stamp = datetime.now().strftime("%Y%m%d_%H%M%S_%f")
        self.id = "_".join([self.software_name, str(self.software_version), time_stamp])
        logging.info(f"Assigned the following ID to this run: {self.id}")

    @staticmethod
    def generate_datapoint(
        intermediate: pd.DataFrame,
        input_format: str,
        user_input: dict,
        scores: ImmunoScores = None,
        r: float = math.nan,
        r_source: str = "fasta",
        max_qvalue_in_input: float = math.nan,
        by_qvalue: Dict[str, Dict[str, Any]] = None,
        n_below_min_length: int = 0,
    ) -> pd.Series:
        """
        Build the datapoint of one submission.

        Parameters
        ----------
        intermediate : pd.DataFrame
            Output of :meth:`ImmunoScores.generate_intermediate`.
        input_format : str
            The software tool (input format).
        user_input : dict
            Workflow parameters; missing keys are left at their defaults.
        scores : ImmunoScores
            The scorer, which carries the dataset and the declared FDR.
        r : float
            Entrapment-to-sample ratio used.
        r_source : str
            "fasta" or "library".
        max_qvalue_in_input : float
            Q-value up to which the submitted file is complete.
        by_qvalue : Dict[str, Dict[str, Any]], optional
            Summary metrics per q-value threshold, from ``ImmunoScores.metrics_by_qvalue``.
        n_below_min_length : int, optional
            Unique accepted peptides left out because they are shorter than ``min_length``.

        Returns
        -------
        pd.Series
            The datapoint as a Series.
        """
        user_input = {k: ("" if v is None else v) for k, v in dict(user_input or {}).items()}
        fields = {f.name for f in dataclasses.fields(ImmunoDatapoint)}
        params = {k: v for k, v in user_input.items() if k in fields and k not in ("id", "results")}
        params["software_name"] = input_format
        params = name_rescored_workflow(params, input_format)
        point = ImmunoDatapoint(
            **params,
            intermediate_hash=hashlib.sha1(intermediate.to_string().encode("utf-8")).hexdigest(),
            comments=user_input.get("comments_for_plotting", ""),
            proteobench_version=proteobench.__version__,
        )
        point.generate_id()

        metrics = ImmunoDatapoint.get_metrics(intermediate, scores=scores, r=r)
        summary = metrics["summary"]
        point.dataset = scores.dataset.name
        point.dataset_version = scores.dataset.version
        point.r, point.r_source = r, r_source
        point.declared_fdr = scores.declared_fdr
        point.max_qvalue_in_input = max_qvalue_in_input
        point.n_below_min_length = int(n_below_min_length)
        point.n_runs = int(summary["n_runs"])
        point.n_sample_peptides = summary.get("n_sample", math.nan)
        point.recall = summary.get("recall", math.nan)
        point.recall_all = summary.get("recall_all", math.nan)
        point.fdp_lower = summary.get("fdp_lower", math.nan)
        point.fdp_upper = summary.get("fdp_upper", math.nan)
        point.category = summary["category"]
        point.results = {**metrics, "by_qvalue": by_qvalue or {}}
        return pd.Series(dataclasses.asdict(point))

    @staticmethod
    def get_metrics(intermediate: pd.DataFrame, scores: ImmunoScores = None, r: float = math.nan) -> Dict[str, Any]:
        """
        Compute the metrics of an intermediate table.

        Parameters
        ----------
        intermediate : pd.DataFrame
            Output of :meth:`ImmunoScores.generate_intermediate`.
        scores : ImmunoScores
            The scorer.
        r : float
            Entrapment-to-sample ratio used.

        Returns
        -------
        Dict[str, Any]
            ``summary`` and ``per_run`` metrics.
        """
        if scores is None:
            raise ValueError("scores must be provided for ImmunoDatapoint.get_metrics")
        return scores.calculate_metrics(intermediate, r)
