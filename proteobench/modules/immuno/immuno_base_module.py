"""
Base module for immunopeptidomics identification benchmarks with a known ground truth.

Workflows are scored on a sample of known content: each accepted peptide is classified against the
ground truth of the module's dataset (see :mod:`proteobench.modules.immuno.dataset`), which gives
recall and the false discovery proportion (FDP) with the combined entrapment estimator.

The dataset is named in ``module_settings.toml`` (``[general] dataset``), so a module can be moved to
another sample by changing that setting, and a new module for another sample is a new settings
folder plus a new dataset folder.
"""

from __future__ import annotations

import os
import tomllib
from functools import cached_property
from pathlib import Path
from tempfile import TemporaryDirectory
from typing import Any, Dict, Optional, Tuple

import pandas as pd
import requests

from proteobench.datapoint.immuno_datapoint import ImmunoDatapoint
from proteobench.exceptions import (
    ConvertStandardFormatError,
    DatapointAppendError,
    DatapointGenerationError,
    ImmunoError,
    IntermediateFormatGenerationError,
    ParseError,
    ParseSettingsError,
)
from proteobench.github.gh import GithubProteobotRepo
from proteobench.io.params.diann import extract_params as extract_params_diann
from proteobench.io.params.sage import extract_params as extract_params_sage
from proteobench.io.parsing.parse_settings import ParseSettingsBuilder
from proteobench.modules.constants import MODULE_SETTINGS_DIRS
from proteobench.modules.entrapment.entrapment_base_module import EntrapmentModule
from proteobench.modules.immuno.dataset import ImmunoDataset
from proteobench.plotting.plot_generator_base import PlotGeneratorBase
from proteobench.plotting.plot_generator_immuno import ImmunoPlotGenerator
from proteobench.score.immunoscores import ImmunoScores

PARAMS_JSON = os.path.join(os.path.dirname(__file__), "..", "..", "io", "params", "json", "Immuno", "immuno_HLAI.json")
LOCAL_RESULTS_ENV = "PROTEOBENCH_IMMUNO_RESULTS_DIR"
DEFAULT_LOCAL_RESULTS = Path.home() / ".cache" / "proteobench" / "immuno" / "results"


def _repo_exists(repo_name: str) -> bool:
    """
    Check whether a public GitHub repository exists, without cloning it.

    Parameters
    ----------
    repo_name : str
        Repository as "owner/name".

    Returns
    -------
    bool
        True when GitHub answers the repository page with HTTP 200.
    """
    try:
        return requests.head(f"https://github.com/{repo_name}", timeout=10).status_code == 200
    except requests.RequestException:
        return False


def _file_name(file: Any) -> str:
    """Name of a path or of an uploaded (Streamlit) file object."""
    return file if isinstance(file, (str, os.PathLike)) else getattr(file, "name", "")


def read_table(file: Any, columns: Optional[set] = None) -> pd.DataFrame:
    """
    Read a tab-separated or parquet table, keeping only ``columns`` when given.

    Parameters
    ----------
    file : str, PathLike or file object
        The table.
    columns : set, optional
        Columns to keep; the others are not loaded.

    Returns
    -------
    pd.DataFrame
        The table.
    """
    name = str(_file_name(file)).lower()
    if name.endswith(".parquet"):
        if not columns:
            return pd.read_parquet(file)
        import pyarrow.parquet as pq

        present = [c for c in pq.read_schema(file).names if c in columns]
        return pd.read_parquet(file, columns=present)
    usecols = (lambda c: c in columns) if columns else None
    sep = "," if name.endswith(".csv") else "\t"
    return pd.read_csv(file, sep=sep, usecols=usecols, low_memory=False)


class ImmunoModule(EntrapmentModule):
    """
    Base module for immunopeptidomics identification benchmarks.

    The GitHub and submission plumbing is shared with the entrapment module; parsing, scoring and
    datapoints are specific to the ground-truth approach.

    Parameters
    ----------
    token : Optional[str]
        The GitHub token.
    proteobench_repo_name : str
        The name of the ProteoBench results repository.
    proteobot_repo_name : str
        The name of the ProteoBot repository for pull requests.
    parse_settings_dir : str
        The directory containing the parse settings and module_settings.toml.
    module_id : str
        The module identifier.
    branch : Optional[str]
        If set, check out this branch of the results repository.
    """

    module_id: str = ""

    EXTRACT_PARAMS_DICT: Dict[str, Any] = {
        "Sage": extract_params_sage,
        "DIA-NN": extract_params_diann,
    }

    def __init__(
        self,
        token: Optional[str],
        proteobench_repo_name: str,
        proteobot_repo_name: str,
        parse_settings_dir: str,
        module_id: str,
        branch: Optional[str] = None,
    ):
        """
        Clone the results repository and load the module settings.

        The entrapment base ``__init__`` is not called: it requires the entrapment mapping file.

        Parameters
        ----------
        token : Optional[str]
            The GitHub token.
        proteobench_repo_name : str
            The name of the ProteoBench results repository.
        proteobot_repo_name : str
            The name of the ProteoBot repository for pull requests.
        parse_settings_dir : str
            The directory containing the parse settings and module_settings.toml.
        module_id : str
            The module identifier.
        branch : Optional[str]
            If set, check out this branch of the results repository.
        """
        self._load_module_settings(parse_settings_dir, module_id)
        self.t_dir = TemporaryDirectory().name
        self.t_dir_pr = TemporaryDirectory().name
        self.github_repo = GithubProteobotRepo(
            token,
            proteobot_repo_name=proteobot_repo_name,
            proteobench_repo_name=proteobench_repo_name,
            clone_dir=self.t_dir,
            clone_dir_pr=self.t_dir_pr,
            branch=branch,
        )
        # The results repositories do not exist yet. Until they do, or when
        # PROTEOBENCH_IMMUNO_RESULTS_DIR is set, results are read from and saved to a local folder.
        self.uses_local_results = LOCAL_RESULTS_ENV in os.environ or not _repo_exists(proteobench_repo_name)
        if self.uses_local_results:
            self.local_results_dir.mkdir(parents=True, exist_ok=True)
            self.github_repo.clone_dir = str(self.local_results_dir)
        else:
            self.github_repo.clone_repo()

    @classmethod
    def offline(cls) -> "ImmunoModule":
        """
        Create the module without cloning the results repository, for local scoring and tests.

        Returns
        -------
        ImmunoModule
            A module that can score submissions; GitHub-dependent methods are unavailable.
        """
        module = cls.__new__(cls)
        module._load_module_settings(MODULE_SETTINGS_DIRS[cls.module_id], cls.module_id)
        return module

    def _load_module_settings(self, parse_settings_dir: str, module_id: str) -> None:
        """Read module_settings.toml."""
        self.parse_settings_dir = parse_settings_dir
        self.module_id = module_id
        self.precursor_column_name = "peptide"
        with open(Path(parse_settings_dir) / "module_settings.toml", "rb") as f:
            settings = tomllib.load(f)
        general = settings["general"]
        dataset = Path(general["dataset"])
        if not dataset.is_absolute() and (Path(parse_settings_dir) / dataset).exists():
            dataset = Path(parse_settings_dir) / dataset
        self.dataset_location = dataset
        self.declared_fdr = float(general.get("declared_fdr", 0.01))
        self.run_mapper: Dict[str, str] = dict(settings.get("run_mapper", {}))
        self.library_sequence_columns = settings.get("library", {}).get("sequence_columns")

    @property
    def local_results_dir(self) -> Path:
        """
        Folder for results while there is no results repository.

        Returns
        -------
        Path
            ``$PROTEOBENCH_IMMUNO_RESULTS_DIR/<module_id>``, by default under ``~/.cache/proteobench/immuno/results``.
        """
        return Path(os.environ.get(LOCAL_RESULTS_ENV, DEFAULT_LOCAL_RESULTS)) / self.module_id

    def save_local_result(self, datapoint: pd.Series, intermediate: pd.DataFrame) -> Path:
        """
        Save a datapoint and its intermediate table to the local results folder.

        Saved datapoints are read back as public results when the module uses local results.

        Parameters
        ----------
        datapoint : pd.Series
            The datapoint.
        intermediate : pd.DataFrame
            Its intermediate table.

        Returns
        -------
        Path
            The datapoint JSON file.
        """
        folder = self.local_results_dir
        folder.mkdir(parents=True, exist_ok=True)
        point = datapoint.drop(labels=["old_new", "Highlight", "Keyword", "color", "hover_text"], errors="ignore")
        point["is_temporary"] = False
        path = folder / f"{point['intermediate_hash']}.json"
        path.write_text(point.to_json(indent=2))
        intermediate.to_csv(folder / f"{point['intermediate_hash']}_intermediate.csv.gz", index=False)
        return path

    def load_intermediate(self, intermediate_hash: str) -> Optional[pd.DataFrame]:
        """
        Load the intermediate table of a locally saved result.

        Parameters
        ----------
        intermediate_hash : str
            Hash of the result.

        Returns
        -------
        Optional[pd.DataFrame]
            The intermediate table, or None when it was not saved locally.
        """
        path = self.local_results_dir / f"{intermediate_hash}_intermediate.csv.gz"
        return pd.read_csv(path, dtype={"run": str}) if path.exists() else None

    @cached_property
    def dataset(self) -> ImmunoDataset:
        """The ground-truth definition named in the module settings."""
        return ImmunoDataset.from_toml(self.dataset_location)

    def is_implemented(self) -> bool:
        """
        Return whether the module is ready for public use.

        Returns
        -------
        bool
            False: the immunopeptidomics modules are in development.
        """
        return False

    # ------------------------------------------------------------------ scoring

    def _map_runs(self, standard: pd.DataFrame) -> pd.DataFrame:
        """
        Replace run names by their labels and reject runs that are not part of the benchmark sample.

        Parameters
        ----------
        standard : pd.DataFrame
            Standardised identifications with a ``Raw file`` column.

        Returns
        -------
        pd.DataFrame
            The same table with run labels in ``Raw file``.
        """
        if not self.run_mapper:
            return standard
        unknown = sorted(set(standard["Raw file"]) - set(self.run_mapper))
        if unknown:
            raise ImmunoError(
                f"{len(unknown)} run(s) in the input are not runs of this module's sample, e.g. {unknown[0]!r}. "
                f"Expected run names (without extension): {sorted(self.run_mapper)}"
            )
        return standard.assign(**{"Raw file": standard["Raw file"].map(self.run_mapper)})

    def library_ratio(self, library_file: Any) -> float:
        """
        Entrapment-to-sample ratio of an uploaded spectral library.

        Parameters
        ----------
        library_file : str, PathLike or file object
            Library as TSV or parquet with a peptide column listed in ``[library] sequence_columns``.

        Returns
        -------
        float
            The ratio r of the library.
        """
        if not self.library_sequence_columns:
            raise ImmunoError("This module does not take a spectral library.")
        if str(_file_name(library_file)).lower().endswith(".speclib"):
            raise ImmunoError("Binary .speclib libraries cannot be read; export the library as TSV or parquet.")
        library = read_table(library_file, set(self.library_sequence_columns) | {"Decoy"})
        present = [c for c in self.library_sequence_columns if c in library.columns]
        if not present:
            raise ImmunoError(f"The library has none of the peptide columns {self.library_sequence_columns}.")
        # DIA-NN libraries may contain decoy entries; they are not part of the searchable space.
        if "Decoy" in library.columns:
            library = library[pd.to_numeric(library["Decoy"], errors="coerce").fillna(0) == 0]
        return self.dataset.library_ratio(library[present[0]].dropna().unique())

    def score(
        self,
        input_file: Any,
        input_format: str,
        user_input: dict,
        input_file_secondary: Any = None,
    ) -> Tuple[pd.DataFrame, pd.Series, pd.DataFrame]:
        """
        Parse, filter and score a submission without touching the results repository.

        Parameters
        ----------
        input_file : str, PathLike or file object
            The workflow output.
        input_format : str
            The software tool, a key of the module's section in parse_settings_files.toml.
        user_input : dict
            Workflow parameters stored in the datapoint.
        input_file_secondary : str, PathLike or file object, optional
            Spectral library of a library search; r is then taken from it instead of the FASTA.

        Returns
        -------
        Tuple[pd.DataFrame, pd.Series, pd.DataFrame]
            The intermediate table (one row per run and peptide), the datapoint and the standardised input.
        """
        try:
            parser = ParseSettingsBuilder(
                parse_settings_dir=self.parse_settings_dir, module_id=self.module_id
            ).build_parser(input_format)
        except Exception as e:
            raise ParseSettingsError(f"No parse settings for {input_format!r} in module {self.module_id}: {e}")

        try:
            input_df = read_table(input_file, parser.needed_columns())
        except Exception as e:
            raise ParseError(f"Could not read the {input_format} file: {e}")

        try:
            standard = parser.convert_to_standard_format(input_df)
        except Exception as e:
            raise ConvertStandardFormatError(f"Error converting to the standard format: {e}")
        coverage = standard.attrs.get("qvalue_coverage")
        standard = self._map_runs(standard)

        targets = standard[(standard["Rank"] == 1) & ~standard["Decoy"]]
        accepted = targets[targets["Q-Value"] <= self.declared_fdr]
        if accepted.empty:
            raise ImmunoError(f"No rank-1 target identifications at a peptide-level q-value <= {self.declared_fdr}.")

        if input_file_secondary is not None:
            r, r_source = self.library_ratio(input_file_secondary), "library"
        else:
            r, r_source = self.dataset.r, "fasta"

        scores = ImmunoScores(self.dataset, declared_fdr=self.declared_fdr)
        try:
            intermediate = scores.generate_intermediate(accepted)
        except Exception as e:
            raise IntermediateFormatGenerationError(f"Error classifying the identifications: {e}")

        try:
            datapoint = ImmunoDatapoint.generate_datapoint(
                intermediate,
                input_format,
                user_input,
                scores=scores,
                r=r,
                r_source=r_source,
                max_qvalue_in_input=float(targets["Q-Value"].max()) if coverage is None else coverage,
                by_qvalue=scores.metrics_by_qvalue(targets, r, max_qvalue=coverage),
            )
        except Exception as e:
            raise DatapointGenerationError(f"Error generating the datapoint: {e}")
        return intermediate, datapoint, standard

    def benchmarking(
        self,
        input_file: Any,
        input_format: str,
        user_input: dict,
        all_datapoints: Optional[pd.DataFrame],
        input_file_secondary: Any = None,
    ) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
        """
        Score a submission and add it to the existing datapoints.

        Parameters
        ----------
        input_file : str, PathLike or file object
            The workflow output.
        input_format : str
            The software tool.
        user_input : dict
            Workflow parameters stored in the datapoint.
        all_datapoints : Optional[pd.DataFrame]
            Existing datapoints; read from the results repository when None.
        input_file_secondary : str, PathLike or file object, optional
            Spectral library of a library search.

        Returns
        -------
        Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]
            The intermediate table, all datapoints including the new one, and the standardised input. The
            submission is the only row with ``old_new == "new"``; when it duplicates an existing result
            (same ``intermediate_hash``) it is not added and no row is "new".
        """
        intermediate, datapoint, standard = self.score(input_file, input_format, user_input, input_file_secondary)
        if isinstance(all_datapoints, pd.DataFrame) and "old_new" in all_datapoints.columns:
            # Only this submission can be "new"; a duplicate of an existing result leaves no "new" row.
            all_datapoints = all_datapoints.assign(old_new="old")
        try:
            all_datapoints = self.add_current_data_point(datapoint, all_datapoints=all_datapoints)
        except Exception as e:
            raise DatapointAppendError(f"Error adding the datapoint: {e}")
        return intermediate, all_datapoints, standard

    def load_params_file(self, input_file: list, input_format: str, json_file: str = PARAMS_JSON):
        """
        Read workflow parameters from the tool's parameter file.

        Parameters
        ----------
        input_file : list
            Paths of the parameter files.
        input_format : str
            The software tool.
        json_file : str
            Parameter definitions of the module.

        Returns
        -------
        ProteoBenchParameters
            The parameters.
        """
        if input_format not in self.EXTRACT_PARAMS_DICT:
            raise ImmunoError(f"No parameter parser for {input_format}; enter the parameters by hand.")
        return super().load_params_file(input_file, input_format, json_file)

    def get_plot_generator(self, y_axis_title: str = None) -> PlotGeneratorBase:
        """
        Get the plot generator of the immunopeptidomics modules.

        Parameters
        ----------
        y_axis_title : str, optional
            Unused; accepted for interface compatibility.

        Returns
        -------
        PlotGeneratorBase
            The plot generator.
        """
        return ImmunoPlotGenerator()
