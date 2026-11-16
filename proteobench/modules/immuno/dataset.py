"""
Ground truth of an immunopeptidomics benchmark sample.

An :class:`ImmunoDataset` is read from a ``dataset.toml`` (see ``datasets/pepyrus_1k/``) and holds
everything that depends on the sample: which library peptides are in it, which are co-present or
absent, which non-library FASTA entries are real background, and the entrapment-to-sample ratio
``r`` of the search space. The modules and the scoring never hard-code any of this, so another
sample only needs another dataset folder.
"""

from __future__ import annotations

import gzip
import io
import re
import tomllib
import urllib.parse
import urllib.request
from collections import defaultdict
from dataclasses import dataclass, field
from functools import cached_property
from pathlib import Path
from typing import Dict, Iterable, Iterator, List, Optional, Set, Tuple

import pandas as pd

from proteobench.exceptions import ImmunoError

DATASETS_DIR = Path(__file__).resolve().parent / "datasets"
CACHE_DIR = Path.home() / ".cache" / "proteobench" / "immuno"

SAMPLE = "sample"
SAMPLE_TRUNCATED = "sample_truncated"
CO_PRESENT = "co_present"
LIBRARY_ABSENT = "library_absent"
LIBRARY_ABSENT_TRUNCATED = "library_absent_truncated"
ENTRAPMENT_PROTEOME = "entrapment_proteome"
FALSE_CLASSES = [LIBRARY_ABSENT, LIBRARY_ABSENT_TRUNCATED, ENTRAPMENT_PROTEOME]
DISPLAY_GROUPS = ["True", "Truncated", "Co-present group", "Background", "False"]

_MODIFICATION = re.compile(r"\[.*?\]|\(.*?\)|\{.*?\}")


def _open_text(path: Path) -> io.TextIOBase:
    """Open a plain or gzip-compressed text file."""
    if path.suffix == ".gz":
        return io.TextIOWrapper(gzip.open(path), encoding="utf-8")
    return open(path, encoding="utf-8")


def read_fasta(path: Path) -> Iterator[Tuple[str, str]]:
    """
    Yield (header without ">", sequence) for every entry of a plain or gzipped FASTA file.

    Parameters
    ----------
    path : Path
        FASTA file.

    Yields
    ------
    Tuple[str, str]
        Header and upper-case sequence.
    """
    header, buf = None, []
    with _open_text(path) as f:
        for line in f:
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(buf).upper()
                header, buf = line[1:].rstrip("\n"), []
            else:
                buf.append(line.strip())
    if header is not None:
        yield header, "".join(buf).upper()


@dataclass
class ImmunoDataset:
    """
    Ground-truth definition of an immunopeptidomics benchmark sample.

    Use :meth:`from_toml` to build one from a dataset folder. Peptide classification and the
    search-space ratio are computed lazily and cached on the instance.

    Attributes
    ----------
    name : str
        Human-readable name, stored in every datapoint.
    version : str
        Version of the definition, stored in every datapoint.
    ground_truth : pd.DataFrame
        Library peptides with columns ``sequence`` and ``group``.
    sample_groups : List[str]
        Library groups that make up the sample.
    co_present_groups : List[str]
        Library groups that are present but not designed into the sample.
    background_fasta : Optional[Path]
        FASTA with the non-library sample-side entries.
    background_categories : List[str]
        FASTA categories counted as background, in the order they are checked.
    entrapment_categories : List[str]
        FASTA categories that are entrapment (only needed to compute ``r``).
    fasta_rules : List[Tuple[re.Pattern, str]]
        Ordered (regex, category) rules for FASTA headers.
    default_category : str
        Category of a header that matches no rule.
    r : float
        Entrapment-to-sample ratio of the search FASTA.
    min_r : float
        Smallest ``r`` for which an upper bound of the FDP is reported.
    min_length, max_length : int
        Peptide length range of the search space.
    collapse_il : bool
        Treat isoleucine and leucine as the same residue.
    exclude_first_residue : List[str]
        Sample peptides starting with these residues are left out of the recall denominator.
    """

    name: str
    version: str
    ground_truth: pd.DataFrame
    sample_groups: List[str]
    co_present_groups: List[str] = field(default_factory=list)
    background_fasta: Optional[Path] = None
    background_categories: List[str] = field(default_factory=list)
    entrapment_categories: List[str] = field(default_factory=list)
    fasta_rules: List[Tuple[re.Pattern, str]] = field(default_factory=list)
    default_category: str = "entrapment"
    r: float = float("nan")
    min_r: float = 1.0
    min_length: int = 8
    max_length: int = 12
    collapse_il: bool = True
    exclude_first_residue: List[str] = field(default_factory=list)
    description: str = ""

    @classmethod
    def from_toml(cls, path: str | Path) -> "ImmunoDataset":
        """
        Load a dataset definition.

        Parameters
        ----------
        path : str or Path
            A ``dataset.toml`` file, or the folder that contains it. A relative folder name is
            looked up in the bundled ``datasets/`` folder first.

        Returns
        -------
        ImmunoDataset
            The dataset definition.
        """
        path = Path(path)
        if not path.exists() and (DATASETS_DIR / path).exists():
            path = DATASETS_DIR / path
        if path.is_dir():
            path = path / "dataset.toml"
        if not path.is_file():
            raise ImmunoError(f"Dataset definition not found: {path}")
        with open(path, "rb") as f:
            cfg = tomllib.load(f)
        base = path.parent

        gt_cfg = cfg["ground_truth"]
        gt = pd.read_csv(_resolve(gt_cfg["file"], base, gt_cfg.get("url")), sep="\t", dtype=str)
        gt = gt.rename(columns={gt_cfg.get("sequence_column", "sequence"): "sequence"})
        gt = gt.rename(columns={gt_cfg.get("group_column", "group"): "group"})[["sequence", "group"]]
        gt["sequence"] = gt["sequence"].str.upper()

        known = set(gt["group"])
        listed = set(gt_cfg["sample_groups"]) | set(gt_cfg.get("co_present_groups", []))
        if listed - known:
            raise ImmunoError(f"Groups not in the ground-truth table: {sorted(listed - known)}")

        background = cfg.get("background", {})
        space = cfg.get("search_space", {})
        peptides = cfg.get("peptides", {})
        rules = [(re.compile(rule["pattern"]), rule["category"]) for rule in cfg.get("fasta_rules", [])]
        return cls(
            name=cfg["dataset"]["name"],
            version=str(cfg["dataset"].get("version", "")),
            description=cfg["dataset"].get("description", ""),
            ground_truth=gt,
            sample_groups=list(gt_cfg["sample_groups"]),
            co_present_groups=list(gt_cfg.get("co_present_groups", [])),
            background_fasta=(
                _resolve(background["fasta"], base, background.get("fasta_url")) if background.get("fasta") else None
            ),
            background_categories=list(background.get("categories", [])),
            entrapment_categories=list(space.get("entrapment_categories", [])),
            fasta_rules=rules,
            default_category=cfg.get("fasta_rules_default", {}).get("category", "entrapment"),
            r=float(space.get("r", float("nan"))),
            min_r=float(space.get("min_r", 1.0)),
            min_length=int(peptides.get("min_length", 8)),
            max_length=int(peptides.get("max_length", 12)),
            collapse_il=bool(peptides.get("collapse_il", True)),
            exclude_first_residue=list(cfg.get("recall", {}).get("exclude_first_residue", [])),
        )

    # ------------------------------------------------------------------ sequences

    def normalise(self, sequence: str) -> str:
        """
        Return the bare sequence used for matching: modifications, charge suffix and non-letters removed.

        Parameters
        ----------
        sequence : str
            Peptide or peptidoform in any common notation (ProForma, UniMod, mass shifts).

        Returns
        -------
        str
            Upper-case residues, with I written as L when ``collapse_il`` is set.
        """
        seq = _MODIFICATION.sub("", str(sequence).split("/")[0])
        seq = re.sub(r"[^A-Za-z]", "", seq).upper()
        return seq.replace("I", "L") if self.collapse_il else seq

    def _subsequences(self, peptides: Iterable[str]) -> Set[str]:
        """All sub-sequences of at least ``min_length`` residues; shorter ones cannot be identified."""
        n = self.min_length
        return {p[i:j] for p in peptides for i in range(len(p)) for j in range(i + n, len(p) + 1)}

    def _windows(self, sequence: str) -> Iterator[str]:
        """All sub-sequences in the length range: the non-specific search space of one entry."""
        for k in range(self.min_length, self.max_length + 1):
            for i in range(len(sequence) - k + 1):
                yield sequence[i : i + k]

    def category(self, header: str) -> str:
        """
        Category of a FASTA entry by its header.

        Parameters
        ----------
        header : str
            FASTA header without the leading ">".

        Returns
        -------
        str
            Category of the first matching rule, or the default category.
        """
        for pattern, category in self.fasta_rules:
            if pattern.search(header):
                return category
        return self.default_category

    # ------------------------------------------------------------------ ground-truth sets

    def _group_peptides(self, groups: Iterable[str]) -> Set[str]:
        gt = self.ground_truth
        return {self.normalise(s) for s in gt.loc[gt["group"].isin(list(groups)), "sequence"]}

    @cached_property
    def sample_peptides(self) -> Set[str]:
        """Normalised peptides of the sample groups."""
        return self._group_peptides(self.sample_groups)

    @cached_property
    def sample_peptides_in_scope(self) -> Set[str]:
        """Sample peptides that count in the recall denominator."""
        gt = self.ground_truth
        excluded = tuple(self.exclude_first_residue)
        seqs = gt.loc[gt["group"].isin(self.sample_groups), "sequence"]
        if excluded:
            seqs = seqs[~seqs.str.startswith(excluded)]
        return {self.normalise(s) for s in seqs}

    @cached_property
    def _sets(self) -> Dict[str, Set[str]]:
        """Disjoint lookup sets for the library classes, built once."""
        library = {self.normalise(s) for s in self.ground_truth["sequence"]}
        sample = self.sample_peptides
        sample_sub = self._subsequences(sample) - sample
        present = self._group_peptides(self.co_present_groups) - sample
        present |= self._subsequences(present) - sample - sample_sub
        absent_sub = self._subsequences(library - sample - present) - library - sample_sub - present
        return {
            SAMPLE: sample,
            SAMPLE_TRUNCATED: sample_sub,
            CO_PRESENT: present,
            LIBRARY_ABSENT: library,
            LIBRARY_ABSENT_TRUNCATED: absent_sub,
        }

    @cached_property
    def _background_index(self) -> Dict[str, Dict[str, List[Tuple[str, int]]]]:
        """Per background category, a k-mer index (k = min_length) of its FASTA entries."""
        k = self.min_length
        index = {c: defaultdict(list) for c in self.background_categories}
        if self.background_fasta is None:
            return index
        for header, seq in read_fasta(self.background_fasta):
            cat = self.category(header)
            if cat not in index:
                continue
            seq = seq.replace("I", "L") if self.collapse_il else seq
            for i in range(max(len(seq) - k + 1, 1)):
                index[cat][seq[i : i + k]].append((seq, i))
        return index

    @property
    def classes(self) -> List[str]:
        """All classes in the order they are checked."""
        return [
            SAMPLE,
            SAMPLE_TRUNCATED,
            CO_PRESENT,
            LIBRARY_ABSENT,
            LIBRARY_ABSENT_TRUNCATED,
            *self.background_categories,
            ENTRAPMENT_PROTEOME,
        ]

    @property
    def display_groups(self) -> Dict[str, List[str]]:
        """Classes grouped as shown in the web interface and the datapoint summary."""
        return {
            "True": [SAMPLE],
            "Truncated": [SAMPLE_TRUNCATED],
            "Co-present group": [CO_PRESENT],
            "Background": list(self.background_categories),
            "False": list(FALSE_CLASSES),
        }

    def classify(self, peptide: str) -> str:
        """
        Assign a normalised peptide its ground-truth class.

        The checks run in this order: sample, truncated sample peptide, co-present group, peptide
        of an absent library group, sub-sequence of one, background categories, and finally the
        entrapment proteome for everything else.

        Parameters
        ----------
        peptide : str
            Normalised peptide (see :meth:`normalise`).

        Returns
        -------
        str
            The class name.
        """
        for name, members in self._sets.items():
            if peptide in members:
                return name
        k = self.min_length
        for cat, index in self._background_index.items():
            if any(seq.startswith(peptide, i) for seq, i in index.get(peptide[:k], ())):
                return cat
        return ENTRAPMENT_PROTEOME

    # ------------------------------------------------------------------ search-space ratio

    def library_ratio(self, peptides: Iterable[str]) -> float:
        """
        Entrapment-to-sample ratio of a spectral library, counted as unique peptides.

        A library search can only report library entries, so its ``r`` comes from the library
        rather than the FASTA. Each library peptide is classified like an identification; the
        false classes are entrapment.

        Parameters
        ----------
        peptides : Iterable[str]
            Library peptides in any notation.

        Returns
        -------
        float
            Number of entrapment peptides divided by the number of the others.
        """
        unique = {self.normalise(p) for p in peptides}
        n_false = sum(self.classify(p) in FALSE_CLASSES for p in unique)
        n_sample = len(unique) - n_false
        return n_false / n_sample if n_sample else float("nan")


def compute_search_space_ratio(dataset: ImmunoDataset, fasta: str | Path) -> float:
    """
    Entrapment-to-sample ratio ``r`` of a search FASTA, counted as unique sub-sequences.

    Sample part: sub-sequences of the sample and co-present groups and of the background
    categories. Entrapment part: sub-sequences of the entrapment categories and the absent library
    groups, minus those that also occur in the sample part (an identification of those is
    classified as sample). Counting unique sequences rather than windows keeps isoforms and
    duplicated entries from inflating the entrapment side. Library entries in the FASTA are skipped;
    the ground-truth table describes them. Sequences are hashed to keep memory low.

    This runs once per dataset, offline; store the result as ``r`` in ``dataset.toml``.

    Parameters
    ----------
    dataset : ImmunoDataset
        The dataset definition.
    fasta : str or Path
        The full search FASTA (plain or gzipped).

    Returns
    -------
    float
        The ratio r.
    """
    sample, entrapment = set(), set()
    for header, seq in read_fasta(Path(fasta)):
        cat = dataset.category(header)
        seq = seq.replace("I", "L") if dataset.collapse_il else seq
        if cat in dataset.background_categories:
            sample.update(hash(w) for w in dataset._windows(seq))
        elif cat in dataset.entrapment_categories:
            entrapment.update(hash(w) for w in dataset._windows(seq))

    in_sample = set(dataset.sample_groups) | set(dataset.co_present_groups)
    for seq, group in zip(dataset.ground_truth["sequence"], dataset.ground_truth["group"]):
        target = sample if group in in_sample else entrapment
        target.update(hash(w) for w in dataset._windows(dataset.normalise(seq)))
    return len(entrapment - sample) / len(sample)


def _resolve(location: str, base: Path, url: Optional[str] = None) -> Path:
    """
    Find a data file of a dataset.

    Parameters
    ----------
    location : str
        Path relative to ``base``, or an http(s) URL.
    base : Path
        Folder of the dataset definition.
    url : str, optional
        Where to download the file when ``location`` is a path that does not exist.

    Returns
    -------
    Path
        The local file, or a downloaded copy in ``CACHE_DIR``.
    """
    if location.startswith(("http://", "https://")):
        url = location
    else:
        local = (base / location).resolve()
        if local.exists() or not url:
            return local
    target = CACHE_DIR / urllib.parse.urlparse(url).path.lstrip("/")
    if not target.exists():
        target.parent.mkdir(parents=True, exist_ok=True)
        tmp = target.with_suffix(target.suffix + ".part")
        try:
            urllib.request.urlretrieve(url, tmp)
        except OSError as e:
            tmp.unlink(missing_ok=True)
            raise ImmunoError(f"Could not download {url}: {e}") from e
        tmp.rename(target)
    return target
