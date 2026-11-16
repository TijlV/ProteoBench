"""Tests of the immunopeptidomics identification modules.

Two datasets are used:

- the bundled Pepyrus 1k dataset, with a subset of a real Sage search of runs 4959 and 4961
  (``data/immuno/sage_1k_subset.tsv``). The expected class counts come from an independent
  implementation of the classification.
- a toy dataset (``data/immuno/toy_dataset``) with other group and column names, another host and
  numbers small enough to check by hand, to make sure nothing is specific to Pepyrus.
"""

import math
import os
from types import SimpleNamespace

import pandas as pd
import pytest

from proteobench.io.parsing.parse_settings import ParseSettingsBuilder
from proteobench.exceptions import ImmunoError
from proteobench.datapoint.immuno_datapoint import name_rescored_workflow
from proteobench.modules.constants import MODULE_SETTINGS_DIRS
from proteobench.modules.immuno.dataset import ImmunoDataset, compute_search_space_ratio
from proteobench.modules.immuno.immuno_HLAI_DDA_timsTOF import DDAImmunoHLAIModuletimsTOF
from proteobench.modules.immuno.immuno_HLAI_DIA_diaPASEF import DIAImmunoHLAIModulediaPASEF
from proteobench.plotting.plot_generator_immuno import ImmunoPlotGenerator, workflow_label
from proteobench.validation import ModuleValidationConfig, validate_submission

TESTDATA_DIR = os.path.join(os.path.dirname(__file__), "data", "immuno")
SAGE_SUBSET = os.path.join(TESTDATA_DIR, "sage_1k_subset.tsv")
TOY_DATASET = os.path.join(TESTDATA_DIR, "toy_dataset")
TOY_CUSTOM = os.path.join(TESTDATA_DIR, "toy_custom.tsv")

# Unique accepted peptides per class in the Sage subset.
SAGE_EXPECTED = {
    "4959": {
        "sample": 25,
        "sample_truncated": 10,
        "co_present": 6,
        "library_absent": 2,
        "library_absent_truncated": 0,
        "construct": 0,
        "host": 3,
        "contaminant": 3,
        "entrapment_proteome": 4,
    },
    "4961": {
        "sample": 25,
        "sample_truncated": 9,
        "co_present": 6,
        "library_absent": 3,
        "library_absent_truncated": 1,
        "construct": 3,
        "host": 3,
        "contaminant": 3,
        "entrapment_proteome": 4,
    },
}


def toy_module(module_cls=DDAImmunoHLAIModuletimsTOF):
    """A module scored against the toy dataset and accepting any run name."""
    module = module_cls.offline()
    module.dataset_location = TOY_DATASET
    module.run_mapper = {}
    return module


@pytest.fixture(scope="module")
def toy():
    return ImmunoDataset.from_toml(TOY_DATASET)


@pytest.fixture(scope="module")
def result():
    return toy_module().score(TOY_CUSTOM, "Custom", {"software_version": "1"})


@pytest.fixture(scope="module")
def pepyrus():
    try:
        return ImmunoDataset.from_toml("pepyrus_1k")
    except ImmunoError as e:
        pytest.skip(f"Pepyrus 1k dataset files not available: {e}")


@pytest.fixture(scope="module")
def sage_result(pepyrus):
    return DDAImmunoHLAIModuletimsTOF.offline().score(SAGE_SUBSET, "Sage", {"software_version": "0.14.6"})


def test_pepyrus_dataset(pepyrus):
    # 1k pool = order1 G4; 68 of its peptides start with P or G and are out of the recall scope.
    assert len(pepyrus.sample_peptides) == 1000
    assert len(pepyrus.sample_peptides_in_scope) == 932
    assert pepyrus.r == pytest.approx(6.47)


class TestToyDataset:
    @pytest.mark.parametrize(
        "peptide,expected",
        [
            ("SAMPLEPEPK", "sample"),
            ("SAMPIEPEPK", "sample"),  # I and L are collapsed
            ("AMPLEPEPK", "sample_truncated"),
            ("CQPRESENTK", "co_present"),
            ("QPRESENTK", "co_present"),  # sub-sequences of co-present peptides are co-present
            ("ASSENTPEPK", "library_absent"),
            ("SSENTPEPK", "library_absent_truncated"),
            ("WWWWYYYY", "host"),
            ("HHHHHHHG", "construct"),
            ("DDDDEEEEK", "entrapment_proteome"),
        ],
    )
    def test_classify(self, toy, peptide, expected):
        assert toy.classify(toy.normalise(peptide)) == expected

    @pytest.mark.parametrize(
        "notation", ["S[Phospho]AMPLEPEPK", "SAMPLEPEPK/2", "S(UniMod:21)AMPLEPEPK", "[+42.0106]-SAMPLEPEPK"]
    )
    def test_normalise_strips_modifications(self, toy, notation):
        assert toy.normalise(notation) == "SAMPLEPEPK"

    def test_search_space_ratio(self, toy):
        # Unique 8-mers. Sample side: sample (3 + 2), co-present (3), host (3), construct (3) = 14.
        # Entrapment side: absent group (3) and the human entry (6, of which 3 are also sample) = 6.
        r = compute_search_space_ratio(toy, os.path.join(TOY_DATASET, "search.fasta"))
        assert r == pytest.approx(6 / 14)


class TestToyModule:
    def test_filtering(self, result):
        intermediate, _, _ = result
        # Decoy, rank-2 and q > 1% rows are dropped; the two notations of SAMPLEPEPK in run2 merge.
        assert len(intermediate) == 8
        run2 = intermediate[intermediate["run"] == "run2"].set_index("peptide")
        assert run2.loc["SAMPLEPEPK", "n_psms"] == 2

    def test_metrics(self, result):
        _, datapoint, _ = result
        run1, run2 = datapoint["results"]["per_run"]["run1"], datapoint["results"]["per_run"]["run2"]
        assert (run1["n_all"], run1["n_false"]) == (6, 2)
        assert run1["fdp_lower"] == pytest.approx(2 / 6)
        assert run1["fdp_upper"] == pytest.approx(2 * (1 + 1 / 4) / 6)
        assert (run1["recall"], run1["recall_all"]) == (1.0, 0.5)
        assert (run2["n_all"], run2["fdp_upper"], run2["recall_all"]) == (2, 0.0, 1.0)
        # The datapoint reports the mean over runs.
        assert datapoint["recall_all"] == pytest.approx(0.75)
        assert datapoint["fdp_upper"] == pytest.approx(2 * 1.25 / 6 / 2)
        assert datapoint["category"] == "invalid"

    def test_metrics_by_qvalue(self, result):
        _, datapoint, _ = result
        by_q = datapoint["results"]["by_qvalue"]
        # The toy file reaches q = 0.5, so 1.0 is left out; the declared FDR equals the headline metrics.
        assert list(by_q) == ["0.001", "0.01", "0.05", "0.1"]
        assert by_q["0.01"]["fdp_upper"] == pytest.approx(datapoint["fdp_upper"])
        assert by_q["0.01"]["recall_all"] == pytest.approx(datapoint["recall_all"])
        assert by_q["0.001"]["n_all"] <= by_q["0.01"]["n_all"] <= by_q["0.1"]["n_all"]
        shown = ImmunoPlotGenerator().plot_main_metric(pd.DataFrame([datapoint]), qvalue_threshold=1.0)
        assert not shown.data  # not reached by the file, so not plotted

    def test_only_the_latest_upload_is_new(self, tmp_path):
        module = toy_module()
        other = tmp_path / "other.tsv"
        pd.read_csv(TOY_CUSTOM, sep="\t").iloc[:-1].to_csv(other, sep="\t", index=False)
        _, points, _ = module.benchmarking(TOY_CUSTOM, "Custom", {}, pd.DataFrame())
        _, points, _ = module.benchmarking(str(other), "Custom", {}, points)
        assert points["old_new"].tolist() == ["old", "new"]
        # A duplicate of an existing result is not added and leaves no "new" row.
        _, points, _ = module.benchmarking(TOY_CUSTOM, "Custom", {}, points)
        assert points["old_new"].tolist() == ["old", "old"]

    def test_library_with_low_r_has_no_upper_bound(self, tmp_path):
        library = tmp_path / "library.tsv"
        pd.DataFrame({"StrippedPeptide": ["SAMPLEPEPK", "GTHEPEPKR", "CQPRESENTK", "ASSENTPEPK"]}).to_csv(
            library, sep="\t", index=False
        )
        _, datapoint, _ = toy_module(DIAImmunoHLAIModulediaPASEF).score(TOY_CUSTOM, "Custom", {}, library)
        assert datapoint["r_source"] == "library"
        assert datapoint["r"] == pytest.approx(1 / 3)
        assert math.isnan(datapoint["fdp_upper"])

    def test_library_decoys_are_ignored(self, tmp_path):
        library = tmp_path / "library.parquet"
        pd.DataFrame(
            {
                "Stripped.Sequence": ["SAMPLEPEPK", "GTHEPEPKR", "CQPRESENTK", "ASSENTPEPK", "KPEPELPMAS"],
                "Decoy": [0, 0, 0, 0, 1],
            }
        ).to_parquet(library)
        assert toy_module(DIAImmunoHLAIModulediaPASEF).library_ratio(str(library)) == pytest.approx(1 / 3)


class TestPepyrusModule:
    def test_sage_class_counts(self, sage_result):
        intermediate, datapoint, _ = sage_result
        counts = intermediate.groupby(["run", "class"]).size().unstack(fill_value=0)
        for run, expected in SAGE_EXPECTED.items():
            assert counts.reindex(columns=list(expected), fill_value=0).loc[run].to_dict() == expected
        assert datapoint["results"]["per_run"]["4959"]["fdp_upper"] == pytest.approx(6 * (1 + 1 / 6.47) / 53)
        assert datapoint["r_source"] == "fasta"
        datapoint.to_json()

    def test_plots(self, sage_result):
        intermediate, datapoint, _ = sage_result
        generator = ImmunoPlotGenerator()
        generator.plot_main_metric(pd.DataFrame([datapoint]))
        plots = generator.generate_in_depth_plots(intermediate)
        assert set(plots) == {p for layout in generator.get_in_depth_plot_layout() for p in layout["plots"]}


@pytest.mark.parametrize(
    "fields,input_format,expected",
    [
        # Upload without parameters: named after the search engine the parse settings expect.
        (
            {"software_name": "MS2Rescore"},
            "MS2Rescore",
            {"software_name": "Sage", "software_version": "", "rescoring": "MS2Rescore"},
        ),
        # Save with the parameters entered: the tool version moves to rescoring_version.
        (
            {
                "software_name": "MS2Rescore",
                "rescoring": "MS2Rescore",
                "software_version": "4.0.2",
                "search_engine_version": "0.14.6",
            },
            "",
            {"software_name": "Sage", "software_version": "0.14.6", "rescoring_version": "4.0.2"},
        ),
        # Not a rescoring tool: unchanged.
        ({"software_name": "Sage", "software_version": "0.14.6"}, "Sage", {"software_name": "Sage"}),
    ],
)
def test_rescored_workflow_is_named_after_search_engine(fields, input_format, expected):
    named = name_rescored_workflow(fields, input_format)
    assert {k: named.get(k) for k in expected} == expected
    assert workflow_label(pd.Series(named)) in ("Sage", "Sage + MS2Rescore")


def test_decoy_column_with_empty_cells():
    parser = ParseSettingsBuilder(
        MODULE_SETTINGS_DIRS["immuno_HLAI_DDA_timsTOF"], "immuno_HLAI_DDA_timsTOF"
    ).build_parser("Custom")
    table = pd.DataFrame(
        {"Raw file": ["a"] * 3, "Sequence": ["SAMPLEPEPK"] * 3, "Q-Value": [0.001] * 3, "Decoy": [0.0, 1.0, None]}
    )
    assert parser.convert_to_standard_format(table)["Decoy"].tolist() == [False, True, False]


def test_diann_qvalue_is_max_of_run_and_global():
    parser = ParseSettingsBuilder(
        MODULE_SETTINGS_DIRS["immuno_HLAI_DIA_diaPASEF"], "immuno_HLAI_DIA_diaPASEF"
    ).build_parser("DIA-NN")
    report = pd.DataFrame(
        {
            "Run": ["a", "a"],
            "Stripped.Sequence": ["SAMPLEPEPK", "CQPRESENTK"],
            "Precursor.Charge": [1, 2],
            "Q.Value": [0.001, 0.02],
            "Global.Q.Value": [0.005, 0.001],  # fallback when Global.Peptidoform.Q.Value is absent
        }
    )
    standard = parser.convert_to_standard_format(report)
    assert standard["Q-Value"].tolist() == [0.005, 0.02]
    # Complete only up to the stricter export cut-off of the two columns.
    assert standard.attrs["qvalue_coverage"] == 0.005
    assert not standard["Decoy"].any()


class TestValidation:
    @pytest.mark.parametrize("module_id", ["immuno_HLAI_DDA_timsTOF", "immuno_HLAI_DIA_diaPASEF"])
    def test_modules_use_immuno_profile(self, module_id):
        config = ModuleValidationConfig.from_parse_settings(MODULE_SETTINGS_DIRS[module_id], module_id, "Custom")
        assert config.validation_profile == "immuno"

    def test_lengths_are_counted_without_modifications(self):
        standard = pd.DataFrame(
            {
                "Raw file": ["a", "a", "a"],
                "Sequence": ["S[Phospho]AMPLEPEPK", "SAMPLEPEPK/2", "SAMPLEPEPKSAMPLEPEPK"],
                "Charge": [1, 2, 2],
            }
        )
        params = SimpleNamespace(
            min_peptide_length=8, max_peptide_length=12, min_precursor_charge=1, max_precursor_charge=3
        )
        report = validate_submission(standard, parameters=params, config=ModuleValidationConfig(), profile="immuno")
        errors = [issue for issue in report.issues if issue.code == "peptide_length_out_of_range"]
        assert len(errors) == 1 and errors[0].examples == ["SAMPLEPEPKSAMPLEPEPK"]
