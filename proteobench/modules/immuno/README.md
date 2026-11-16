# Immunopeptidomics modules

Identification benchmarks for HLA class I immunopeptidomics on samples of known content. Each
accepted peptide is classified against a ground truth, which gives the number of designed pool
peptides identified and the false discovery proportion (FDP) measured with the combined entrapment
estimator (Wen et al., *Nat. Methods* 2025).

Status: in development (`is_implemented()` returns False). The web pages run locally; there are no
results repositories yet.

| Module | `module_id` | Tools |
|---|---|---|
| `DDAImmunoHLAIModuletimsTOF` | `immuno_HLAI_DDA_timsTOF` | Sage, MS2Rescore, Custom |
| `DIAImmunoHLAIModulediaPASEF` | `immuno_HLAI_DIA_diaPASEF` | DIA-NN, Custom |

## How it fits together

```
module_settings.toml          dataset = "pepyrus_1k", declared_fdr, run_mapper, [library] (DIA)
parse_settings_<tool>.toml    tool columns -> Raw file, Sequence, Q-Value, Rank, Decoy, Charge
        │
        ▼
ImmunoModule.score()          read -> ParseSettingsImmuno -> run labels -> rank 1, targets, length >= min_length,
                              q <= declared_fdr
        │
        ▼
ImmunoScores                  classify unique peptides per run (ImmunoDataset.classify), recall, FDP bounds
ImmunoDatapoint               mean over runs + per-run results, dataset name and version, r and its source
```

- `dataset.py`: `ImmunoDataset` (ground truth, classification, library `r`) and
  `compute_search_space_ratio()` (FASTA `r`, run once offline).
- `immuno_base_module.py`: `ImmunoModule`. `score()` needs no GitHub access;
  `benchmarking()` also adds the datapoint to the existing ones. Submission plumbing (`clone_pr`,
  hashes) is inherited from the entrapment module.
- Settings: `proteobench/io/parsing/io_parse_settings/immuno/HLAI/{DDA/timsTOF,DIA/diaPASEF}/`.
- Score, datapoint, plots: `proteobench/score/immunoscores.py`,
  `proteobench/datapoint/immuno_datapoint.py`, `proteobench/plotting/plot_generator_immuno.py`.

## Classes

Checked in this order; the first match wins. Sequences are compared without modifications, with
I and L collapsed.

| Class | Meaning | Group | False |
|---|---|---|---|
| `sample` | peptide of a sample group | True | no |
| `sample_truncated` | sub-sequence of a sample peptide (mostly in-source fragments) | Truncated | no |
| `co_present` | peptide or sub-sequence of a co-present group | Co-present group | no |
| `library_absent` | peptide of an absent library group | False | yes |
| `library_absent_truncated` | sub-sequence of one | False | yes |
| one per background category | found in the background FASTA (e.g. `construct`, `host`, `contaminant`) | Background | no |
| `entrapment_proteome` | anything else (the human proteome in the Pepyrus search) | False | yes |

FDP lower = N_false / N_all; upper = N_false × (1 + 1/r) / N_all. The upper bound is not reported
when r is below the dataset's `min_r` (an experimental library built from the sample itself holds
almost no entrapment). The category (valid, inconclusive, invalid) follows the entrapment module.

## Q-value thresholds

The headline metrics use the benchmark rule (q ≤ `declared_fdr`). Each datapoint also stores the
summary metrics at q ≤ 0.001, 0.01, 0.05, 0.1 and 1.0 in `results["by_qvalue"]`
(`ImmunoScores.metrics_by_qvalue`), for the q-value selector of the main plot, as in the entrapment
module. A threshold is stored only when the submitted file covers it: its q-values reach at least 90%
of the threshold, taking the strictest q-value column (for DIA-NN the run-level `--qvalue` cut-off,
not the global q-value, which goes up to 1). Thresholds without identifications are left out.

## Using another dataset

Everything sample-specific is in a dataset folder (`datasets/pepyrus_1k/`):

- `dataset.toml`: sample and co-present groups, peptide length range, recall exclusions,
  background categories, FASTA header rules, and `r`;
- the ground-truth table (one row per library peptide: sequence and group);
- the background FASTA (host, construct, spike-ins, contaminants; not the entrapment proteome).

To score another sample, for example a purpose-made Pepyrus pool:

1. Copy `datasets/pepyrus_1k/` to `datasets/<name>/` and replace the two data files. Data files
   are not stored in the repository: give each one a download location on the ProteoBench server
   (`url` for the ground truth, `fasta_url` for the background FASTA). A file that is in the
   dataset folder is used as is; otherwise it is downloaded once to `~/.cache/proteobench/immuno/`.
   Tests that need the dataset are skipped when its files cannot be found or downloaded.
2. Edit `dataset.toml`. Column and group names are free, and so are the background categories,
   as long as `fasta_rules` assigns them.
3. Compute `r` for the search FASTA you distribute with the module and write it into
   `dataset.toml`:

   ```python
   from proteobench.modules.immuno.dataset import ImmunoDataset, compute_search_space_ratio
   compute_search_space_ratio(ImmunoDataset.from_toml("<name>"), "search.fasta")
   ```

4. Either point an existing module at it (`dataset = "<name>"` in its `module_settings.toml`, and
   update `run_mapper`), or add a module for it: a settings folder, an entry in
   `MODULE_SETTINGS_DIRS` (`modules/constants.py`), `parse_settings_files.toml` and
   `MODULE_TO_CLASS` (`io/parsing/parse_settings.py`), and a subclass of `ImmunoModule` that sets
   `module_id`.

Datapoints store `dataset` and `dataset_version`, so results scored against different samples
can be told apart.

## Scoring a result locally

```python
from proteobench.modules.immuno.immuno_HLAI_DDA_timsTOF import DDAImmunoHLAIModuletimsTOF

module = DDAImmunoHLAIModuletimsTOF.offline()  # no GitHub clone
intermediate, datapoint, standard = module.score("results.sage.tsv", "Sage", {"software_version": "0.14.6"})
datapoint["results"]["per_run"]
```

For a DIA-NN spectral-library search, pass the library (TSV or parquet) as
`input_file_secondary`; `r` is then computed from the library.

## Web pages

`webinterface/pages/15_Immuno_HLAI_DDA_timsTOF.py` and `16_Immuno_HLAI_DIA_diaPASEF.py`, with
`webinterface/pages/base_pages/immuno.py` (`ImmunoUIObjects`) and the variables in
`webinterface/pages/pages_variables/Immuno/`. Run with `cd webinterface && streamlit run Home.py`.

While the results repositories `Results_immuno_HLAI_{DDA_timsTOF,DIA_diaPASEF}` do not exist, the
module reads and writes results in a local folder: `~/.cache/proteobench/immuno/results/<module_id>/`,
or `$PROTEOBENCH_IMMUNO_RESULTS_DIR/<module_id>/` when that variable is set (which also forces local
mode). The *Submit New Results* tab saves the uploaded result there (datapoint JSON plus the
classified peptides), after which it shows as a public result.

## Validation and resubmission

- Validation profile `immuno` (`proteobench/validation/profiles.py`), declared in both
  `module_settings.toml` files: charge range, peptide length (counted without modifications),
  mass tolerances, PSM FDR and software identity. The protein-ID check (no protein FASTA; peptides
  are checked against the ground truth) and the tryptic enzyme check are left out. It only runs
  once submission goes through the standard Submit tab; the local-results tab does not call it.
- `resubmit_datapoints.py` lists both modules in `REPO_MODULE_REGISTRY`. It only finds datapoints
  once the results repositories exist.

## Not done yet

- The results repositories and pull-request submission.
- Parsers for MSFragger + Percolator, MSBooster and MuMDIA; Percolator needs two input files.
- `MODULE_CLASSES` in `utils/server_io.py` (scripted submission) and the module documentation page.
