"""Streamlit interface of the immunopeptidomics identification modules."""

import json
import math
from typing import Any, Dict, Optional

import pandas as pd
import streamlit as st

from proteobench.datapoint.immuno_datapoint import name_rescored_workflow

from .base import BaseUIModule
from .tabs import (
    tab1_view_public_results,
    tab2_upload_results,
    tab3_view_single_result,
    tab4_view_public_and_new_results,
)

FDP_BOUNDS = {"Upper bound": "upper", "Lower bound": "lower"}
GROUP_COLUMNS = ["True", "Truncated", "Co-present group", "Background", "False"]
QVALUE_THRESHOLDS = ["0.001", "0.01", "0.05", "0.1", "1.0"]


def _pct(value: Any) -> str:
    """Format a fraction as a percentage, or 'not estimable' when missing."""
    try:
        value = float(value)
    except (TypeError, ValueError):
        return "not estimable"
    return "not estimable" if math.isnan(value) else f"{value:.2%}"


def _count(value: Any) -> str:
    """Format a peptide count (a mean over runs gets one decimal), or 'n/a' when missing."""
    try:
        value = float(value)
    except (TypeError, ValueError):
        return "n/a"
    if math.isnan(value):
        return "n/a"
    return f"{value:,.0f}" if value.is_integer() else f"{value:,.1f}"


def _results_dict(datapoint: pd.Series) -> Dict[str, Any]:
    """Return the nested results of a datapoint, whether it was read from JSON or generated in this session."""
    results = datapoint.get("results")
    if isinstance(results, str):
        results = json.loads(results)
    return results or {}


class ImmunoUIObjects(BaseUIModule):
    """
    Streamlit UI of the immunopeptidomics identification modules.

    Parameters
    ----------
    variables : dataclass
        The page variables (session-state keys, texts, settings paths).
    ionmodule : ImmunoModule
        The module.
    parsesettingsbuilder : ParseSettingsBuilder
        The parse settings builder.
    page_name : str
        The sidebar label of the page.
    """

    def __init__(self, variables, ionmodule, parsesettingsbuilder, page_name: str = "/") -> None:
        super().__init__(
            variables=variables, ionmodule=ionmodule, parsesettingsbuilder=parsesettingsbuilder, page_name=page_name
        )
        st.session_state[self.variables.submit] = False

    # ------------------------------------------------------------------ shared pieces

    def _key(self, name: str) -> str:
        """Session-state key unique to this page."""
        return f"{self.variables.prefix_params}{name}"

    def _local_results_note(self) -> None:
        """Tell the user where results come from while there is no results repository."""
        if getattr(self.ionmodule, "uses_local_results", False):
            st.info(
                f"There is no public results repository for this module yet. Results are read from the local folder "
                f"`{self.ionmodule.local_results_dir}`; save a result there in the *Submit New Results* tab.",
                icon="ℹ️",
            )

    def _intro(self) -> None:
        """Expander with the module description, the ground truth and the metric definitions."""
        dataset = self.ionmodule.dataset
        with st.expander("Quick Intro", expanded=False):
            st.markdown(open(self.variables.description_module_md, encoding="utf-8").read())
            st.markdown(
                f"**Ground truth:** {dataset.name} (version {dataset.version}). {len(dataset.sample_peptides):,} "
                "designed pool peptides."
            )
            st.markdown("##### Ground-truth groups")
            st.markdown(
                "- **True**: peptide of the sample.\n"
                "- **Truncated**: sub-sequence of a sample peptide, mostly N-terminal; usually in-source fragments.\n"
                "- **Co-present group**: peptide of a library group that is in the sample although not in the "
                "design. Never counted as entrapment.\n"
                f"- **Background**: {', '.join(dataset.background_categories)} entries. Real, not false.\n"
                "- **False**: peptide (or sub-sequence) of an absent library group, or any other peptide (the "
                "human proteome)."
            )
            st.markdown("##### FDP estimates")
            st.latex(r"\mathrm{FDP}_{\mathrm{lower}} = \frac{N_{\mathrm{false}}}{N_{\mathrm{all}}}")
            st.latex(
                r"\mathrm{FDP}_{\mathrm{upper}} = \frac{N_{\mathrm{false}}\,(1 + 1/r)}{N_{\mathrm{all}}}",
            )
            st.markdown(
                f"*r* is the ratio of entrapment to sample search space: {dataset.r:g} for the benchmark FASTA "
                f"(unique {dataset.min_length}- to {dataset.max_length}-mers). For a spectral-library search it is "
                f"computed from the uploaded library; below {dataset.min_r:g} the upper bound is not estimable. "
                "Workflows are **valid** when the upper bound is at most the declared FDR, **invalid** when the "
                "lower bound exceeds it, and **inconclusive** in between."
            )

    def _plot_options(self, submitted: bool) -> Dict[str, Any]:
        """Render the plot options and return them as plot_main_metric keyword arguments."""
        suffix = "submitted" if submitted else "main"
        selectbox_uuid = self.variables.selectbox_id_submitted_uuid if submitted else self.variables.selectbox_id_uuid
        tab1_view_public_results.initialize_main_selectbox(selectbox_id_uuid=selectbox_uuid, default_value="None")

        def render_bound():
            return st.radio(
                "FDP on the x axis", list(FDP_BOUNDS), horizontal=True, key=self._key(f"fdp_bound_{suffix}")
            )

        def render_threshold():
            return st.selectbox(
                "Q-value threshold",
                QVALUE_THRESHOLDS,
                index=QVALUE_THRESHOLDS.index("0.01"),
                help="Pool peptides identified and FDP of every workflow when its identifications are accepted at "
                "this peptide-level q-value. The benchmark rule is 0.01. Workflows whose result file does not "
                "reach the threshold (the tool already filtered at a stricter q-value) are not shown.",
                key=self._key(f"qvalue_threshold_{suffix}"),
            )

        def render_label():
            tab1_view_public_results.generate_main_selectbox(self.variables, selectbox_id_uuid=selectbox_uuid)

        def render_colorblind():
            return tab1_view_public_results.display_colorblindmode_selector(self.variables, use_submitted=submitted)

        results = self.render_plot_options_expander(
            filter_callbacks=[],
            selector_callbacks=[render_bound, render_threshold, render_label, render_colorblind],
            selector_cols_spec=[1, 1, 1, 1],
            expander_key=self._key(f"plot_options_{suffix}"),
        )
        threshold = float(results[1])
        new = self._new_datapoint() if submitted else None
        if new is not None and f"{threshold:g}" not in _results_dict(new).get("by_qvalue", {}):
            st.warning(
                f"The uploaded result has a largest q-value of {float(new['max_qvalue_in_input']):.4f}, so it is not "
                f"shown at a q-value threshold of {threshold:g}. Export the identifications at a higher q-value to "
                "compare it at this threshold.",
                icon="⚠️",
            )
        return {
            "fdp_bound": FDP_BOUNDS[results[0]],
            "y_metric": "n_sample_peptides",
            "qvalue_threshold": threshold,
            "colorblind_mode": bool(results[3]),
            "label": st.session_state.get(st.session_state.get(selectbox_uuid, ""), "None"),
            "alpha_warning": getattr(self.variables, "alpha_warning", False),
            "beta_warning": getattr(self.variables, "beta_warning", False),
            "archived_warning": getattr(self.variables, "archived_warning", False),
        }

    def _fdp_interval_plot(self, data: pd.DataFrame, plot_params: Dict[str, Any], suffix: str) -> None:
        """FDP interval (forest) plot of every shown workflow, laid out as in the entrapment module."""
        st.subheader("FDP Interval Plot")
        st.markdown(
            "Each row is one workflow. The thick horizontal bar spans the estimated FDP from the **lower bound** "
            "(left) to the **upper bound** (right); open circles mark both endpoints. Bar colour indicates the "
            "validity category at the q-value threshold: **green** = valid, **orange** = inconclusive, "
            "**red** = invalid. A grey open circle is a workflow whose upper bound is not estimable, drawn at its "
            "lower bound. With *Declared FDR*, the diamond marker (◆) shows each workflow's declared FDR; a star "
            "(★) marks the uploaded result. The bars on the right give the pool peptides identified."
        )
        options = QVALUE_THRESHOLDS + ["Declared FDR"]
        main_threshold = float(plot_params.get("qvalue_threshold", 0.01))
        default = next((o for o in QVALUE_THRESHOLDS if float(o) == main_threshold), "0.01")
        col_sort, col_thresh = st.columns([2, 1])
        with col_sort:
            sort_dir = st.radio(
                "Sort by pool peptides identified",
                ["Ascending ↑", "Descending ↓"],
                index=1,
                horizontal=True,
                key=self._key(f"forest_sort_{suffix}"),
            )
        with col_thresh:
            threshold_str = st.selectbox(
                "Q-value threshold",
                options,
                index=options.index(default) if default in options else options.index("0.01"),
                help="FDP bounds and pool peptides of every workflow at this peptide-level q-value. *Declared FDR* "
                "shows the headline metrics at each workflow's declared FDR.",
                key=self._key(f"forest_threshold_{suffix}"),
            )
        threshold = None if threshold_str == "Declared FDR" else float(threshold_str)
        try:
            fig = self.ionmodule.get_plot_generator().plot_fdp_interval(
                data, qvalue_threshold=threshold, sort_ascending=sort_dir.startswith("Ascending")
            )
        except Exception as e:  # keep the main plot and table usable
            st.error(f"Could not render forest plot: {e}", icon="🚨")
            return
        st.plotly_chart(fig, key=self._key(f"fdp_interval_{suffix}"))

    def _all_datapoints(self) -> pd.DataFrame:
        """Public (or locally saved) datapoints plus the one uploaded in this session, if any."""
        for key in (self.variables.all_datapoints_submitted, self.variables.all_datapoints):
            df = st.session_state.get(key)
            if isinstance(df, pd.DataFrame) and not df.empty:
                return df
        return pd.DataFrame()

    def _new_datapoint(self) -> Optional[pd.Series]:
        """The datapoint uploaded in this session, if any."""
        df = st.session_state.get(self.variables.all_datapoints_submitted)
        if not isinstance(df, pd.DataFrame) or df.empty or "old_new" not in df.columns:
            return None
        new = df[df["old_new"] == "new"]
        return None if new.empty else new.iloc[-1]

    def _intermediate(self, datapoint: pd.Series) -> Optional[pd.DataFrame]:
        """Intermediate table of a datapoint: from this session, the local results folder or server storage."""
        new = self._new_datapoint()
        if new is not None and new["intermediate_hash"] == datapoint["intermediate_hash"]:
            return st.session_state.get(self.variables.result_perf)
        local = self.ionmodule.load_intermediate(datapoint["intermediate_hash"])
        if local is not None or getattr(self.ionmodule, "uses_local_results", False):
            return local
        return tab3_view_single_result.load_public_performance_data(datapoint["intermediate_hash"])

    def _pool_lengths(self) -> Optional[pd.Series]:
        """Lengths of the designed pool peptides of the module's dataset."""
        try:
            return pd.Series([len(p) for p in self.ionmodule.dataset.sample_peptides])
        except Exception:  # the plots are drawn without the pool bars
            return None

    def _length_note(self, datapoint: pd.Series) -> None:
        """Note on the accepted peptides that were left out because they are too short."""
        try:
            n = int(datapoint.get("n_below_min_length", 0) or 0)
        except (TypeError, ValueError):
            return
        if n:
            k = self.ionmodule.dataset.min_length
            st.info(
                f"{n:,} accepted peptide{'s' if n != 1 else ''} shorter than {k} residues "
                f"{'were' if n != 1 else 'was'} left out. The benchmark covers peptides of at least {k} "
                "residues, so these count neither as true nor as false.",
                icon="ℹ️",
            )

    def _summary_metrics(self, datapoint: pd.Series) -> None:
        """Headline metrics of one datapoint."""
        cols = st.columns(5)
        cols[0].metric(
            "Pool peptides identified",
            _count(datapoint.get("n_sample_peptides")),
            help="Unique designed pool peptides accepted per run (mean over runs).",
        )
        cols[1].metric("FDP lower bound", _pct(datapoint.get("fdp_lower")))
        cols[2].metric("FDP upper bound", _pct(datapoint.get("fdp_upper")))
        cols[3].metric("Category", str(datapoint.get("category", "")))
        cols[4].metric("r", f"{float(datapoint.get('r', math.nan)):.3g}", help=f"from the {datapoint.get('r_source')}")

    def _per_run_table(self, datapoint: pd.Series) -> None:
        """Per-run counts and FDP of one datapoint."""
        per_run = pd.DataFrame(_results_dict(datapoint).get("per_run", {})).T
        if per_run.empty:
            return
        groups = {f"n_group_{g}": g for g in GROUP_COLUMNS}
        table = per_run[[c for c in groups if c in per_run.columns]].rename(columns=groups).astype(int)
        table["All"] = per_run["n_all"].astype(int)
        table["FDP lower"] = per_run["fdp_lower"].astype(float).map(_pct)
        table["FDP upper"] = per_run["fdp_upper"].astype(float).map(_pct)
        table.index.name = "Run"
        st.dataframe(table, width="stretch")

    def _datapoint_selector(self, label: str, key: str, index: int = 0) -> Optional[pd.Series]:
        """Select one datapoint by its ID; the uploaded one comes first."""
        df = self._all_datapoints()
        if df.empty:
            return None
        df = df.sort_values("old_new", ascending=True, key=lambda s: s.ne("new")) if "old_new" in df else df
        labels = [
            f"{row['id']} (uploaded)" if row.get("old_new") == "new" else str(row["id"]) for _, row in df.iterrows()
        ]
        choice = st.selectbox(
            label, range(len(df)), index=min(index, len(df) - 1), format_func=labels.__getitem__, key=key
        )
        return df.iloc[choice]

    # ------------------------------------------------------------------ tabs

    @st.fragment
    def display_all_data_results_main(self) -> None:
        """Tab 1: public results."""
        tab1_view_public_results.initialize_main_data_points(variables=self.variables, ionmodule=self.ionmodule)
        self._intro()
        self._local_results_note()
        st.subheader("Pool peptides identified and FDP of all workflows")
        st.markdown(
            "Each point is one workflow: pool peptides identified (mean over the submitted runs) against the "
            "estimated FDP. The dashed line marks the selected q-value threshold (the benchmark rule is 0.01). "
            "An open marker means the upper bound is not estimable and the point is drawn at the lower bound."
        )
        plot_params = self._plot_options(submitted=False)
        data = tab1_view_public_results.load_and_filter_main_data(self.variables, self.ionmodule, use_slider=False)
        if data.empty:
            st.info("No results available yet.", icon="ℹ️")
            return
        tab1_view_public_results.render_main_plot(
            self.variables,
            self.ionmodule,
            plot_params,
            data,
            render_forest_plot=lambda: self._fdp_interval_plot(data, plot_params, suffix="main"),
        )
        tab1_view_public_results.render_results_table(self.variables, data)

    def display_submission_form(self) -> None:
        """Tab 2: upload and score a result file."""
        tab2_upload_results.show_software_selector_and_alphadia_info(
            variables=self.variables, parsesettingsbuilder=self.parsesettingsbuilder, user_input=self.user_input
        )
        with st.form(key=self._key("main_form")):
            self.user_input["input_csv"] = st.file_uploader(
                "Software tool result file", help="The result file of the selected software tool."
            )
            self.user_input["input_csv_secondary"] = None
            if getattr(self.ionmodule, "library_sequence_columns", None):
                self.user_input["input_csv_secondary"] = st.file_uploader(
                    "Spectral library (only for library searches)",
                    type=["tsv", "parquet", "csv"],
                    help="r is computed from the library instead of the FASTA.",
                )
            self.user_input["submission_keyword"] = st.text_input(
                "Submission keyword (optional)",
                help="Shown with this result in the tables during this session only.",
                max_chars=50,
            )
            tab2_upload_results.generate_additional_parameters_fields(
                variables=self.variables, user_input=self.user_input
            )
            submit_button = st.form_submit_button("Parse and bench", help=self.variables.texts.Help.parse_button)

        if submit_button:
            with st.spinner("Classifying identifications against the ground truth..."):
                success = tab2_upload_results.process_submission_form(
                    variables=self.variables, ionmodule=self.ionmodule, user_input=self.user_input
                )
            if not success:
                return
            new = self._new_datapoint()
            if new is None:
                # add_current_data_point skips a result whose identifications equal an existing one.
                st.info("This result is identical to a result that is already in the results.", icon="ℹ️")
                return
            st.subheader("Result")
            self._summary_metrics(new)
            self._length_note(new)
            self._per_run_table(new)

    def display_indepth_plots(self) -> None:
        """Tab 3: one result in detail."""
        if self._all_datapoints().empty:
            st.info("Upload a result or save one to the local results folder first.", icon="ℹ️")
            return
        datapoint = self._datapoint_selector("Select a result", key=self._key("indepth_selector"))
        if datapoint is None:
            return
        self._summary_metrics(datapoint)
        self._length_note(datapoint)
        st.subheader("Per run")
        self._per_run_table(datapoint)

        intermediate = self._intermediate(datapoint)
        if intermediate is None:
            st.warning("The identifications of this result are not available, so no plots can be drawn.", icon="⚠️")
            return
        generator = self.ionmodule.get_plot_generator()
        plots = generator.generate_in_depth_plots(intermediate, pool_lengths=self._pool_lengths())
        tab3_view_single_result.display_plots_with_layout(plots, generator, self.variables, str(datapoint["id"]))
        st.download_button(
            "Download classified peptides",
            data=intermediate.to_csv(index=False).encode("utf-8"),
            file_name=f"{datapoint['id']}_classified_peptides.csv",
            mime="text/csv",
            key=self._key("download_intermediate"),
            icon=":material/download:",
        )

    @st.fragment
    def display_all_data_results_submitted(self) -> None:
        """Tab 4: public results together with the uploaded one."""
        tab4_view_public_and_new_results.initialize_submitted_data_points(
            variables=self.variables, ionmodule=self.ionmodule
        )
        st.subheader("Pool peptides identified and FDP, including the uploaded result")
        plot_params = self._plot_options(submitted=True)
        data = tab4_view_public_and_new_results.load_and_filter_submitted_data(self.variables, self.ionmodule)
        if data.empty:
            st.info("No results available yet.", icon="ℹ️")
            return
        tab4_view_public_and_new_results.render_submitted_main_plot(
            self.variables,
            self.ionmodule,
            plot_params,
            data,
            render_forest_plot=lambda: self._fdp_interval_plot(data, plot_params, suffix="submitted"),
        )
        tab4_view_public_and_new_results.render_submitted_results_table(data, self.variables)

    def display_workflow_comparison(self) -> None:
        """Tab 5: compare the peptides of two results."""
        st.subheader("Compare two results")
        if len(self._all_datapoints()) < 2:
            st.info("At least two results are needed.", icon="ℹ️")
            return
        col_a, col_b = st.columns(2)
        with col_a:
            a = self._datapoint_selector("Result A", key=self._key("compare_a"), index=0)
        with col_b:
            b = self._datapoint_selector("Result B", key=self._key("compare_b"), index=1)
        if a is None or b is None or a["intermediate_hash"] == b["intermediate_hash"]:
            st.info("Select two different results.", icon="ℹ️")
            return

        fields = ["software_name", "software_version", "search_engine"]
        fields += ["rescoring", "rescoring_version", "rescoring_model", "n_runs", "n_sample_peptides"]
        fields += ["fdp_lower", "fdp_upper", "category", "r", "r_source"]
        st.dataframe(pd.DataFrame({"A": a.reindex(fields), "B": b.reindex(fields)}).astype(str), width="stretch")

        inter_a, inter_b = self._intermediate(a), self._intermediate(b)
        if inter_a is None or inter_b is None:
            st.warning("The identifications of one of the results are not available.", icon="⚠️")
            return
        st.markdown("Peptides found in at least one run, per ground-truth group.")
        rows = []
        for group in GROUP_COLUMNS:
            pa = set(inter_a.loc[inter_a["group"] == group, "peptide"])
            pb = set(inter_b.loc[inter_b["group"] == group, "peptide"])
            rows.append({"Group": group, "Only A": len(pa - pb), "Both": len(pa & pb), "Only B": len(pb - pa)})
        st.dataframe(pd.DataFrame(rows).set_index("Group"), width="stretch")

        st.subheader("Length and charge distributions")
        st.markdown(
            "Number of identified pool peptides per length (left) and of their precursors per charge state "
            "(right), per run and averaged over the runs. The grey bars give the number of designed pool peptides "
            "of each length."
        )
        results = {f"A: {a['id']}": inter_a, f"B: {b['id']}": inter_b}
        fig = self.ionmodule.get_plot_generator().plot_length_charge_comparison(
            results, pool_lengths=self._pool_lengths()
        )
        st.plotly_chart(fig, key=self._key("compare_length_charge"))

        false_a = set(inter_a.loc[inter_a["group"] == "False", "peptide"])
        false_b = set(inter_b.loc[inter_b["group"] == "False", "peptide"])
        if false_a & false_b:
            with st.expander(f"False peptides reported by both ({len(false_a & false_b)})"):
                st.write(", ".join(sorted(false_a & false_b)))

    def display_public_submission_ui(self) -> None:
        """Tab 6: save the uploaded result (no public results repository yet)."""
        st.markdown(open(self.variables.description_submission_md, encoding="utf-8").read())
        new = self._new_datapoint()
        if new is None:
            st.info("Upload a result in the *Upload New Results* tab first.", icon="ℹ️")
            return
        st.markdown(f"Uploaded result: `{new['id']}`")
        self._summary_metrics(new)

        with open(self.variables.additional_params_json, encoding="utf-8") as f:
            fields = json.load(f)
        with st.form(key=self._key("save_form")):
            st.markdown("Workflow parameters, stored with the result:")
            values = {}
            cols = st.columns(2)
            for i, (name, spec) in enumerate(fields.items()):
                if name == "software_name":
                    continue
                with cols[i % 2]:
                    if spec.get("type") == "checkbox":
                        values[name] = st.checkbox(spec["label"], value=bool(spec.get("value", False)))
                    else:
                        values[name] = st.text_input(spec["label"], placeholder=str(spec.get("placeholder", "")))
            values["comments"] = st.text_area("Comments")
            save = st.form_submit_button(
                "Save to local results", disabled=not getattr(self.ionmodule, "uses_local_results", False)
            )
        if not getattr(self.ionmodule, "uses_local_results", False):
            st.warning("This module has a results repository; pull-request submission is not implemented yet.")
            return
        if save:
            point = new.copy()
            for name, value in values.items():
                if value not in ("", None):
                    point[name] = value
            point = pd.Series(name_rescored_workflow(point.to_dict()), name=point.name)
            path = self.ionmodule.save_local_result(point, st.session_state[self.variables.result_perf])
            # Reload the public results on the next run so the saved one is included.
            for key in (self.variables.all_datapoints, self.variables.all_datapoints_submitted):
                st.session_state.pop(key, None)
            st.success(f"Saved to `{path}`. It now shows up as a public result.", icon="✅")
