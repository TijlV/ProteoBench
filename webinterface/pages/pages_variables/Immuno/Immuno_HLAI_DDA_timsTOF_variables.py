"""
Variables for the immunopeptidomics module: HLA class I identification, DDA on a timsTOF.
"""

import os
from dataclasses import dataclass, field
from typing import List, Type

from pages.texts.generic_texts import WebpageTexts


@dataclass
class VariablesImmunoHLAIDDAtimsTOF:
    """
    Variables for the HLA class I identification module, DDA on a timsTOF.
    """

    all_datapoints: str = "all_datapoints_immuno_hlai_dda_timstof"
    all_datapoints_submission: str = "all_datapoints_submission_immuno_hlai_dda_timstof"
    input_df_submission: str = "input_df_submission_immuno_hlai_dda_timstof"
    result_performance_submission: str = "result_performance_submission_immuno_hlai_dda_timstof"
    submit: str = "submit_immuno_hlai_dda_timstof"
    fig_metric: str = "fig_metric_immuno_hlai_dda_timstof"
    fig_metric_submitted: str = "fig_metric_submitted_immuno_hlai_dda_timstof"
    fig_prefix: str = "fig_immuno_hlai_dda_timstof_"
    result_perf: str = "result_perf_immuno_hlai_dda_timstof"
    meta_data: str = "meta_data_immuno_hlai_dda_timstof"
    input_df: str = "input_df_immuno_hlai_dda_timstof"
    meta_file_uploader_uuid: str = "meta_file_uploader_uuid_immuno_hlai_dda_timstof"
    comments_submission_uuid: str = "comments_submission_uuid_immuno_hlai_dda_timstof"
    check_submission_uuid: str = "check_submission_uuid_immuno_hlai_dda_timstof"
    meta_data_text: str = "comments_for_submission_immuno_hlai_dda_timstof"
    check_submission: str = "check_submission_immuno_hlai_dda_timstof"
    button_submission_uuid: str = "button_submission_uuid_immuno_hlai_dda_timstof"
    df_head: str = "df_head_immuno_hlai_dda_timstof"

    all_datapoints_submitted: str = "all_datapoints_submitted_immuno_hlai_dda_timstof"
    highlight_list_submitted: List[str] = field(default_factory=list)
    selectbox_id_submitted_uuid: str = "selectbox_id_submitted_immuno_hlai_dda_timstof"
    selectbox_id_uuid: str = "selectbox_id_immuno_hlai_dda_timstof"
    colorblind_mode_selector_uuid: str = "colorblind_mode_selector_immuno_hlai_dda_timstof"
    colorblind_mode_selector_submitted_uuid: str = "colorblind_mode_selector_submitted_immuno_hlai_dda_timstof"
    download_selector_id_uuid: str = "download_selector_id_immuno_hlai_dda_timstof"
    dataset_selector_id_uuid: str = "dataset_selector_id_immuno_hlai_dda_timstof"
    placeholder_dataset_selection_container: str = "placeholder_dataset_selection_container_immuno_hlai_dda_timstof"
    placeholder_downloads_container: str = "placeholder_downloads_container_immuno_hlai_dda_timstof"
    metric_plot_labels: List[str] = field(
        default_factory=lambda: [
            "None",
            "software_version",
            "search_engine",
            "rescoring",
            "predictors_library",
            "precursor_mass_tolerance",
            "fragment_mass_tolerance",
            "max_mods",
            "r_source",
            "category",
        ]
    )

    highlight_list: List[str] = field(default_factory=list)
    first_new_plot: bool = True
    alpha_warning: bool = True
    beta_warning: bool = False
    archived_warning: bool = False
    enable_pmultiqc: bool = False
    github_link_pr: str = "github.com/Proteobot/Results_immuno_HLAI_DDA_timsTOF.git"

    # Sidebar metadata
    sidebar_label: str = "Immuno HLA-I DDA timsTOF"
    homepage_title: str = "HLA-I Identification on a timsTOF"
    graphical_abstract: str = ""
    documentation_description: str = (
        "Benchmark the pool peptides identified and the false discovery proportion of HLA class I "
        "identification workflows on timsTOF DDA runs of a synthetic peptide pool of known composition."
    )
    sidebar_path: str = "/Immuno_HLAI_DDA_timsTOF"
    sidebar_category: str = "DDA"
    keywords: List[str] = field(
        default_factory=lambda: [
            "DDA",
            "immunopeptidomics",
            "HLA",
            "MHC",
            "timsTOF",
            "peptide",
            "identification",
            "FDR",
            "entrapment",
            "Pepyrus",
        ]
    )

    description_module_md: str = "pages/markdown_files/Immuno/HLAI/introduction.md"
    description_input_file_md: str = "pages/markdown_files/Immuno/HLAI/DDA/input_file_description.md"
    description_table_md: str = "pages/markdown_files/Immuno/HLAI/table_description.md"
    description_slider_md: str = "pages/markdown_files/Immuno/HLAI/table_description.md"
    description_results_md: str = "pages/markdown_files/Immuno/HLAI/result_description.md"
    description_submission_md: str = "pages/markdown_files/Immuno/HLAI/submit_description.md"

    parse_settings_dir: str = "../proteobench/io/parsing/io_parse_settings/immuno/HLAI/DDA/timsTOF"

    texts: Type[WebpageTexts] = WebpageTexts
    # No module documentation page yet; see proteobench/modules/immuno/README.md.
    doc_url: str = "https://proteobench.readthedocs.io/en/latest/"
    # The benchmark runs on MassIVE, or an archive of them under $PROTEOBENCH_IMMUNO_DOWNLOAD_URL when set.
    raw_data_url: str = (
        f"{os.environ['PROTEOBENCH_IMMUNO_DOWNLOAD_URL'].rstrip('/')}/all_data_Immuno_HLAI_DDA_timsTOF.tar"
        if os.environ.get("PROTEOBENCH_IMMUNO_DOWNLOAD_URL")
        else "https://massive.ucsd.edu/ProteoSAFe/dataset.jsp?accession=MSV000096555"
    )

    additional_params_json: str = "../proteobench/io/params/json/Immuno/immuno_HLAI.json"
    title: str = "HLA class I identification - DDA timsTOF (alpha)"
    prefix_params: str = "immuno_hlai_dda_timstof_"
    params_json_dict: str = "params_json_dict_immuno_hlai_dda_timstof"
    params_file_dict: str = "params_file_dict_immuno_hlai_dda_timstof"
