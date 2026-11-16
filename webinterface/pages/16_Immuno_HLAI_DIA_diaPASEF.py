"""
Streamlit UI for the HLA class I immunopeptidomics identification module, diaPASEF.
"""

import logging

from pages.base import BaseStreamlitUI
from pages.base_pages.immuno import ImmunoUIObjects
from pages.pages_variables.Immuno.Immuno_HLAI_DIA_diaPASEF_variables import VariablesImmunoHLAIDIAdiaPASEF
from pages.texts.generic_texts import WebpageTexts

from proteobench.io.parsing.parse_settings import ParseSettingsBuilder
from proteobench.modules.immuno.immuno_HLAI_DIA_diaPASEF import DIAImmunoHLAIModulediaPASEF

if __name__ == "__main__":
    logging.basicConfig(format="%(levelname)s:%(message)s", level=logging.INFO)
    st_ui = BaseStreamlitUI(
        variables=VariablesImmunoHLAIDIAdiaPASEF(),
        texts=WebpageTexts,
        ionmodule=DIAImmunoHLAIModulediaPASEF,
        parsesettingsbuilder=ParseSettingsBuilder,
        uiobjects=ImmunoUIObjects,
        page_name="Immuno HLA-I DIA diaPASEF",
    )
    st_ui.main_page()
