"""
Streamlit UI for the HLA class I immunopeptidomics identification module, DDA on a timsTOF.
"""

import logging

from pages.base import BaseStreamlitUI
from pages.base_pages.immuno import ImmunoUIObjects
from pages.pages_variables.Immuno.Immuno_HLAI_DDA_timsTOF_variables import VariablesImmunoHLAIDDAtimsTOF
from pages.texts.generic_texts import WebpageTexts

from proteobench.io.parsing.parse_settings import ParseSettingsBuilder
from proteobench.modules.immuno.immuno_HLAI_DDA_timsTOF import DDAImmunoHLAIModuletimsTOF

if __name__ == "__main__":
    logging.basicConfig(format="%(levelname)s:%(message)s", level=logging.INFO)
    st_ui = BaseStreamlitUI(
        variables=VariablesImmunoHLAIDDAtimsTOF(),
        texts=WebpageTexts,
        ionmodule=DDAImmunoHLAIModuletimsTOF,
        parsesettingsbuilder=ParseSettingsBuilder,
        uiobjects=ImmunoUIObjects,
        page_name="Immuno HLA-I DDA timsTOF",
    )
    st_ui.main_page()
