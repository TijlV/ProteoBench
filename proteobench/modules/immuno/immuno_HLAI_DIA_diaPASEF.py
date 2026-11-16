"""
HLA class I immunopeptidomics identification module, diaPASEF.
"""

from __future__ import annotations

from typing import Optional

from proteobench.modules.constants import MODULE_SETTINGS_DIRS
from proteobench.modules.immuno.immuno_base_module import ImmunoModule


class DIAImmunoHLAIModulediaPASEF(ImmunoModule):
    """
    Identification of HLA class I peptides from diaPASEF runs of a sample with a known ground truth.

    For spectral-library searches, upload the library as the secondary file so that r is taken from
    the library.

    Parameters
    ----------
    token : str
        GitHub token for the user.
    proteobot_repo_name : str, optional
        Repository for pull requests, by default "Proteobot/Results_immuno_HLAI_DIA_diaPASEF".
    proteobench_repo_name : str, optional
        Repository with the public results, by default "Proteobench/Results_immuno_HLAI_DIA_diaPASEF".
    branch : Optional[str], optional
        Branch of the results repository to check out.

    Attributes
    ----------
    module_id : str
        Module identifier for configuration.
    """

    module_id: str = "immuno_HLAI_DIA_diaPASEF"

    def __init__(
        self,
        token: str,
        proteobot_repo_name: str = "Proteobot/Results_immuno_HLAI_DIA_diaPASEF",
        proteobench_repo_name: str = "Proteobench/Results_immuno_HLAI_DIA_diaPASEF",
        branch: Optional[str] = None,
    ):
        """
        Initialize the module.

        Parameters
        ----------
        token : str
            GitHub token for the user.
        proteobot_repo_name : str, optional
            Repository for pull requests.
        proteobench_repo_name : str, optional
            Repository with the public results.
        branch : Optional[str], optional
            Branch of the results repository to check out.
        """
        super().__init__(
            token,
            proteobot_repo_name=proteobot_repo_name,
            proteobench_repo_name=proteobench_repo_name,
            parse_settings_dir=MODULE_SETTINGS_DIRS[self.module_id],
            module_id=self.module_id,
            branch=branch,
        )
