"""Emission model package."""

from .cleaning import clean_and_reconstruct, reconstruct_with_svd, subtract_background
from .configs import (
    param_config_CYP92C5_gene_expression,
    param_config_DMNT,
    param_config_DMNT_single_os_adjusted,
    param_config_DMNT_triple_damage,
    param_config_HAC_hires,
    param_config_hexa_hires,
    param_config_hexo_hires,
    param_config_IGL_gene_expression,
    param_config_indole,
    param_config_mono,
    param_config_sesq,
    param_config_TMTT,
    param_config_TPS10_gene_expression,
    param_config_TPS2_gene_expression,
)
from .fitting import process_single_group, process_single_group_three, r_squared, run_de_then_slsqp
from .io import save_results_to_csv
from .model import duration, model_gamma, shape, total_integral_model_gamma

__all__ = [
    "clean_and_reconstruct",
    "reconstruct_with_svd",
    "subtract_background",
    "param_config_CYP92C5_gene_expression",
    "param_config_DMNT",
    "param_config_DMNT_single_os_adjusted",
    "param_config_DMNT_triple_damage",
    "param_config_HAC_hires",
    "param_config_hexa_hires",
    "param_config_hexo_hires",
    "param_config_IGL_gene_expression",
    "param_config_indole",
    "param_config_mono",
    "param_config_sesq",
    "param_config_TMTT",
    "param_config_TPS10_gene_expression",
    "param_config_TPS2_gene_expression",
    "process_single_group",
    "process_single_group_three",
    "r_squared",
    "run_de_then_slsqp",
    "save_results_to_csv",
    "duration",
    "model_gamma",
    "shape",
    "total_integral_model_gamma",
]
