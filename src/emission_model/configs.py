"""Parameter configuration definitions extracted from the notebook workflow."""

import copy

# Base single-curve configs reused across multiple analyses.
param_config_DMNT = {
    # 'R_peak': {'min': 5, 'max':  60 },
    # 't_peak': {'min': 1.5, 'max':  3.5 },
    "t_onset": {"min": 0, "max": 10},
    "t_mean": {"min": 0, "max": 20},
}

param_config_indole = {
    # 'R_peak': {'min': 30, 'max':  90 },
    # 't_peak': {'min': 2, 'max':  4 },
    "t_onset": {"min": 0, "max": 10},
    "t_mean": {"min": 1, "max": 20},
}

param_config_TMTT = {
    # 'R_peak': {'min': 2, 'max':  10 },
    # 't_peak': {'min': 6, 'max':  10 },
    "t_onset": {"min": 0, "max": 10},
    "t_mean": {"min": 1, "max": 20},
}

param_config_sesq = {
    # 'R_peak': {'min': 10, 'max':  50 },
    # 't_peak': {'min': 6, 'max':  9 },
    "t_onset": {"min": 0, "max": 10},
    "t_mean": {"min": 1, "max": 20},
}

param_config_mono = {
    # 'R_peak': {'min': 25, 'max':  150 },
    # 't_peak': {'min': 3.1, 'max':  8 },
    "t_onset": {"min": 0, "max": 10},
    "t_mean": {"min": 1, "max": 20},
}

# Triple-damage config aligned with the "New" notebook.
param_config_DMNT_triple_damage = {
    # Curve 1
    "R_peak": {"min": 5, "max": 150},
    "t_peak": {"min": 1, "max": 4},
    "t_onset": {"min": 0, "max": 10},
    "t_mean": {"min": 1, "max": 20},
    # Curve 2
    "R_peak2": {"min": 5, "max": 150},
    "t_peak2": {"min": 1, "max": 4},
    "t_onset2": {"min": 0, "max": 10},
    "t_mean2": {"min": 1, "max": 20},
    # Curve 3
    "R_peak3": {"min": 5, "max": 150},
    "t_peak3": {"min": 1, "max": 4},
    "t_onset3": {"min": 0, "max": 10},
    "t_mean3": {"min": 1, "max": 20},
}

# GLV high-resolution configs.
param_config_HAC_hires = {
    # 'R_peak': {'min': 10, 'max':  100 },
    # 't_peak': {'min': 4, 'max':  15 },
    "t_onset": {"min": 0, "max": 20},
    "t_mean": {"min": 1, "max": 40},
}

param_config_hexo_hires = {
    # 'R_peak': {'min': 60, 'max':  110 },
    # 't_peak': {'min': 4, 'max':  8 },
    "t_onset": {"min": 0.0, "max": 20},
    "t_mean": {"min": 1, "max": 40},
}

param_config_hexa_hires = {
    # 'R_peak': {'min': 5000, 'max':  8000 },
    # 't_peak': {'min': 2, 'max':  6 },
    "t_onset": {"min": 0.0, "max": 1},
    "t_mean": {"min": 4, "max": 15},
}

# Gene-expression configs.
param_config_gene_expression_common = {
    # 'R_peak': {'min': 0.02, 'max':  0.06 },
    # 't_peak': {'min': 0.21, 'max':  0.6 },
    "t_onset": {"min": 0, "max": 2},
    "t_mean": {"min": 0.1, "max": 5},
}

param_config_CYP92C5_gene_expression = copy.deepcopy(param_config_gene_expression_common)
param_config_IGL_gene_expression = copy.deepcopy(param_config_gene_expression_common)

param_config_TPS10_gene_expression = {
    # 'R_peak': {'min': 3, 'max':  6 },
    # 't_peak': {'min': 2, 'max':  5 },
    "t_onset": {"min": 0, "max": 2},
    "t_mean": {"min": 0.1, "max": 5},
}

param_config_TPS2_gene_expression = {
    # 'R_peak': {'min': 0.4, 'max':  1.8 },
    # 't_peak': {'min': 0.8, 'max':  2 },
    "t_onset": {"min": 0, "max": 2},
    "t_mean": {"min": 0.1, "max": 5},
}

# DMNT configs used by the single-OS sensitivity analyses.
param_config_DMNT_single_os_adjusted = {
    # 'R_peak': {'min': 1, 'max':  600 },
    # 't_peak': {'min': 1, 'max':  15 },
    "t_onset": {"min": 0, "max": 10},
    "t_mean": {"min": 1, "max": 20},
    # 't_mean': {'min': 1, 'max': 18 }
}


def get_default_param_configs():
    """Return deep-copied default per-compound parameter configs."""
    return {
        "DMNT": copy.deepcopy(param_config_DMNT),
        "indole": copy.deepcopy(param_config_indole),
        "TMTT": copy.deepcopy(param_config_TMTT),
        "sesq": copy.deepcopy(param_config_sesq),
        "mono": copy.deepcopy(param_config_mono),
    }
