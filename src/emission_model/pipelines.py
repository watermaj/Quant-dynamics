"""High-level analysis orchestration extracted from the notebook workflow."""

from pathlib import Path
import copy

import pandas as pd

from .cleaning import clean_and_reconstruct
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
from .fitting import process_single_group, process_single_group_three
from .io import save_results_to_csv


def run_dmnt_dose(data_dir="data", output_dir="results/Submitted_results"):
    data_dose = clean_and_reconstruct(pd.read_csv(Path(data_dir) / "singledose_sub.csv"), norm="Leaf3")

    param_config = copy.deepcopy(param_config_DMNT)

    comp = "DMNT"
    res = []
    data = data_dose
    for type_ in sorted(data.Type.unique()):
        if type_ not in ["a"]:
            print(type_)
            individual_results = process_single_group(data, type_, comp, param_config, lambda_prior=1)
            res.append(individual_results)

    save_results_to_csv(res, "Results_DMNT_dose", output_dir=output_dir)


def run_dmnt_time(data_dir="data", output_dir="results/Submitted_results"):
    data_cir = clean_and_reconstruct(pd.read_csv(Path(data_dir) / "time_of_day_sub.csv"), norm="Leaf3")

    param_config = copy.deepcopy(param_config_DMNT)

    comp = "DMNT"
    res = []
    data = data_cir
    for type_ in sorted(data.Type.unique()):
        if type_ not in ["a"]:
            print(type_)

            # Adjust onset/mean bounds for treatment-time shift
            data_group = data[(data["Type"] == type_) & (data["comp"] == comp)]
            time_shift = data_group["d1_time"].iloc[0] - data["d1_time"].min()

            param_local = copy.deepcopy(param_config)

            for key in ["t_onset", "t_mean"]:
                param_local[key]["min"] = param_config[key]["min"] + time_shift
                param_local[key]["max"] = param_config[key]["max"] + time_shift

            individual_results = process_single_group(data, type_, comp, param_local, lambda_prior=1)
            res.append(individual_results)

    save_results_to_csv(res, "Results_DMNT_time", output_dir=output_dir)


def run_dmnt_leaves(data_dir="data", output_dir="results/Submitted_results"):
    leaves = clean_and_reconstruct(pd.read_csv(Path(data_dir) / "leaves_sub.csv"), norm=None)

    param_config = copy.deepcopy(param_config_DMNT)

    comp = "DMNT"
    res = []
    data = leaves
    for type_ in sorted(data.Type.unique()):
        if type_ in ["l2", "l3", "l4"]:
            print(type_)
            individual_results = process_single_group(data, type_, comp, param_config, lambda_prior=2)
            res.append(individual_results)

    save_results_to_csv(res, "Results_DMNT_leaves", output_dir=output_dir)


def run_dmnt_genotype(data_dir="data", output_dir="results/Submitted_results"):
    dirty_data = pd.read_csv(Path(data_dir) / "genotype_sub.csv")

    # Remove genotype specific backgrounds
    CML287 = dirty_data[dirty_data["gt"] == "CML287"]
    CML287 = clean_and_reconstruct(CML287, norm="Leaf3")
    MO17 = dirty_data[dirty_data["gt"] == "MO17"]
    MO17 = clean_and_reconstruct(MO17, norm="Leaf3")
    NC300 = dirty_data[dirty_data["gt"] == "NC3000"]
    NC300 = clean_and_reconstruct(NC300, norm="Leaf3")
    # Rejoin into a single dataframe
    geno_single = pd.concat([CML287, MO17, NC300], ignore_index=True)

    param_config = copy.deepcopy(param_config_DMNT)

    comp = "DMNT"
    res = []
    data = geno_single
    for type_ in sorted(data.Type.unique()):
        if type_ in ["CML287", "MO17", "NC3000"]:
            print(type_)
            individual_results = process_single_group(data, type_, comp, param_config, lambda_prior=1)
            res.append(individual_results)

    save_results_to_csv(res, "Results_DMNT_Genotype", output_dir=output_dir)


def run_single_os_compounds(data_dir="data", output_dir="results/Submitted_results"):
    os_with = clean_and_reconstruct(pd.read_csv(Path(data_dir) / "os_with_sub.csv"), norm="Leaf3")
    os_without = clean_and_reconstruct(pd.read_csv(Path(data_dir) / "os_without_sub.csv"), norm="Leaf3")
    os_single = pd.concat([os_with, os_without], ignore_index=True)

    param_config = copy.deepcopy(param_config_DMNT)
    comp = "DMNT"
    res = []
    data = os_single
    for type_ in sorted(data.Type.unique()):
        if type_ in ["wo", "w"]:
            print(type_)
            individual_results = process_single_group(data, type_, comp, param_config, lambda_prior=1)
            res.append(individual_results)
    save_results_to_csv(res, "Results_DMNT_Single_OS", output_dir=output_dir)

    param_config = copy.deepcopy(param_config_indole)
    comp = "indole"
    res = []
    data = os_single
    for type_ in sorted(data.Type.unique()):
        if type_ in ["wo", "w"]:
            print(type_)
            print()
            individual_results = process_single_group(data, type_, comp, param_config, lambda_prior=1)
            res.append(individual_results)
    save_results_to_csv(res, "Results_indole_Single_OS", output_dir=output_dir)

    param_config = copy.deepcopy(param_config_TMTT)
    comp = "TMTT"
    res = []
    data = os_single
    for type_ in sorted(data.Type.unique()):
        if type_ in ["wo", "w"]:
            print(type_)
            individual_results = process_single_group(data, type_, comp, param_config, lambda_prior=1)
            res.append(individual_results)
    save_results_to_csv(res, "Results_TMTT_Single_OS", output_dir=output_dir)

    param_config = copy.deepcopy(param_config_sesq)
    comp = "sesq"
    res = []
    data = os_single
    for type_ in sorted(data.Type.unique()):
        if type_ in ["wo", "w"]:
            print(type_)
            individual_results = process_single_group(data, type_, comp, param_config, lambda_prior=1)
            res.append(individual_results)
    save_results_to_csv(res, "Results_sesq_Single_OS", output_dir=output_dir)

    param_config = copy.deepcopy(param_config_mono)
    comp = "mono"
    res = []
    data = os_single
    for type_ in sorted(data.Type.unique()):
        if type_ in ["wo", "w"]:
            print(type_)
            individual_results = process_single_group(data, type_, comp, param_config, lambda_prior=1)
            res.append(individual_results)
    save_results_to_csv(res, "Results_mono_Single_OS", output_dir=output_dir)


def run_triple_damage(data_dir="data", output_dir="results/Submitted_results"):
    triple = clean_and_reconstruct(pd.read_csv(Path(data_dir) / "triple_sub.csv"), norm="Leaf3")

    param_config = copy.deepcopy(param_config_DMNT_triple_damage)
    comp = "DMNT"
    res = []
    data = triple
    for type_ in sorted(data.Type.unique()):
        if type_ not in ["a"]:
            print(type_)
            individual_results = process_single_group_three(data, type_, comp, param_config, lambda_prior=0)
            res.append(individual_results)

    array_cols = [
        "time",
        "emission",
        "prior",
        "fitted_emission1",
        "fitted_emission2",
        "fitted_emission3",
        "fitted_emission",
    ]
    save_results_to_csv(res, "Results_DMNT_trip", array_cols=array_cols, output_dir=output_dir)


def run_single_herbivore(data_dir="data", output_dir="results/Submitted_results"):
    single_herb = clean_and_reconstruct(pd.read_csv(Path(data_dir) / "herbreal_sub.csv"), norm=None)

    param_config = copy.deepcopy(param_config_DMNT)
    comp = "DMNT"
    res = []
    data = single_herb
    for type_ in sorted(data.Type.unique()):
        if type_ in ["b"]:
            print(type_)
            individual_results = process_single_group(data, type_, comp, param_config, lambda_prior=1)
            res.append(individual_results)
    save_results_to_csv(res, "Results_DMNT_single_herb", output_dir=output_dir)

    param_config = copy.deepcopy(param_config_indole)
    comp = "indole"
    res = []
    data = single_herb
    for type_ in sorted(data.Type.unique()):
        if type_ in ["b"]:
            print(type_)
            individual_results = process_single_group(data, type_, comp, param_config, lambda_prior=1)
            res.append(individual_results)
    save_results_to_csv(res, "Results_indole_single_herb", output_dir=output_dir)

    param_config = copy.deepcopy(param_config_TMTT)
    comp = "TMTT"
    res = []
    data = single_herb
    for type_ in sorted(data.Type.unique()):
        if type_ in ["b"]:
            print(type_)
            individual_results = process_single_group(data, type_, comp, param_config, lambda_prior=1)
            res.append(individual_results)
    save_results_to_csv(res, "Results_TMTT_single_herb", output_dir=output_dir)

    param_config = copy.deepcopy(param_config_sesq)
    comp = "sesq"
    res = []
    data = single_herb
    for type_ in sorted(data.Type.unique()):
        if type_ in ["b"]:
            print(type_)
            individual_results = process_single_group(data, type_, comp, param_config, lambda_prior=1)
            res.append(individual_results)
    save_results_to_csv(res, "Results_sesq_single_herb", output_dir=output_dir)

    param_config = copy.deepcopy(param_config_mono)
    comp = "mono"
    res = []
    data = single_herb
    for type_ in sorted(data.Type.unique()):
        if type_ in ["b"]:
            print(type_)
            individual_results = process_single_group(data, type_, comp, param_config, lambda_prior=1)
            res.append(individual_results)
    save_results_to_csv(res, "Results_mono_single_herb", output_dir=output_dir)


def run_glv_hires(data_dir="data", output_dir="results/Submitted_results"):
    dirty_data = pd.read_csv(Path(data_dir) / "glvkin_sub.csv")

    # Remove an obvious background from hexo
    dirty_data.loc[dirty_data["comp"] == "hexo", "Emission"] -= 5
    data_glvkin = clean_and_reconstruct(dirty_data, background=False, norm=None)

    param_config = copy.deepcopy(param_config_HAC_hires)

    comp = "HAC"
    res = []
    data = data_glvkin
    for type_ in sorted(data.Type.unique()):
        if type_ in ["b"]:
            print(type_)
            individual_results = process_single_group(data, type_, comp, param_config, lambda_prior=1)
            res.append(individual_results)
    save_results_to_csv(res, "Results_HAC_HiRes", output_dir=output_dir)

    param_config = copy.deepcopy(param_config_hexo_hires)

    comp = "hexo"
    res = []
    data = data_glvkin
    for type_ in sorted(data.Type.unique()):
        if type_ in ["b"]:
            print(type_)
            individual_results = process_single_group(data, type_, comp, param_config, lambda_prior=0)
            res.append(individual_results)
    save_results_to_csv(res, "Results_hexo_HiRes", output_dir=output_dir)

    param_config = copy.deepcopy(param_config_hexa_hires)

    comp = "hexa"
    res = []
    data = data_glvkin
    for type_ in sorted(data.Type.unique()):
        if type_ in ["b"]:
            print(type_)
            individual_results = process_single_group(data, type_, comp, param_config, lambda_prior=0)
            res.append(individual_results)
    save_results_to_csv(res, "Results_hexa_HiRes", output_dir=output_dir)


def run_gene_expression(data_dir="data", output_dir="results/Submitted_results"):
    gene = clean_and_reconstruct(pd.read_csv(Path(data_dir) / "genexpression_real_sub.csv"), norm=None)

    param_config = copy.deepcopy(param_config_CYP92C5_gene_expression)

    comp = "CYP92C5"
    res = []
    data = gene
    for type_ in sorted(data.Type.unique()):
        if type_ in ["b"]:
            print(type_)
            individual_results = process_single_group(data, type_, comp, param_config, lambda_prior=1)
            res.append(individual_results)
    save_results_to_csv(res, "Results_CYP92C5_GeneExpression", output_dir=output_dir)

    param_config = copy.deepcopy(param_config_IGL_gene_expression)

    comp = "IGL"
    res = []
    data = gene
    for type_ in sorted(data.Type.unique()):
        if type_ in ["b"]:
            print(type_)
            individual_results = process_single_group(data, type_, comp, param_config, lambda_prior=1)
            res.append(individual_results)
    save_results_to_csv(res, "Results_IGL_GeneExpression", output_dir=output_dir)

    param_config = copy.deepcopy(param_config_TPS10_gene_expression)

    comp = "TPS10"
    res = []
    data = gene
    for type_ in sorted(data.Type.unique()):
        if type_ in ["b"]:
            print(type_)
            individual_results = process_single_group(data, type_, comp, param_config, lambda_prior=1)
            res.append(individual_results)
    save_results_to_csv(res, "Results_TPS10_GeneExpression", output_dir=output_dir)

    param_config = copy.deepcopy(param_config_TPS2_gene_expression)

    comp = "TPS2"
    res = []
    data = gene
    for type_ in sorted(data.Type.unique()):
        if type_ in ["b"]:
            print(type_)
            individual_results = process_single_group(data, type_, comp, param_config, lambda_prior=1)
            res.append(individual_results)
    save_results_to_csv(res, "Results_TPS2_GeneExpression", output_dir=output_dir)


def run_os_time_cut(data_dir="data", output_dir="results/Submitted_results"):
    os_with = clean_and_reconstruct(pd.read_csv(Path(data_dir) / "os_with_sub.csv"), norm="Leaf3")

    new_dfs = []
    next_channel = os_with["Channel_number"].max() + 1

    filtered = os_with[os_with["Type"] != "a"]

    for perc in range(10, 89, 1):
        for ch, df_ch in filtered.groupby("Channel_number"):
            # Calculate how many rows to keep
            keep_n = int(len(df_ch) * (100 - perc) / 100)

            # Take only the first part of the data (truncate the end)
            df_trunc = df_ch.sort_values("time").iloc[:keep_n].copy()

            # Assign new channel number
            df_trunc["Channel_number"] = next_channel
            next_channel += 1

            # Modify intensity1 and Type
            df_trunc["intensity1"] = df_trunc["intensity1"] + perc
            df_trunc["Type"] = df_trunc["Type"] + str(perc)

            new_dfs.append(df_trunc)

    # Combine original + truncated copies
    os_shortened = pd.concat([filtered] + new_dfs, ignore_index=True)

    param_config = copy.deepcopy(param_config_indole)
    comp = "indole"
    res = []
    data = os_shortened
    for type_ in sorted(data.Type.unique()):
        if type_ not in ["a"]:
            print(type_)
            try:
                individual_results = process_single_group(data, type_, comp, param_config, lambda_prior=1)
                res.append(individual_results)
            except Exception:
                continue
    save_results_to_csv(res, "Results_indole_Single_OS_time_cut_1percent", output_dir=output_dir)

    param_config = copy.deepcopy(param_config_DMNT_single_os_adjusted)

    comp = "DMNT"
    res = []
    data = os_shortened
    for type_ in sorted(data.Type.unique()):
        if type_ not in ["a"]:
            print(type_)
            try:
                individual_results = process_single_group(data, type_, comp, param_config, lambda_prior=1)
                res.append(individual_results)
            except Exception:
                continue
    save_results_to_csv(res, "Results_DMNT_Single_OS_time_cut_1percent", output_dir=output_dir)

    param_config = copy.deepcopy(param_config_TMTT)
    comp = "TMTT"
    res = []
    data = os_shortened
    for type_ in sorted(data.Type.unique()):
        if type_ not in ["a"]:
            print(type_)
            try:
                individual_results = process_single_group(data, type_, comp, param_config, lambda_prior=1)
                res.append(individual_results)
            except Exception:
                continue
    save_results_to_csv(res, "Results_TMTT_Single_OS_time_cut_1percent", output_dir=output_dir)


def run_os_resolution_adjusted(data_dir="data", output_dir="results/Submitted_results"):
    os_with = clean_and_reconstruct(pd.read_csv(Path(data_dir) / "os_with_sub.csv"), norm="Leaf3")

    new_dfs = []
    next_channel = os_with["Channel_number"].max() + 1

    filtered = os_with[os_with["Type"] != "a"]

    for perc in range(10, 90, 5):  # 10%, 15%, ... 85% removed
        for ch, df_ch in filtered.groupby("Channel_number"):
            df_ch = df_ch.sort_values("time")

            # Calculate how many rows to keep
            keep_frac = (100 - perc) / 100
            step = int(round(1 / keep_frac))  # e.g. 50% keep -> step=2

            # Keep every nth row to simulate lower resolution
            df_down = df_ch.iloc[::step].copy()

            # Assign new channel number
            df_down["Channel_number"] = next_channel
            next_channel += 1

            # Modify intensity and Type
            df_down["intensity1"] = df_down["intensity1"] + perc
            df_down["Type"] = df_down["Type"] + str(perc)

            new_dfs.append(df_down)

    # Combine
    os_downsampled = pd.concat([filtered] + new_dfs, ignore_index=True)

    param_config = copy.deepcopy(param_config_indole)
    comp = "indole"
    res = []
    data = os_downsampled
    for type_ in sorted(data.Type.unique()):
        if type_ not in ["a"]:
            print(type_)
            try:
                individual_results = process_single_group(data, type_, comp, param_config, lambda_prior=1)
                res.append(individual_results)
            except Exception:
                continue
    save_results_to_csv(res, "Results_indole_Single_OS_resolution_adjusted", output_dir=output_dir)

    param_config = copy.deepcopy(param_config_DMNT_single_os_adjusted)

    comp = "DMNT"
    res = []
    data = os_downsampled
    for type_ in sorted(data.Type.unique()):
        if type_ not in ["r"]:
            print(type_)
            try:
                individual_results = process_single_group(data, type_, comp, param_config, lambda_prior=1)
                res.append(individual_results)
            except Exception:
                continue
    save_results_to_csv(res, "Results_DMNT_Single_OS_resolution_adjusted", output_dir=output_dir)

    param_config = copy.deepcopy(param_config_TMTT)
    comp = "TMTT"
    res = []
    data = os_downsampled
    for type_ in sorted(data.Type.unique()):
        if type_ not in ["a"]:
            print(type_)
            try:
                individual_results = process_single_group(data, type_, comp, param_config, lambda_prior=1)
                res.append(individual_results)
            except Exception:
                continue
    save_results_to_csv(res, "Results_TMTT_Single_OS_resolution_adjusted", output_dir=output_dir)


PIPELINE_STEPS = {
    "dmnt_dose": run_dmnt_dose,
    "dmnt_time": run_dmnt_time,
    "dmnt_leaves": run_dmnt_leaves,
    "dmnt_genotype": run_dmnt_genotype,
    "single_os": run_single_os_compounds,
    "triple_damage": run_triple_damage,
    "single_herbivore": run_single_herbivore,
    "glv_hires": run_glv_hires,
    "gene_expression": run_gene_expression,
    "os_time_cut": run_os_time_cut,
    "os_resolution_adjusted": run_os_resolution_adjusted,
}


def run_selected_analyses(selected=None, data_dir="data", output_dir="results/Submitted_results"):
    if selected is None:
        selected = list(PIPELINE_STEPS.keys())

    for name in selected:
        if name not in PIPELINE_STEPS:
            raise ValueError(f"Unknown analysis step: {name}")
        print(f"\n=== Running {name} ===")
        PIPELINE_STEPS[name](data_dir=data_dir, output_dir=output_dir)
