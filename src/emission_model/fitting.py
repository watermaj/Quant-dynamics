"""Fitting utilities for single and triple response models."""

import copy

import numpy as np
from scipy.optimize import differential_evolution, minimize

from .model import duration, model_gamma, shape, total_integral_model_gamma


def run_de_then_slsqp(objective, bounds, constraints):
    """
    Run differential evolution followed by SLSQP for robust optimization.

    Args:
        objective (callable): Objective function to minimize.
        bounds (iterable): Parameter bounds.
        constraints (iterable): Constraints for the SLSQP stage.

    Returns:
        tuple: Optimized parameters and success flag.
    """

    try:
        # --- Stage 1: Global optimization ---
        de_result = differential_evolution(
            objective, bounds=bounds, strategy="best1bin", maxiter=50000, polish=False, seed=42
        )

        # --- Stage 2: Local refinement ---
        local_result = minimize(
            objective,
            de_result.x,
            method="SLSQP",
            bounds=bounds,
            constraints=constraints,
            options={"maxiter": 50000},
        )

        return local_result.x if local_result.success else de_result.x, True

    except Exception as e:
        print(f"Optimization error: {e}")
        return None, False


def r_squared(emission, fit):
    """
    Compute the coefficient of determination between data and fit.

    Args:
        emission (array): Observed values.
        fit (array): Fitted values.

    Returns:
        float: R squared value.
    """

    ss_tot = np.sum((emission - emission.mean()) ** 2)
    ss_res = np.sum((emission - fit) ** 2)
    r2 = 1 - ss_res / ss_tot if ss_tot > 0 else 0
    return r2


def process_single_group(data, type_value, comp_value, param_config, r_quared_limit=0.0, lambda_prior=1):
    """
    Fit a single gamma emission curve for one type-compound group.

    Args:
        data (DataFrame): Input dataset containing emission channels.
        type_value (str): Type label to filter.
        comp_value (str): Compound label to filter.
        param_config (dict): Parameter bounds.
        r_quared_limit (float): Minimum R^2 threshold for accepting fits.
        lambda_prior (float): Weight of the prior-curve penalty.

    Returns:
        list: Fitting results and metadata for accepted channels.
    """

    # Filter data for the specified type and comp.
    data_group = data[(data["Type"] == type_value) & (data["comp"] == comp_value)]
    results = []

    def fit_emission_curve(time, emission, prior, param_config):
        param_config = copy.deepcopy(param_config)
        t_peak_data = time[np.argmax(emission)]
        em_max = emission.max()
        mask = time >= 0

        if "t_peak" not in param_config:
            param_config["t_peak"] = {"min": 0.8 * t_peak_data, "max": 1.2 * t_peak_data}

        if "R_peak" not in param_config:
            param_config["R_peak"] = {"min": 0.8 * em_max, "max": 1.2 * em_max}

        param_names = ["R_peak", "t_peak", "t_onset", "t_mean"]
        bounds_list = [(param_config[n]["min"], param_config[n]["max"]) for n in param_names]

        def objective(x):
            params = dict(zip(param_names, x))
            model_pred = model_gamma(time[mask], **params)
            loss_prior = np.sum((model_pred - prior[mask]) ** 2)
            return np.sum((model_pred - emission[mask]) ** 2) + lambda_prior * loss_prior

        constraints = [
            {"type": "ineq", "fun": lambda x: (x[3] - x[1]) - 0.1 * (x[3] - x[2])},
            {"type": "ineq", "fun": lambda x: x[1] - x[2] - 0.1},
            {"type": "ineq", "fun": lambda x: x[3] - x[1] - 0.1},
        ]

        result_x, success = run_de_then_slsqp(objective, bounds_list, constraints)
        if not success:
            return None

        # --- Postprocessing ---
        fitted_params = dict(zip(param_names, result_x))
        fitted_emission = model_gamma(time, **fitted_params)

        integral = total_integral_model_gamma(**fitted_params)
        shape_val = shape(fitted_params["t_onset"], fitted_params["t_peak"], fitted_params["t_mean"])
        duration_val = duration(fitted_params["t_onset"], fitted_params["t_peak"], fitted_params["t_mean"])

        return {
            **{f"{k}_fit": v for k, v in fitted_params.items()},
            "duration": duration_val,
            "shape": shape_val,
            "fitted_emission": fitted_emission,
            "prior": prior,
            "time": time,
            "emission": emission,
            "integral": integral,
        }

    # --- Loop over channels ---
    for channel in data_group["Channel_number"].unique():
        data_channel = data_group[data_group["Channel_number"] == channel]
        time = data_channel["time"].values
        emission = data_channel["Emission_RAW"].values
        prior = data_channel["Emission"].values

        try:
            fit_results = fit_emission_curve(time, emission, prior, param_config)
        except Exception as e:
            print(f"Fitting failed for channel {channel}: {e}")
            continue

        if fit_results is None:
            print(f"Fitting failed for channel {channel}: optimizer returned no result.")
            continue

        # R^2
        r_2 = r_squared(emission, fit_results["fitted_emission"])
        print(f"Channel {channel}, Type {type_value}, R^2 = {r_2:.3f}")
        if r_2 < r_quared_limit:
            print(f"Channel {channel} not accepted (R^2 below threshold).")
            continue

        if np.isnan(fit_results["integral"]):
            print("Curve fit failed (NaN integral).")
            continue

        fit_results["Type"] = type_value
        fit_results["comp"] = comp_value
        fit_results["Channel_number"] = channel
        fit_results["intensity1"] = data_channel["intensity1"].iloc[0]
        fit_results["Totalbio"] = data_channel["Totalbio"].iloc[0]
        fit_results["Leaf3"] = data_channel["Leaf3"].iloc[0]
        fit_results["d1_time"] = data_channel["d1_time"].iloc[0]

        results.append(fit_results)

    return results


def process_single_group_three(data, type_value, comp_value, param_config, r_squared_limit=0.0, lambda_prior=1):
    """
    Fit three gamma emission curves for a specific type and compound group.

    Args:
        data (DataFrame): Input dataset containing time-series emission data.
        type_value (str): Type label to filter.
        comp_value (str): Compound label to filter.
        param_config (dict): Parameter configuration with bounds.
        r_squared_limit (float): Minimum acceptable R^2 for accepting a fit.
        lambda_prior (float): Weight of the prior-curve loss term.

    Returns:
        list: Fitted results and metadata for all accepted channels.
    """

    data_group = data[(data["Type"] == type_value) & (data["comp"] == comp_value)]
    results = []

    def fit_emission_curve_three(time, emission_raw, emission_prior, param_config, delta_t1, delta_t2):
        cfg = copy.deepcopy(param_config)
        mask = time >= 0

        # ------------------------------------------------------------
        # Apply damage-time shifts to second and third gamma curves
        # ------------------------------------------------------------
        # Curve 2
        for name in ("t_onset2", "t_peak2", "t_mean2"):
            cfg[name]["min"] += delta_t1
            cfg[name]["max"] += delta_t1

        # Curve 3
        for name in ("t_onset3", "t_peak3", "t_mean3"):
            cfg[name]["min"] += delta_t1 + delta_t2
            cfg[name]["max"] += delta_t1 + delta_t2

        # ------------------------------------------------------------
        # Parameter vector definition
        # ------------------------------------------------------------
        param_names = [
            "R_peak",
            "t_peak",
            "t_onset",
            "t_mean",
            "R_peak2",
            "t_peak2",
            "t_onset2",
            "t_mean2",
            "R_peak3",
            "t_peak3",
            "t_onset3",
            "t_mean3",
        ]

        bounds_list = []
        for pname in param_names:
            entry = cfg[pname]
            bounds_list.append((entry["min"], entry["max"]))

        # ------------------------------------------------------------
        # Objective
        # ------------------------------------------------------------
        def objective(parameters):
            (E1, tp1, ton1, tm1, E2, tp2, ton2, tm2, E3, tp3, ton3, tm3) = parameters

            model1 = model_gamma(time[mask], E1, tp1, ton1, tm1)
            model2 = model_gamma(time[mask], E2, tp2, ton2, tm2)
            model3 = model_gamma(time[mask], E3, tp3, ton3, tm3)
            model_total = model1 + model2 + model3

            loss_fit = np.sum((model_total - emission_raw[mask]) ** 2)
            loss_prior = np.sum((model_total - emission_prior[mask]) ** 2)

            return loss_fit + lambda_prior * loss_prior

        # ------------------------------------------------------------
        # Constraints: monotonic order of onset, peak, mean
        # ------------------------------------------------------------
        constraints = [
            {"type": "ineq", "fun": lambda x: x[6] - x[2]},  # t_onset2 >= t_onset1
            {"type": "ineq", "fun": lambda x: x[10] - x[6]},  # t_onset3 >= t_onset2
            {"type": "ineq", "fun": lambda x: x[5] - x[1]},  # t_peak2 >= t_peak1
            {"type": "ineq", "fun": lambda x: x[9] - x[5]},  # t_peak3 >= t_peak2
            {"type": "ineq", "fun": lambda x: x[7] - x[3]},  # t_mean2 >= t_mean1
            {"type": "ineq", "fun": lambda x: x[11] - x[7]},  # t_mean3 >= t_mean2
        ]

        # ------------------------------------------------------------
        # Optimization
        # ------------------------------------------------------------
        best_parameters, success = run_de_then_slsqp(objective, bounds_list, constraints)
        if not success:
            return None

        # ------------------------------------------------------------
        # Unpack fitted parameters
        # ------------------------------------------------------------
        fit_parameters = {p: best_parameters[i] for i, p in enumerate(param_names)}

        # Model curves
        model_curve_1 = model_gamma(
            time,
            fit_parameters["R_peak"],
            fit_parameters["t_peak"],
            fit_parameters["t_onset"],
            fit_parameters["t_mean"],
        )
        model_curve_2 = model_gamma(
            time,
            fit_parameters["R_peak2"],
            fit_parameters["t_peak2"],
            fit_parameters["t_onset2"],
            fit_parameters["t_mean2"],
        )
        model_curve_3 = model_gamma(
            time,
            fit_parameters["R_peak3"],
            fit_parameters["t_peak3"],
            fit_parameters["t_onset3"],
            fit_parameters["t_mean3"],
        )

        # Integrals
        integral1 = total_integral_model_gamma(
            fit_parameters["R_peak"],
            fit_parameters["t_peak"],
            fit_parameters["t_onset"],
            fit_parameters["t_mean"],
        )
        integral2 = total_integral_model_gamma(
            fit_parameters["R_peak2"],
            fit_parameters["t_peak2"],
            fit_parameters["t_onset2"],
            fit_parameters["t_mean2"],
        )
        integral3 = total_integral_model_gamma(
            fit_parameters["R_peak3"],
            fit_parameters["t_peak3"],
            fit_parameters["t_onset3"],
            fit_parameters["t_mean3"],
        )

        # Durations
        duration1 = duration(
            fit_parameters["t_onset"], fit_parameters["t_peak"], fit_parameters["t_mean"]
        )
        duration2 = duration(
            fit_parameters["t_onset2"], fit_parameters["t_peak2"], fit_parameters["t_mean2"]
        )
        duration3 = duration(
            fit_parameters["t_onset3"], fit_parameters["t_peak3"], fit_parameters["t_mean3"]
        )

        # Shapes
        shape1 = shape(fit_parameters["t_onset"], fit_parameters["t_peak"], fit_parameters["t_mean"])
        shape2 = shape(fit_parameters["t_onset2"], fit_parameters["t_peak2"], fit_parameters["t_mean2"])
        shape3 = shape(fit_parameters["t_onset3"], fit_parameters["t_peak3"], fit_parameters["t_mean3"])

        return {
            **{f"{k}_fit": v for k, v in fit_parameters.items()},
            "fitted_emission1": model_curve_1,
            "fitted_emission2": model_curve_2,
            "fitted_emission3": model_curve_3,
            "fitted_emission": model_curve_1 + model_curve_2 + model_curve_3,
            "prior": emission_prior,
            "time": time,
            "emission": emission_raw,
            "integral1": integral1,
            "integral2": integral2,
            "integral3": integral3,
            "duration1": duration1,
            "duration2": duration2,
            "duration3": duration3,
            "shape1": shape1,
            "shape2": shape2,
            "shape3": shape3,
        }

    # ------------------------------------------------------------
    # Loop over channels
    # ------------------------------------------------------------
    for channel in data_group["Channel_number"].unique():
        channel_data = data_group[data_group["Channel_number"] == channel]

        time = channel_data["time"].values
        emission_raw = channel_data["Emission_RAW"].values
        emission_prior = channel_data["Emission"].values
        delta_t1 = channel_data["dt1"].iloc[0]
        delta_t2 = channel_data["dt2"].iloc[0]

        try:
            fit = fit_emission_curve_three(
                time, emission_raw, emission_prior, param_config, delta_t1, delta_t2
            )

        except Exception as e:
            print(f"Fitting failed for channel {channel}: {e}")
            continue

        if fit is None:
            print(f"Fitting failed for channel {channel}: optimizer returned no result.")
            continue

        # R^2
        r2 = r_squared(emission_raw, fit["fitted_emission"])
        print(f"Channel {channel}, Type {type_value}, R^2 = {r2:.3f}")
        if r2 < r_squared_limit:
            print(f"Channel {channel} not accepted (R^2 below threshold).")
            continue

        # Integral sanity check
        if any(np.isnan([fit["integral1"], fit["integral2"], fit["integral3"]])):
            print("Curve fit failed (NaN integral).")
            continue

        # ------------------------------------------------------------
        # Attach metadata
        # ------------------------------------------------------------
        fit.update(
            {
                "Type": type_value,
                "comp": comp_value,
                "Channel_number": channel,
                "intensity1": channel_data["intensity1"].iloc[0],
                "intensity2": channel_data["intensity2"].iloc[0],
                "intensity3": channel_data["intensity3"].iloc[0],
                "Totalbio": channel_data["Totalbio"].iloc[0],
                "Leaf3": channel_data["Leaf3"].iloc[0],
                "d1_time": channel_data["d1_time"].iloc[0],
                "d2_time": channel_data["d2_time"].iloc[0],
                "d3_time": channel_data["d3_time"].iloc[0],
                "dt1": delta_t1,
                "dt2": delta_t2,
            }
        )

        results.append(fit)

    return results
