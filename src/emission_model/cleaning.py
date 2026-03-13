"""Data cleaning and reconstruction helpers extracted from the notebook."""

import numpy as np
import pandas as pd
from scipy.optimize import minimize
from sklearn.preprocessing import StandardScaler


def subtract_background(data, type="a"):
    """
    Subtract background channels from emission data using per-comp scaling.

    Args:
        data (DataFrame): Input dataframe containing emission measurements.
            Must include an `Emission_RAW` column (created during reconstruction).
        type (str): Label identifying background rows to subtract.

    Returns:
        DataFrame: Data with background-subtracted emission values.
    """

    processed_frames = []
    grouped_data = data.groupby("comp")

    for comp, group_data in grouped_data:
        background_data = group_data[group_data["Type"] == type]

        if not background_data.empty:
            # Pivot background
            bg_emission = background_data.pivot(index="time", columns="Channel_number", values="Emission")
            bg_mean = np.mean(bg_emission, axis=1)

            bg_raw = background_data.pivot(index="time", columns="Channel_number", values="Emission_RAW")
            bg_raw_mean = np.mean(bg_raw, axis=1)

            non_bg_data = group_data[group_data["Type"] != type].copy()
            non_bg_emission = non_bg_data.pivot(index="time", columns="Channel_number", values="Emission")
            non_bg_raw = non_bg_data.pivot(index="time", columns="Channel_number", values="Emission_RAW")

            adjusted_em = pd.DataFrame(index=non_bg_emission.index, columns=non_bg_emission.columns)
            adjusted_raw = pd.DataFrame(index=non_bg_raw.index, columns=non_bg_raw.columns)

            for ch in non_bg_emission.columns:
                e_curve = non_bg_emission[ch]
                bg_curve = bg_mean

                def residuals(scale, e, b):
                    idx = e.index[e.index <= 0]
                    return np.nansum((e.loc[idx] - b.loc[idx] * scale) ** 2)

                result = minimize(residuals, 1.0, args=(e_curve, bg_curve), method="Nelder-Mead")
                scale = result.x[0]

                adjusted_em[ch] = e_curve - bg_curve * scale

                raw_curve = non_bg_raw[ch]
                adjusted_raw[ch] = raw_curve - bg_raw_mean * scale

            adj_df = adjusted_em.reset_index().melt(
                id_vars="time", var_name="Channel_number", value_name="Emission"
            )
            meta = non_bg_data.drop(columns=["Emission", "Emission_RAW"])
            adj_df = adj_df.merge(meta, on=["time", "Channel_number"], how="left")

            raw_df = adjusted_raw.reset_index().melt(
                id_vars="time", var_name="Channel_number", value_name="Emission_RAW"
            )
            adj_df = adj_df.merge(raw_df, on=["time", "Channel_number"], how="left")

            adjusted_group = pd.concat([adj_df, background_data], ignore_index=True)
        else:
            adjusted_group = group_data

        processed_frames.append(adjusted_group)

    return pd.concat(processed_frames, ignore_index=True)


def reconstruct_with_svd(emission_data):
    """
    Reconstruct emission signals using the strongest SVD component.

    Args:
        emission_data (DataFrame): Long-format emission data with time and channel columns.

    Returns:
        DataFrame: Reconstructed emissions merged with original raw values.
    """
    # Reshape data: Time as rows, Channels as columns
    emission_matrix = emission_data.pivot(index="time", columns="Channel_number", values="Emission")

    # Standardize the data
    scaler = StandardScaler()
    standardized_matrix = scaler.fit_transform(emission_matrix.fillna(0))

    # Perform SVD decomposition
    U, S, Vt = np.linalg.svd(standardized_matrix, full_matrices=False)

    # Reconstruct with the strongest singular component
    strongest_component = np.outer(U[:, 0], S[0] * Vt[0, :])
    reconstructed_matrix = scaler.inverse_transform(strongest_component)

    # Convert back to long-format DataFrame
    reconstructed_df = pd.DataFrame(
        reconstructed_matrix, index=emission_matrix.index, columns=emission_matrix.columns
    )
    reconstructed_long = reconstructed_df.reset_index().melt(
        id_vars="time", var_name="Channel_number", value_name="Emission"
    )

    # Merge with original data to get Emission_RAW
    merged = emission_data.copy()
    merged = merged.rename(columns={"Emission": "Emission_RAW"})
    reconstructed_long = reconstructed_long.merge(merged, on=["time", "Channel_number"])

    return reconstructed_long


def clean_and_reconstruct(data, background=True, norm=None):
    """
    Clean emission data, reconstruct signals with SVD, and optionally subtract background.

    Args:
        data (DataFrame): Input emission dataset.
        background (bool): Whether to subtract background after reconstruction.
        norm (str or None): Column name used to normalize emissions.

    Returns:
        DataFrame: Cleaned and reconstructed emission data.
    """

    processed_frames = []

    if norm:
        data["Emission"] = data["Emission"] / data[norm]

    grouped_data = data.groupby(["Type", "comp"])

    for (type_group, comp_group), group_data in grouped_data:
        reconstructed_emissions = reconstruct_with_svd(group_data[["time", "Channel_number", "Emission"]])

        reconstructed_group = group_data.drop(columns="Emission").merge(
            reconstructed_emissions, on=["time", "Channel_number"], how="left"
        )

        processed_frames.append(reconstructed_group)

    cleaned_data = pd.concat(processed_frames, ignore_index=True)

    if "Unnamed: 0" in cleaned_data.columns:
        cleaned_data = cleaned_data.drop(columns=["Unnamed: 0"])

    if background:
        final_data = subtract_background(cleaned_data)
    else:
        final_data = cleaned_data

    return final_data
