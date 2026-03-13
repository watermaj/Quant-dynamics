"""Input/output helpers for datasets and result files."""

from pathlib import Path

import pandas as pd


def save_results_to_csv(results, name, array_cols=None, output_dir="results/Submitted_results"):
    """
    Save flattened fit results to CSV, exploding array columns into rows.

    Args:
        results (list): Nested list of fit-result dictionaries.
        name (str): Output filename (without extension).
        array_cols (list or None): Columns containing arrays to expand.
        output_dir (str): Output directory for result files.

    Returns:
        None: Writes a CSV file.
    """
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)

    path = output_path / f"{name}.csv"

    if array_cols is None:
        array_cols = ["time", "emission", "prior", "fitted_emission"]
    flat = [d for group in results for d in group]
    df = pd.DataFrame(flat)
    for col in array_cols:
        if col in df.columns:
            df[col] = df[col].apply(lambda arr: arr.tolist())
    df = df.explode(array_cols, ignore_index=True)
    df.to_csv(path, index=False, encoding="utf-8-sig")
    print(f"Saved {len(df)} rows to {path}")
