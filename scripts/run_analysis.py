"""Command-line runner for the extracted emission analysis pipeline."""

from __future__ import annotations

import argparse
from pathlib import Path
import sys
import warnings

# Allow running as `python scripts/run_analysis.py` without package install.
ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from emission_model.pipelines import PIPELINE_STEPS, run_selected_analyses

# Match notebook behavior (it suppresses warnings globally).
warnings.filterwarnings("ignore")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Run extracted analyses. "
            "By default runs all analyses that produce files in results/Submitted_results."
        )
    )
    parser.add_argument(
        "--analysis",
        action="append",
        dest="analyses",
        choices=sorted(PIPELINE_STEPS.keys()),
        help="Analysis step to run (repeat flag to run multiple). Default: all.",
    )
    parser.add_argument(
        "--config",
        default="configs/analyses.yaml",
        help=(
            "Optional YAML config listing analysis step names (simple '- step_name' list). "
            "Used when --analysis is not provided."
        ),
    )
    parser.add_argument("--data-dir", default="data", help="Directory containing *_sub.csv files.")
    parser.add_argument(
        "--output-dir",
        default="results/Submitted_results",
        help="Directory where result CSV files are written.",
    )
    return parser.parse_args()


def load_analysis_list(config_path: Path) -> list[str]:
    path = config_path
    if not path.exists():
        return []

    selected = []
    for raw in path.read_text(encoding="utf-8").splitlines():
        line = raw.strip()
        if not line or line.startswith("#"):
            continue
        if line.startswith("- "):
            selected.append(line[2:].strip())
    return selected


def main() -> int:
    args = parse_args()

    config_path = Path(args.config)
    if not config_path.is_absolute():
        config_path = ROOT / config_path

    data_dir = Path(args.data_dir)
    if not data_dir.is_absolute():
        data_dir = ROOT / data_dir

    output_dir = Path(args.output_dir)
    if not output_dir.is_absolute():
        output_dir = ROOT / output_dir

    selected = args.analyses
    if selected is None:
        from_config = load_analysis_list(config_path)
        selected = from_config if from_config else None

    run_selected_analyses(
        selected=selected,
        data_dir=str(data_dir),
        output_dir=str(output_dir),
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
