# Quantitative Modelling of Biological Response Dynamics

This repository contains the data, analysis code, and reference outputs for the paper:

**"Quantitative modelling of biological response dynamics reveals novel patterns in plant volatile signalling"**

Authors: Jamie M. Waterman*<sup>&dagger;1,2</sup>, Gareth J. Moore<sup>&dagger;3</sup>, Loren K. Amdahl-Culleton<sup>3</sup>, Sara Hoefer<sup>2</sup>, Matthias Erb<sup>2</sup>

Affiliations:
<sup>1</sup> Discipline of Botany, School of Natural Sciences, Trinity College Dublin, Dublin, Ireland<br>
<sup>2</sup> Institute of Plant Sciences, University of Bern, Bern, Switzerland<br>
<sup>3</sup> Independent

*Corresponding author. Email: watermaj@tcd.ie<br>
<sup>&dagger;</sup> Authors contributed equally to this work

The paper itself is included in this repository as:

- [Quantitative modelling of biological response dynamics reveals novel patterns in plant volatile signalling.pdf](./Quantitative%20modelling%20of%20biological%20response%20dynamics%20reveals%20novel%20patterns%20in%20plant%20volatile%20signalling.pdf)

## What This Project Is

This project studies how plant volatile emissions change over time after damage or treatment. The analysis fits mathematical response curves to measured time-series data and produces the submission result tables used in the paper.

In practical terms, this repository is set up so that someone can:

1. use the provided input data,
2. run the analysis from the command line,
3. regenerate the result files,
4. confirm that the regenerated files exactly match the reference outputs.

## Where To Start

If you are new to the repository, the simplest path is:

1. install the Python dependencies,
2. run the main analysis script,
3. compare the generated outputs with the reference outputs.

You do not need to use the notebook to reproduce the submitted results.

## Repository Layout

- `data/`
  Submission input datasets (`*_sub.csv`)
- `results/Submitted_results/`
  Generated output CSV files
- `tests/regression/reference/`
  Reference output CSV files used for exact comparison
- `scripts/run_analysis.py`
  Main script that runs the full analysis pipeline
- `src/emission_model/`
  Reusable analysis code extracted from the original workflow
- `tests/regression/compare_results.py`
  Exact file comparison tool

## Quick Setup

On Windows:

```bash
python -m venv .venv
.\.venv\Scripts\activate
pip install -r requirements.txt
```

## Run The Full Analysis

From the repository root, run:

```bash
python scripts/run_analysis.py
```

This will regenerate the result files in:

`results/Submitted_results/`

## Check That The Results Match

To confirm that the outputs are exactly the expected ones, run:

```bash
python tests/regression/compare_results.py "results/Submitted_results" "tests/regression/reference"
```

Expected result:

`No differences found.`

## Run One Part Of The Analysis

If you only want to run one analysis block, you can do that too. For example:

```bash
python scripts/run_analysis.py --analysis dmnt_dose
```

## For Readers Of The Paper

If your main goal is to understand the scientific context, start with the PDF paper. If your main goal is to reproduce the submitted computational outputs, start with the script runner and the regression comparison command above.

## Troubleshooting

- Run commands from the repository root.
- Make sure the same Python environment is used for both installing dependencies and running the script.
- If you rerun the pipeline, rerun the regression comparison as well.
