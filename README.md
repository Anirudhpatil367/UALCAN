# UALCAN Research Scripts Repository

This repository is a centralized location for the Python scripts created for UALCAN-related gene-expression and cancer-data analysis workflows.

For every Python script added to this repository, two small sample input files will also be included. These sample files will demonstrate the required input structure and allow users to understand and test the corresponding code without downloading the complete research dataset.

The repository will contain:

- Python scripts for data processing, metadata preparation, expression analysis, correlation analysis, and related workflows.
- Two representative sample input files for every Python script.
- Documentation describing each script, its required inputs, and its generated outputs.
- A shared Python dependency file.
- Tests and example commands where applicable.

Large datasets, complete generated results, confidential data, and patient-identifiable information will not be stored in this repository.

## Repository tree

```text
ualcan/
├── README.md
├── requirements.txt
├── .gitignore
│
├── scripts/
│   ├── cancer_vs_cancer.py
│   ├── GSE201284_all_feature_correlations.py
│   ├── GSE120741_all_feature_correlations.py
│   └── additional_scripts_will_be_added_here.py
│
├── samples/
│   ├── cancer_vs_cancer/
│   │   ├── GSE183019_series_matrix_sample.txt
│   │   └── GSE183019_processed_TPM_sample.txt
│   │
│   ├── GSE201284_all_feature_correlations/
│   │   ├── GSE201284_processed_TPM_sample.txt
│   │   └── GSE201284_sample_metadata_sample.tsv
│   │
│   ├── GSE120741_all_feature_correlations/
│   │   ├── GSE120741_Porto_ge_table_sample.txt
│   │   └── GSE120741_series_matrix_sample.txt
│   │
│   └── additional_script_name/
│       ├── sample_input_1.txt
│       └── sample_input_2.txt
│
├── docs/
│   ├── scripts.md
│   ├── datasets.md
│   ├── input_formats.md
│   └── output_formats.md
│
├── tests/
│   ├── test_cancer_vs_cancer.py
│   ├── test_GSE201284_correlations.py
│   └── test_GSE120741_correlations.py
│
├── outputs/
│   └── .gitkeep
│
└── logs/
    └── .gitkeep
```

## Current scripts

### 1. `cancer_vs_cancer.py`

This script processes the GSE183019 dataset. It maps the expression-table sample names to GEO GSM identifiers, identifies cancer samples, exports gene-expression values, and calculates pairwise Pearson correlations between genes.

Two sample inputs are provided:

```text
samples/cancer_vs_cancer/GSE183019_series_matrix_sample.txt
samples/cancer_vs_cancer/GSE183019_processed_TPM_sample.txt
```

### 2. `GSE201284_all_feature_correlations.py`

This script performs all-gene correlation analysis for GSE201284 using clinical and sample metadata. It supports subgroup analysis for disease state, PSA, Gleason score, race, pathologic stage, age, five-year disease status, sample type, and source name.

Two sample inputs are provided:

```text
samples/GSE201284_all_feature_correlations/GSE201284_processed_TPM_sample.txt
samples/GSE201284_all_feature_correlations/GSE201284_sample_metadata_sample.tsv
```

### 3. `GSE120741_all_feature_correlations.py`

This script performs all-gene correlation analysis for GSE120741 using ERG-expression and recurrence-status groups. It uses the normalized Porto gene-expression table together with the GEO series-matrix metadata.

Two sample inputs are provided:

```text
samples/GSE120741_all_feature_correlations/GSE120741_Porto_ge_table_sample.txt
samples/GSE120741_all_feature_correlations/GSE120741_series_matrix_sample.txt
```

## Rule for adding a new script

Every new analysis should include the following items:

1. One Python file inside `scripts/`.
2. One matching folder inside `samples/`.
3. Exactly two small, representative sample input files inside that sample folder.
4. A description of the script in `docs/scripts.md`.
5. Input-column and file-format information in `docs/input_formats.md`.
6. Output-file information in `docs/output_formats.md`.
7. A test inside `tests/` when practical.

Use the following naming pattern:

```text
scripts/<script_name>.py
samples/<script_name>/sample_input_1.<extension>
samples/<script_name>/sample_input_2.<extension>
tests/test_<script_name>.py
```

## Installation

Clone the repository and enter the project directory:

```bash
git clone <repository-url>
cd ualcan
```

Create and activate a virtual environment:

```bash
python3 -m venv .venv
source .venv/bin/activate
```

Install the required Python packages:

```bash
python -m pip install --upgrade pip
pip install -r requirements.txt
```

The initial `requirements.txt` should contain:

```text
numpy
pandas
scipy
pytest
```

## Running a script

Run a script from the repository root:

```bash
python scripts/<script_name>.py --help
```

Example:

```bash
python scripts/GSE201284_all_feature_correlations.py \
  --tpm samples/GSE201284_all_feature_correlations/GSE201284_processed_TPM_sample.txt \
  --metadata samples/GSE201284_all_feature_correlations/GSE201284_sample_metadata_sample.tsv \
  --output-root outputs/GSE201284_sample_results \
  --max-genes 10
```

The command used for each script should be documented in `docs/scripts.md`.

## Sample-file requirements

The two sample files provided for each script should:

- Contain only a small number of representative rows and columns.
- Preserve the same column names and general structure as the complete input files.
- Be sufficient to demonstrate how the script reads and processes its inputs.
- Use de-identified or publicly available data.
- Avoid confidential information and patient identifiers.
- Include a clear header where the original format supports one.

Sample files are examples only. They are not substitutes for the complete datasets used in the full analysis.

## Files that must not be committed

Do not upload:

- Complete GEO datasets or other large raw-data files.
- Complete TPM or gene-expression matrices.
- Generated correlation and expression JSON collections.
- HPC job logs or temporary files.
- Python virtual environments.
- API keys, passwords, access tokens, or configuration secrets.
- Patient-identifiable or confidential research information.

The `.gitignore` file should exclude these files and directories while allowing the small example files under `samples/`.

## Output storage

Generated results should be written to `outputs/`. This directory is retained in the repository using `.gitkeep`, but its generated contents should be ignored by Git.

Runtime and HPC logs should be written to `logs/`. Log contents should also be excluded from Git.

## Project status

This repository will grow as additional Python workflows are created. Each new script should follow the same structure so that the code, two associated sample inputs, documentation, and tests remain easy to locate.

## Use and validation

These scripts are intended for research and reproducible data-analysis workflows. All results must be validated before they are used in publications, downstream databases, or clinical interpretation.

## License

No open-source license is currently provided. Approval should be obtained from the supervising laboratory or principal investigator before the repository is made public or the scripts are authorized for external reuse.
