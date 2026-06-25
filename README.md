# Tumour Growth Analysis

Reusable analysis script for mouse tumour-volume Excel workbooks.

The script reads dated `Tumor Volumes*.xlsx` files, extracts tumour measurements by worksheet, mouse/cage ID, tumour model, and tumour site, then generates:

- tumour-growth plots by cage and worksheet
- free-scale, same-linear-scale, and log-scale plots
- audit CSVs for extracted rows and plotted values
- current PDX projections for active tumours above a configurable volume threshold

An example workbook is included in `examples/`.

## Quick Start

```bash
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
python plot_tumour_growth_by_cage.py --config config.yml
```

Outputs are written to the directory configured in `config.yml`.

## Configuration

The analysis is not tied to a hard-coded data directory. Edit `config.yml`:

```yaml
input_dir: examples
output_dir: outputs
pattern: "Tumor Volumes*.xlsx"
aggregate: sum
include_zero: false
scale: all
min_latest_volume_for_projection: 200
projection_thresholds_mm3: [500, 1000]
projection_fit_windows: [5, 4, 3]
min_projection_r2: 0.5
```

Command-line arguments can override the config:

```bash
python plot_tumour_growth_by_cage.py --config config.yml /path/to/excel/files --output-dir /path/to/output
```

## Workbook Extraction

The script scans each workbook for sheets with `Mouse ID` and `Volume (mm3)` headers. It extracts:

- date from the workbook filename, for example `2026.06.22`
- worksheet name, for example `Xenograft` or `Organoid Study`
- mouse ID and cage ID
- tumour model from the comments column, for example `SA501 Breast`
- tumour site, for example `L`, `R`, `B`, or `N`
- tumour volume in mm3
- source workbook and Excel row number for auditability

Rows with blank mouse IDs inherit the most recent valid mouse ID, matching the template structure used in the source workbooks. Rows marked as passage notes are kept as notes but do not overwrite the mouse ID.

## Growth Plots

The main growth plots group traces by cage, tumour model, and worksheet. Trace labels include mouse ID and tumour site.

Generated plot variants:

- free y-axis
- same linear y-axis
- log y-axis

Date axes are formatted with limited, readable ticks. Plot filenames include both the generation date and the latest measurement date, for example:

```text
tumour_growth_by_cage_sum_generated_2026-06-23_through_2026-06-22.png
```

## PDX Projection Logic

Current PDX projection includes only traces that:

- have a tumour model beginning with `SA`
- were measured on the latest measurement date in the dataset
- have latest volume greater than `min_latest_volume_for_projection`

For each trace, the script fits recent exponential growth using log-linear regression. It tests recent windows from `projection_fit_windows`, such as 5, 4, and 3 points. It selects the best projectable fit with `R2 >= min_projection_r2`; if none pass, it keeps the best projectable fit and labels it as weak in the CSV.

The projection CSV includes:

- latest measured volume and date
- selected fit window
- daily log growth rate
- doubling time
- R2
- projected or observed dates for each threshold

The projection plot shows observed data as solid lines and projected data as dashed lines. It labels only the final threshold date per trace, in bold, to reduce overlap.

## Audit Outputs

The extracted-row CSV is the source-of-truth audit file. It records the source workbook, worksheet, row number, mouse ID, cage, tumour model, tumour site, and extracted volume for each point.

Use this file to verify any plotted point against the original Excel source.
