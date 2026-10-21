#!/usr/bin/env python3
"""
Read Tumor Volumes Excel files and plot tumour growth per mouse by cage.

Default input folder is this script's folder. Outputs are written to a
`tumour_growth_plots` subfolder unless --output-dir is supplied.
"""

from __future__ import annotations

import argparse
import math
import os
import re
import tempfile
from datetime import date
from pathlib import Path
from typing import Any

plot_cache_dir = Path(tempfile.gettempdir()) / f"pdx-growth-cache-{os.getuid()}"
os.environ.setdefault("MPLCONFIGDIR", str(plot_cache_dir / "matplotlib"))
os.environ.setdefault("XDG_CACHE_HOME", str(plot_cache_dir))

import matplotlib
import matplotlib.dates as mdates

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


DATE_RE = re.compile(r"(\d{4})\.(\d{2})\.(\d{2})")
PASSAGE_RE = re.compile(r"^passage\b", re.IGNORECASE)
MODEL_RE = re.compile(r"\b(SA\d+[A-Z]?)\b(?:\s+([A-Za-z]+))?", re.IGNORECASE)
FIT_WINDOWS = (5, 4, 3)
MIN_PROJECTION_R2 = 0.5
MIN_LATEST_VOLUME_FOR_PROJECTION = 200.0
PROJECTION_THRESHOLDS = (500, 1000)


def parse_args() -> argparse.Namespace:
    script_dir = Path(__file__).resolve().parent
    parser = argparse.ArgumentParser(
        description="Plot tumour growth per mouse, faceted by cage."
    )
    parser.add_argument(
        "--config",
        type=Path,
        default=script_dir / "config.yml",
        help="YAML-like config file with input/output paths and analysis settings.",
    )
    parser.add_argument(
        "input_dir",
        nargs="?",
        type=Path,
        default=None,
        help="Folder containing Tumor Volumes*.xlsx files. Overrides config input_dir.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=None,
        help="Output folder. Overrides config output_dir.",
    )
    parser.add_argument(
        "--pattern",
        default=None,
        help="Excel filename glob pattern to read.",
    )
    parser.add_argument(
        "--aggregate",
        choices=["sum", "max", "mean"],
        default=None,
        help="How to combine multiple tumour sites for the same mouse/date.",
    )
    parser.add_argument(
        "--include-zero",
        action="store_true",
        help="Include rows with recorded volume 0. Default skips zero-volume rows.",
    )
    parser.add_argument(
        "--scale",
        choices=["free", "same", "log", "all"],
        default=None,
        help="Y-axis scale to plot. Default writes free, same-linear, and log versions.",
    )
    return parser.parse_args()


def parse_scalar_config(value: str) -> Any:
    value = value.strip()
    if not value:
        return ""
    if value[0] in {"'", '"'} and value[-1:] == value[0]:
        return value[1:-1]
    lowered = value.lower()
    if lowered in {"true", "false"}:
        return lowered == "true"
    if value.startswith("[") and value.endswith("]"):
        items = [item.strip() for item in value[1:-1].split(",") if item.strip()]
        return [parse_scalar_config(item) for item in items]
    try:
        if "." in value:
            return float(value)
        return int(value)
    except ValueError:
        return value


def load_config(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {}
    config: dict[str, Any] = {}
    current_list_key: str | None = None
    for raw_line in path.read_text(encoding="utf-8").splitlines():
        line = raw_line.split("#", 1)[0].rstrip()
        if not line.strip():
            continue
        stripped = line.strip()
        if stripped.startswith("- ") and current_list_key:
            config.setdefault(current_list_key, []).append(parse_scalar_config(stripped[2:]))
            continue
        if ":" not in stripped:
            continue
        key, value = stripped.split(":", 1)
        key = key.strip()
        value = value.strip()
        current_list_key = key if value == "" else None
        config[key] = [] if value == "" else parse_scalar_config(value)
    return config


def config_path(value: Any, base_dir: Path) -> Path | None:
    if not value:
        return None
    path = Path(str(value)).expanduser()
    return path if path.is_absolute() else (base_dir / path).resolve()


def as_text(value: Any) -> str:
    if pd.isna(value):
        return ""
    return str(value).strip()


def parse_date_from_filename(path: Path) -> pd.Timestamp | None:
    match = DATE_RE.search(path.name)
    if not match:
        return None
    year, month, day = map(int, match.groups())
    return pd.Timestamp(year=year, month=month, day=day)


def normalize_cage(mouse_id: str) -> str:
    cage = mouse_id.strip().replace("*", "")
    # Parenthesised identifiers distinguish physical cages, for example
    # 251220(1) and 251220(3). Only a trailing mouse suffix is removed.
    cage = re.sub(r"-\d+\s*$", "", cage)
    return cage.strip() or mouse_id


def extract_tumour_model(comment: Any) -> str:
    text = as_text(comment)
    match = MODEL_RE.search(text)
    if not match:
        return ""
    model = match.group(1).upper()
    tissue = match.group(2) or ""
    if tissue:
        return f"{model} {tissue.capitalize()}"
    return model


def looks_like_mouse_id(value: Any) -> bool:
    text = as_text(value)
    if not text or PASSAGE_RE.match(text):
        return False
    if text.lower().startswith("cage "):
        return False
    if "tumor" in text.lower() or "tumour" in text.lower():
        return False
    return bool(re.search(r"\d", text))


def numeric(value: Any) -> float | None:
    number = pd.to_numeric(value, errors="coerce")
    if pd.isna(number):
        return None
    return float(number)


def configure_date_axis(ax: plt.Axes, max_ticks: int = 6) -> None:
    locator = mdates.AutoDateLocator(minticks=3, maxticks=max_ticks)
    formatter = mdates.ConciseDateFormatter(locator)
    formatter.formats = ["%Y", "%b", "%d", "%H:%M", "%H:%M", "%S.%f"]
    formatter.zero_formats = ["", "%Y", "%b", "%b %d", "%H:%M", "%H:%M"]
    formatter.offset_formats = ["", "", "", "", "", ""]
    ax.xaxis.set_major_locator(locator)
    ax.xaxis.set_major_formatter(formatter)
    ax.tick_params(axis="x", labelrotation=30, labelsize=8)
    for label in ax.get_xticklabels():
        label.set_horizontalalignment("right")


def annotation_offset(index: int, threshold: int) -> tuple[int, int]:
    x_offsets = [4, 14, 24, -18, -30, 36, -42, 48]
    y_base = 11 if threshold == 500 else -17
    y_step = 8 if threshold == 500 else -8
    return x_offsets[index % len(x_offsets)], y_base + (index // len(x_offsets)) * y_step


def threshold_label_offset(series_index: int, threshold: int) -> tuple[int, int]:
    offsets_500 = [(16, 18), (-24, 18), (30, -22), (-36, -22), (44, 2), (-48, 2), (20, 34), (-24, -38)]
    offsets_1000 = [(18, -22), (-26, -24), (34, 20), (-40, 20), (50, -4), (-54, -4), (22, -40), (-26, 36)]
    offsets = offsets_500 if threshold == 500 else offsets_1000
    return offsets[series_index % len(offsets)]


def display_date_label(value: pd.Timestamp) -> str:
    timestamp = pd.Timestamp(value)
    return f"{timestamp.day} {timestamp.strftime('%b')}"


def fit_recent_exponential(trace: pd.DataFrame) -> dict[str, Any]:
    candidates: list[dict[str, Any]] = []
    positive = trace[trace["tumour_volume_mm3"] > 0].copy()
    for window in FIT_WINDOWS:
        if len(positive) < window:
            continue
        fit = positive.tail(window).copy()
        fit_start = fit["date"].min()
        x = (fit["date"] - fit_start).dt.days.to_numpy(dtype=float)
        y = np.log(fit["tumour_volume_mm3"].to_numpy(dtype=float))
        growth_rate, intercept = np.polyfit(x, y, 1)
        fitted = intercept + growth_rate * x
        residual = float(np.sum((y - fitted) ** 2))
        total = float(np.sum((y - y.mean()) ** 2))
        r_squared = 1 - residual / total if total > 0 else np.nan
        fit_status = "projectable" if growth_rate > 0 else "non-growing fit"
        candidates.append(
            {
                "fit": fit,
                "fit_window": window,
                "daily_log_growth_rate": growth_rate,
                "fit_intercept": intercept,
                "fit_r_squared": r_squared,
                "fit_status": fit_status,
            }
        )

    projectable = [c for c in candidates if c["fit_status"] == "projectable"]
    acceptable = [
        c for c in projectable
        if pd.notna(c["fit_r_squared"]) and c["fit_r_squared"] >= MIN_PROJECTION_R2
    ]
    if acceptable:
        return sorted(acceptable, key=lambda c: (c["fit_r_squared"], c["fit_window"]), reverse=True)[0]
    if projectable:
        best = sorted(projectable, key=lambda c: c["fit_r_squared"] if pd.notna(c["fit_r_squared"]) else -1, reverse=True)[0]
        best = best.copy()
        best["fit_status"] = "weak projectable fit"
        return best
    if candidates:
        return candidates[0]
    return {
        "fit": positive,
        "fit_window": len(positive),
        "daily_log_growth_rate": np.nan,
        "fit_intercept": np.nan,
        "fit_r_squared": np.nan,
        "fit_status": "insufficient data",
    }


def find_header(df: pd.DataFrame) -> tuple[int, int, int, int | None, int | None] | None:
    for row_idx, row in df.iterrows():
        values = [as_text(v) for v in row.tolist()]
        mouse_cols = [i for i, v in enumerate(values) if v == "Mouse ID"]
        volume_cols = [i for i, v in enumerate(values) if "Volume" in v]
        if not mouse_cols or not volume_cols:
            continue

        mouse_col = mouse_cols[0]
        volume_col = volume_cols[0]
        length_col = None
        comment_col = None
        for col_idx, value in enumerate(values):
            if value == "Comments":
                comment_col = col_idx
                break
        next_row = df.iloc[row_idx + 1].tolist() if row_idx + 1 < len(df) else []
        for col_idx, value in enumerate(next_row):
            if "length" in as_text(value).lower():
                length_col = col_idx
                break

        # Some sheets have a blank tumour-site column between mouse and length.
        site_col = None
        if length_col is not None and length_col - mouse_col > 1:
            site_col = mouse_col + 1

        return row_idx, mouse_col, volume_col, site_col, comment_col
    return None


def extract_sheet(
    path: Path,
    sheet_name: str,
    date: pd.Timestamp,
    include_zero: bool,
) -> list[dict[str, Any]]:
    df = pd.read_excel(path, sheet_name=sheet_name, header=None)
    header = find_header(df)
    if header is None:
        return []

    header_row, mouse_col, volume_col, site_col, comment_col = header
    current_mouse = ""
    current_model = ""
    rows: list[dict[str, Any]] = []

    for row_idx in range(header_row + 1, len(df)):
        raw_mouse = df.iat[row_idx, mouse_col] if mouse_col < df.shape[1] else None
        comment = df.iat[row_idx, comment_col] if comment_col is not None and comment_col < df.shape[1] else ""
        if looks_like_mouse_id(raw_mouse):
            current_mouse = as_text(raw_mouse)
            current_model = ""

        row_model = extract_tumour_model(comment)
        if row_model:
            current_model = row_model

        if not current_mouse:
            continue

        volume = numeric(df.iat[row_idx, volume_col])
        if volume is None:
            continue
        if volume == 0 and not include_zero:
            continue

        site = ""
        if site_col is not None and site_col < df.shape[1]:
            site = as_text(df.iat[row_idx, site_col])

        note = as_text(raw_mouse) if PASSAGE_RE.match(as_text(raw_mouse)) else ""
        rows.append(
            {
                "date": date.date().isoformat(),
                "source_file": path.name,
                "sheet": sheet_name,
                "excel_row": row_idx + 1,
                "mouse_id": current_mouse,
                "cage": normalize_cage(current_mouse),
                "tumour_model": current_model,
                "tumour_site": site,
                "tumour_volume_mm3": volume,
                "row_note": note,
            }
        )

    return rows


def read_all(input_dir: Path, pattern: str, include_zero: bool) -> pd.DataFrame:
    records: list[dict[str, Any]] = []
    for path in sorted(input_dir.glob(pattern)):
        if path.name.startswith("~$"):
            continue
        date = parse_date_from_filename(path)
        if date is None:
            print(f"Skipping {path.name}: no YYYY.MM.DD date in filename")
            continue
        for sheet_name in pd.ExcelFile(path).sheet_names:
            records.extend(extract_sheet(path, sheet_name, date, include_zero))

    if not records:
        return pd.DataFrame()

    data = pd.DataFrame(records)
    data["date"] = pd.to_datetime(data["date"])
    return data.sort_values(["cage", "mouse_id", "date", "sheet", "excel_row"])


def aggregate_for_plot(data: pd.DataFrame, method: str) -> pd.DataFrame:
    plot_data = data.copy()
    plot_data["tumour_site"] = plot_data["tumour_site"].replace("", "site")
    group_models = (
        plot_data[plot_data["tumour_model"].astype(str) != ""]
        .groupby(["cage", "sheet"])["tumour_model"]
        .agg(lambda values: ", ".join(sorted(set(values))))
        .to_dict()
    )
    plot_data["group_model"] = plot_data.apply(
        lambda row: group_models.get((row["cage"], row["sheet"]), ""),
        axis=1,
    )
    plot_data["plot_group"] = plot_data.apply(
        lambda row: (
            f"{row['cage']} {row['group_model']} ({row['sheet']})"
            if row["group_model"]
            else f"{row['cage']} ({row['sheet']})"
        ),
        axis=1,
    )
    plot_data["plot_label"] = plot_data["mouse_id"].astype(str) + " " + plot_data["tumour_site"].astype(str)

    grouped = plot_data.groupby(["date", "plot_group", "plot_label"], as_index=False)[
        "tumour_volume_mm3"
    ]
    if method == "sum":
        aggregated = grouped.sum()
    elif method == "max":
        aggregated = grouped.max()
    else:
        aggregated = grouped.mean()

    return aggregated.sort_values(["plot_group", "plot_label", "date"])


def plot_by_cage(
    plot_data: pd.DataFrame,
    output_png: Path,
    output_pdf: Path,
    y_mode: str = "free",
    generated_date: str = "",
    latest_measurement_date: str = "",
) -> None:
    plot_groups = sorted(plot_data["plot_group"].dropna().unique())
    n_groups = len(plot_groups)
    n_cols = min(3, max(1, n_groups))
    n_rows = math.ceil(n_groups / n_cols)
    shared_y = y_mode in {"same", "log"}
    max_volume = plot_data["tumour_volume_mm3"].max()

    fig, axes = plt.subplots(
        n_rows,
        n_cols,
        figsize=(5.2 * n_cols, 3.9 * n_rows),
        sharex=False,
        sharey=shared_y,
        squeeze=False,
    )

    for ax, plot_group in zip(axes.flat, plot_groups):
        group_data = plot_data[plot_data["plot_group"] == plot_group]
        for plot_label, mouse_data in group_data.groupby("plot_label"):
            mouse_data = mouse_data.sort_values("date")
            ax.plot(
                mouse_data["date"],
                mouse_data["tumour_volume_mm3"],
                marker="o",
                linewidth=1.8,
                label=plot_label,
            )

        ax.set_title(f"Cage {plot_group}", fontsize=11)
        ax.set_ylabel("Tumour volume (mm3)")
        ax.grid(True, alpha=0.25)
        ax.legend(fontsize=8, frameon=False)
        configure_date_axis(ax)
        if y_mode == "same":
            ax.set_ylim(0, max_volume * 1.08)
        elif y_mode == "log":
            ax.set_yscale("log")
            ax.set_ylim(max(1, plot_data["tumour_volume_mm3"].min() * 0.85), max_volume * 1.15)

    for ax in axes.flat[n_groups:]:
        ax.axis("off")

    for ax in list(axes.flat)[:n_groups]:
        ax.set_xlabel("Date", fontsize=9)

    scale_label = {
        "free": "free y-axis",
        "same": "same linear y-axis",
        "log": "log y-axis",
    }[y_mode]
    fig.suptitle(f"Tumour Growth per Mouse by Cage ({scale_label})", fontsize=15, y=0.995)
    fig.text(
        0.5,
        0.975,
        f"Generated {generated_date} | Data through {latest_measurement_date}",
        ha="center",
        va="top",
        fontsize=10,
    )
    fig.tight_layout(rect=(0, 0, 1, 0.955))
    fig.savefig(output_png, dpi=300, bbox_inches="tight")
    fig.savefig(output_pdf, bbox_inches="tight")
    plt.close(fig)


def build_pdx_projections(data: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    pdx = data[data["tumour_model"].astype(str).str.match(r"^SA\d", na=False)].copy()
    analysis_date = pdx["date"].max()
    pdx["tumour_site"] = pdx["tumour_site"].replace("", "site")
    group_cols = ["cage", "sheet", "tumour_model", "mouse_id", "tumour_site"]
    observed_rows: list[pd.DataFrame] = []
    summaries: list[dict[str, Any]] = []

    for keys, trace in pdx.groupby(group_cols):
        cage, sheet, model, mouse_id, site = keys
        trace = (
            trace.groupby("date", as_index=False)["tumour_volume_mm3"]
            .mean()
            .sort_values("date")
        )
        latest = trace.iloc[-1]
        if latest["date"] != analysis_date or latest["tumour_volume_mm3"] <= MIN_LATEST_VOLUME_FOR_PROJECTION:
            continue

        trace = trace.copy()
        trace["cage"] = cage
        trace["sheet"] = sheet
        trace["tumour_model"] = model
        trace["mouse_id"] = mouse_id
        trace["tumour_site"] = site
        trace["trace_label"] = f"{mouse_id} {site}"
        observed_rows.append(trace)

        fit_result = fit_recent_exponential(trace)
        fit = fit_result["fit"]
        growth_rate = fit_result["daily_log_growth_rate"]
        r_squared = fit_result["fit_r_squared"]
        fit_status = fit_result["fit_status"]
        doubling_days = math.log(2) / growth_rate if pd.notna(growth_rate) and growth_rate > 0 else np.nan

        summary: dict[str, Any] = {
            "cage": cage,
            "worksheet": sheet,
            "tumour_model": model,
            "mouse_id": mouse_id,
            "tumour_site": site,
            "latest_measurement_date": latest["date"].date().isoformat(),
            "latest_volume_mm3": latest["tumour_volume_mm3"],
            "fit_measurements": len(fit),
            "fit_window_selected": fit_result["fit_window"],
            "daily_log_growth_rate": growth_rate,
            "doubling_time_days": doubling_days,
            "fit_r_squared": r_squared,
            "fit_status": fit_status,
        }

        for threshold in PROJECTION_THRESHOLDS:
            crossing = trace[trace["tumour_volume_mm3"] >= threshold]
            if not crossing.empty:
                crossing_date = crossing.iloc[0]["date"]
                status = "observed"
            elif fit_status in {"projectable", "weak projectable fit"}:
                days = math.log(threshold / latest["tumour_volume_mm3"]) / growth_rate
                crossing_date = latest["date"] + pd.to_timedelta(max(0, days), unit="D")
                status = "projected"
            else:
                crossing_date = pd.NaT
                status = "not projected"
            summary[f"date_{threshold}_mm3"] = (
                crossing_date.date().isoformat() if pd.notna(crossing_date) else ""
            )
            summary[f"date_{threshold}_status"] = status

        summaries.append(summary)

    observed = pd.concat(observed_rows, ignore_index=True) if observed_rows else pd.DataFrame()
    return observed, pd.DataFrame(summaries)


def plot_pdx_projections(
    observed: pd.DataFrame,
    summary: pd.DataFrame,
    output_png: Path,
    output_pdf: Path,
    generated_date: str,
    latest_measurement_date: str,
) -> None:
    if observed.empty or summary.empty:
        return

    groups = sorted(observed[["cage", "sheet", "tumour_model"]].drop_duplicates().itertuples(index=False, name=None))
    n_cols = min(3, len(groups))
    n_rows = math.ceil(len(groups) / n_cols)
    fig, axes = plt.subplots(
        n_rows,
        n_cols,
        figsize=(7.0 * n_cols, 5.0 * n_rows + 0.8),
        sharex=False,
        sharey=True,
        squeeze=False,
    )
    y_max = max(1300, float(observed["tumour_volume_mm3"].max()) * 1.18)

    colors = plt.get_cmap("tab10")
    for ax, (cage, sheet, model) in zip(axes.flat, groups):
        group_observed = observed[
            (observed["cage"] == cage)
            & (observed["sheet"] == sheet)
            & (observed["tumour_model"] == model)
        ]
        group_summary = summary[
            (summary["cage"] == cage)
            & (summary["worksheet"] == sheet)
            & (summary["tumour_model"] == model)
        ]

        for color_index, (label, trace) in enumerate(group_observed.groupby("trace_label")):
            color = colors(color_index % 10)
            trace = trace.sort_values("date")
            ax.plot(
                trace["date"],
                trace["tumour_volume_mm3"],
                color=color,
                marker="o",
                linewidth=1.8,
                label=label,
            )
            row = group_summary[
                (group_summary["mouse_id"].astype(str) + " " + group_summary["tumour_site"].astype(str)) == label
            ].iloc[0]
            if row["fit_status"] in {"projectable", "weak projectable fit"} and row["latest_volume_mm3"] < 1000:
                latest_date = trace.iloc[-1]["date"]
                latest_volume = trace.iloc[-1]["tumour_volume_mm3"]
                projected_dates = [latest_date]
                for threshold in PROJECTION_THRESHOLDS:
                    value = row[f"date_{threshold}_mm3"]
                    if row[f"date_{threshold}_status"] == "projected" and value:
                        projected_dates.append(pd.Timestamp(value))
                end_date = max(projected_dates)
                projection_dates = pd.date_range(latest_date, end_date, periods=40)
                days = (projection_dates - latest_date).days.to_numpy(dtype=float)
                projection_values = latest_volume * np.exp(row["daily_log_growth_rate"] * days)
                projection_values = np.minimum(projection_values, max(PROJECTION_THRESHOLDS))
                ax.plot(projection_dates, projection_values, color=color, linestyle="--", linewidth=1.6)

            threshold_markers = {500: "x", 1000: "D"}
            final_label_threshold = max(
                [threshold for threshold in PROJECTION_THRESHOLDS if row.get(f"date_{threshold}_mm3", "")],
                default=max(PROJECTION_THRESHOLDS),
            )
            for threshold in PROJECTION_THRESHOLDS:
                marker = threshold_markers.get(threshold, "o")
                value = row[f"date_{threshold}_mm3"]
                if not value:
                    continue
                threshold_date = pd.Timestamp(value)
                ax.scatter(threshold_date, threshold, color=color, marker=marker, s=45, zorder=5)
                if threshold == final_label_threshold:
                    xytext = threshold_label_offset(color_index, threshold)
                    ax.annotate(
                        display_date_label(threshold_date),
                        (threshold_date, threshold),
                        xytext=xytext,
                        textcoords="offset points",
                        fontsize=8,
                        fontweight="bold",
                        color=color,
                        ha="left" if xytext[0] >= 0 else "right",
                        va="bottom" if xytext[1] >= 0 else "top",
                        arrowprops={"arrowstyle": "-", "color": color, "lw": 0.5, "alpha": 0.65},
                    )

        ax.axhline(500, color="#555555", linestyle=":", linewidth=1)
        ax.axhline(1000, color="#222222", linestyle=":", linewidth=1)
        ax.set_title(f"Cage {cage} {model} ({sheet})", fontsize=11)
        ax.set_ylabel("Tumour volume (mm3)")
        ax.set_ylim(0, y_max)
        ax.grid(True, alpha=0.25)
        configure_date_axis(ax)

    for ax in axes.flat[len(groups):]:
        ax.axis("off")
    for ax in list(axes.flat)[:len(groups)]:
        ax.set_xlabel("Date", fontsize=9)

    fig.suptitle("Current PDX Growth Projection for Volume >200 mm3", fontsize=15, y=0.98)
    fig.text(
        0.5,
        0.925,
        f"Measured on {latest_measurement_date} | Solid = observed | Dashed = projected | Generated {generated_date}",
        ha="center",
        va="top",
        fontsize=10,
    )
    fig.tight_layout(rect=(0, 0, 1, 0.87))
    fig.savefig(output_png, dpi=300, bbox_inches="tight")
    fig.savefig(output_pdf, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    global FIT_WINDOWS, MIN_PROJECTION_R2, MIN_LATEST_VOLUME_FOR_PROJECTION, PROJECTION_THRESHOLDS

    args = parse_args()
    script_dir = Path(__file__).resolve().parent
    config_file = args.config.expanduser().resolve()
    config = load_config(config_file)
    config_base = config_file.parent

    FIT_WINDOWS = tuple(int(value) for value in config.get("projection_fit_windows", FIT_WINDOWS))
    MIN_PROJECTION_R2 = float(config.get("min_projection_r2", MIN_PROJECTION_R2))
    MIN_LATEST_VOLUME_FOR_PROJECTION = float(
        config.get("min_latest_volume_for_projection", MIN_LATEST_VOLUME_FOR_PROJECTION)
    )
    PROJECTION_THRESHOLDS = tuple(int(value) for value in config.get("projection_thresholds_mm3", PROJECTION_THRESHOLDS))

    input_dir = (
        args.input_dir.expanduser().resolve()
        if args.input_dir
        else config_path(config.get("input_dir"), config_base) or script_dir
    )
    pattern = args.pattern or str(config.get("pattern", "Tumor Volumes*.xlsx"))
    aggregate = args.aggregate or str(config.get("aggregate", "sum"))
    scale = args.scale or str(config.get("scale", "all"))
    include_zero = bool(args.include_zero or config.get("include_zero", False))
    output_dir = (
        args.output_dir.expanduser().resolve()
        if args.output_dir
        else config_path(config.get("output_dir"), config_base) or input_dir / "tumour_growth_plots"
    )
    output_dir.mkdir(parents=True, exist_ok=True)

    data = read_all(input_dir, pattern, include_zero)
    if data.empty:
        raise SystemExit(f"No tumour-volume records found in {input_dir}")

    generated_date = date.today().isoformat()
    latest_measurement_date = data["date"].max().date().isoformat()
    date_tag = f"generated_{generated_date}_through_{latest_measurement_date}"

    data["analysis_generated_date"] = generated_date
    data["latest_measurement_date"] = latest_measurement_date
    extracted_csv = output_dir / f"tumour_growth_extracted_rows_{date_tag}.csv"
    plot_csv = output_dir / f"tumour_growth_per_mouse_{aggregate}_{date_tag}.csv"
    data.to_csv(extracted_csv, index=False)
    plot_data = aggregate_for_plot(data, aggregate)
    plot_data["analysis_generated_date"] = generated_date
    plot_data["latest_measurement_date"] = latest_measurement_date
    plot_data.to_csv(plot_csv, index=False)
    plot_modes = ["free", "same", "log"] if scale == "all" else [scale]
    written_plots: list[Path] = []
    for mode in plot_modes:
        suffix = "" if mode == "free" else f"_{mode}"
        mode_png = output_dir / f"tumour_growth_by_cage_{aggregate}{suffix}_{date_tag}.png"
        mode_pdf = output_dir / f"tumour_growth_by_cage_{aggregate}{suffix}_{date_tag}.pdf"
        plot_by_cage(
            plot_data,
            mode_png,
            mode_pdf,
            y_mode=mode,
            generated_date=generated_date,
            latest_measurement_date=latest_measurement_date,
        )
        written_plots.extend([mode_png, mode_pdf])

    pdx_observed, pdx_summary = build_pdx_projections(data)
    pdx_csv = output_dir / f"pdx_growth_projection_over_200_{date_tag}.csv"
    pdx_png = output_dir / f"pdx_growth_projection_over_200_{date_tag}.png"
    pdx_pdf = output_dir / f"pdx_growth_projection_over_200_{date_tag}.pdf"
    if not pdx_summary.empty:
        pdx_summary["analysis_generated_date"] = generated_date
        pdx_summary["data_through_date"] = latest_measurement_date
        pdx_summary.to_csv(pdx_csv, index=False)
        plot_pdx_projections(
            pdx_observed,
            pdx_summary,
            pdx_png,
            pdx_pdf,
            generated_date,
            latest_measurement_date,
        )

    print(f"Read {data['source_file'].nunique()} Excel files")
    print(f"Extracted {len(data)} tumour-site rows")
    print(
        f"Plotted {plot_data['plot_label'].nunique()} mouse/site traces "
        f"across {plot_data['plot_group'].nunique()} cage worksheet groups"
    )
    print(f"Wrote {extracted_csv}")
    print(f"Wrote {plot_csv}")
    for path in written_plots:
        print(f"Wrote {path}")
    if not pdx_summary.empty:
        print(f"Projected {len(pdx_summary)} PDX mouse/site traces with latest volume >200 mm3")
        print(f"Wrote {pdx_csv}")
        print(f"Wrote {pdx_png}")
        print(f"Wrote {pdx_pdf}")


if __name__ == "__main__":
    main()
