from __future__ import annotations

import argparse
import math
import re
import warnings
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

try:
    import pyreadstat
except ImportError:
    pyreadstat = None


SUPPORTED = {".csv", ".tsv", ".txt", ".xls", ".xlsx", ".sav", ".zsav", ".dta", ".sas7bdat", ".xpt"}
REPORT_FILES = {"audit_outputs/district_counts.csv", "audit_outputs/completeness_tables.csv"}
SPECIAL_RE = re.compile(r"(?:^|[_ ])(?:888|999|997|998|99)(?:$|[_ ])", re.I)


def clean_name(value: Any) -> str:
    return re.sub(r"[^a-z0-9]+", "", str(value).strip().lower())


def display_value(value: Any) -> str:
    if pd.isna(value):
        return "<MISSING>"
    if isinstance(value, float) and value.is_integer():
        return str(int(value))
    return str(value)


def read_table(path: Path) -> tuple[pd.DataFrame, dict[str, str], dict[str, dict[Any, Any]], str]:
    suffix = path.suffix.lower()
    labels: dict[str, str] = {}
    value_labels: dict[str, dict[Any, Any]] = {}
    notes = ""
    if suffix == ".csv" or suffix in {".tsv", ".txt"}:
        header_row = 0
        if path.name.lower() == "nationalsurveydata.csv":
            header_row = 4
            notes = "Parsed with the verified five-line BGS preamble; the second row contains measurement units."
        df = pd.read_csv(path, header=header_row, sep="\t" if suffix == ".tsv" else ",", low_memory=False)
        if path.name.lower() == "nationalsurveydata.csv" and len(df) and str(df.iloc[0, 0]).strip() == "":
            df = df.iloc[1:].reset_index(drop=True)
        return df, labels, value_labels, notes
    if suffix in {".xls", ".xlsx"}:
        df = pd.read_excel(path)
        return df, labels, value_labels, "Excel labels are represented by column names; no source variable-label metadata was available."
    if pyreadstat is None:
        raise RuntimeError("pyreadstat is required for SPSS, Stata, SAS, and XPT files.")
    if suffix in {".sav", ".zsav"}:
        df, meta = pyreadstat.read_sav(path, apply_value_formats=False)
    elif suffix == ".dta":
        df, meta = pyreadstat.read_dta(path, apply_value_formats=False)
    elif suffix == ".sas7bdat":
        df, meta = pyreadstat.read_sas7bdat(path, apply_value_formats=False)
    elif suffix == ".xpt":
        df, meta = pyreadstat.read_xport(path, apply_value_formats=False)
    else:
        raise ValueError(f"Unsupported format: {path}")
    labels = {str(k): str(v) for k, v in getattr(meta, "column_names_to_labels", {}).items()}
    value_labels = {str(k): dict(v) for k, v in getattr(meta, "variable_value_labels", {}).items()}
    return df, labels, value_labels, "Raw coded values retained; SPSS/Stata value labels are printed alongside frequencies."


def resolve_column(df: pd.DataFrame, labels: dict[str, str], names: list[str], label_terms: list[str] | None = None) -> str | None:
    by_clean = {clean_name(c): str(c) for c in df.columns}
    for name in names:
        if clean_name(name) in by_clean:
            return by_clean[clean_name(name)]
    if label_terms:
        for column in df.columns:
            label = labels.get(str(column), "").lower()
            if all(term.lower() in label for term in label_terms):
                return str(column)
    return None


def numeric_series(series: pd.Series) -> pd.Series:
    return pd.to_numeric(series, errors="coerce")


def is_missing(series: pd.Series) -> pd.Series:
    return series.isna() | series.astype("string").str.strip().eq("")


def frequency_text(series: pd.Series, value_labels: dict[Any, Any] | None = None, limit: int = 200) -> str:
    counts = series.map(display_value).value_counts(dropna=False).head(limit)
    lines = []
    for value, count in counts.items():
        label = ""
        if value_labels:
            for raw, value_label in value_labels.items():
                if display_value(raw) == value:
                    label = f" [{value_label}]"
                    break
        lines.append(f"{value}: {int(count)}{label}")
    if len(series.map(display_value).value_counts(dropna=False)) > limit:
        lines.append(f"... {len(series.map(display_value).value_counts(dropna=False)) - limit} more values omitted")
    return "; ".join(lines)


def markdown_table(frame: pd.DataFrame) -> str:
    if frame.empty:
        return "No rows."
    columns = [str(column) for column in frame.columns]
    rows = [[str(value).replace("|", "\\|").replace("\n", " ") for value in row] for row in frame.fillna("").itertuples(index=False, name=None)]
    widths = [max(len(column), *(len(row[index]) for row in rows)) for index, column in enumerate(columns)]
    header = "| " + " | ".join(column.ljust(widths[index]) for index, column in enumerate(columns)) + " |"
    separator = "| " + " | ".join("-" * widths[index] for index in range(len(columns))) + " |"
    body = ["| " + " | ".join(row[index].ljust(widths[index]) for index in range(len(columns))) + " |" for row in rows]
    return "\n".join([header, separator, *body])


def column_table(df: pd.DataFrame, labels: dict[str, str]) -> str:
    rows = []
    for column in df.columns:
        rows.append({"column": str(column), "label": labels.get(str(column), ""), "dtype": str(df[column].dtype)})
    return markdown_table(pd.DataFrame(rows))


def variable_audit(df: pd.DataFrame, labels: dict[str, str], value_labels: dict[str, dict[Any, Any]], columns: list[str]) -> str:
    rows = []
    for column in columns:
        if column not in df:
            continue
        series = df[column]
        missing = int(is_missing(series).sum())
        rows.append({
            "variable": column,
            "label": labels.get(column, ""),
            "dtype": str(series.dtype),
            "missing_n": missing,
            "missing_pct": round(100 * missing / len(df), 3) if len(df) else 0,
            "frequencies (raw; labels in brackets)": frequency_text(series, value_labels.get(column)),
        })
    return markdown_table(pd.DataFrame(rows)) if rows else "No matching variables found."


def parse_as(series: pd.Series) -> tuple[pd.Series, pd.Series]:
    text = series.astype("string").str.strip()
    censored = text.str.startswith("<")
    numeric = pd.to_numeric(text.str.replace(r"^[<>]=?\s*", "", regex=True), errors="coerce")
    return numeric, censored.fillna(False)


def bgs_audit(path: Path, df: pd.DataFrame, labels: dict[str, str]) -> tuple[str, pd.DataFrame]:
    def col(*names: str) -> str | None:
        return resolve_column(df, labels, list(names))

    as_col = col("As", "arsenic", "as_ug_l")
    district_col = col("DISTRICT", "district")
    id_col = col("SAMPLE_ID", "sample_id", "well_id", "sample")
    lat_col = col("LAT_DEG", "latitude", "lat")
    lon_col = col("LONG_DEG", "longitude", "lon", "lng")
    if not as_col or not district_col:
        return "BGS columns As or district were not identified.", pd.DataFrame()
    as_numeric, censored = parse_as(df[as_col])
    missing = is_missing(df[as_col])
    district = df[district_col].astype("string").fillna("<MISSING>")
    records = []
    for name, group_index in district.groupby(district).groups.items():
        values = as_numeric.loc[group_index]
        valid = values.dropna()
        records.append({
            "district": name, "total_rows": len(group_index), "valid_numeric_As": int(valid.notna().sum()),
            "censored_As": int(censored.loc[group_index].sum()), "missing_As": int(missing.loc[group_index].sum()),
            "zeros": int((valid == 0).sum()), "negatives": int((valid < 0).sum()),
            "min": valid.min(), "Q1": valid.quantile(.25), "median": valid.median(), "mean": valid.mean(),
            "SD": valid.std(), "Q3": valid.quantile(.75), "P90": valid.quantile(.90), "P95": valid.quantile(.95),
            "P99": valid.quantile(.99), "max": valid.max(),
        })
    table = pd.DataFrame(records).sort_values(["valid_numeric_As", "district"], ascending=[False, True])
    lines = ["### Identified columns", f"As: `{as_col}`; district: `{district_col}`; sample ID: `{id_col}`; coordinates: `{lat_col}`, `{lon_col}.", "", "### District arsenic audit", markdown_table(table)]
    lines.append("\nDistricts meeting usable numeric As thresholds (no validity judgment is made):")
    for threshold in [20, 30, 40, 50, 75, 100]:
        names = table.loc[table["valid_numeric_As"] >= threshold, "district"].astype(str).tolist()
        lines.append(f"- n >= {threshold}: {', '.join(names) if names else 'none'}")
    duplicate_lines = []
    if id_col:
        duplicate_lines.append(f"- Duplicate sample IDs: {int(df[id_col].duplicated(keep=False).sum())} rows across {int(df[id_col].duplicated().sum())} repeated IDs.")
    if lat_col and lon_col:
        duplicate_lines.append(f"- Duplicate coordinates: {int(df.duplicated([lat_col, lon_col], keep=False).sum())} rows across {int(df.duplicated([lat_col, lon_col]).sum())} repeated coordinate pairs.")
    lines.extend(["", "### Duplicate identifiers and coordinates", *duplicate_lines])
    lines.extend(["", "### Detection-limit observations", f"Rows with preserved censored As strings (for example `<6`): {int(censored.sum())}", f"Rows with numeric As after parsing the numeric component for audit only: {int(as_numeric.notna().sum())}"])
    return "\n".join(lines), table


def survey_audit(name: str, df: pd.DataFrame, labels: dict[str, str], value_labels: dict[str, dict[Any, Any]], kind: str) -> tuple[str, pd.DataFrame]:
    if kind == "steps":
        specs = {
            "weight": ("wstep2", ["weight"]), "age": ("age", ["age"]), "sex": ("sex", ["sex"]),
            "height_cm": ("m11", ["height", "length"]), "weight_kg": ("m12", ["weight", "kilogram"]),
        }
    else:
        specs = {
            "weight": ("chweight", ["children", "sample weight"]), "age_months": ("CAGE", ["age", "months"]),
            "sex": ("HL4", ["sex"]), "district": ("HH7A", ["district"]), "child_weight_kg": ("AN8", ["weight", "kilograms"]),
            "height_cm": ("AN11", ["length", "height"]),
        }
    resolved = {key: resolve_column(df, labels, [candidate], terms) for key, (candidate, terms) in specs.items()}
    present = [c for c in resolved.values() if c]
    lines = [f"### Identified {name} variables", ", ".join(f"{key}: `{value}`" for key, value in resolved.items() if value) or "No expected variables found.", "", "### Missingness and raw frequencies", variable_audit(df, labels, value_labels, present)]
    for key, column in resolved.items():
        if column and key in {"age", "age_months", "sex", "weight", "height_cm", "weight_kg", "child_weight_kg", "district"}:
            special = df[column].astype("string").str.contains(r"888|999|997|998", na=False)
            lines.append(f"- Raw suspicious/special-code pattern counts in `{column}`: {int(special.sum())}; values are retained and require later review.")
    required_keys = ["weight_kg", "age", "sex", "weight"] if kind == "steps" else ["child_weight_kg", "height_cm", "age_months", "sex", "district", "weight"]
    required = [resolved.get(key) for key in required_keys]
    required = [c for c in required if c]
    complete = ~pd.concat([is_missing(df[c]) for c in required], axis=1).any(axis=1) if required else pd.Series(False, index=df.index)
    lines.append(f"\nValid complete combinations of {' + '.join(required)}: {int(complete.sum())} of {len(df)}")
    completeness = pd.DataFrame([{"dataset": name, "combination": "+".join(required), "total_rows": len(df), "complete_rows": int(complete.sum()), "complete_pct": round(100 * complete.mean(), 3) if len(df) else 0}])
    if kind == "mics" and resolved.get("age_months"):
        age = numeric_series(df[resolved["age_months"]])
        usable = complete & age.notna()
        for group_name, groups in [("age_month", age), ("completed_year", np.floor(age / 12)), ("sex", df[resolved["sex"]] if resolved.get("sex") else pd.Series(index=df.index)), ("district", df[resolved["district"]] if resolved.get("district") else pd.Series(index=df.index))]:
            if groups.empty:
                continue
            counts = groups[usable].map(display_value).value_counts().rename_axis("category").reset_index(name="usable_complete_rows")
            counts.insert(0, "dataset", "MICS")
            counts.insert(1, "breakdown", group_name)
            completeness = pd.concat([completeness, counts], ignore_index=True, sort=False)
            lines.extend([f"\nUsable complete rows by {group_name}", markdown_table(counts)])
    return "\n".join(lines), completeness


def main() -> None:
    parser = argparse.ArgumentParser(description="Audit arsenic, STEPS, and MICS datasets without cleaning or imputing data.")
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parent)
    args = parser.parse_args()
    root = args.root.resolve()
    output_dir = root / "audit_outputs"
    output_dir.mkdir(exist_ok=True)
    paths = sorted(p for p in root.rglob("*") if p.is_file() and p.suffix.lower() in SUPPORTED and p.name != Path(__file__).name and "audit_outputs" not in p.parts)
    sections = ["# Dataset Audit Report", "", f"Generated: {pd.Timestamp.now().isoformat()}", "", "This report is descriptive only. No source values were removed, replaced, imputed, transformed in place, or judged invalid. Numeric As parsing is used only to summarize censored and numeric observations separately."]
    district_tables = []
    completeness_tables = []
    loaded: list[tuple[Path, pd.DataFrame, dict[str, str], dict[str, dict[Any, Any]], str]] = []
    for path in paths:
        try:
            df, labels, value_labels, notes = read_table(path)
            loaded.append((path, df, labels, value_labels, notes))
        except Exception as exc:
            sections.extend([f"\n## {path.relative_to(root)}", f"Could not read file: `{type(exc).__name__}: {exc}`"])
    for path, df, labels, value_labels, notes in loaded:
        relative = path.relative_to(root).as_posix()
        sections.extend([f"\n## {relative}", f"Shape: {df.shape[0]} rows x {df.shape[1]} columns.", notes, "", "### Complete column-name, label, and dtype table", column_table(df, labels)])
        if path.name.lower() == "bgd2018.csv" or resolve_column(df, labels, ["pid"]) and resolve_column(df, labels, ["wstep2"]):
            text, completeness = survey_audit("STEPS 2018", df, labels, value_labels, "steps")
            sections.append(text)
            completeness_tables.append(completeness)
        elif path.name.lower() == "nationalsurveydata.csv" or resolve_column(df, labels, ["As", "arsenic"]):
            text, table = bgs_audit(path, df, labels)
            sections.append(text)
            if not table.empty:
                table.insert(0, "source_file", relative)
                district_tables.append(table)
        elif path.name.lower() == "ch.sav" or all(resolve_column(df, labels, names, terms) for names, terms in [(["AN8"], ["weight", "kilograms"]), (["AN11"], ["length", "height"]), (["CAGE"], ["age", "months"]), (["HL4"], ["sex"]), (["HH7A"], ["district"]), (["chweight"], ["sample", "weight"])]):
            text, table = survey_audit("MICS", df, labels, value_labels, "mics")
            sections.append(text)
            completeness_tables.append(table)
        else:
            relevant = [str(c) for c in df.columns if any(term in clean_name(c) for term in ["district", "arsenic", "weight", "height", "age", "sex", "psu", "stratum"])]
            sections.extend(["### Relevant-variable audit", variable_audit(df, labels, value_labels, relevant)])
    sections.extend(["\n## Decision Inputs", "", "Use the district numeric-As counts, censored-As counts, missingness, duplicate identifiers/coordinates, and quantiles above to choose up to 10 districts. For distribution fitting, decide separately how censored observations are handled and whether district counts meet the desired threshold; this script does not make that decision.", "", "For adult BW preprocessing, compare STEPS complete weight+age+sex+survey-weight counts, raw special-code frequencies, and measurement ranges. For child BW preprocessing, compare MICS complete child-weight+height+age-months+sex+district+survey-weight counts by age month, completed year, sex, and district. Extreme values are flagged for later review, never automatically rejected."])
    report = "\n".join(sections) + "\n"
    md_path = root / "dataset_audit_report.md"
    txt_path = root / "dataset_audit_report.txt"
    md_path.write_text(report, encoding="utf-8")
    txt_path.write_text(report, encoding="utf-8")
    district_path = output_dir / "district_counts.csv"
    completeness_path = output_dir / "completeness_tables.csv"
    pd.concat(district_tables, ignore_index=True).to_csv(district_path, index=False) if district_tables else pd.DataFrame().to_csv(district_path, index=False)
    pd.concat(completeness_tables, ignore_index=True).to_csv(completeness_path, index=False) if completeness_tables else pd.DataFrame().to_csv(completeness_path, index=False)
    print(f"Generated: {md_path}")
    print(f"Generated: {txt_path}")
    print(f"Generated: {district_path}")
    print(f"Generated: {completeness_path}")


if __name__ == "__main__":
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        main()