import csv
import json
import re
from pathlib import Path
from typing import Any


def _resolve(path: str, root: Path) -> Path:
    return root / path


def _first_float(text: str) -> float:
    match = re.search(r"[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?", text)
    if match is None:
        raise AssertionError("No numeric value found in output text.")
    return float(match.group(0))


def _validate_exists(spec: dict[str, Any], out_dir: Path, _expected_dir: Path) -> None:
    path = _resolve(spec["path"], out_dir)
    assert path.exists(), f"missing expected output: {path}"


def _validate_nonempty(spec: dict[str, Any], out_dir: Path, _expected_dir: Path) -> None:
    path = _resolve(spec["path"], out_dir)
    assert path.exists(), f"missing expected output: {path}"
    assert path.stat().st_size > 0, f"empty expected output: {path}"


def _validate_contains_text(spec: dict[str, Any], out_dir: Path, _expected_dir: Path) -> None:
    path = _resolve(spec["path"], out_dir)
    needle = spec["text"]
    contents = path.read_text(encoding="utf-8")
    assert needle in contents, f"'{needle}' not found in {path}"


def _validate_exact_text(spec: dict[str, Any], out_dir: Path, expected_dir: Path) -> None:
    actual_path = _resolve(spec["path"], out_dir)
    expected_path = _resolve(spec["expected_path"], expected_dir)
    actual = actual_path.read_text(encoding="utf-8")
    expected = expected_path.read_text(encoding="utf-8")
    if spec.get("strip", True):
        actual = actual.strip()
        expected = expected.strip()
    assert actual == expected, f"text mismatch for {actual_path} against {expected_path}"


def _validate_json_equals(spec: dict[str, Any], out_dir: Path, expected_dir: Path) -> None:
    actual_path = _resolve(spec["path"], out_dir)
    expected_path = _resolve(spec["expected_path"], expected_dir)
    actual = json.loads(actual_path.read_text(encoding="utf-8"))
    expected = json.loads(expected_path.read_text(encoding="utf-8"))
    assert actual == expected, f"JSON mismatch for {actual_path} against {expected_path}"


def _validate_csv_columns(spec: dict[str, Any], out_dir: Path, _expected_dir: Path) -> None:
    path = _resolve(spec["path"], out_dir)
    expected_columns = spec["columns"]
    with path.open("r", encoding="utf-8", newline="") as fh:
        reader = csv.reader(fh)
        header = next(reader)
    assert header == expected_columns, f"CSV header mismatch for {path}: {header} != {expected_columns}"


def _validate_numeric_tolerance(spec: dict[str, Any], out_dir: Path, expected_dir: Path) -> None:
    actual_path = _resolve(spec["path"], out_dir)
    actual_text = actual_path.read_text(encoding="utf-8")
    actual_value = _first_float(actual_text)

    if "expected_path" in spec:
        expected_text = _resolve(spec["expected_path"], expected_dir).read_text(encoding="utf-8")
        expected_value = _first_float(expected_text)
    else:
        expected_value = float(spec["expected"])

    atol = float(spec.get("atol", 0.0))
    rtol = float(spec.get("rtol", 0.0))
    tolerance = atol + abs(expected_value) * rtol
    diff = abs(actual_value - expected_value)
    assert (
        diff <= tolerance
    ), f"numeric mismatch for {actual_path}: actual={actual_value}, expected={expected_value}, tolerance={tolerance}"


VALIDATORS = {
    "exists": _validate_exists,
    "nonempty": _validate_nonempty,
    "contains_text": _validate_contains_text,
    "exact_text": _validate_exact_text,
    "json_equals": _validate_json_equals,
    "csv_columns": _validate_csv_columns,
    "numeric_tolerance": _validate_numeric_tolerance,
}


def run_validations(validations: list[dict[str, Any]], out_dir: Path, expected_dir: Path) -> None:
    for spec in validations:
        kind = spec["kind"]
        if kind not in VALIDATORS:
            raise ValueError(f"Unknown validator kind: {kind}")
        VALIDATORS[kind](spec, out_dir, expected_dir)
