from pathlib import Path
from typing import Iterable, Optional, Sequence, Union

PathLike = Union[str, Path]


def normalize_path(path: Optional[PathLike]) -> Optional[Path]:
    """Return a Path object for a path-like value, or None."""
    if path is None:
        return None
    return Path(path)


def normalize_paths(paths: Optional[Sequence[PathLike]]) -> list[Path]:
    """Return a list of Path objects from a sequence of path-like values."""
    if paths is None:
        return []
    return [Path(p) for p in paths]


def require_exactly_one(**kwargs: object) -> str:
    """
    Validate that exactly one named option is provided.

    Returns the key name that was provided.
    """
    provided = [name for name, value in kwargs.items() if value is not None]
    if len(provided) != 1:
        options = ", ".join(kwargs.keys())
        raise ValueError(f"Provide exactly one of: {options}")
    return provided[0]


def ensure_output_file_parent(output_file: PathLike) -> Path:
    """Create parent directory for an output file and return normalized path."""
    output = Path(output_file)
    output.parent.mkdir(parents=True, exist_ok=True)
    return output


def ensure_output_dir(output_dir: PathLike) -> Path:
    """Create an output directory and return normalized path."""
    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)
    return out


def ensure_output_prefix_parent(output_prefix: PathLike) -> Path:
    """Create parent directory for an output prefix and return normalized path."""
    prefix = Path(output_prefix)
    prefix.parent.mkdir(parents=True, exist_ok=True)
    return prefix


def validate_suffix(path: PathLike, allowed_suffixes: Iterable[str], label: str = "path") -> Path:
    """Ensure a path uses one of the allowed suffixes."""
    p = Path(path)
    allowed = tuple(allowed_suffixes)
    if p.suffix not in allowed:
        allowed_str = ", ".join(allowed)
        raise ValueError(f"{label} must have one of suffixes: {allowed_str}")
    return p
