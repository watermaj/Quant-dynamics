from __future__ import annotations

import argparse
import hashlib
from pathlib import Path
import difflib

ROOT = Path(__file__).resolve().parents[2]


def resolve_path(raw: str) -> Path:
    path = Path(raw)
    if path.is_absolute():
        return path
    return ROOT / path


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def read_text_lines(path: Path) -> list[str] | None:
    for enc in ("utf-8-sig", "utf-8", "latin-1"):
        try:
            with path.open("r", encoding=enc, newline="") as f:
                return f.readlines()
        except UnicodeDecodeError:
            continue
    return None


def first_text_diff(a: Path, b: Path) -> str:
    a_lines = read_text_lines(a)
    b_lines = read_text_lines(b)
    if a_lines is None or b_lines is None:
        return "binary/content differs"

    for idx, (la, lb) in enumerate(zip(a_lines, b_lines), start=1):
        if la != lb:
            a_visible = la.rstrip("\r\n")
            b_visible = lb.rstrip("\r\n")
            if a_visible == b_visible:
                return (
                    f"line {idx} differs only by invisible characters\n"
                    f"  A repr: {la!r}\n"
                    f"  B repr: {lb!r}"
                )
            return (
                f"line {idx}\n"
                f"  A: {a_visible}\n"
                f"  B: {b_visible}"
            )

    if len(a_lines) != len(b_lines):
        return f"line count differs (A={len(a_lines)}, B={len(b_lines)})"

    # Fallback: if hashes differ but no line-level difference found for some reason.
    diff = difflib.unified_diff(a_lines, b_lines, n=1)
    preview = "".join(list(diff)[:8]).strip()
    return preview if preview else "content differs"


def compare_files(a: Path, b: Path) -> int:
    if not a.exists():
        print(f"Missing A: {a}")
        return 1
    if not b.exists():
        print(f"Missing B: {b}")
        return 1
    if a.is_dir() or b.is_dir():
        print("Both paths must be files for file mode.")
        return 2

    if sha256(a) == sha256(b):
        print("Files are identical.")
        return 0

    print("Files differ:")
    print(first_text_diff(a, b))
    return 1


def compare_dirs(a_dir: Path, b_dir: Path) -> int:
    if not a_dir.exists() or not a_dir.is_dir():
        print(f"Missing or invalid directory A: {a_dir}")
        return 2
    if not b_dir.exists() or not b_dir.is_dir():
        print(f"Missing or invalid directory B: {b_dir}")
        return 2

    a_files = {p.relative_to(a_dir) for p in a_dir.rglob("*") if p.is_file()}
    b_files = {p.relative_to(b_dir) for p in b_dir.rglob("*") if p.is_file()}

    only_a = sorted(a_files - b_files)
    only_b = sorted(b_files - a_files)
    common = sorted(a_files & b_files)

    if only_a:
        print("Only in A:")
        for rel in only_a:
            print(f"  {rel}")
    if only_b:
        print("Only in B:")
        for rel in only_b:
            print(f"  {rel}")

    diff_count = 0
    for rel in common:
        pa = a_dir / rel
        pb = b_dir / rel
        if sha256(pa) != sha256(pb):
            diff_count += 1
            print(f"\nDifferent: {rel}")
            print(first_text_diff(pa, pb))

    if not only_a and not only_b and diff_count == 0:
        print("No differences found.")
        return 0

    total = len(only_a) + len(only_b) + diff_count
    print(f"\nTotal differences: {total}")
    return 1


def main() -> int:
    parser = argparse.ArgumentParser(
        description=(
            "Compare two files or two directories and print where content differs. "
            "Defaults to results/Submitted_results vs tests/regression/reference."
        )
    )
    parser.add_argument("path_a", nargs="?", default="results/Submitted_results")
    parser.add_argument("path_b", nargs="?", default="tests/regression/reference")
    args = parser.parse_args()

    a = resolve_path(args.path_a)
    b = resolve_path(args.path_b)

    if a.is_file() or b.is_file():
        return compare_files(a, b)
    return compare_dirs(a, b)


if __name__ == "__main__":
    raise SystemExit(main())
