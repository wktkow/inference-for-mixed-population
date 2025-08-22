#!/usr/bin/env python3

"""
Splice PDFs in a directory to a maximum number of pages and write new files
with a "-splice" suffix next to the originals. Non-destructive: inputs remain.

Usage (default directory is the repo's pdf_source):
  python3 pdf_splice.py --dir /path/to/pdf_source --max-pages 210 --suffix -splice
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import List

try:
    from pypdf import PdfReader, PdfWriter
except Exception as import_error:  # noqa: F841
    print("Missing dependency 'pypdf'. Install with: python3 -m pip install --user pypdf", file=sys.stderr)
    sys.exit(1)


def find_pdfs(directory: Path) -> List[Path]:
    """Return a sorted list of PDF files in the given directory."""
    return sorted(directory.glob("*.pdf"))


def splice_pdf(input_path: Path, output_path: Path, max_pages: int) -> int:
    """Write a copy of input_path truncated to max_pages pages into output_path.

    Returns the total number of pages in the input PDF.
    """
    reader = PdfReader(str(input_path))
    writer = PdfWriter()

    total_pages = len(reader.pages)
    limit_pages = min(max_pages, total_pages)

    for page_index in range(limit_pages):
        writer.add_page(reader.pages[page_index])

    with output_path.open("wb") as output_file:
        writer.write(output_file)

    # Some PdfWriter implementations expose close(); ignore if unavailable
    try:  # pragma: no cover - best-effort cleanup
        writer.close()  # type: ignore[attr-defined]
    except Exception:
        pass

    return total_pages


def main() -> int:
    parser = argparse.ArgumentParser(description="Splice PDFs to a maximum number of pages.")
    default_dir = Path(__file__).resolve().parent / "pdf_source"
    parser.add_argument(
        "--dir",
        dest="directory",
        type=Path,
        default=default_dir,
        help="Directory containing PDFs to splice (default: repo's pdf_source)",
    )
    parser.add_argument(
        "--max-pages",
        dest="max_pages",
        type=int,
        default=210,
        help="Maximum number of pages to keep (default: 210)",
    )
    parser.add_argument(
        "--suffix",
        dest="suffix",
        type=str,
        default="-splice",
        help="Suffix to append before .pdf for outputs (default: -splice)",
    )

    args = parser.parse_args()

    directory: Path = args.directory
    if not directory.exists():
        print(f"Directory not found: {directory}", file=sys.stderr)
        return 2

    pdf_paths = find_pdfs(directory)
    if not pdf_paths:
        print(f"No PDFs found in {directory}")
        return 0

    for pdf_path in pdf_paths:
        output_path = pdf_path.with_name(pdf_path.stem + args.suffix + pdf_path.suffix)
        try:
            total_pages = splice_pdf(pdf_path, output_path, args.max_pages)
            kept_pages = min(args.max_pages, total_pages)
            print(f"Wrote {output_path} ({kept_pages}/{total_pages} pages)")
        except Exception as exc:  # noqa: BLE001
            print(f"Failed to process {pdf_path}: {exc}", file=sys.stderr)

    return 0


if __name__ == "__main__":
    sys.exit(main())


