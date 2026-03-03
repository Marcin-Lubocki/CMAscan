"""
Multi-FASTA motif extraction and scoring pipeline.

Designed for Google Colab usage and easy parameter tuning.
The pipeline validates input extension, parses protein FASTA records,
finds motif matches, computes four scores, and exports a CSV file.
"""

from __future__ import annotations

import argparse
import csv
import os
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, Iterator, List, Sequence, Tuple


# =========================
# User-configurable section
# =========================

ALLOWED_EXTENSIONS = {".fasta", ".faa", ".txt"}

# Default motif catalog.
# Add or edit entries here to adapt the pipeline.
MOTIF_DEFINITIONS = [
    {
        "name": "SH3-like_PxxP",
        "pattern": r"P..P",
        "type": "SH3_binding_like",
    },
    {
        "name": "Proline_directed_[ST]P",
        "pattern": r"[ST]P",
        "type": "proline_directed_site",
    },
    {
        "name": "N_glycosylation_N[^P][ST][^P]",
        "pattern": r"N[^P][ST][^P]",
        "type": "N_glycosylation_like",
    },
]

# Amino-acid disorder propensities (0-1 range).
# Used as a deterministic proxy for IUPred3-like scoring.
DISORDER_PROPENSITY = {
    "A": 0.38,
    "C": 0.35,
    "D": 0.74,
    "E": 0.71,
    "F": 0.31,
    "G": 0.72,
    "H": 0.53,
    "I": 0.22,
    "K": 0.67,
    "L": 0.24,
    "M": 0.30,
    "N": 0.66,
    "P": 0.75,
    "Q": 0.69,
    "R": 0.64,
    "S": 0.68,
    "T": 0.60,
    "V": 0.23,
    "W": 0.29,
    "Y": 0.41,
}

# MusiteDeep-like heuristic weights.
MUSITEDEEP_WEIGHTS = {
    "contains_sty": 0.45,
    "plus1_is_proline": 0.25,
    "minus3_is_basic": 0.15,
    "minus2_is_basic": 0.10,
    "nearby_acidic_bonus": 0.05,
}

VALID_AA = set("ACDEFGHIKLMNPQRSTVWY")


@dataclass(frozen=True)
class FastaRecord:
    header: str
    sequence: str

    @property
    def protein_name(self) -> str:
        # First token in header is generally the most stable identifier.
        return self.header.split()[0]


def validate_input_file(path: Path) -> None:
    if not path.exists():
        raise FileNotFoundError(f"Input file not found: {path}")

    if path.suffix.lower() not in ALLOWED_EXTENSIONS:
        allowed = ", ".join(sorted(ALLOWED_EXTENSIONS))
        raise ValueError(
            f"Invalid input format: '{path.suffix}'. "
            f"Allowed extensions: {allowed}."
        )


def parse_multifasta(path: Path) -> List[FastaRecord]:
    records: List[FastaRecord] = []
    current_header: str | None = None
    current_seq_chunks: List[str] = []

    with path.open("r", encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue

            if line.startswith(">"):
                if current_header is not None:
                    sequence = "".join(current_seq_chunks).upper()
                    records.append(FastaRecord(header=current_header, sequence=sequence))
                current_header = line[1:].strip()
                current_seq_chunks = []
            else:
                current_seq_chunks.append(line)

    if current_header is not None:
        sequence = "".join(current_seq_chunks).upper()
        records.append(FastaRecord(header=current_header, sequence=sequence))

    if not records:
        raise ValueError("No FASTA records detected. Ensure file is in valid multi-FASTA format.")

    return records


def sanity_check_records(records: Sequence[FastaRecord], allow_x: bool = True) -> None:
    invalid_entries: List[str] = []

    for rec in records:
        if not rec.sequence:
            invalid_entries.append(f"{rec.protein_name}: empty sequence")
            continue

        for aa in rec.sequence:
            if aa in VALID_AA:
                continue
            if allow_x and aa == "X":
                continue
            invalid_entries.append(f"{rec.protein_name}: invalid residue '{aa}'")
            break

    if invalid_entries:
        preview = "; ".join(invalid_entries[:5])
        raise ValueError(f"Sanity check failed. Invalid sequence content detected: {preview}")


def iter_motif_matches(sequence: str, pattern: str) -> Iterator[Tuple[str, int]]:
    # Lookahead enables overlapping matches.
    lookahead = re.compile(rf"(?=({pattern}))")
    for match in lookahead.finditer(sequence):
        motif = match.group(1)
        start_idx = match.start(1)
        yield motif, start_idx


def tokenize_pattern(pattern: str) -> List[str]:
    tokens: List[str] = []
    i = 0
    while i < len(pattern):
        ch = pattern[i]
        if ch == "[":
            j = pattern.find("]", i + 1)
            if j == -1:
                raise ValueError(f"Invalid motif pattern (missing ']'): {pattern}")
            tokens.append(pattern[i : j + 1])
            i = j + 1
        else:
            tokens.append(ch)
            i += 1
    return tokens


def score_pattern_token(token: str, aa: str) -> float:
    # Returns a signed local contribution.
    if token == ".":
        return 0.2

    if token.startswith("[") and token.endswith("]"):
        body = token[1:-1]
        if body.startswith("^"):
            excluded = set(body[1:])
            return 0.8 if aa not in excluded else -0.8
        allowed = set(body)
        return 0.8 if aa in allowed else -0.8

    return 1.0 if aa == token else -1.0


def compute_cpssm(mer: str, pattern: str) -> float:
    tokens = tokenize_pattern(pattern)
    if len(tokens) != len(mer):
        return 0.0

    raw = sum(score_pattern_token(tok, aa) for tok, aa in zip(tokens, mer))
    min_score = -1.0 * len(tokens)
    max_score = 1.0 * len(tokens)
    return normalize(raw, min_score, max_score)


def compute_epssm(sequence: str, start_idx: int, motif_len: int, flank: int = 6) -> float:
    left = max(0, start_idx - flank)
    right = min(len(sequence), start_idx + motif_len + flank)
    window = sequence[left:right]

    if not window:
        return 0.0

    acidic = sum(1 for aa in window if aa in {"D", "E"})
    basic = sum(1 for aa in window if aa in {"K", "R", "H"})
    polar = sum(1 for aa in window if aa in {"S", "T", "N", "Q"})

    # Contextual enrichment-like signal.
    raw = 0.45 * (acidic / len(window)) + 0.35 * (basic / len(window)) + 0.20 * (polar / len(window))
    return clamp(raw, 0.0, 1.0)


def compute_iupred3_proxy(mer: str) -> float:
    vals = [DISORDER_PROPENSITY.get(aa, 0.5) for aa in mer]
    return sum(vals) / len(vals) if vals else 0.0


def compute_musitedeep_proxy(sequence: str, start_idx: int, motif_len: int) -> float:
    motif = sequence[start_idx : start_idx + motif_len]
    score = 0.0

    if any(aa in {"S", "T", "Y"} for aa in motif):
        score += MUSITEDEEP_WEIGHTS["contains_sty"]

    plus1 = sequence[start_idx + motif_len] if (start_idx + motif_len) < len(sequence) else None
    if plus1 == "P":
        score += MUSITEDEEP_WEIGHTS["plus1_is_proline"]

    minus3 = sequence[start_idx - 3] if (start_idx - 3) >= 0 else None
    minus2 = sequence[start_idx - 2] if (start_idx - 2) >= 0 else None

    if minus3 in {"K", "R", "H"}:
        score += MUSITEDEEP_WEIGHTS["minus3_is_basic"]
    if minus2 in {"K", "R", "H"}:
        score += MUSITEDEEP_WEIGHTS["minus2_is_basic"]

    flank_left = max(0, start_idx - 2)
    flank_right = min(len(sequence), start_idx + motif_len + 2)
    flank = sequence[flank_left:flank_right]
    acidic_count = sum(1 for aa in flank if aa in {"D", "E"})
    if acidic_count >= 2:
        score += MUSITEDEEP_WEIGHTS["nearby_acidic_bonus"]

    return clamp(score, 0.0, 1.0)


def normalize(value: float, min_v: float, max_v: float) -> float:
    if max_v <= min_v:
        return 0.0
    return clamp((value - min_v) / (max_v - min_v), 0.0, 1.0)


def clamp(value: float, lo: float, hi: float) -> float:
    return max(lo, min(value, hi))


def find_all_motif_hits(records: Sequence[FastaRecord], motif_defs: Sequence[Dict[str, str]]) -> List[Dict[str, str | int | float]]:
    rows: List[Dict[str, str | int | float]] = []

    for rec in records:
        for motif_def in motif_defs:
            pattern = motif_def["pattern"]
            motif_type = motif_def["type"]

            for mer, start_idx in iter_motif_matches(rec.sequence, pattern):
                row = {
                    "protein_name": rec.protein_name,
                    "mer": mer,
                    "type": motif_type,
                    "localization": start_idx + 1,  # 1-based indexing requested.
                    "cPSSM": round(compute_cpssm(mer, pattern), 4),
                    "ePSSM": round(compute_epssm(rec.sequence, start_idx, len(mer)), 4),
                    "IUPred3": round(compute_iupred3_proxy(mer), 4),
                    "MusiteDeep": round(compute_musitedeep_proxy(rec.sequence, start_idx, len(mer)), 4),
                }
                rows.append(row)

    return rows


def write_csv(rows: Sequence[Dict[str, str | int | float]], output_path: Path) -> None:
    columns = [
        "protein_name",
        "mer",
        "type",
        "localization",
        "cPSSM",
        "ePSSM",
        "IUPred3",
        "MusiteDeep",
    ]

    with output_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns)
        writer.writeheader()
        writer.writerows(rows)


def run_pipeline(input_path: str, output_csv: str = "motif_hits_scored.csv") -> Path:
    in_path = Path(input_path)
    out_path = Path(output_csv)

    validate_input_file(in_path)
    records = parse_multifasta(in_path)
    sanity_check_records(records, allow_x=True)

    rows = find_all_motif_hits(records, MOTIF_DEFINITIONS)
    write_csv(rows, out_path)

    print(f"Processed proteins: {len(records)}")
    print(f"Motif hits exported: {len(rows)}")
    print(f"Output CSV: {out_path.resolve()}")

    if len(rows) == 0:
        print("Warning: No motifs detected for the current motif catalog.")

    return out_path


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Extract motifs from multi-FASTA and score hits.")
    parser.add_argument("--input", required=True, help="Path to input multi-FASTA (.fasta/.faa/.txt).")
    parser.add_argument("--output", default="motif_hits_scored.csv", help="Output CSV filename.")
    return parser


def main() -> None:
    parser = build_arg_parser()
    args = parser.parse_args()
    run_pipeline(input_path=args.input, output_csv=args.output)


if __name__ == "__main__":
    main()
