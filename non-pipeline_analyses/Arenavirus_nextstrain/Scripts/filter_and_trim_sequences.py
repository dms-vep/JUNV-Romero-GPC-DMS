"""
filter_and_trim_sequences.py
----------------------------
Snakemake script that:
  1. Aligns each query sequence (forward + reverse complement) against a
     reference using pairwise2 (Biopython) or parasail (faster, optional).
  2. Keeps sequences that meet a percent-identity threshold.
  3. Extracts only the portion of each sequence that aligns with the reference
     (trims overhanging ends).
  4. Writes a filtered .fasta and a filtered .tsv metadata file.

Expected Snakemake rule (example)
----------------------------------
rule filter_and_trim:
    input:
        sequences = "data/sequences.fasta",
        metadata  = "data/metadata.tsv",
        reference = "data/reference.fasta",
    output:
        sequences = "results/filtered.fasta",
        metadata  = "results/filtered_metadata.tsv",
    params:
        identity_threshold = 0.90,   # 0–1 float
        metadata_id_col    = "strain" # column in the TSV that matches FASTA IDs
    script:
        "scripts/filter_and_trim_sequences.py"

Dependencies
------------
    pip install biopython

Optional faster aligner (recommended for large datasets):
    pip install parasail
"""

import sys
import csv
import logging
from pathlib import Path

# ---------------------------------------------------------------------------
# Logging
# ---------------------------------------------------------------------------
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s  %(levelname)-8s  %(message)s",
    datefmt="%H:%M:%S",
    stream=sys.stderr,
)
log = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Try to import alignment backends
# ---------------------------------------------------------------------------
try:
    import parasail  # fast SIMD Smith–Waterman
    _HAS_PARASAIL = True
    log.info("Using parasail for alignments (fast)")
except ImportError:
    _HAS_PARASAIL = False
    log.info("parasail not found – falling back to Biopython pairwise2 (slower)")

from Bio import SeqIO
from Bio.Seq import Seq


# ---------------------------------------------------------------------------
# Alignment helpers
# ---------------------------------------------------------------------------

def _align_parasail(query_str: str, ref_str: str):
    """
    Semi-global alignment (query ends free) via parasail.
    Returns (identity_fraction, aligned_query_subsequence).
    """
    matrix = parasail.blosum62  # nucleotide; NUC44 is also fine
    matrix = parasail.matrix_create("ACGTN", 2, -1)
    result = parasail.sg_trace_striped_32(query_str, ref_str, 10, 1, matrix)

    cigar = result.cigar
    ref_beg = result.cigar.beg_ref
    ref_end = result.end_ref
    query_beg = result.cigar.beg_query
    query_end = result.end_query

    # Count matches in the aligned region
    matches = result.matches
    length  = result.length  # alignment length (including gaps)

    identity = matches / length if length > 0 else 0.0

    # Extract query bases that cover the reference span
    trimmed_query = query_str[query_beg : query_end + 1]
    return identity, trimmed_query


def _align_biopython(query_str: str, ref_str: str):
    """
    Semi-global alignment via Biopython PairwiseAligner.
    Returns (identity_fraction, aligned_query_subsequence).
    """
    from Bio import Align

    aligner = Align.PairwiseAligner()
    aligner.mode = "global"

    # Score settings favouring nucleotide comparison
    aligner.match_score    =  2
    aligner.mismatch_score = -1
    aligner.open_gap_score    = -2
    aligner.extend_gap_score  = -0.5

    # Make query-end gaps free (semi-global: query can extend past reference)
    aligner.query_left_open_gap_score    = 0
    aligner.query_left_extend_gap_score  = 0
    aligner.query_right_open_gap_score   = 0
    aligner.query_right_extend_gap_score = 0

    alignments = aligner.align(ref_str, query_str)
    if not alignments:
        return 0.0, ""

    best = alignments[0]

    # Alignment coordinates: best.aligned is ((ref_blocks), (query_blocks))
    ref_blocks   = best.aligned[0]
    query_blocks = best.aligned[1]

    if len(ref_blocks) == 0:
        return 0.0, ""

    # Reference span covered
    ref_start  = ref_blocks[0][0]
    ref_end    = ref_blocks[-1][1]

    # Corresponding query span
    query_start = query_blocks[0][0]
    query_end   = query_blocks[-1][1]

    # Count identical positions
    matches = sum(
        1
        for rb, qb in zip(ref_blocks, query_blocks)
        for ri, qi in zip(range(rb[0], rb[1]), range(qb[0], qb[1]))
        if ref_str[ri] == query_str[qi]
    )

    alignment_len = ref_end - ref_start  # length of covered reference region
    identity = matches / alignment_len if alignment_len > 0 else 0.0
    trimmed_query = query_str[query_start:query_end]

    return identity, trimmed_query


def align(query_str: str, ref_str: str):
    """Dispatch to the best available aligner."""
    if _HAS_PARASAIL:
        return _align_parasail(query_str, ref_str)
    return _align_biopython(query_str, ref_str)


def best_orientation(query_seq: Seq, ref_str: str):
    """
    Align forward and reverse-complement; return the better alignment.
    Returns (identity, trimmed_seq_string, orientation).
    orientation is '+' or '-'.
    """
    fwd_str = str(query_seq)
    rev_str = str(query_seq.reverse_complement())

    fwd_identity, fwd_trimmed = align(fwd_str, ref_str)
    rev_identity, rev_trimmed = align(rev_str, ref_str)

    if fwd_identity >= rev_identity:
        return fwd_identity, fwd_trimmed, "+"
    return rev_identity, rev_trimmed, "-"


# ---------------------------------------------------------------------------
# Main logic
# ---------------------------------------------------------------------------

def main():
    # ------------------------------------------------------------------
    # Pull parameters from Snakemake object (or fall back to CLI for testing)
    # ------------------------------------------------------------------
    try:
        input_fasta     = str(snakemake.input.sequences)   # noqa: F821
        input_metadata  = str(snakemake.input.metadata)
        input_reference = str(snakemake.input.reference)
        output_fasta    = str(snakemake.output.sequences)
        output_metadata = str(snakemake.output.metadata)
        threshold       = float(snakemake.params.get("identity_threshold", 0.90))
        id_col          = str(snakemake.params.get("metadata_id_col", "strain"))
    except NameError:
        # Fallback: allow running directly for testing
        import argparse
        p = argparse.ArgumentParser()
        p.add_argument("--sequences",          required=True)
        p.add_argument("--metadata",           required=True)
        p.add_argument("--reference",          required=True)
        p.add_argument("--output-sequences",   required=True)
        p.add_argument("--output-metadata",    required=True)
        p.add_argument("--identity-threshold", type=float, default=0.90)
        p.add_argument("--metadata-id-col",    default="strain")
        args = p.parse_args()
        input_fasta     = args.sequences
        input_metadata  = args.metadata
        input_reference = args.reference
        output_fasta    = args.output_sequences
        output_metadata = args.output_metadata
        threshold       = args.identity_threshold
        id_col          = args.metadata_id_col

    log.info("Input sequences : %s", input_fasta)
    log.info("Input metadata  : %s", input_metadata)
    log.info("Reference       : %s", input_reference)
    log.info("Identity cutoff : %.1f%%", threshold * 100)

    # ------------------------------------------------------------------
    # Load reference (first record only)
    # ------------------------------------------------------------------
    ref_records = list(SeqIO.parse(input_reference, "fasta"))
    if not ref_records:
        raise ValueError(f"No sequences found in reference file: {input_reference}")
    if len(ref_records) > 1:
        log.warning(
            "Reference file contains %d sequences; using only the first (%s).",
            len(ref_records),
            ref_records[0].id,
        )
    ref_seq = ref_records[0]
    ref_str = str(ref_seq.seq).upper()
    log.info("Reference: %s  (%d bp)", ref_seq.id, len(ref_str))

    # ------------------------------------------------------------------
    # Load metadata into a dict keyed by sequence ID
    # ------------------------------------------------------------------
    metadata_rows: dict[str, dict] = {}
    fieldnames: list[str] = []

    with open(input_metadata, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError("Metadata TSV appears to be empty or has no header.")
        fieldnames = list(reader.fieldnames)
        if id_col not in fieldnames:
            raise ValueError(
                f"ID column '{id_col}' not found in metadata. "
                f"Available columns: {fieldnames}"
            )
        for row in reader:
            metadata_rows[row[id_col]] = dict(row)

    log.info("Loaded %d metadata rows", len(metadata_rows))

    # ------------------------------------------------------------------
    # Process sequences
    # ------------------------------------------------------------------
    Path(output_fasta).parent.mkdir(parents=True, exist_ok=True)
    Path(output_metadata).parent.mkdir(parents=True, exist_ok=True)

    kept_ids: list[str] = []
    total = 0

    with open(output_fasta, "w") as fasta_out:
        for record in SeqIO.parse(input_fasta, "fasta"):
            total += 1
            seq_id = record.id

            identity, trimmed, orientation = best_orientation(
                record.seq.upper(), ref_str
            )

            if identity < threshold:
                log.debug(
                    "SKIP  %-40s  identity=%.2f%%  (< %.1f%%)",
                    seq_id, identity * 100, threshold * 100,
                )
                continue

            log.info(
                "KEEP  %-40s  identity=%.2f%%  orientation=%s  trimmed_len=%d",
                seq_id, identity * 100, orientation, len(trimmed),
            )

            if not trimmed:
                log.warning("Trimmed sequence is empty for %s – skipping.", seq_id)
                continue

            # Write trimmed sequence; append orientation tag to description
            fasta_out.write(f">{seq_id} orientation={orientation}\n")
            # Wrap at 60 characters
            for i in range(0, len(trimmed), 60):
                fasta_out.write(trimmed[i:i+60] + "\n")

            kept_ids.append(seq_id)

    log.info("Kept %d / %d sequences (%.1f%%)", len(kept_ids), total,
             100 * len(kept_ids) / total if total else 0)

    # ------------------------------------------------------------------
    # Write filtered metadata
    # ------------------------------------------------------------------
    kept_set = set(kept_ids)
    missing_meta = kept_set - set(metadata_rows.keys())
    if missing_meta:
        log.warning(
            "%d kept sequence(s) have no metadata row: %s",
            len(missing_meta),
            ", ".join(sorted(missing_meta)[:10]),
        )

    with open(output_metadata, "w", newline="") as meta_out:
        writer = csv.DictWriter(meta_out, fieldnames=fieldnames, delimiter="\t",
                                extrasaction="ignore")
        writer.writeheader()
        for seq_id in kept_ids:
            if seq_id in metadata_rows:
                writer.writerow(metadata_rows[seq_id])

    log.info("Done. Filtered FASTA → %s", output_fasta)
    log.info("      Filtered TSV   → %s", output_metadata)


if __name__ == "__main__" or "snakemake" in dir():
    main()
