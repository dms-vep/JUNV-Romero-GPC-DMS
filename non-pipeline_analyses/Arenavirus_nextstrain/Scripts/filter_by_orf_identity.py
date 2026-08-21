"""
Snakemake rule script: filter_by_orf_identity.py

Reads a multi-sequence FASTA and its TSV metadata, extracts ORFs from each
sequence using EMBOSS getorf, aligns every candidate ORF to the reference ORF,
and selects the best-matching one (Options B + D):

  Option B – test ALL candidate ORFs and keep the one with the highest
             percent similarity to the reference, rather than just the longest.
  Option D – use getorf -find 3 (ATG-initiated ORFs that end at a stop codon)
             with a configurable -minsize to suppress spuriously short ORFs.

If the best similarity meets the identity threshold the original nucleotide
sequence and its metadata row are written to the output files.

Expected Snakemake rule (example):
--------------------------------------
rule filter_by_orf_identity:
    input:
        fasta    = "data/sequences.fasta",
        metadata = "data/metadata.tsv",
        ref      = "data/reference.fasta",
    params:
        identity_threshold = 0.90,  # 0-1 float, e.g. 90 %
        # Minimum ORF size in *nucleotides* passed to getorf -minsize.
        # Set to ~50-80 % of your expected CDS length to exclude tiny spurious
        # ORFs while still allowing partial sequences.  Default: 300 nt.
        min_orf_nt = 300,
    output:
        fasta    = "results/filtered.fasta",
        metadata = "results/filtered_metadata.tsv",
    script:
        "scripts/filter_by_orf_identity.py"
--------------------------------------
"""

import csv
import os
import re
import subprocess
import sys
import tempfile
from pathlib import Path

# ---------------------------------------------------------------------------
# Snakemake bindings (populated automatically when run as a snakemake script)
# ---------------------------------------------------------------------------
fasta_in:     str   = snakemake.input.fasta
metadata_in:  str   = snakemake.input.metadata
reference_in: str   = snakemake.input.ref
id_threshold: float = float(snakemake.params.identity_threshold)
min_orf_nt:   int   = int(getattr(snakemake.params, "min_orf_nt", 300))
fasta_out:    str   = snakemake.output.fasta
metadata_out: str   = snakemake.output.metadata

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def parse_fasta(path: str) -> dict[str, str]:
    """Return {header: sequence} for every record in *path*."""
    records: dict[str, str] = {}
    header = None
    parts: list[str] = []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith(">"):
                if header is not None:
                    records[header] = "".join(parts)
                header = line[1:]
                parts = []
            else:
                parts.append(line)
    if header is not None:
        records[header] = "".join(parts)
    return records


def write_fasta(path: str, records: dict[str, str]) -> None:
    """Write *records* ({header: sequence}) to *path* in FASTA format."""
    with open(path, "w") as fh:
        for header, seq in records.items():
            fh.write(f">{header}\n")
            for i in range(0, len(seq), 60):
                fh.write(seq[i:i + 60] + "\n")


def write_single_fasta(path: str, header: str, seq: str) -> None:
    """Write a single FASTA record to *path*."""
    write_fasta(path, {header: seq})


def run_getorf_all(seq: str, seq_id: str, tmpdir: str, min_nt: int) -> dict[str, str]:
    """
    Run EMBOSS getorf on *seq* (nucleotide) with:
        -find 3    : ATG-initiated ORFs that end at a stop codon only
                     (Option D — avoids spurious between-stop translations)
        -minsize   : minimum ORF size in nucleotides (Option D)

    Both strands are searched (getorf default).

    Returns a dict {orf_header: aa_sequence} for ALL reported ORFs.
    Returns an empty dict if getorf fails or finds nothing.
    """
    in_fa  = os.path.join(tmpdir, f"{seq_id}_in.fa")
    out_fa = os.path.join(tmpdir, f"{seq_id}_orfs.fa")

    write_single_fasta(in_fa, seq_id, seq)

    cmd = [
        "getorf",
        "-sequence", in_fa,
        "-outseq",   out_fa,
        "-minsize",  str(min_nt),
        "-find",     "3",        # ATG … stop only  (Option D)
        "-circular", "false",
        "-auto",
    ]
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        print(
            f"[WARNING] getorf failed for {seq_id}:\n{result.stderr.strip()}",
            file=sys.stderr,
        )
        return {}

    if not os.path.exists(out_fa):
        return {}

    return parse_fasta(out_fa)


def run_needle(query_seq: str, ref_seq: str, tmpdir: str, tag: str) -> str | None:
    """
    Run EMBOSS needle (global pairwise alignment) between *query_seq* and
    *ref_seq* (both amino-acid).  Returns the raw alignment text or None.
    """
    q_fa  = os.path.join(tmpdir, f"{tag}_query.fa")
    r_fa  = os.path.join(tmpdir, f"{tag}_ref.fa")
    aln_f = os.path.join(tmpdir, f"{tag}_aln.txt")

    write_single_fasta(q_fa, "query", query_seq)
    write_single_fasta(r_fa, "ref",   ref_seq)

    cmd = [
        "needle",
        "-asequence", q_fa,
        "-bsequence",  r_fa,
        "-outfile",    aln_f,
        "-gapopen",    "10",
        "-gapextend",  "0.5",
        "-aformat3",   "pair",
        "-auto",
    ]
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        print(
            f"[WARNING] needle failed for {tag}:\n{result.stderr.strip()}",
            file=sys.stderr,
        )
        return None

    with open(aln_f) as fh:
        return fh.read()


def parse_percent_similarity(alignment_text: str, ref_length: int) -> float:
    """
    Parse an EMBOSS needle 'pair' format alignment and compute:

        percent_similarity = matching_positions / ref_length

    Only columns where both residues are non-gap and identical are counted.
    """
    query_blocks: list[str] = []
    ref_blocks:   list[str] = []

    for line in alignment_text.splitlines():
        m = re.match(r"^(query|ref)\s+\d+\s+([A-Za-z\-]+)\s+\d+", line)
        if m:
            if m.group(1) == "query":
                query_blocks.append(m.group(2))
            else:
                ref_blocks.append(m.group(2))

    aligned_query = "".join(query_blocks)
    aligned_ref   = "".join(ref_blocks)

    if len(aligned_query) != len(aligned_ref):
        raise ValueError(
            f"Aligned sequences have different lengths "
            f"({len(aligned_query)} vs {len(aligned_ref)})"
        )

    matches = sum(
        1
        for q, r in zip(aligned_query, aligned_ref)
        if q != "-" and r != "-" and q.upper() == r.upper()
    )

    if ref_length == 0:
        return 0.0
    return matches / ref_length


def best_orf_similarity(
    orfs: dict[str, str],
    ref_orf: str,
    ref_orf_len: int,
    seq_id: str,
    tmpdir: str,
) -> tuple[float, str | None]:
    """
    Option B: align every candidate ORF in *orfs* to *ref_orf* and return
    (best_similarity, best_orf_header).  Returns (0.0, None) if *orfs* is empty.
    """
    best_pct    = 0.0
    best_header = None

    for orf_idx, (orf_header, orf_aa) in enumerate(orfs.items()):
        # Use a short stable tag to avoid filesystem issues with long headers
        tag = f"{seq_id}_orf{orf_idx}"
        aln_text = run_needle(orf_aa, ref_orf, tmpdir, tag)
        if aln_text is None:
            continue

        try:
            pct = parse_percent_similarity(aln_text, ref_orf_len)
        except ValueError as exc:
            print(f"  [WARNING] alignment parse error for {tag}: {exc}", file=sys.stderr)
            continue

        if pct > best_pct:
            best_pct    = pct
            best_header = orf_header

    return best_pct, best_header


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> None:
    # --- Load inputs --------------------------------------------------------
    print(f"[INFO] Reading sequences from : {fasta_in}")
    sequences = parse_fasta(fasta_in)
    print(f"[INFO] {len(sequences)} sequences loaded.")

    print(f"[INFO] Reading metadata from  : {metadata_in}")
    with open(metadata_in, newline="") as fh:
        reader     = csv.DictReader(fh, delimiter="\t")
        fieldnames = reader.fieldnames or []
        metadata_rows: list[dict] = list(reader)

    if not fieldnames:
        raise ValueError("Metadata TSV appears to be empty or has no header.")
    id_field   = fieldnames[0]
    meta_by_id = {row[id_field]: row for row in metadata_rows}

    print(f"[INFO] Reading reference from  : {reference_in}")
    ref_records = parse_fasta(reference_in)
    if not ref_records:
        raise ValueError("Reference FASTA is empty.")
    ref_header, ref_seq = next(iter(ref_records.items()))
    print(f"[INFO] Reference '{ref_header}', length {len(ref_seq)} bp")
    print(f"[INFO] Identity threshold      : {id_threshold}")
    print(f"[INFO] Min ORF size (nt)       : {min_orf_nt}")

    # --- Translate reference ORF (also uses -find 3 / -minsize) -------------
    with tempfile.TemporaryDirectory() as tmpdir:
        print("[INFO] Extracting reference ORF(s) …")
        ref_orfs = run_getorf_all(ref_seq, "reference", tmpdir, min_orf_nt)

        if not ref_orfs:
            raise RuntimeError(
                f"getorf (-find 3, -minsize {min_orf_nt}) found no ORF in the "
                "reference sequence.  Try lowering min_orf_nt or check that the "
                "reference contains a complete CDS."
            )

        # For the reference, take the longest ATG-initiated ORF
        ref_orf_header, ref_orf = max(ref_orfs.items(), key=lambda kv: len(kv[1]))
        ref_orf_len = len(ref_orf)
        print(
            f"[INFO] Reference ORF selected : '{ref_orf_header}' "
            f"({ref_orf_len} aa, from {len(ref_orfs)} candidate(s))"
        )

        # --- Process each query sequence ------------------------------------
        passing_seqs: dict[str, str] = {}
        passing_meta: list[dict]     = []
        total = len(sequences)

        for idx, (header, seq) in enumerate(sequences.items(), 1):
            seq_id = header.split()[0]
            print(f"[INFO] ({idx}/{total}) {seq_id}")

            # Option D: -find 3 (ATG…stop), -minsize min_orf_nt
            orfs = run_getorf_all(seq, seq_id, tmpdir, min_orf_nt)

            if not orfs:
                print(
                    f"  → SKIP  (no ATG-initiated ORF ≥ {min_orf_nt} nt found)"
                )
                continue

            print(f"  → {len(orfs)} candidate ORF(s) found; aligning all to reference …")

            # Option B: pick the ORF with the best similarity, not the longest
            best_pct, best_header = best_orf_similarity(
                orfs, ref_orf, ref_orf_len, seq_id, tmpdir
            )

            print(
                f"  → best similarity = {best_pct:.3f}  "
                f"(ORF: '{best_header}', threshold: {id_threshold})"
            )

            if best_pct >= id_threshold:
                passing_seqs[header] = seq
                if seq_id in meta_by_id:
                    passing_meta.append(meta_by_id[seq_id])
                else:
                    print(
                        f"  [WARNING] No metadata row found for '{seq_id}'",
                        file=sys.stderr,
                    )

    # --- Write outputs ------------------------------------------------------
    print(f"\n[INFO] {len(passing_seqs)}/{total} sequences passed the threshold.")

    Path(fasta_out).parent.mkdir(parents=True, exist_ok=True)
    Path(metadata_out).parent.mkdir(parents=True, exist_ok=True)

    write_fasta(fasta_out, passing_seqs)
    print(f"[INFO] Filtered FASTA written to    : {fasta_out}")

    with open(metadata_out, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(passing_meta)
    print(f"[INFO] Filtered metadata written to : {metadata_out}")


if __name__ == "__main__":
    main()
