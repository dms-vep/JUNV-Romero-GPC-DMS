# Description:
# Python script to download all fasta sequences
# specified by a list of accession numbers as well
# as metadata

# Author:
# Caleb Carr
# (Refactored for efficiency)

# Imports
import os
import re
import time
import datetime
import pandas as pd
from urllib.error import HTTPError
from Bio import Entrez, AlignIO, SeqIO
from Bio.Seq import Seq

Entrez.email = "JUNVexample@JUNVexample.com"     # Always tell NCBI who you are
ENTREZ_BATCH_SIZE = 500                           # Max records per NCBI request

# Strain names (case-insensitive substring match) that indicate a
# Candid#1-related sequence.  "Candid" is intentionally included so that
# "Candid#1" and plain "Candid" both match without needing to list them
# separately.
CANDID_RELATED_STRAINS = {
    "candid",
    "candid#1",
    "xj13",
    "xj44",
    "xj17",
    "xj34",
    "xj39",
    "xj48",
    "xj-44",
    "xjcl3",
}


# ---------------------------------------------------------------------------
# Entrez helpers
# ---------------------------------------------------------------------------

def entrez_fetch_with_retry(accession_ids, retries=3, delay=5):
    """
    Fetch a batch of GenBank records from NCBI with exponential back-off
    retry logic to handle transient HTTP errors.
    """
    for attempt in range(retries):
        try:
            handle = Entrez.efetch(
                db="nucleotide",
                id=",".join(accession_ids),
                rettype="gb",
                retmode="text",
            )
            return handle
        except HTTPError as e:
            if attempt < retries - 1:
                wait = delay * (2 ** attempt)   # exponential back-off
                time.sleep(wait)
            else:
                raise e


def fetch_records_in_batches(accession_ids):
    """
    Yield Biopython SeqRecord objects fetched from NCBI in batches of
    ENTREZ_BATCH_SIZE, which is far more efficient than one request per
    accession.
    """
    for i in range(0, len(accession_ids), ENTREZ_BATCH_SIZE):
        batch = accession_ids[i : i + ENTREZ_BATCH_SIZE]
        handle = entrez_fetch_with_retry(batch)
        records = list(SeqIO.parse(handle, "genbank"))
        handle.close()
        yield from records


# ---------------------------------------------------------------------------
# Shared parsing / filtering helpers
# ---------------------------------------------------------------------------

def percent_ids(seq1, seq2):
    """
    Calculate percent similarity between two aligned sequences,
    ignoring gap positions in seq1.
    """
    assert len(seq1) == len(seq2), "Aligned sequences do not have same length!"
    num_gaps   = seq1.count("-")
    num_bases  = len(seq1) - num_gaps
    matching   = sum(
        1 for i in range(len(seq1))
        if seq1[i] == seq2[i] and seq1[i] != "-"
    )
    return matching / num_bases


def reverse_complement(sequence):
    """Return the reverse complement of a DNA sequence using Biopython."""
    return str(Seq(sequence).reverse_complement())


def align_sequences_and_remove_low_matches(sequence):
    """
    Align the sequence, its reverse, and its reverse complement to the
    reference genome and return the orientation with the highest percent
    identity along with that identity score.
    """
    S_segment_reference = str(snakemake.params.reference_genome)
    temp_fasta_file     = str(snakemake.params.temp_fasta_file)
    temp_alignment_file = str(snakemake.params.temp_alignment_file)

    candidates = {
        "not_adjusted":      sequence,
        "reverse":           sequence[::-1],
        "reverse_complement": reverse_complement(sequence),
    }

    # Write all three orientations to a temp FASTA file
    with open(temp_fasta_file, "w") as fh:
        for name, seq in candidates.items():
            fh.write(f">{name}\n{seq}\n")

    # Align against reference
    os.system(
        f"mafft --6merpair --keeplength --quiet "
        f"--addfragments {temp_fasta_file} {S_segment_reference} "
        f"> {temp_alignment_file}"
    )

    alignment  = AlignIO.read(temp_alignment_file, "fasta")
    ref_seq    = alignment[0].seq

    best_name  = None
    best_seq   = None
    best_pid   = -1.0

    for record in alignment[1:]:   # skip the reference (index 0)
        pid = percent_ids(ref_seq, record.seq)
        if pid > best_pid:
            best_pid  = pid
            best_name = record.id
            best_seq  = candidates[record.id]

    return best_name, best_seq, best_pid


def format_date(raw_date, log_file):
    """
    Parse a GenBank collection_date string into YYYY-MM-DD format,
    using 'XX' for unknown month / day components.
    Returns '?' on failure.
    """
    parts = raw_date.split("-")

    try:
        if len(parts) == 1:
            year_str = raw_date.split("/")[0]
            return datetime.datetime.strptime(year_str, "%Y").strftime("%Y") + "-XX-XX"

        if len(parts) == 2:
            for fmt in ("%b-%Y", "%Y-%m"):
                try:
                    return datetime.datetime.strptime(raw_date, fmt).strftime("%Y-%m") + "-XX"
                except ValueError:
                    log_file.write(f"ERROR! {raw_date} not recognised as {fmt}\n")

        if len(parts) == 3:
            for fmt in ("%d-%b-%Y", "%Y-%m-%d"):
                try:
                    return datetime.datetime.strptime(raw_date, fmt).strftime("%Y-%m-%d")
                except ValueError:
                    log_file.write(f"ERROR! {raw_date} not recognised as {fmt}\n")

    except Exception:
        pass

    log_file.write(f"ERROR! Could not parse date: {raw_date}\n")
    return "?"


def country_extraction(location):
    """
    Normalise a GenBank geo_loc_name / country string into a
    (region, country) tuple, mapping the country to a continent region
    using the lists defined in the Snakemake config.
    """
    location = location.lower().split(":")[0].strip()

    special = {
        "viet nam":                       "Vietnam",
        "usa":                            "USA",
        "zaire":                          "Democratic Republic of the Congo",
        "democratic republic of the congo": "Democratic Republic of the Congo",
        "cote d'ivoire":                  "Cote d'Ivoire",
        "bosnia and herzegovina":         "Bosnia and Herzegovina",
        "ussr":                           "Russia",
    }

    extracted_country = special.get(
        location,
        " ".join(word.capitalize() for word in location.split()),
    )

    region_map = {
        "Asia":          snakemake.params.asia,
        "Oceania":       snakemake.params.oceania,
        "Africa":        snakemake.params.africa,
        "Europe":        snakemake.params.europe,
        "South America": snakemake.params.south_america,
        "North America": snakemake.params.north_america,
    }

    region = next(
        (r for r, countries in region_map.items() if extracted_country in countries),
        "?",
    )

    return region, extracted_country


def virus_grouper(virus):
    """Map an organism name to its configured virus group, or 'UNKNOWN'."""
    virus_dict = snakemake.params.virus_grouper
    for grouped, virus_list in virus_dict.items():
        if virus.upper() in virus_list:
            return grouped
    return "UNKNOWN"


def host_grouper(host):
    """Map a host description to its configured host group, or 'UNKNOWN'."""
    host_dict = snakemake.params.host_grouper
    for grouped, host_list in host_dict.items():
        if host.lower() in host_list:
            return grouped
    return "UNKNOWN"


METADATA_HEADER = "\t".join([
    "strain", "virus", "gene", "host", "accession", "date",
    "location", "region", "country", "database", "authors",
    "url", "title", "journal", "paper_url", "candid_related",
]) + "\n"


# ---------------------------------------------------------------------------
# Shared record parser (SeqRecord → metadata dict + sequence string)
# ---------------------------------------------------------------------------

def extract_record_fields(record, desired_gene, log_file):
    """
    Given a Biopython SeqRecord (parsed from GenBank format), extract all
    metadata fields and the nucleotide sequence.  Returns a dict.
    """
    accession = record.id
    organism  = record.annotations.get("organism", "?")
    date      = "?"

    # Sequence length from the record itself
    length = len(record.seq)

    # NCBI submission / update date as backup
    raw_backup = record.annotations.get("date", "")
    backup_date = format_date(raw_backup, log_file) if raw_backup else "XXXX-XX-XX"

    # GenBank URL
    genbank_base_url = "https://www.ncbi.nlm.nih.gov/nuccore/"
    url = genbank_base_url + accession

    # Database cross-reference (e.g. BioProject)
    db_xrefs = record.dbxrefs
    database = ", ".join(db_xrefs) if db_xrefs else "?"

    # Publication metadata (first reference only)
    authors   = "?"
    title     = "?"
    journal   = "?"
    paper_url = "?"
    if record.annotations.get("references"):
        ref = record.annotations["references"][0]
        raw_authors = ref.authors or ""
        authors     = raw_authors.split(",")[0].strip() + " et al" if raw_authors else "?"
        title       = ref.title   or "?"
        journal     = ref.journal or "?"
        if ref.pubmed_id:
            paper_url = f"https://pubmed.ncbi.nlm.nih.gov/{ref.pubmed_id}"

    # Source feature: strain, host, location, collection date, CDS
    strain   = "MISSING"
    host     = "?"
    location = "?"
    region   = "?"
    country  = "?"
    gene     = "?"
    CDS      = ""

    for feature in record.features:
        if feature.type == "source":
            isolate = (
                feature.qualifiers.get("isolate")
                or feature.qualifiers.get("strain")
            )
            if isolate:
                strain = isolate[0]
                for ch in ["@", " "]:
                    strain = strain.replace(ch, "-")
                for ch in list("()[]{}|#><;"):
                    strain = strain.replace(ch, "")

            raw_host = (
                feature.qualifiers.get("host")
                or feature.qualifiers.get("lab_host")
            )
            if raw_host:
                raw_host_str = raw_host[0].split(";")[0]
                grouped      = host_grouper(raw_host_str)
                if grouped != "UNKNOWN":
                    host = grouped
                else:
                    host = raw_host_str   # preserve raw value rather than losing it
                    log_file.write(f"{raw_host_str} host not known!\n")

            raw_loc = (
                feature.qualifiers.get("geo_loc_name")
                or feature.qualifiers.get("country")
            )
            if raw_loc:
                location             = raw_loc[0]
                region, country      = country_extraction(location)

            raw_date = feature.qualifiers.get("collection_date")
            if raw_date:
                date = format_date(raw_date[0], log_file)

        # CDS feature — look for desired gene if required
        if feature.type == "CDS" and gene != desired_gene:
            product = (
                feature.qualifiers.get("product")
                or feature.qualifiers.get("gene")
            )
            if product:
                product_name = product[0].upper()
                gene_names   = snakemake.params.gene_names

                if desired_gene is None or product_name in gene_names:
                    gene = desired_gene or product_name
                    log_file.write(f"{product_name} is {gene}!\n")

                    # Check for partial coordinates (<, >) and log a warning
                    loc_str = str(feature.location)
                    if "<" in loc_str or ">" in loc_str:
                        log_file.write(
                            f"WARNING: {product_name} CDS has partial coordinates "
                            f"({loc_str}); sequence may be incomplete at one or "
                            f"both ends.\n"
                        )

                    # Use Biopython's feature.extract() directly -- this correctly
                    # handles all coordinate styles (simple, join, complement,
                    # partial </>), is 0-based internally, and avoids the
                    # off-by-one errors that arise from re-parsing the location
                    # string manually in parse_CDS.
                    CDS = str(feature.extract(record.seq))

                    # Apply /codon_start offset if present.
                    # codon_start=1 means no offset (default), =2 skip 1 base,
                    # =3 skip 2 bases -- i.e. the offset is always codon_start - 1.
                    raw_codon_start = feature.qualifiers.get("codon_start", ["1"])
                    codon_start = int(raw_codon_start[0])
                    if codon_start not in (1, 2, 3):
                        log_file.write(
                            f"WARNING: {product_name} has unexpected codon_start "
                            f"value ({codon_start}); defaulting to 1 (no offset).\n"
                        )
                        codon_start = 1
                    offset = codon_start - 1
                    if offset:
                        log_file.write(
                            f"{product_name} codon_start={codon_start}; "
                            f"trimming {offset} base(s) from CDS start.\n"
                        )
                        CDS = CDS[offset:]
                else:
                    log_file.write(
                        f"ERROR: {product_name} product not known or not {desired_gene}\n"
                    )

    nucleotide_sequence = str(record.seq)

    return {
        "strain":              strain,
        "virus":               organism,
        "gene":                gene,
        "host":                host,
        "accession":           accession,
        "date":                date,
        "backup_date":         backup_date,
        "location":            location,
        "region":              region,
        "country":             country,
        "database":            database,
        "authors":             authors,
        "url":                 url,
        "title":               title,
        "journal":             journal,
        "paper_url":           paper_url,
        "CDS":                 CDS,
        "length":              length,
        "nucleotide_sequence": nucleotide_sequence,
    }


# ---------------------------------------------------------------------------
# Shared filtering + writing step
# ---------------------------------------------------------------------------

def filter_and_write_record(
    fields,
    accession_id,
    desired_gene,
    length_threshold,
    remove_duplicates,
    seen_sequences,
    accessions_to_include,
    output_fasta_file,
    output_metadata_file,
    log_file,
):
    """
    Apply all quality / relevance filters to an extracted record dict and,
    if it passes, write it to the FASTA and metadata output files.
    Returns True if the record was retained, False if it was removed.
    """

    virus  = fields["virus"]
    length = fields["length"]
    gene   = fields["gene"]
    CDS    = fields["CDS"]
    nucleotide_sequence = fields["nucleotide_sequence"]

    # Strip version suffix (e.g. "AF123456.1" -> "AF123456") for include-list
    # lookups, matching the original script's behaviour.
    accession_id_unversioned = accession_id.split(".")[0]
    force_include = accession_id_unversioned in accessions_to_include

    # --- virus filter ---
    grouped_virus = virus_grouper(virus)
    if grouped_virus == "UNKNOWN" and not force_include:
        log_file.write(
            f"{accession_id} excluded because not desired virus ({virus})!\n"
        )
        return False
    fields["virus"] = grouped_virus if grouped_virus != "UNKNOWN" else virus

    # --- length filter (bypassed for force-included accessions) ---
    lo, hi = length_threshold
    if not force_include and not (lo <= length <= hi):
        log_file.write(
            f"{accession_id} excluded because {length} not within length limits!\n"
        )
        return False

    # --- gene / CDS extraction ---
    if desired_gene is not None:
        if gene != desired_gene and not force_include:
            log_file.write(
                f"{accession_id} excluded because {gene} is not "
                f"the desired gene ({desired_gene})!\n"
            )
            return False

        if CDS:
            # CDS is already the extracted sequence string from feature.extract()
            # in extract_record_fields(); no further coordinate parsing needed.
            log_file.write(f"CDS found and extracted ({len(CDS)} bp)\n")
            nucleotide_sequence = CDS
        else:
            # Fall back to alignment-based orientation detection
            best_name, best_seq, pct = align_sequences_and_remove_low_matches(
                nucleotide_sequence
            )
            if best_name != "not_adjusted":
                fields["strain"] += f"_{best_name}"
            if pct >= snakemake.params.min_id_cutoff:
                nucleotide_sequence = best_seq
                log_file.write(
                    f"{accession_id} has no {desired_gene} CDS but "
                    f"{best_name} found to match by {pct}\n"
                )
            elif force_include:
                log_file.write(
                    f"{accession_id} has low percent match ({pct}) but "
                    f"force-included; using best orientation ({best_name}).\n"
                )
                nucleotide_sequence = best_seq
            else:
                log_file.write(
                    f"{accession_id} excluded because low percent match ({pct})!\n"
                )
                return False

    fields["nucleotide_sequence"] = nucleotide_sequence

    # --- duplicate filter ---
    if remove_duplicates == "Yes":
        if (
            nucleotide_sequence in seen_sequences
            and not force_include
        ):
            log_file.write(
                f"{accession_id} excluded because duplicate sequence!\n"
            )
            return False
        seen_sequences.append(nucleotide_sequence)

    # --- missing date fallback ---
    if fields["date"] == "?":
        fields["date"] = "XXXX-XX-XX"

    # --- ambiguous base filter ---
    if not nucleotide_sequence:
        log_file.write(
            f"{accession_id} excluded because sequence is empty "
            f"(CDS extraction may have failed)!\n"
        )
        return False
    n_frac = nucleotide_sequence.upper().count("N") / len(nucleotide_sequence)
    if n_frac > snakemake.params.max_frac_N:
        if force_include:
            log_file.write(
                f"{accession_id} has high N fraction ({n_frac:.2%}) but "
                f"force-included; retaining.\n"
            )
        else:
            log_file.write(
                f"{accession_id} excluded because of high N count!\n"
            )
            return False

    # --- build final strain name ---
    strain    = fields["strain"]
    accession = fields["accession"]
    date      = fields["date"]

    if strain != "MISSING":
        strain = f"{strain}_{accession}_{date}"
    else:
        strain = f"{accession}_{date}"

    for old, new in [("/", "_"), (". ", "-"), (" ", "_"), (".", "-")]:
        strain = strain.replace(old, new)

    assert strain != "MISSING", "Virus strain name is missing"
    fields["strain"] = strain

    # --- write outputs ---
    # Determine Candid#1-related status by checking whether any of the
    # known Candid-related strain identifiers appear in the strain name.
    # Comparison is case-insensitive so "CANDID", "Candid", etc. all match.
    strain_lower   = fields["strain"].lower()
    candid_related = "yes" if any(
        s in strain_lower for s in CANDID_RELATED_STRAINS
    ) else "no"

    meta_line = "\t".join([
        fields["strain"],
        fields["virus"],
        fields["gene"],
        fields["host"],
        fields["accession"],
        fields["date"],
        fields["location"],
        fields["region"],
        fields["country"],
        fields["database"],
        fields["authors"],
        fields["url"],
        fields["title"],
        fields["journal"],
        fields["paper_url"],
        candid_related,
    ]) + "\n"

    output_metadata_file.write(meta_line)
    output_fasta_file.write(f">{strain}\n{nucleotide_sequence}\n")

    return True


# ---------------------------------------------------------------------------
# Main processing functions
# ---------------------------------------------------------------------------

def read_and_process_accession_list(
    input_file_name,
    output_fasta_file_name,
    output_metadata_file_name,
    length_threshold,
    remove_duplicates,
    desired_gene,
    output_log_file,
):
    """
    Read a list of accession IDs, fetch their GenBank records from NCBI in
    efficient batches, parse each record, apply filters, and write passing
    records to FASTA + metadata outputs.
    """
    accessions_to_exclude = snakemake.params.accesstions_to_exclude
    accessions_to_include = snakemake.params.accessions_to_include

    with open(input_file_name, "r") as fh:
        all_accessions = [line.split()[0] for line in fh if line.strip()]

    with (
        open(output_fasta_file_name,    "a+") as out_fasta,
        open(output_metadata_file_name, "a+") as out_meta,
        open(output_log_file,           "a+") as log,
    ):
        out_meta.write(METADATA_HEADER)

        total    = 0
        retained = 0
        seen_sequences = []

        # Separate excluded accessions before batching
        to_fetch = []
        for acc in all_accessions:
            total += 1
            if acc[:-2] in accessions_to_exclude:
                log.write(f"{acc} excluded based on config file!\n")
            else:
                to_fetch.append(acc)

        # Fetch all remaining accessions in batches
        for record in fetch_records_in_batches(to_fetch):
            acc_id = record.id
            log.write(f"\nProcessing {acc_id}\n")

            fields = extract_record_fields(record, desired_gene, log)
            kept   = filter_and_write_record(
                fields,
                acc_id,
                desired_gene,
                length_threshold,
                remove_duplicates,
                seen_sequences,
                accessions_to_include,
                out_fasta,
                out_meta,
                log,
            )
            if kept:
                retained += 1

        removed = total - retained - (len(all_accessions) - len(to_fetch))
        log.write(
            f"\nA total of {total} were processed and {retained} were retained!\n"
        )


def read_and_process_local_files(
    directory,
    output_fasta_file_name,
    output_metadata_file_name,
    length_threshold,
    remove_duplicates,
    desired_gene,
    output_log_file,
):
    """
    Parse local GenBank files in a directory, apply the same filters as the
    NCBI fetcher, and append passing records to the shared FASTA + metadata
    outputs.
    """
    accessions_to_include = snakemake.params.accessions_to_include

    with (
        open(output_fasta_file_name,    "a+") as out_fasta,
        open(output_metadata_file_name, "a+") as out_meta,
        open(output_log_file,           "a+") as log,
    ):
        total    = 0
        retained = 0
        seen_sequences = []

        for entry in os.scandir(directory):
            if not entry.is_file():
                continue

            log.write(f"\nProcessing {entry.path}\n")
            total += 1

            with open(entry.path, "r") as fh:
                records = list(SeqIO.parse(fh, "genbank"))

            if not records:
                log.write(f"WARNING: No GenBank records found in {entry.path}\n")
                continue

            # Local files are expected to contain a single record each
            record = records[0]
            acc_id = str(entry).split("/")[-1].split(".")[0]

            fields = extract_record_fields(record, desired_gene, log)
            kept   = filter_and_write_record(
                fields,
                acc_id,
                desired_gene,
                length_threshold,
                remove_duplicates,
                seen_sequences,
                accessions_to_include,
                out_fasta,
                out_meta,
                log,
            )
            if kept:
                retained += 1

        log.write(
            f"\nA total of {total} were processed and {retained} were retained!\n"
        )


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

def main():
    """Main method"""

    # Inputs
    list_of_accessions = str(snakemake.input.accession_list)
    local_file_dir     = str(snakemake.params.local_file_dir)

    # Params
    length_threshold = (
        int(str(snakemake.params.genome_size_threshold_lower)),
        int(str(snakemake.params.genome_size_threshold_upper)),
    )
    desired_gene     = str(snakemake.params.desired_gene) or None
    remove_duplicates = str(snakemake.params.remove_duplicates)

    # Outputs
    fasta_output    = str(snakemake.output.fasta_sequences)
    metadata_output = str(snakemake.output.metadata)
    output_log_file = str(snakemake.log)

    # Fetch + process NCBI accessions
    read_and_process_accession_list(
        list_of_accessions,
        fasta_output,
        metadata_output,
        length_threshold,
        remove_duplicates,
        desired_gene,
        output_log_file,
    )

    # Optionally process local GenBank files
    if local_file_dir != "no files":
        read_and_process_local_files(
            local_file_dir,
            fasta_output,
            metadata_output,
            length_threshold,
            remove_duplicates,
            desired_gene,
            output_log_file,
        )


if __name__ == "__main__":
    main()
