# scripts/build_tx2gene_from_gff.py
# Purpose:
#   Parse a JGI GFF3 (Phaglo1 GeneCatalog) and build a transcript->gene mapping table.
#   We keep it robust and informative:
#     - Primary columns used by the workflow: transcript_id, gene_id
#     - Extra helpful columns (if present in the GFF attributes): transcript_num_id, gene_num_id, protein_id, transcript_name, gene_name
#
# How it maps to example GFF lines:
#   scaffold_1 ... gene ...   ID=gene_1;Name=jgi.p|Phaglo1|1;proteinId=1;transcriptId=1
#   scaffold_1 ... mRNA ...   ID=mRNA_1;Name=jgi.p|Phaglo1|1;Parent=gene_1;proteinId=1;transcriptId=1
#   -> transcript_id = "mRNA_1"
#      gene_id       = "gene_1"
#      transcript_num_id = "1"      (if present)
#      gene_num_id       = "1"      (rare; if present under 'gene_id' or similar)
#      protein_id        = "1"
#      transcript_name   = "jgi.p|Phaglo1|1"
#      gene_name         = "jgi.p|Phaglo1|1"
#
# Notes:
#   - We treat features with type 'mRNA' or 'transcript' as transcripts.
#   - We prefer attributes['ID'] for the transcript_id and attributes['Parent'] for the gene_id,
#     because those are the canonical GFF links.
#   - We also capture numeric IDs (e.g., 'transcriptId=1') because JGI FASTA headers sometimes
#     reference those; having them in the table can help if we need to reconcile later.
#
# Inputs (from Snakemake):
#   snakemake.input[0]  -> path to *.gff3(.gz)
# Outputs:
#   snakemake.output[0] -> path to tx2gene.tsv
#
# The output is tab-delimited with at least:
#   transcript_id  gene_id
# and may include extra columns if present.
# The script is meant to be run as part of the snakemake workflow.

import os
import gzip

import pandas as pd

# ---- Get IO from Snakemake (Snakemake injects `snakemake` in the script namespace) ----
gff_path = snakemake.input[0]
out_path = snakemake.output[0]

# ---- Helper functions ----
def open_maybe_gzip(path, mode="rt"):
    """Open plain text or gzipped file transparently."""
    if str(path).endswith(".gz"):
        return gzip.open(path, mode)
    return open(path, mode, encoding="utf-8")

def parse_gff_attributes(attr_field):
    """
    Parse the 9th GFF column (attributes) into a dict.
    Attributes look like: key1=value1;key2=value2;...
    """
    out = {}
    for kv in attr_field.strip().split(";"):
        if not kv:
            continue
        if "=" in kv:
            k, v = kv.split("=", 1)
            out[k] = v
    return out

# ---- Parse the GFF and extract transcript -> gene pairs ----
records = []
with open_maybe_gzip(gff_path, "rt") as fh:
    for line in fh:
        if not line or line.startswith("#"):
            continue
        parts = line.rstrip("\n").split("\t")
        if len(parts) < 9:
            continue

        seqid, source, feature, start, end, score, strand, phase, attrs_raw = parts
        feature_lower = feature.lower()
        if feature_lower not in ("mrna", "transcript"):
            # only consider transcript-like features
            continue

        attrs = parse_gff_attributes(attrs_raw)

        # Primary IDs for mapping (GFF3 best practice):
        #   - transcript: attributes['ID']
        #   - gene:       attributes['Parent']
        t_id = attrs.get("ID")
        g_id = attrs.get("Parent") or attrs.get("gene_id")  # some GFFs may use 'gene_id'

        # If either is missing, skip this row (cannot build a valid mapping)
        if not t_id or not g_id:
            continue

        # Nice-to-have extras (if present in this GFF):
        t_num = attrs.get("transcriptId")      # numeric transcript id in many JGI GFFs
        g_num = attrs.get("gene_id")           # rarely present; keep if available
        prot  = attrs.get("proteinId")
        t_nm  = attrs.get("Name")
        # There isn't a separate gene-name attribute here, but keep the transcript Name as a proxy
        g_nm  = None

        records.append(
            {
                "transcript_id": t_id,
                "gene_id": g_id,
                "transcript_num_id": t_num,
                "gene_num_id": g_num,
                "protein_id": prot,
                "transcript_name": t_nm,
                "gene_name": g_nm,
            }
        )

# ---- Build DataFrame, drop duplicates, and write ----
df = pd.DataFrame.from_records(records)

# Keep at least the two required columns
if df.empty:
    raise RuntimeError(
        f"No transcript->gene pairs were found in {gff_path}. "
        "Check that the GFF has 'mRNA' or 'transcript' features with ID and Parent attributes."
    )

# De-duplicate on the primary keys
df = df.drop_duplicates(subset=["transcript_id", "gene_id"]).reset_index(drop=True)

# Ensure output directory exists
os.makedirs(os.path.dirname(out_path), exist_ok=True)

# Write as TSV
df.to_csv(out_path, sep="\t", index=False)