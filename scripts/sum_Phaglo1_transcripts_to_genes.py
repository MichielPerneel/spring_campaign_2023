import os, re
import pandas as pd

tx2gene_path = snakemake.input[0]
quant_paths  = list(snakemake.input[1:])
tpm_out      = snakemake.output[0]
num_out      = snakemake.output[1]

tx2gene = pd.read_csv(tx2gene_path, sep="\t", dtype=str)
tx2gene.columns = [c.strip() for c in tx2gene.columns]
if "transcript_id" not in tx2gene.columns or "gene_id" not in tx2gene.columns:
    raise RuntimeError("tx2gene must contain 'transcript_id' and 'gene_id'")
has_tnum = "transcript_num_id" in tx2gene.columns

def sample_from_quant_path(p):
    return os.path.basename(os.path.dirname(p))

def jgi_pipe3(name: str):
    """For JGI headers like jgi|Phaglo1|10000|Phglo.0020s0046.1, return '10000'."""
    if not isinstance(name, str):
        return None
    tok = name.split()[0]  # just in case
    parts = tok.split("|")
    if len(parts) >= 3:
        return parts[2]
    return None

def name_variants_basic(n):
    """Other reasonable variants to try against transcript_id (string IDs)."""
    if not isinstance(n, str):
        return []
    v = []
    n0 = n
    v.append(n0)
    n1 = n0.split()[0]
    v.append(n1)
    # last pipe token
    if "|" in n1:
        v.append(n1.split("|")[-1])
    # strip trailing .version
    v.append(re.sub(r"\.\d+$", "", n1))
    # drop mRNA_/transcript_ prefix
    v.append(re.sub(r"^(mRNA_|transcript_)", "", n1))
    # dedupe
    seen, out = set(), []
    for x in v:
        if x and x not in seen:
            seen.add(x); out.append(x)
    return out

tpm_tables = []
num_tables = []

for qf in quant_paths:
    sample = sample_from_quant_path(qf)
    qs = pd.read_csv(qf, sep="\t", usecols=["Name", "TPM", "NumReads"])
    qs["Name"] = qs["Name"].astype(str)

    # Build candidate columns
    qs["cand_raw"]   = qs["Name"]
    qs["cand_token"] = qs["Name"].str.split().str[0]
    qs["cand_pipe_last"] = qs["cand_token"].apply(lambda x: x.split("|")[-1] if isinstance(x,str) else x)
    qs["cand_novers"] = qs["cand_token"].str.replace(r"\.\d+$", "", regex=True)
    qs["cand_nopref"] = qs["cand_pipe_last"].str.replace(r"^(mRNA_|transcript_)", "", regex=True)

    # JGI-specific: 3rd pipe field => numeric transcriptId; match to transcript_num_id
    qs["cand_pipe3num"] = qs["cand_token"].apply(jgi_pipe3)

    attempts = []
    best = None
    # Try to match to transcript_id first (string IDs like mRNA_123)
    for col in ["cand_raw","cand_token","cand_pipe_last","cand_novers","cand_nopref"]:
        m = qs.merge(tx2gene, left_on=col, right_on="transcript_id", how="inner")
        attempts.append((col, "transcript_id", len(m)))
        if best is None or len(m) > len(best[2]):
            best = (col, "transcript_id", m)

    # Try to match to transcript_num_id using pipe3 numeric (JGI)
    if has_tnum:
        m_num = qs.merge(tx2gene, left_on="cand_pipe3num", right_on="transcript_num_id", how="inner")
        attempts.append(("cand_pipe3num", "transcript_num_id", len(m_num)))
        if best is None or len(m_num) > len(best[2]):
            best = ("cand_pipe3num", "transcript_num_id", m_num)

    if best is None or best[2].empty:
        raise RuntimeError(
            f"No overlap between quant.sf ({qf}) and tx2gene ({tx2gene_path}). "
            f"Tried keys: {attempts}"
        )

    g = best[2]
    tpm_gene = g.groupby("gene_id", as_index=True)["TPM"].sum().to_frame(name=sample)
    num_gene = g.groupby("gene_id", as_index=True)["NumReads"].sum().to_frame(name=sample)
    tpm_tables.append(tpm_gene)
    num_tables.append(num_gene)

tpm_wide = pd.concat(tpm_tables, axis=1).fillna(0).sort_index()
num_wide = pd.concat(num_tables, axis=1).fillna(0).sort_index()

os.makedirs(os.path.dirname(tpm_out), exist_ok=True)
os.makedirs(os.path.dirname(num_out), exist_ok=True)
tpm_wide.to_csv(tpm_out, sep="\t", index=True)
num_wide.to_csv(num_out, sep="\t", index=True)