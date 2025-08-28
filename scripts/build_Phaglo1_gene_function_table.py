# scripts/build_Phaglo1_gene_function_table.py
# Purpose
#   Join JGI functional annotations (KEGG/GO/InterPro/KOG/SignalP + our added dbCAN annotations)
#   at the protein level, map to gene_id via protein->transcript->gene, and aggregate per gene.
#
# Inputs (from Snakemake):
#   snakemake.input['tx2gene']  -> TSV with columns: transcript_id, gene_id
#   snakemake.input['prot_tx']  -> protein<tab>transcript map (JGI), headered or not
#   snakemake.input['kegg']     -> JGI KEGG tab.gz
#   snakemake.input['go']       -> JGI GO tab.gz
#   snakemake.input['ipr']      -> JGI InterPro tab.gz
#   snakemake.input['kog']      -> JGI KOG tab.gz
#   snakemake.input['sigp']     -> JGI SignalP tab.gz
#   snakemake.input['dbcan']    -> dbCAN overview.txt (or empty placeholder if it did not run)
#
# Output:
#   snakemake.output['func']    -> gene_functions.tsv (per gene aggregates)
#
# Notes:
#   - We normalize all tables and keep only columns we need.
#   - We read with comment='#' to ignore commented headers.
#   - We’re robust to small format changes by probing likely column names.

import os
import re
import pandas as pd

def read_table(path, **kwargs):
    kw = dict(sep="\t", comment="#", dtype=str)
    kw.update(kwargs)
    return pd.read_csv(path, **kw)

def ensure_cols(df, lower=True):
    """Lower-case columns for robust access."""
    df = df.copy()
    if lower:
        df.columns = [c.strip().lower() for c in df.columns]
    return df

def pick_col(df, candidates, required=False, label=""):
    """Pick the first existing column among candidates."""
    for c in candidates:
        if c in df.columns:
            return c
    if required:
        raise ValueError(f"Required column not found for {label}. Tried: {candidates}")
    return None

# Define Regex pattern for JGI header
jgi_id_pat = re.compile(r"jgi\|[^|]*\|(\d+)\|")

def normalize_protein_id_to_numeric(x: str) -> str:
    """Convert 'jgi|Species|12345|...' to '12345'; else return original."""
    if not isinstance(x, str):
        return x
    m = jgi_id_pat.search(x)
    return m.group(1) if m else x

def semijoin(series):
    return ";".join(sorted(set([s for s in series.dropna().astype(str) if s and s != r"\N"])))

# ------------------------ Load core transcript to gene maps ------------------------

tx2gene = read_table(snakemake.input['tx2gene'])
tx2gene = ensure_cols(tx2gene)

# required columns
if not {'transcript_id', 'gene_id'} <= set(tx2gene.columns):
    raise ValueError("tx2gene must have columns: transcript_id, gene_id")

# keep helper columns if present
has_tnum = 'transcript_num_id' in tx2gene.columns
if has_tnum:
    # normalize numeric-like strings
    tx2gene['transcript_num_id'] = (
        tx2gene['transcript_num_id']
        .astype(str)
        .str.replace(r'\.0$', '', regex=True)
        .str.strip()
    )

# Normalize primary keys to strings
tx2gene['transcript_id'] = tx2gene['transcript_id'].astype(str).str.strip()
tx2gene['gene_id']       = tx2gene['gene_id'].astype(str).str.strip()

# --- protein <-> transcript map (JGI) ---
prot_tx = read_table(snakemake.input['prot_tx'], header=0)
prot_tx = ensure_cols(prot_tx)

# Detect columns
pcol = pick_col(prot_tx, ["proteinid", "protein_id", "protein"], required=True,  label="protein-transcript map")
tcol = pick_col(prot_tx, ["transcriptid", "transcript_id", "transcript"], required=True, label="protein-transcript map")

prot_tx = prot_tx.rename(columns={pcol: "protein_id", tcol: "transcript_id"})
prot_tx = prot_tx[["protein_id", "transcript_id"]].dropna()

# Normalize types (treat everything as string)
prot_tx['protein_id']    = prot_tx['protein_id'].astype(str).str.strip()
prot_tx['transcript_id'] = prot_tx['transcript_id'].astype(str).str.replace(r'\.0$', '', regex=True).str.strip()

# ---- try to map to gene via transcript_id (string) first ----
map1 = prot_tx.merge(tx2gene[["transcript_id", "gene_id"]], on="transcript_id", how="left")
n1 = map1['gene_id'].notna().sum()

# ---- if poor overlap, try numeric transcript id column ----
use_map = map1
if (n1 == 0 or (n1 / len(prot_tx) < 0.05)) and has_tnum:
    tx2mini = tx2gene[['transcript_num_id', 'gene_id']].dropna().copy()
    tx2mini['transcript_num_id'] = tx2mini['transcript_num_id'].astype(str).str.strip()
    map2 = prot_tx.merge(tx2mini, left_on="transcript_id", right_on="transcript_num_id", how="left")
    n2 = map2['gene_id'].notna().sum()
    if n2 > n1:
        use_map = map2.drop(columns=['transcript_num_id'])
    # If still terrible, we’ll proceed but aggregation will just drop NAs.

# Final clean map
prot_tx = use_map.drop_duplicates(subset=["protein_id"]).reset_index(drop=True)

# ------------------------ Load functional annotation tables ------------------------
# KEGG
def load_kegg(path):
    # JGI header has '#proteinId ...'
    df = pd.read_csv(path, sep="\t", comment=None, header=0, dtype=str)
    df.columns = [c.lstrip("#").strip() for c in df.columns]
    df = ensure_cols(df)

    pcol = pick_col(df, ["proteinid", "protein_id", "protein"], required=True, label="KEGG")
    df = df.rename(columns={pcol: "protein_id"})
    df["protein_id"] = df["protein_id"].astype(str).str.strip()

    # definition / pathway columns if present
    defcol   = pick_col(df, ["definition", "def", "desc"], required=False, label="KEGG")
    pathcol  = pick_col(df, ["pathway", "pathway_name"], required=False, label="KEGG")
    eccol    = pick_col(df, ["ecnum", "ec_num", "ec"], required=False, label="KEGG")

    # EC: if no explicit column, regex from row text (support both EC:1.2.3.4 and 1.2.3.4)
    if eccol:
        ec_join = df[eccol].fillna("").astype(str).str.replace(r"\s+", "", regex=True)
    else:
        text = df.astype(str).agg(" ".join, axis=1)
        ec_from_text = text.str.findall(r"\b(?:EC:)?\d+\.\d+\.\d+\.\d+\b") \
                           .apply(lambda xs: sorted(set([x.replace("EC:", "") for x in xs])) if isinstance(xs, list) else [])
        ec_join = ec_from_text.apply(lambda xs: ";".join(xs))

    out = pd.DataFrame({
        "protein_id": df["protein_id"],
        "EC": ec_join,
        "KEGG_def": df[defcol].astype(str) if defcol else "",
        "KEGG_pathway": df[pathcol].astype(str) if pathcol else "",
    })
    return out.drop_duplicates()

kegg = load_kegg(snakemake.input['kegg'])

# GO
def load_go(path):
    df = pd.read_csv(path, sep="\t", comment=None, header=0, dtype=str)
    df.columns = [c.lstrip("#").strip() for c in df.columns]
    df = ensure_cols(df)

    pcol   = pick_col(df, ["proteinid", "protein_id", "protein"])
    tcol   = pick_col(df, ["transcriptid", "transcript_id", "transcript"])
    gocol  = pick_col(df, ["goacc", "go_acc", "goid", "go_id", "go"], required=False, label="GO")
    gname  = pick_col(df, ["goname", "go_name", "name"], required=False, label="GO")
    gtype  = pick_col(df, ["gotermtype", "go_type", "type"], required=False, label="GO")

    if gocol is None:
        text = df.astype(str).agg(" ".join, axis=1)
        df["go_extracted"] = text.str.findall(r"\bGO:\d{7}\b").apply(lambda lst: sorted(set(lst)) if isinstance(lst, list) else [])
        gocol = "go_extracted"

    if df[gocol].apply(lambda x: isinstance(x, list)).any():
        df = df.explode(gocol, ignore_index=True)
    df[gocol] = df[gocol].fillna("").astype(str)

    if pcol:
        df = df.rename(columns={pcol: "protein_id"})
    elif tcol:
        df = df.rename(columns={tcol: "transcript_id"})
        df = df.merge(prot_tx[["protein_id", "transcript_id"]].drop_duplicates(), on="transcript_id", how="left")
    else:
        raise ValueError("GO table must have proteinId or transcriptId")

    # Keep accession + name + type
    out = pd.DataFrame({
        "protein_id": df["protein_id"].astype(str),
        "GO": df[gocol].astype(str),
        "GO_name": df[gname].astype(str) if gname else "",
        "GO_type": df[gtype].astype(str) if gtype else "",
    })
    return out.dropna(subset=["protein_id"]).drop_duplicates()

go = load_go(snakemake.input['go'])

# InterPro
def load_ipr(path):
    # header is '#proteinId ...'
    df = pd.read_csv(path, sep="\t", comment=None, header=0, dtype=str)
    df.columns = [c.lstrip("#").strip() for c in df.columns]
    df = ensure_cols(df)
    pcol = pick_col(df, ["proteinid", "protein_id", "protein"], required=True, label="InterPro")
    df = df.rename(columns={pcol: "protein_id"})
    text = df.astype(str).agg(" ".join, axis=1)
    ipr_ids = text.str.findall(r"IPR\d{6}").apply(lambda lst: sorted(set(lst)) if isinstance(lst, list) else [])
    out = pd.DataFrame({
        "protein_id": df["protein_id"].astype(str),
        "IPR": ipr_ids.apply(lambda xs: ";".join(xs)),
    })
    return out.drop_duplicates()

ipr = load_ipr(snakemake.input['ipr'])

# KOG: transcriptId  proteinId  kogid  kogdefline  kogClass  kogGroup
def load_kog(path):
    df = pd.read_csv(path, sep="\t", comment=None, header=0, dtype=str)
    df.columns = [c.lstrip("#").strip() for c in df.columns]
    df = ensure_cols(df)
    pcol = pick_col(df, ["proteinid", "protein_id", "protein"], required=False, label="KOG")
    tcol = pick_col(df, ["transcriptid", "transcript_id", "transcript"], required=False, label="KOG")
    if pcol:
        df = df.rename(columns={pcol: "protein_id"})
    elif tcol:
        df = df.rename(columns={tcol: "transcript_id"})
        df = df.merge(prot_tx[["protein_id", "transcript_id"]].drop_duplicates(), on="transcript_id", how="left")
    else:
        raise ValueError("KOG table must have proteinId or transcriptId")
    kogcol   = pick_col(df, ["kogid", "kog", "kog_id"], required=True, label="KOG")
    kogdef   = pick_col(df, ["kogdefline", "kog_defline", "definition", "desc"], required=False, label="KOG")
    out = pd.DataFrame({
        "protein_id": df["protein_id"].astype(str),
        "KOG": df[kogcol].astype(str),
        "KOG_desc": df[kogdef].astype(str) if kogdef else "",
    })
    return out.drop_duplicates()

kog = load_kog(snakemake.input['kog'])

# SignalP
def load_sigp(path):
    # file has a free-text comment line, then a '#proteinid ...' header line
    # easiest: skip all '#' lines, then assign expected column names
    names = ["proteinid", "nn_cutpos", "neuro_net_vote", "hmm_cutpos", "hmm_signalpep_probability"]
    df = pd.read_csv(path, sep="\t", comment="#", header=None, names=names, dtype=str)
    df = ensure_cols(df)
    df = df.rename(columns={"proteinid": "protein_id"})
    def has_sp(row):
        # heuristic: either prob>0.5 or votes>0 or any ‘yes’ strings
        try:
            prob = float(row.get("hmm_signalpep_probability", 0) or 0)
        except Exception:
            prob = 0.0
        votes = str(row.get("neuro_net_vote", "")).strip()
        txt = " ".join(map(str, row.values)).lower()
        return (prob >= 0.5) or (votes not in ["", "0", "0.0"]) or ("yes" in txt) or ("signalp" in txt and "no" not in txt)
    df["signalp"] = df.apply(has_sp, axis=1)
    return df[["protein_id", "signalp"]].drop_duplicates()

sigp = load_sigp(snakemake.input['sigp'])

# dbCAN
def load_dbcan(dbcan_stub_path):
    """
    Robustly parse dbCAN outputs into a tidy table with columns:
      protein_id, CAZy

    Works whether snakemake passes .../dbcan/overview.txt (file)
    or the directory itself (.../dbcan/). We ignore overview.* and
    read the per-tool raw outputs if present:
      - dbCAN_hmm_results.tsv
      - dbCANsub_hmm_results.tsv
      - diamond.out

    Strategy:
      - Extract CAZy family names (GH/GT/PL/CE/CBM/AA + number)
      - Union across methods by default (see knobs below to change)
    """
    import os
    import re
    import pandas as pd

    if dbcan_stub_path is None:
        return pd.DataFrame(columns=["protein_id", "CAZy"])

    # Accept either a file (overview.*) or the directory
    dbdir = dbcan_stub_path
    if os.path.isfile(dbcan_stub_path):
        dbdir = os.path.dirname(dbcan_stub_path)

    # Known outputs
    f_hmm     = os.path.join(dbdir, "dbCAN_hmm_results.tsv")
    f_dbsub   = os.path.join(dbdir, "dbCANsub_hmm_results.tsv")
    f_diamond = os.path.join(dbdir, "diamond.out")

    frames = []

    # Helper: tidy a CAZy string to e.g. "GH16", "GT77", "CBM10"
    # (take the first family-like token if multiple; ignore subfamily suffix)
    fam_pat = re.compile(r"\b(?:GH|GT|PL|CE|CBM|AA)\d{1,3}\b", flags=re.I)

    def normalize_fams(text):
        if not isinstance(text, str):
            text = " ".join(map(str, text)) if not pd.isna(text) else ""
        hits = sorted(set(m.upper() for m in fam_pat.findall(text)))
        return hits

    # --- 1) dbCAN HMM results ---
    if os.path.exists(f_hmm) and os.path.getsize(f_hmm) > 0:
        df = pd.read_csv(f_hmm, sep="\t", comment="#")
        # Typical columns: "HMM Name", "Target Name", "Coverage", etc.
        cols = {c.lower(): c for c in df.columns}
        hmm_name  = cols.get("hmm name") or cols.get("hmm_name")
        target    = cols.get("target name") or cols.get("target_name")
        if hmm_name and target:
            fams = df[hmm_name].astype(str).str.replace(".hmm$", "", regex=True).apply(normalize_fams)
            # explode to long form
            tmp = pd.DataFrame({"protein_id": df[target].astype(str), "CAZy": fams})
            tmp["protein_id"] = tmp["protein_id"].astype(str).str.strip().apply(normalize_protein_id_to_numeric)
            tmp = tmp.explode("CAZy").dropna()
            tmp["CAZy"] = tmp["CAZy"].astype(str)
            tmp["method"] = "HMM"
            frames.append(tmp[["protein_id", "CAZy", "method"]])

    # --- 2) dbCAN-sub HMM results ---
    if os.path.exists(f_dbsub) and os.path.getsize(f_dbsub) > 0:
        df = pd.read_csv(f_dbsub, sep="\t", comment="#")
        # Same column names as HMM; HMM Name often like "GH5_sub1"
        cols = {c.lower(): c for c in df.columns}
        hmm_name  = cols.get("hmm name") or cols.get("hmm_name")
        target    = cols.get("target name") or cols.get("target_name")
        if hmm_name and target:
            # Extract the *main* family from a subfamily label
            # e.g., "GH5_sub1" -> "GH5"
            main_fam = (
                df[hmm_name]
                .astype(str)
                .str.replace(".hmm$", "", regex=True)
                .str.replace(r"[_\-].*$", "", regex=True)
            )
            fams = main_fam.apply(normalize_fams)
            tmp = pd.DataFrame({"protein_id": df[target].astype(str), "CAZy": fams})
            tmp["protein_id"] = tmp["protein_id"].astype(str).str.strip().apply(normalize_protein_id_to_numeric)
            tmp = tmp.explode("CAZy").dropna()
            tmp["CAZy"] = tmp["CAZy"].astype(str)
            tmp["method"] = "dbCANsub"
            frames.append(tmp[["protein_id", "CAZy", "method"]])

    # --- 3) DIAMOND results ---
    if os.path.exists(f_diamond) and os.path.getsize(f_diamond) > 0:
        # dbCAN typically writes BLAST-like 12-column output without header.
        # We read as strings and scan each row for family tokens (robust).
        df = pd.read_csv(f_diamond, sep="\t", header=None, dtype=str)
        # join row text, extract families
        row_txt = df.astype(str).agg(" ".join, axis=1)
        fams = row_txt.apply(normalize_fams)
        # qseqid is usually column 0
        qids = df[0].astype(str) if 0 in df.columns else pd.Series([], dtype=str)
        tmp = pd.DataFrame({"protein_id": qids, "CAZy": fams})
        tmp["protein_id"] = tmp["protein_id"].astype(str).str.strip().apply(normalize_protein_id_to_numeric)
        tmp = tmp.explode("CAZy").dropna()
        tmp["CAZy"] = tmp["CAZy"].astype(str)
        tmp["method"] = "DIAMOND"
        frames.append(tmp[["protein_id", "CAZy", "method"]])

    if not frames:
        # nothing found — return empty
        return pd.DataFrame(columns=["protein_id", "CAZy"])

    allhits = pd.concat(frames, ignore_index=True).drop_duplicates()

    # ----------------------------
    # Tune the function behaviour
    # ----------------------------

    # (A) Default: UNION of all families per protein (robust/sensitive)
    keep = (
        allhits[["protein_id", "CAZy"]]
        .drop_duplicates()
        .reset_index(drop=True)
    )

    # (B) Consensus: uncomment to require ≥2 methods support the same family
    # support = (allhits.groupby(["protein_id", "CAZy"])["method"]
    #                   .nunique().reset_index(name="n_methods"))
    # keep = support.loc[support["n_methods"] >= 2, ["protein_id", "CAZy"]]

    # (C) Priority: keep HMM first, else dbCANsub, else DIAMOND
    # method_rank = {"HMM": 1, "dbCANsub": 2, "DIAMOND": 3}
    # allhits["rank"] = allhits["method"].map(method_rank).fillna(99)
    # keep = (allhits.sort_values(["protein_id", "CAZy", "rank"])
    #               .drop_duplicates(["protein_id", "CAZy"])
    #               .drop_duplicates(["protein_id"])   # one family per protein
    #               [["protein_id", "CAZy"]])

    # --- Make readable CAZy descriptions from fam-substrate-mapping.tsv ---
    cazy_map = None
    map_path = os.path.join(snakemake.config.get("dbcan_db_dir", ""), "fam-substrate-mapping.tsv")
    if os.path.exists(map_path):
        # Let pandas auto-detect tabs vs spaces; then normalize headers.
        m = pd.read_csv(map_path, sep=None, engine="python", dtype=str).fillna("")
        m.columns = [c.strip() for c in m.columns]          # trim whitespace
        # expected columns: Substrate_high_level, Substrate_curated, Family, Name, EC_Number

        # normalize family tag to match our CAZy (e.g., GH5, CBM10)
        if "Family" not in m.columns:
            print(f"[dbCAN] fam-substrate-mapping.tsv columns: {list(m.columns)}")
            # bail gracefully if it's not the file we expect
            cazy_map = None
        else:
            m["CAZy"] = (
                m["Family"].astype(str)
                .str.upper()
                .str.replace(r"\.hmm$", "", regex=True)
                .str.strip()
            )

            # e.g. "beta-agarase [agarose] (EC 3.2.1.81)"
            def mk_desc(row):
                name = (row.get("Name", "") or "").strip().strip('"')
                sub  = (row.get("Substrate_curated", "") or "").strip()
                ec   = (row.get("EC_Number", "") or "").strip()
                parts = []
                if name:
                    parts.append(name)
                if sub:
                    parts.append(f"[{sub}]")
                if ec:
                    parts.append(f"(EC {ec})")
                return " ".join(parts).strip()

            m["CAZy_desc"] = m.apply(mk_desc, axis=1)
            cazy_map = m[["CAZy", "CAZy_desc"]].drop_duplicates()

    if cazy_map is not None:
        keep = keep.merge(cazy_map, on="CAZy", how="left")

    # --- Logging summary ---
    try:
        n_prots = keep["protein_id"].nunique()
        n_fams  = keep["CAZy"].nunique()
        by_method = (
            allhits.groupby("method")["protein_id"]
            .nunique()
            .to_dict()
        )
        print(f"[dbCAN] Annotated {n_prots} proteins with {n_fams} unique CAZy families.")
        print(f"[dbCAN] Breakdown by method: {by_method}")
    except Exception as e:
        print(f"[dbCAN] Summary logging failed: {e}")

    return keep

dbcan = load_dbcan(snakemake.input.get('dbcan'))

# ------------------------ Merge at protein level ------------------------
prot_func = prot_tx[["protein_id", "gene_id"]].dropna().drop_duplicates()
prot_func["protein_id"] = prot_func["protein_id"].astype(str).str.strip()

# Make signalp a proper boolean to avoid weirdness later
if "signalp" in prot_func.columns:
    prot_func["signalp"] = prot_func["signalp"].fillna(False).astype(bool)
else:
    prot_func["signalp"] = False

for name, piece in [("KEGG", kegg), ("GO", go), ("IPR", ipr), ("KOG", kog), ("SIGP", sigp), ("dbCAN", dbcan)]:
    if piece is not None and not piece.empty:
        piece = piece.copy()
        piece["protein_id"] = piece["protein_id"].astype(str).str.strip()
        prot_func = prot_func.merge(piece, on="protein_id", how="left")

if "CAZy_desc" not in prot_func.columns:
    prot_func["CAZy_desc"] = ""

# ------------------------ Aggregate at gene level ------------------------
def set_join(series, split_pat=None):
    vals = set()
    for x in series.dropna().astype(str):
        x = x.strip()
        if not x or x == r"\N":
            continue
        if split_pat:
            vals.update([y for y in re.split(split_pat, x) if y])
        else:
            vals.update([y for y in x.split(";") if y])
    return ";".join(sorted(vals))

agg = prot_func.groupby("gene_id").agg({
    "protein_id": "count",
    "EC":      lambda s: set_join(s, split_pat=r"[;,\s]+"),
    "KEGG_def":      semijoin,
    "KEGG_pathway":  semijoin,
    "GO":      lambda s: set_join(s, split_pat=r"[;,\s]+"),
    "GO_name": semijoin,
    "GO_type": semijoin,
    "IPR":     lambda s: set_join(s, split_pat=r"[;,\s]+"),
    "KOG":     lambda s: set_join(s, split_pat=r"[;,\s]+"),
    "KOG_desc": semijoin,
    "CAZy":    lambda s: set_join(s, split_pat=r"[;,\s]+"),
    # Only present if mapping file found; harmless if missing
    "CAZy_desc": semijoin if "CAZy_desc" in prot_func.columns else lambda s: "",
    "signalp": "any",
}).reset_index().rename(columns={"protein_id": "n_proteins"})

# Ensure output dir exists and write
os.makedirs(os.path.dirname(snakemake.output['func']), exist_ok=True)
agg.to_csv(snakemake.output['func'], sep="\t", index=False)

# Log final summary
print("[SUMMARY] prot–tx rows:", len(prot_tx), "unique proteins:", prot_tx["protein_id"].nunique())
for label, df in [("KEGG", kegg), ("GO", go), ("IPR", ipr), ("KOG", kog), ("SIGP", sigp), ("dbCAN", dbcan)]:
    if df is not None and not df.empty:
        print(f"[SUMMARY] {label}: proteins={df['protein_id'].nunique()} rows={len(df)}")