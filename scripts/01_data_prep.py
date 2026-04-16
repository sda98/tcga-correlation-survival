#!/usr/bin/env python3
"""
01_data_prep.py

Downloads and preprocesses TCGA pan-cancer expression and survival data
from UCSC Xena for downstream correlation and survival analyses.

Pipeline:
    1. Download gene expression (TCGA RSEM gene TPM, log2(TPM+0.001))
       and survival metadata if not already cached locally.
    2. Convert expression values from log2(TPM+0.001) to log2(TPM+1)
       for more conventional downstream handling.
    3. Trim TCGA sample ID suffixes (e.g., TCGA-A1-A0SB-01A-11R-A144-07
       → TCGA-A1-A0SB-01A-11R-A144).
    4. Map Ensembl gene IDs to HUGO gene symbols via the MyGene.info API,
       deduplicating by keeping the most recent Ensembl ID per symbol.
    5. Clip negative expression values to zero.
    6. Save cleaned expression matrix and survival table to the results
       directory for use by 02_correlation.py and 03_survival.py.

Outputs:
    - results/expression_clean.tsv
    - results/survival_clean.tsv
"""

import os
import urllib.request
import pandas as pd
import numpy as np
import mygene


# ============================================================
# Configuration
# ============================================================

DATA_DIR = "data"
RESULTS_DIR = "results"

# UCSC Xena download URLs
EXPRESSION_URL = (
    "https://toil-xena-hub.s3.us-east-1.amazonaws.com/download/tcga_RSEM_gene_tpm.gz"
)
SURVIVAL_URL = (
    "https://tcga-pancan-atlas-hub.s3.us-east-1.amazonaws.com/download/Survival_SupplementalTable_S1_20171025_xena_sp"
)

# Raw file paths
EXPRESSION_GZ = os.path.join(DATA_DIR, "tcga_RSEM_gene_tpm.gz")
SURVIVAL_RAW = os.path.join(DATA_DIR, "Survival_SupplementalTable_S1_20171025_xena_sp")

# Clean output paths
EXPRESSION_CLEAN = os.path.join(RESULTS_DIR, "expression_clean.tsv")
SURVIVAL_CLEAN = os.path.join(RESULTS_DIR, "survival_clean.tsv")


# ============================================================
# Step 0: Create directories
# ============================================================

os.makedirs(DATA_DIR, exist_ok=True)
os.makedirs(RESULTS_DIR, exist_ok=True)


# ============================================================
# Step 1: Download data
# ============================================================

def download_file(url, dest):
    """Download a file if it doesn't already exist."""
    if os.path.exists(dest):
        print(f"  Already exists: {dest}")
        return
    print(f"  Downloading: {url}")
    urllib.request.urlretrieve(url, dest)
    print(f"  Saved to: {dest}")


print("Step 1: Downloading data files...")
download_file(EXPRESSION_URL, EXPRESSION_GZ)
download_file(SURVIVAL_URL, SURVIVAL_RAW)


# ============================================================
# Step 2: Load data into DataFrames
# ============================================================

print("Step 2: Loading data...")

# Expression matrix: genes (rows) x samples (columns)
chunks = []
for chunk in pd.read_csv(EXPRESSION_GZ, sep="\t", index_col=0, chunksize=5000):
    chunks.append(chunk.astype("float32"))
tpm = pd.concat(chunks)
del chunks
print(f"  Expression matrix: {tpm.shape[0]} genes x {tpm.shape[1]} samples")

# Survival/clinical data: samples (rows) x clinical variables (columns)
survival = pd.read_csv(SURVIVAL_RAW, sep="\t")
# Rename 'sample' column to 'Sample_ID' 
if "sample" in survival.columns:
    survival = survival.rename(columns={"sample": "Sample_ID"})
print(f"  Survival data: {survival.shape[0]} samples x {survival.shape[1]} columns")


# ============================================================
# Step 3: Convert log2(TPM + 0.001) to log2(TPM + 1)
# The raw Xena data stores values as log2(TPM + 0.001).
# We reverse that transform, then reapply with +1 instead.
# ============================================================

print("Step 3: Converting log2(TPM+0.001) to log2(TPM+1)...")
tpm = 2 ** tpm        
tpm -= 0.001          
tpm += 1              
tpm = np.log2(tpm)    


# ============================================================
# Step 4: Trim sample ID suffixes
# TCGA barcodes look like: TCGA-A1-A0SB-01A-11R-A144-07
# This strips the last hyphen-separated segment: -07
# Result: TCGA-A1-A0SB-01A-11R-A144
# ============================================================

print("Step 4: Trimming sample ID suffixes...")
tpm.columns = tpm.columns.str.replace(r"-[^-]+$", "", regex=True)


# ============================================================
# Step 5: Map Ensembl IDs to HUGO gene symbols
# ============================================================

print("Step 5: Mapping Ensembl IDs to HUGO gene symbols...")

# Remove version suffix from Ensembl IDs: ENSG00000141510.11 -> ENSG00000141510
ensg_ids = tpm.index.str.replace(r"\.\d+$", "", regex=True).tolist()

# Query MyGene.info for HUGO symbols
mg = mygene.MyGeneInfo()
results = mg.querymany(
    ensg_ids,
    scopes="ensembl.gene",     # search by Ensembl gene ID
    fields="symbol",           # Return HUGO symbol
    species="human",
    verbose=False,
)

# Build a DataFrame of Ensembl ID -> HUGO symbol mappings
mapping_records = []
for r in results:
    if "symbol" in r and r["symbol"]:
        mapping_records.append({
            "ensembl_gene_id": r["query"],
            "hgnc_symbol": r["symbol"],
        })
mapping_df = pd.DataFrame(mapping_records)
print(f"  MyGene.info returned {len(mapping_df)} mappings")


# ============================================================
# Step 5b: Deduplicate — keep largest ENSG ID per HUGO symbol
# ============================================================

# Convert ENSG prefix to numeric for sorting
mapping_df["ensg_numeric"] = (
    mapping_df["ensembl_gene_id"].str.replace("ENSG", "", regex=False).astype(int)
)

# Sort by symbol (ascending), then by ENSG numeric (descending)
mapping_df = mapping_df.sort_values(
    ["hgnc_symbol", "ensg_numeric"], ascending=[True, False]
)

# Keep only the first (largest ENSG) per HUGO symbol
mapping_df = mapping_df.drop_duplicates(subset="hgnc_symbol", keep="first")

# Drop empty symbols (some genes have no HUGO name)
mapping_df = mapping_df[mapping_df["hgnc_symbol"] != ""]

# Build lookup dict: Ensembl ID -> HUGO symbol
ensembl_to_symbol = dict(
    zip(mapping_df["ensembl_gene_id"], mapping_df["hgnc_symbol"])
)
print(f"  Unique HUGO symbols after dedup: {len(ensembl_to_symbol)}")


# ============================================================
# Step 6: Filter expression matrix and rename rows to HUGO symbols
# ============================================================

print("Step 6: Filtering and renaming genes...")

# Replace versioned index with clean Ensembl IDs
tpm.index = ensg_ids

# Keep only rows that have a mapping
mask = tpm.index.isin(ensembl_to_symbol.keys())
tpm = tpm[mask].copy()

# Rename index from Ensembl IDs to HUGO symbols
tpm.index = tpm.index.map(ensembl_to_symbol)
tpm.index.name = "gene"
print(f"  Final expression matrix: {tpm.shape[0]} genes x {tpm.shape[1]} samples")


# ============================================================
# Step 7: Clip negative values to zero
# ============================================================

print("Step 7: Clipping negative values to 0...")
tpm = tpm.clip(lower=0)


# ============================================================
# Step 8: Save cleaned outputs
# ============================================================

print("Step 8: Saving cleaned data...")
tpm.to_csv(EXPRESSION_CLEAN, sep="\t")
survival.to_csv(SURVIVAL_CLEAN, sep="\t", index=False)
print(f"  Expression saved to: {EXPRESSION_CLEAN}")
print(f"  Survival saved to:   {SURVIVAL_CLEAN}")
print("Done.")
