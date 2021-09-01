# Remove unwanted genes, names or terms that are uninformative of would pollute the dataset.
# cmdoret, 20210825

import re
import numpy as np
import pandas as pd
import pathlib

DATA_DIR = pathlib.Path("../data/train")

# Change filtering parameters here
MIN_SEQLEN, MAX_SEQLEN = 10, 10000
MIN_NTERMS, MAX_NTERMS = 2, 100000
MIN_TERMFREQ, MAX_TERMFREQ = 10, 1000000
MIN_NAMELEN, MAX_NAMELEN = 5, 100
BLACKLIST = [
    r'.*hypothetical.*',
    r'.*[0-9]+-[0-9]+$',
    r'^[0-9A-Z]+$',
    'Protein X',
    'Negative factor',
    'Uncharacterized protein',
    'Predicted protein',
    'Putative secreted protein',
    'Negative factor',
    '3\'ORF'
]

def load_tbl(path):
    """Helper function to load tables with specific options"""
    df = pd.read_csv(path, sep="\t", compression="gzip")
    return df


features, names, sequences, terms = map(
    lambda fname: load_tbl(DATA_DIR / f"{fname}.tsv.gz"),
    ["features", "names", "sequences", "terms"],
)

### GENE-LEVEL FILTERS
# Genes whose a.a. seqs are too short/long
bad_seqlen = (
    sequences.seq.apply(len)
    < MIN_SEQLEN | sequences.seq.str.len()
    >= MAX_SEQLEN
)

# Genes whose names are all too short/long
namelen = names.name.str.len()
bad_namelen = (
    names.loc[MIN_NAMELEN < namelen <= MAX_NAMELEN, :]
    .groupby("id")
    .apply(lambda g: g[g["name"].str.len() == g["name"].str.len().max()])
    .reset_index(drop=True)
    .merge(features, by="id", how="right")["name"]
    .isnull()
)

# Genes with foo few/many terms
good_nterms = (
    terms.groupby("id").size().apply(lambda x: MIN_NTERMS < x <= MAX_NTERMS)
)
bad_nterms = ~features["id"].isin(good_nterms.index[good_nterms].values)

### TERM-LEVEL FILTERS
keep_terms = (
    terms.groupby("term")
    .size()
    .apply(lambda x: MIN_TERMFREQ < x <= MAX_TERMFREQ)
)

### NAME-LEVEL FILTERS
