"""Contains utility functions for the SRC package."""

import os
import logging
from dataclasses import dataclass, field
import pandas as pd
import pysam


@dataclass
class PRSConfig:
    """Common configuration class for PRS scoring."""

    bcf: str
    col: str = field(default="GT")
    build: str = "hg38"
    impute: bool = False
    refbcf: str = None
    parallel: bool = False
    ntasks: int = 1
    batch_size: int = 1

    def __post_init__(self):
        """Validate the 'col' parameter."""
        if self.col not in {"GT", "GP"}:
            raise ValueError(
                f"Invalid value for 'col': {self.col}.  Must be 'GT' or 'GP'."
            )


def configure_logging():
    """Configure logging with both file and stream handlers."""
    log_file = f"{__name__}.log"
    log_path = os.path.join(os.path.dirname(__file__), "..", log_file)

    handlers = [logging.StreamHandler()]

    try:
        file_handler = logging.FileHandler(log_path)
        handlers.insert(0, file_handler)
    except (OSError, IOError) as e:
        print(
            f"Warning: Could not write to log file '{log_path}'. Using console logging only. ({e})"
        )

    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
        handlers=handlers,
    )


def get_samples(var_obj):
    """Extract and return sample names from a bcf object."""
    return list(var_obj.header.samples)


def check_bed_type(bed):
    if isinstance(bed, pd.DataFrame):
        return bed
    if not os.path.isfile(bed):
        raise FileNotFoundError(f"'{bed}' not found.")

    df = pd.read_csv(bed, sep=r"\s+", engine="python")
    if df.shape[1] < 3:
        raise InvalidBedFormatError(f"Invalid BED format in '{bed}'.")
    return df


def determine_bcf_type(bcf):
    """Determine if input is a mapping file or a single BCF/VCF."""
    bcf = os.path.abspath(bcf)

    if os.path.splitext(bcf)[1].lower() == ".txt":
        mapping_dir = os.path.dirname(bcf)
        df = pd.read_csv(bcf, sep=r"\s+", header=None, engine="python")
        if df.shape[1] < 2:
            raise ValueError(f"Invalid mapping file: {bcf}")

        out = {}
        for contig, path in df.iloc[:, :2].values:
            path = os.path.expanduser(path)
            if not os.path.isabs(path):
                path = os.path.join(mapping_dir, path)
            path = os.path.abspath(path)
            if not os.path.exists(path):
                raise FileNotFoundError(f"{path} not found (from mapping file: {bcf})")
            out[contig] = path
        return out

    logging.info(f"Processing single BCF/VCF file: {bcf}")
    return {"all": bcf}


def normalize_bed_contigs(snplist, bcf_file):
    """Modify BED contigs to match the prefix style of bcf contigs."""
    if isinstance(bcf_file, str) and bcf_file.endswith(".txt"):
        mapping_dir = os.path.dirname(bcf_file)
        df = pd.read_csv(bcf_file, sep=r"\s+", header=None, engine="python")
        if df.shape[1] < 2:
            raise ValueError(f"Invalid mapping file: {bcf_file}")
        first_path = str(df.iloc[0, 1])
        if not os.path.isabs(first_path):
            first_path = os.path.join(mapping_dir, first_path)
        bcf_file = first_path

    elif isinstance(bcf_file, dict):
        bcf_file = next(iter(bcf_file.values()))

    with pysam.VariantFile(bcf_file, "r") as bcf:
        bcf_contigs = set(bcf.header.contigs)

    bcf_has_chr_prefix = any(contig.startswith("chr") for contig in bcf_contigs)

    if bcf_has_chr_prefix:
        snplist["contig_id"] = snplist["contig_id"].apply(
            lambda x: f"chr{x}" if not str(x).startswith("chr") else x
        )
    else:
        snplist["contig_id"] = snplist["contig_id"].apply(
            lambda x: str(x).lstrip("chr")
        )

    return snplist


def load_meta_data(path):
    """Load the JSON configuration file."""
    import json

    logging.info(f"Loading metadata from {path}")
    with open(path, "r") as f:
        return json.load(f)


def fetch_db(db_path, table):
    """Fetch data from the SQLite database."""
    import sqlite3

    with sqlite3.connect(db_path) as conn:
        df = pd.read_sql_query(f"SELECT * FROM {table}", conn)
    logging.info(f"Loaded {len(df)} SNP rows from table '{table}' in {db_path}")
    return df


class InvalidBedFormatError(Exception):
    """Exception raised for invalid BED file format."""

    pass
