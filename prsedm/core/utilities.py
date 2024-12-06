"""Contains utility functions for the SRC package."""
import os
import logging
from importlib.resources import files
from dataclasses import dataclass, field
import pandas as pd
import pysam


@dataclass
class PRSConfig:
    """Common configuration class for PRS scoring."""
    bcf: str  # Mandatory, fixed bcf path
    col: str = field(default="GT")  # Optional, default to "GT"
    build: str = "hg38"  # Optional, default to "hg38"
    impute: bool = False  # Optional, default to False
    refbcf: str = None  # Optional, can be None
    parallel: bool = False  # Optional, default to False
    ntasks: int = os.cpu_count()  # Optional, defaults to CPU count
    batch_size: int = 1  # Optional, default to 1

    def __post_init__(self):
        """Validate the 'col' parameter."""
        if self.col not in {"GT", "GP"}:
            raise ValueError(
                f"Invalid value for 'col': {self.col}. "
                " Must be 'GT' or 'GP'."
            )


def configure_logging():
    log_file = "f{__name__}.log"
    """Configure logging with both file and stream handlers."""
    log_path = os.path.join(os.path.dirname(__file__), "..", log_file)
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
        handlers=[logging.FileHandler(log_path), logging.StreamHandler()]
    )


def get_samples(var_obj):
    """Extract and return sample names from a bcf object."""
    return list(var_obj.header.samples)


def check_bed_type(bed):
    """Parse a BED file or DataFrame and return it as a DataFrame."""
    if isinstance(bed, pd.DataFrame):
        return bed
    if not os.path.isfile(bed):
        raise FileNotFoundError(f"'{bed}' not found.")
    try:
        df = pd.read_whitespace(bed)
        if df.shape[1] < 3:
            raise InvalidBedFormatError(f"Invalid BED format in '{bed}'.")
        return df
    except Exception as e:
        raise ValueError(f"Error reading '{bed}': {e}")


def read_bcf_mapping(file_path):
    """Read bcf-to-contig mapping from a text file."""
    bcf_files = {}
    with open(file_path, 'r') as f:
        for line in f:
            filename, contig = line.strip().split()
            bcf_files[contig] = filename
    return bcf_files


def determine_bcf_type(bcf):
    """Determine if the bcf input is a mapping text file or a single bcf file."""
    try:
        with open(bcf, 'r') as f:
            f.read(1024)  # Test if it's a text file
        logging.info(f"Reading bcf mapping from {bcf}")
        mapping_dir = os.path.dirname(bcf)
        df = read_whitespace(bcf)
        return {contig: os.path.join(mapping_dir, path)
                for contig, path in df.values}
    except (UnicodeDecodeError, OSError):
        logging.info(f"Processing a single bcf or binary file: {bcf}")
        return {"all": bcf}


def check_index_bcf(bcf_path):
    """Ensure the bcf/vcf file is indexed."""
    index_path = f"{bcf_path}.tbi" if bcf_path.endswith(
        ".bcf.gz") else f"{bcf_path}.csi"
    if not os.path.isfile(index_path):
        logging.info(f"Indexing {bcf_path}...")
        try:
            preset = "bcf" if bcf_path.endswith(".bcf.gz") else None
            pysam.tabix_index(
                bcf_path,
                preset=preset,
                force=True,
                csi=bcf_path.endswith(".bcf"))
        except Exception as e:
            raise RuntimeError(f"Failed to index {bcf_path}: {e}")


def normalize_bed_contigs(snplist, bcf_file):
    """Modify BED contigs to match the prefix style of bcf contigs."""
    # Handle unloaded txt list, loaded map, or single file
    if isinstance(bcf_file, str) and bcf_file.endswith('.txt'):
        with open(bcf_file, 'r') as f:
            bcf_file = f.readline().strip()  # Use the first bcf file listed in the text file
    # Handle a dictionary of bcf files by selecting the first one
    elif isinstance(bcf_file, dict):
        bcf_file = next(iter(bcf_file.values()))

    with pysam.VariantFile(bcf_file, 'r') as bcf:
        bcf_contigs = set(bcf.header.contigs)
    bcf_has_chr_prefix = any(contig.startswith('chr')
                             for contig in bcf_contigs)
    if bcf_has_chr_prefix:
        snplist['contig_id'] = snplist['contig_id'].apply(
            lambda x: f"chr{x}" if not str(x).startswith('chr') else x
        )
    else:
        snplist['contig_id'] = snplist['contig_id'].apply(
            lambda x: x.lstrip('chr'))
    return snplist


def read_whitespace(file):
    """Read a whitespace-separated file into a DataFrame."""
    return pd.read_csv(file, header=None, sep=r'\t|\s{4}', engine='python')


def save_csv_plain(df, file):
    """Save a DataFrame to a tab-separated file without headers or indices."""
    df.to_csv(file, header=False, index=False, sep="\t")


def load_meta_data(path):
    """Load the JSON configuration file."""
    import json
    logging.info(f"Loading metadata from {path}")
    with open(path, 'r') as f:
        return json.load(f)


def fetch_db(db_path, table):
    """Fetch data from the SQLite database."""
    import sqlite3
    return pd.read_sql_query(
        f"SELECT * FROM {table}", sqlite3.connect(db_path))


def get_snp_db(score_name):
    """Return filtered SNP data from the database."""
    ext_path = files(__name__.split('.')[0]) / 'extensions'
    meta = load_meta_data(ext_path / 'JSON' / 'prs_meta.json')
    db_path = ext_path / 'SQL' / 'variants.db'
    tables = [
        meta[score_name].get(t) for t in (
            "db_table",
            "db_dq") if meta[score_name].get(t)]
    df = pd.concat([fetch_db(db_path, t) for t in tables], ignore_index=True)
    cols = [c for c in df.columns if c.startswith(
        "position")] + ['rsid', 'contig_id', 'effect_allele']
    return df[cols]


class InvalidBedFormatError(Exception):
    """Exception raised for invalid BED file format."""
    pass


class VariantNotFoundError(Exception):
    """Exception raised when a variant cannot be found in the bcf or bcf."""
    pass


class FileReadError(Exception):
    """Exception raised for errors in reading files."""
    pass


class InvalidbcfFormatError(Exception):
    """Exception raised for invalid bcf file format."""
    pass
