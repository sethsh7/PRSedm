#!/usr/bin/env python3
"""
Locate or download the PRS-EDM DM SQL database (variants.db) and metadata
(prs_meta.json).

Robust for:
- Normal laptops/desktops
- Unix clusters where ~/.local/share may be restricted

Public API:
- get_dm_sql()         -> returns path to variants.db
- get_dm_meta()        -> returns path to prs_meta.json
- get_snp_db(score)    -> returns SNP rows for a score from variants.db
"""

import os
import sys
import urllib.request
import shutil
import logging
import pandas as pd

from prsedm.core.utilities import load_meta_data, fetch_db

ZENODO_BASE_URL = "https://zenodo.org/records/19373192/files/"
ZENODO_SQL_URL = f"{ZENODO_BASE_URL}/variants.db?download=1"
ZENODO_META_URL = f"{ZENODO_BASE_URL}/prs_meta.json?download=1"

logger = logging.getLogger(__name__)


def _user_data_dir():
    """
    Return a platform-appropriate intended user data directory for prsedm.
    """
    home = os.path.expanduser("~")

    if sys.platform.startswith("win"):
        appdata = os.environ.get("APPDATA", home)
        return os.path.join(appdata, "prsedm")

    elif sys.platform == "darwin":
        return os.path.join(home, "Library", "Application Support", "prsedm")

    else:
        base = os.environ.get("XDG_DATA_HOME", os.path.join(home, ".local", "share"))
        return os.path.join(base, "prsedm")


def _fallback_tmp_dir():
    """
    Fallback location if user data dir is not writable (common on clusters).
    """
    user = os.environ.get("USER", "unknown")
    tmp_root = os.environ.get("TMPDIR", "/tmp")
    return os.path.join(tmp_root, f"prsedm-{user}")


def _ensure_dir_writable(path):
    """
    Create path and verify it is writable.
    Return path if OK, otherwise None.
    """
    try:
        os.makedirs(path, exist_ok=True)
        test_file = os.path.join(path, ".write_test")
        with open(test_file, "w") as f:
            f.write("ok")
        os.remove(test_file)
        return path
    except Exception as e:
        logger.warning("Directory %s not writable: %s", path, e)
        return None


def _get_data_dir():
    """
    Return a writable PRSEDM data directory.
    """
    data_dir_candidate = _user_data_dir()
    data_dir = _ensure_dir_writable(data_dir_candidate)

    if data_dir is None:
        tmp_dir_candidate = _fallback_tmp_dir()
        data_dir = _ensure_dir_writable(tmp_dir_candidate)

    if data_dir is None:
        raise RuntimeError(
            "Could not find a writable directory for PRSEDM data. "
            "Set PRSEDM_SQL_PATH / PRSEDM_META_PATH to existing local files."
        )

    return data_dir


def _get_or_download_file(env_path_var, env_url_var, default_url, filename, label):
    """
    Resolve a local file path, downloading if needed.

    Search order:
      1. Explicit override path via env_path_var
      2. Writable PRSEDM data dir
         - Use existing local file if present
         - Otherwise download from env_url_var or default_url
    """
    env_path = os.environ.get(env_path_var)
    if env_path:
        if os.path.isfile(env_path):
            logger.info("Using %s from %s: %s", label, env_path_var, env_path)
            return env_path
        else:
            logger.warning(
                "%s is set but file does not exist: %s", env_path_var, env_path
            )

    data_dir = _get_data_dir()
    local_path = os.path.join(data_dir, filename)

    if os.path.isfile(local_path):
        logger.info("Using %s from data dir: %s", label, local_path)
        return local_path

    url = os.environ.get(env_url_var, default_url)
    logger.info("Downloading %s from %s", label, url)
    logger.info("Saving to %s", local_path)

    try:
        with urllib.request.urlopen(url) as r, open(local_path, "wb") as out_f:
            shutil.copyfileobj(r, out_f)
    except Exception as e:
        logger.error("Failed to download %s: %s", label, e, exc_info=True)
        raise RuntimeError(
            f"Could not download {filename}. "
            f"Set {env_path_var} to a local file or fix the download URL."
        )

    logger.info("Download complete for %s.", label)
    return local_path


def get_dm_sql():
    """
    Return the path to the DM SQL database (variants.db).

    Search order:
      1. PRSEDM_SQL_PATH
      2. PRSEDM data dir
      3. Download from PRSEDM_SQL_URL or built-in Zenodo URL
    """
    return _get_or_download_file(
        env_path_var="PRSEDM_SQL_PATH",
        env_url_var="PRSEDM_SQL_URL",
        default_url=ZENODO_SQL_URL,
        filename="variants.db",
        label="SQL database",
    )


def get_dm_meta():
    """
    Return the path to the PRS metadata JSON (prs_meta.json).

    Search order:
      1. PRSEDM_META_PATH
      2. PRSEDM data dir
      3. Download from PRSEDM_META_URL or built-in Zenodo URL
    """
    return _get_or_download_file(
        env_path_var="PRSEDM_META_PATH",
        env_url_var="PRSEDM_META_URL",
        default_url=ZENODO_META_URL,
        filename="prs_meta.json",
        label="PRS metadata JSON",
    )


def get_snp_db(score_name):
    """
    Return filtered SNP data for `score_name` from variants.db,
    independent of where the database lives.
    """
    meta_path = get_dm_meta()
    meta = load_meta_data(meta_path)

    if score_name not in meta:
        raise KeyError(f"Score '{score_name}' not found in prs_meta.json")

    db_path = get_dm_sql()

    tables = [
        meta[score_name].get(t)
        for t in ("db_table", "db_dq")
        if meta[score_name].get(t)
    ]

    if not tables:
        raise ValueError(
            f"No db_table/db_dq defined for score '{score_name}' in prs_meta.json"
        )

    try:
        df = pd.concat([fetch_db(db_path, t) for t in tables], ignore_index=True)
    except Exception as e:
        raise RuntimeError(
            f"Failed to load SNP data for '{score_name}' from database: {db_path}. "
            f"Tables attempted: {tables}. Underlying error: {e}"
        )

    return df


if __name__ == "__main__":
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s - %(levelname)s - %(message)s",
    )

    print("SQL:", get_dm_sql())
    print("META:", get_dm_meta())
