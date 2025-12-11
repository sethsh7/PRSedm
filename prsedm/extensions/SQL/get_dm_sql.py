#!/usr/bin/env python3
"""
Locate or download the PRS-EDM DM SQL database (variants.db).

Robust for:
- Normal laptops/desktops
- Unix clusters where ~/.local/share may be restricted
"""

import os
import sys
import urllib.request
import shutil
import logging

ZENODO_SQL_URL = "https://zenodo.org/records/17903390/files/variants.db?download=1"

logger = logging.getLogger(__name__)


def _user_data_dir():
	"""
	Return a platform-appropriate *intended* user data directory for prsedm, e.g.:

	- Linux/Unix:  $XDG_DATA_HOME/prsedm or ~/.local/share/prsedm
	- macOS:       ~/Library/Application Support/prsedm
	- Windows:     %APPDATA%\\prsedm

	We will still check that it's writable before using it.
	"""
	home = os.path.expanduser("~")

	if sys.platform.startswith("win"):
		appdata = os.environ.get("APPDATA", home)
		return os.path.join(appdata, "prsedm")

	elif sys.platform == "darwin":
		return os.path.join(home, "Library", "Application Support", "prsedm")

	else:
		# Linux / Unix (respect XDG if set)
		base = os.environ.get("XDG_DATA_HOME", os.path.join(home, ".local", "share"))
		return os.path.join(base, "prsedm")


def _fallback_tmp_dir():
	"""
	Fallback location if user data dir is not writable, typically for clusters.
	"""
	user = os.environ.get("USER", "unknown")
	tmp_root = os.environ.get("TMPDIR", "/tmp")
	return os.path.join(tmp_root, f"prsedm-{user}")


def _ensure_dir_writable(path):
	"""
	Try to create `path` and verify it is writable.
	Return the path if OK, otherwise return None.
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


def get_dm_sql():
	"""
	Return the path to the DM SQL database (variants.db).

	Search order:
	  1. PRSEDM_SQL_PATH (if set and exists)
	  2. Dev repo copy: ../../PBC_for_SQL/variants.db
	  3. User data dir (XDG / ~/.local/share / macOS AppSupport):
	     - If not present, download to this location.
	     - If not writable, fall back to /tmp/prsedm-$USER/
	  4. If all fails, raise RuntimeError.

	No environment variables are required for normal use.
	"""

	# 1) Explicit override (optional)
	env_path = os.environ.get("PRSEDM_SQL_PATH")
	if env_path:
		if os.path.isfile(env_path):
			logger.info("Using SQL database from PRSEDM_SQL_PATH: %s", env_path)
			return env_path
		else:
			logger.warning(
				"PRSEDM_SQL_PATH is set but file does not exist: %s",
				env_path,
			)

	# 2) Development database: ../../PBC_for_SQL/variants.db (for you)
	this_dir = os.path.dirname(os.path.abspath(__file__))
	repo_root = os.path.abspath(os.path.join(this_dir, "..", ".."))
	dev_db = os.path.join(repo_root, "PBC_for_SQL", "variants.db")
	if os.path.isfile(dev_db):
		logger.info("Using SQL database from development folder: %s", dev_db)
		return dev_db

	# 3) User data directory (preferred persistent location)
	data_dir_candidate = _user_data_dir()
	data_dir = _ensure_dir_writable(data_dir_candidate)

	# 3b) If that fails (e.g. weird cluster home), fall back to /tmp/prsedm-$USER
	if data_dir is None:
		tmp_dir_candidate = _fallback_tmp_dir()
		data_dir = _ensure_dir_writable(tmp_dir_candidate)

	if data_dir is None:
		# We genuinely have nowhere we can write
		raise RuntimeError(
			"Could not find a writable directory for PRSEDM data. "
			"Set PRSEDM_SQL_PATH to an existing variants.db, or "
			"set PRSEDM_SQL_URL and ensure a writable TMPDIR."
		)

	user_db = os.path.join(data_dir, "variants.db")

	# If already present, use it
	if os.path.isfile(user_db):
		logger.info("Using SQL database from data dir: %s", user_db)
		return user_db

	# Not present anywhere → download into data_dir
	url = os.environ.get("PRSEDM_SQL_URL", ZENODO_SQL_URL)
	logger.info("Downloading SQL database from %s", url)
	logger.info("Saving to %s", user_db)

	try:
		with urllib.request.urlopen(url) as r, open(user_db, "wb") as out_f:
			shutil.copyfileobj(r, out_f)
	except Exception as e:
		logger.error("Failed to download SQL database: %s", e, exc_info=True)
		raise RuntimeError(
			"Could not download variants.db. "
			"Set PRSEDM_SQL_PATH to a local file or fix the download URL."
		)

	logger.info("Download complete.")
	return user_db


if __name__ == "__main__":
	logging.basicConfig(
		level=logging.INFO,
		format="%(asctime)s - %(levelname)s - %(message)s",
	)
	path = get_dm_sql()
	print(path)