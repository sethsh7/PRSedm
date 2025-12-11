"""Main score_dm extension for predefined diabetes PRS."""
import argparse
import logging
import os
from dataclasses import asdict
from importlib.resources import files

import pandas as pd

from .grouped_scoring import score_grouped
from .hla_int_grs import score_int_hla
from .SQL.get_dm_sql import get_dm_sql
from ..core.score_bcf import score_bcf
from ..core.utilities import (
	configure_logging,
	fetch_db,
	load_meta_data,
	PRSConfig,
	normalize_bed_contigs,
)


configure_logging()


def process_flag(flag, meta, db_path, config, norm=False, full=False):
	"""Process a flag and add its total optionally normalized score as a new column."""
	logging.info(f"Processing flag: {flag}")
	score = fetch_db(db_path, meta[flag]["db_table"])
	method = meta[flag].get("method", "additive")
	if "stream" in meta[flag]:
		full = not meta[flag]["stream"]

	if method == "grouped":
		logging.info(f"Generating grouped scores for {flag}")
		result, group_total_cols = score_grouped(bed=score, config=config, full=full)
		result[f"{flag}_total"] = result[group_total_cols].sum(axis=1)

		if not full:
			# compact mode: just group totals + total
			result = result[group_total_cols + [f"{flag}_total"]]

	elif method == "hla_int":
		logging.info(f"Generating HLA interaction scores for {flag}")
		dq, int_df, rank = (
			fetch_db(db_path, meta[flag][k]) for k in ("db_dq", "db_int", "db_rank")
		)
		dq = normalize_bed_contigs(dq, config.bcf)
		result = score_int_hla(score, dq, int_df, rank, config, full=full, flag=flag)

	else:
		logging.info(f"Generating additive scores for {flag}")

		# full=False  → stream=True  (PRS sum only, memory efficient)
		# full=True   → stream=False (variant matrix + sum)
		stream = not full

		result = score_bcf(
			bed=score,
			stream=stream,
			**asdict(config),
		)

		if stream:
			# PRS-only output; score_bcf returns a column named 'sum'
			result.rename(columns={"sum": f"{flag}_total"}, inplace=True)

		else:
			# score_bcf returned variant columns + 'sum'
			rename_map = {}
			for col in result.columns:
				if col == "sum":
					rename_map[col] = f"{flag}_total"
				else:
					rename_map[col] = f"{flag}_{col}"
			result.rename(columns=rename_map, inplace=True)

			if not full:
				# Should not happen normally since stream = not full
				result = result[[f"{flag}_total"]]

	# Normalization block (unchanged)
	if norm:
		flag_total_col = f"{flag}_total"
		norm_col = f"{flag}_norm"
		min_value = meta.get(flag, {}).get("min", result[flag_total_col].min())
		max_value = meta.get(flag, {}).get("max", result[flag_total_col].max())
		result[norm_col] = (result[flag_total_col] - min_value) / (
			max_value - min_value
		)
		logging.info(
			f"Added normalized column '{norm_col}' using min={min_value} and max={max_value}."
		)

	logging.info(f"Processed flag '{flag}' with total column '{flag}_total'.")

	return result


def gen_dm(
	vcf,
	col,
	scores,
	build="hg38",
	impute=False,
	refbcf=None,
	norm=False,
	parallel=False,
	ntasks=os.cpu_count(),
	batch_size=1,
	full=False,
):
	"""Generate DM-related PRS scores using a configuration object."""
	if impute and not refbcf:
		raise ValueError("Reference directory is required if imputation is enabled.")

	logging.info("Starting PRS score generation...")

	# Store the parameters in PRSConfig
	config = PRSConfig(
		bcf=vcf,
		col=col,
		build=build,
		impute=impute,
		refbcf=refbcf,
		parallel=parallel,
		ntasks=ntasks,
		batch_size=batch_size,
	)

	# Resolve database and metadata paths
	db_path = get_dm_sql()
	ext_path = files(__name__.split(".")[0]) / "extensions"
	meta_path = ext_path / "JSON" / "prs_meta.json"
	meta = load_meta_data(meta_path)

	# Process each flag and collect outputs
	outputs = []
	for flag in scores.split(","):
		if flag not in meta:
			logging.warning(f"Flag '{flag}' not found in metadata. Skipping.")
			continue
		outputs.append(process_flag(flag, meta, db_path, config, norm, full))

	# Concatenate the results or return an empty DataFrame
	result = pd.concat(outputs, axis=1) if outputs else pd.DataFrame()
	logging.info("PRS score generation completed.")
	return result


def main():
	"""Main CLI entry point."""
	parser = argparse.ArgumentParser(
		description="Calculate predefined DM-related PRS."
	)
	parser.add_argument(
		"--vcf",
		help="Path to a vcf(.gz)/bcf file, or text file mapping BCF files to contigs.",
	)
	parser.add_argument(
		"--col",
		default="GT",
		choices=["GT", "GP"],
		help="Genotype column to score (GT=WGS or GP=Imputed array).",
	)
	parser.add_argument(
		"--scores",
		help="Comma-separated list of scores (e.g., 'PRS1,PRS2').",
	)
	parser.add_argument(
		"--build",
		default="hg38",
		choices=["hg19", "hg38"],
		help="Genome build to use (hg19/hg38).",
	)
	parser.add_argument(
		"--impute",
		action="store_true",
		help="Enable imputation.",
	)
	parser.add_argument(
		"--refvcf",
		default=None,
		help="Reference directory (required if imputation is enabled).",
	)
	parser.add_argument(
		"--norm",
		action="store_true",
		help="Perform fixed MinMax normalization (not recommended if impute is off).",
	)
	parser.add_argument(
		"--parallel",
		action="store_true",
		help="Enable parallel processing",
	)
	parser.add_argument(
		"--ntasks",
		type=int,
		default=os.cpu_count(),
		help="Number of tasks to use (default=max available).",
	)
	parser.add_argument(
		"--batch-size",
		type=int,
		default=1,
		help="Number of variants per batch (default=1).",
	)
	parser.add_argument(
		"--output",
		default="results.csv",
		help="Path to save the output file (default: 'results.csv').",
	)
	parser.add_argument(
		"--full",
		action="store_true",
		help="Include individual variant scores with PRS name prepended.",
	)
	parser.add_argument(
		"--getsql",
		action="store_true",
		help="Download or locate the PRS SQL database (variants.db) and exit.",
	)

	args = parser.parse_args()
	logging.info("Parsed command line arguments.")

	# Handle SQL-only mode *before* enforcing other required args
	if args.getsql:
		try:
			db_path = get_dm_sql()
			print(db_path)
			logging.info(f"SQL database available at: {db_path}")
		except Exception:
			logging.error("Failed to locate or download SQL database.", exc_info=True)
			raise
		return

	# Enforce required args only when actually scoring
	if args.impute and not args.refvcf:
		parser.error("--refvcf is required if imputation is enabled.")

	if not args.vcf:
		parser.error("--vcf is required unless --getsql is used.")

	if not args.scores:
		parser.error("--scores is required unless --getsql is used.")

	try:
		# Call the gen_dm function with parsed CLI arguments
		result = gen_dm(
			vcf=args.vcf,
			col=args.col,
			scores=args.scores,
			build=args.build,
			impute=args.impute,
			refbcf=args.refvcf,
			norm=args.norm,
			parallel=args.parallel,
			ntasks=args.ntasks,
			batch_size=args.batch_size,
			full=args.full,
		)
		logging.info("Generated PRS scores successfully.")

		# Save the result to the specified output file
		result.to_csv(args.output, index=True)
		logging.info(f"Results saved to {args.output}")

	except Exception as e:
		logging.error(f"Error during PRS generation: {e}", exc_info=True)
		raise  # Re-raise the exception for further handling if needed


if __name__ == "__main__":
	main()