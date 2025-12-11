# score_bcf.py
import os
import logging
import pandas as pd
import pysam
from .utilities import (
	configure_logging,
	check_bed_type,
	get_samples,
	determine_bcf_type,
	normalize_bed_contigs,
)
from .bcf_parallel import process_batch, process_batches_parallel


# Configure logging
configure_logging()


def score_bcf(
		bcf,
		bed,
		col="GT",
		build="hg38",
		impute=False,
		refbcf=None,
		parallel=False,
		ntasks=os.cpu_count(),
		batch_size=1,
		stream=True):
	"""
	Score variants from BCF files with optional imputation.

	Parameters
	----------
	stream : bool
		If True (default, recommended for large PRS):
			- Returns a single per-sample PRS column ('sum').
			- Does NOT keep SNP-level DataFrames in memory at once.
		If False:
			- Returns full IID × SNP matrix (one column per variant)
			  PLUS a 'sum' column.

	parallel : bool
		Whether to use joblib parallelization.
	"""
	logging.info("Starting PRS scoring.")
	logging.info(f"stream={stream}, parallel={parallel}")

	# Prepare BCF files, sample data, and SNP list
	bcf_files = determine_bcf_type(bcf)
	samples = get_samples(
		pysam.VariantFile(
			next(iter(bcf_files.values())),
			'r'
		)
	)

	bed_df = check_bed_type(bed)
	snplist = normalize_bed_contigs(
		bed_df.rename(columns={f'position_{build}': 'position'}),
		bcf_files
	)

	print(f"SNP rows loaded from SQL: {len(bed_df)} rows")
	print(f"After contig normalization: {len(snplist)} rows")

	# Create batches for processing
	batches = [
		snplist[i:i + batch_size]
		for i in range(0, len(snplist), batch_size)
	]

	total_batched = sum(len(b) for b in batches)
	print(f"Number of SNP batches to process: {len(batches)}")
	print(f"Total SNPs across all batches: {total_batched}")

	# Use parallel processing if enabled
	if parallel:
		results = process_batches_parallel(
			batches,
			bcf_files,
			samples,
			col,
			impute,
			refbcf,
			ntasks,
			stream=stream,  # key: use generator when stream=True
		)
	else:
		# SERIAL: generator so we don't materialize all results at once
		results = (
			process_batch(batch, bcf_files, samples, col, impute, refbcf)
			for batch in batches
		)

	# Aggregate results
	var_out_list = [] if not stream else None  # only collect SNPs if stream=False
	total_genotyped, total_imputed = 0, 0

	# Streaming PRS total (per-sample)
	prs = pd.Series(
		0.0,
		index=pd.Index(samples, name='IID'),
		dtype='float64'
	)

	# Works whether results is a list or a generator
	for batch_results, genotyped, imputed in results:
		total_genotyped += genotyped
		total_imputed += imputed

		if not batch_results:
			continue

		# Update streaming PRS from each SNP column
		for df in batch_results:
			# df is (IID × 1)
			col_series = df.iloc[:, 0]
			# Align to prs index just in case
			col_series = col_series.reindex(prs.index)
			prs = prs.add(col_series, fill_value=0.0)

		if not stream:
			# Keep SNP-level results for full output mode
			var_out_list.extend(batch_results)
		# If stream=True, batch_results falls out of scope after this loop
		# iteration and its DataFrames can be garbage-collected.

	# Build final output
	if stream:
		# Memory-efficient mode: PRS-only
		score_out = prs.to_frame(name='sum')
		score_out.index.name = 'IID'
	else:
		# Full IID × SNP matrix + sum column
		if var_out_list:
			score_out = pd.concat(var_out_list, axis=1)
		else:
			score_out = pd.DataFrame(index=pd.Index(samples, name='IID'))
		score_out.index.name = 'IID'
		score_out['sum'] = prs

	logging.info(
		f"Completed with {total_genotyped} genotyped and {total_imputed} imputed variants."
	)
	logging.info(f"Final output shape: {score_out.shape}")
	return score_out