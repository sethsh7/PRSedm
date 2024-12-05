# score_bcf.py
import os
import argparse
import logging
import pandas as pd
from dataclasses import asdict
import pysam
from .utilities import (
	configure_logging,
	check_bed_type,
	get_samples,
	determine_bcf_type,
	normalize_bed_contigs,
	PRSConfig,
)
from .bcf_parallel import process_batch, process_batches_parallel

# Configure logging
configure_logging()

def score_bcf(bcf, bed, col="GT", build="hg38", impute=False, refbcf=None, 
              parallel=False, ntasks=os.cpu_count(), batch_size=1):
	"""Score variants from BCF files with optional imputation."""
	logging.info("Starting PRS scoring.")

	# Prepare BCF files, sample data, and SNP list
	bcf_files = determine_bcf_type(bcf)
	samples = get_samples(pysam.VariantFile(next(iter(bcf_files.values())), 'r'))
	snplist = normalize_bed_contigs(
		check_bed_type(bed).rename(columns={f'position_{build}': 'position'}), bcf_files
	)

	# Create batches for processing
	batches = [snplist[i:i + batch_size] for i in range(0, len(snplist), batch_size)]

	# Use parallel processing if enabled
	results = (
		process_batches_parallel(batches, bcf_files, samples, col, impute, refbcf, ntasks)
		if parallel else
		[process_batch(batch, bcf_files, samples, col, impute, refbcf) for batch in batches]
	)

	# Aggregate results
	var_out_list, total_genotyped, total_imputed = [], 0, 0
	for batch_results, genotyped, imputed in results:
		var_out_list.extend(batch_results)
		total_genotyped += genotyped
		total_imputed += imputed

	# Combine all batch results into a DataFrame
	score_out = pd.concat(var_out_list, axis=1) if var_out_list else pd.DataFrame()
	score_out.index.name = 'IID'

	logging.info(f"Completed with {total_genotyped} genotyped and {total_imputed} imputed variants.")
	return score_out