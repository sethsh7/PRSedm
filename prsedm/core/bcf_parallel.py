# bcf_parallel.py
import os
import logging
from joblib import Parallel, delayed, parallel_backend
import pysam
from .scoring import score_geno, impute_score_ref
from .variant_processing import fetch_variant, geno_to_df

def process_batch(batch, bcf_files, samples, col, impute, refbcf):
	"""Process a batch of variants and return results."""
	batch_results, genotyped, imputed = [], 0, 0

	contig = batch.iloc[0]['contig_id']
	bcf_file = (
		bcf_files.get(contig) or 
		bcf_files.get(f"chr{contig}") or 
		bcf_files.get(contig.lstrip('chr')) or 
		bcf_files.get('all')
	)

	if not bcf_file or not os.path.isfile(bcf_file):
		logging.warning(f"Skipping batch: BCF file not found for contig {contig}")
		return [], 0, 0  # Skip this batch

	try:
		with pysam.VariantFile(bcf_file, 'r') as var_obj:
			for _, variant in batch.iterrows():
				logging.info(f"Processing variant {variant['contig_id']}:{variant['position']}...")
				geno_vcf = fetch_variant(var_obj, variant, samples)
				if geno_vcf:
					geno_df = geno_to_df(geno_vcf, samples)
					batch_results.append(score_geno(geno_df, variant, col))
					genotyped += 1
				elif impute and refbcf:
					logging.info(f"Imputing {variant['contig_id']}:{variant['position']}...")
					batch_results.append(impute_score_ref(samples, variant, refbcf))
					imputed += 1
	except Exception as e:
		logging.error(f"Error processing contig {contig}: {e}")

	return batch_results, genotyped, imputed

def process_batches_parallel(batches, bcf_files, samples, col, impute, refbcf, ntasks):
	"""Run batch processing in parallel."""
	with parallel_backend('loky', n_jobs=ntasks):
		results = Parallel()(
			delayed(process_batch)(batch, bcf_files, samples, col, impute, refbcf) 
			for batch in batches
		)
	return results