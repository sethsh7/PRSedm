# score_bcf.py
import os
import logging
import numpy as np
import pandas as pd
import pysam
from collections import deque
from joblib import Parallel, delayed
from .utilities import (
    configure_logging,
    check_bed_type,
    get_samples,
    determine_bcf_type,
    normalize_bed_contigs,
)
from .bcf_parallel import process_batch

# Configure logging
configure_logging()


def make_batches(snplist, batch_size=5000):
    """
    Split SNP list into contig-specific batches (batch_size only) with minimal logging,
    then interleave batches from different contigs to spread I/O across files.
    """
    contig_batches = []

    # First, create batches per contig
    for contig, contig_df in snplist.groupby("contig_id", sort=False):
        contig_df = contig_df.sort_values("position").reset_index(drop=True)
        start_idx = 0
        n = len(contig_df)
        batch_num = 1
        batches_for_contig = []

        while start_idx < n:
            end_idx = min(start_idx + batch_size, n)
            batch = contig_df.iloc[start_idx:end_idx].copy()
            batches_for_contig.append(batch)
            batch_num += 1
            start_idx = end_idx

        contig_batches.append(deque(batches_for_contig))

    # Interleave batches from different contigs
    interleaved = []
    while any(contig_batches):
        for dq in contig_batches:
            if dq:
                interleaved.append(dq.popleft())

    return interleaved


def score_bcf(
    bcf,
    bed,
    col="GT",
    build="hg38",
    estimate=None,
    ntasks=1,
    batch_size=5000,
    variant_log_path=None,
    prs_name=None,
    full=False,
):
    """
    Score variants from BCF files with optional estimation of missing variants.

    Parameters
    ----------
    full : bool
            If True:
                    - Write the full variant matrix to disk per batch
                    - Still maintains PRS sum
    """
    impute = bool(estimate)
    refbcf = estimate
    parallel = ntasks > 1

    logging.info("Starting PRS scoring.")
    logging.info(f"full={full}, parallel={parallel}, batch_size={batch_size}")

    # Prepare BCF files, sample data, and SNP list
    bcf_files = determine_bcf_type(bcf)

    with pysam.VariantFile(next(iter(bcf_files.values())), "r") as var_obj:
        samples = get_samples(var_obj)

    bed_df = check_bed_type(bed)
    snplist = normalize_bed_contigs(
        bed_df.rename(columns={f"position_{build}": "position"}), bcf_files
    )

    # Create interleaved contig-aware batches
    batches = make_batches(snplist, batch_size=batch_size)
    total_batched = sum(len(b) for b in batches)
    logging.info(f"Number of SNP batches to process: {len(batches)}")
    logging.info(f"Total SNPs across all batches: {total_batched}")

    # Prepare variant log file
    log_handle = None
    if variant_log_path:
        need_header = (
            not os.path.exists(variant_log_path)
            or os.path.getsize(variant_log_path) == 0
        )
        log_handle = open(variant_log_path, "a")
        if need_header:
            header = (
                (
                    "PRS\tcontig_id\tposition\tref\talt\teffect_allele\t"
                    "status\tvariant_af\tvariant_r2\n"
                )
                if prs_name
                else (
                    "contig_id\tposition\tref\talt\teffect_allele\tstatus\tvariant_af\tvariant_r2\n"
                )
            )
            log_handle.write(header)

    # Prepare variant matrix output file if full=True
    variant_matrix_path = None
    if full:
        variant_matrix_path = (
            f"{prs_name}_variant_matrix.tsv" if prs_name else "variant_matrix.tsv"
        )
        matrix_header_written = False
        variant_matrix_file = open(variant_matrix_path, "w")

    # Run parallel or sequential batches
    if parallel:
        logging.info(f"Attempting to use {ntasks} cores for parallel batch processing.")

        # Wrap process_batch with core logging
        def logged_process_batch(
            batch, bcf_files, samples, col, impute, refbcf, batch_idx
        ):
            import threading

            core_id = threading.get_ident()  # thread/process ID
            contigs_in_batch = batch["contig_id"].unique()
            logging.info(
                f"Core {core_id} starting batch {batch_idx} for contigs: {contigs_in_batch}"
            )
            result = process_batch(batch, bcf_files, samples, col, impute, refbcf)
            logging.info(f"Core {core_id} finished batch {batch_idx}")
            return result

        results = Parallel(n_jobs=ntasks, backend="loky", verbose=10)(
            delayed(logged_process_batch)(
                batch, bcf_files, samples, col, impute, refbcf, idx
            )
            for idx, batch in enumerate(batches, 1)
        )
    else:
        results = (
            process_batch(batch, bcf_files, samples, col, impute, refbcf)
            for batch in batches
        )

    # Aggregate results
    var_out_list = None if full else []
    var_names = None if full else []
    total_genotyped, total_imputed, total_missing = 0, 0, 0
    prs = np.zeros(len(samples), dtype=float)

    for (
        batch_results,
        batch_names,
        batch_log_rows,
        genotyped,
        imputed,
        missing,
    ) in results:
        total_genotyped += genotyped
        total_imputed += imputed
        total_missing += missing

        # Write variant log rows
        if log_handle and batch_log_rows:
            for row in batch_log_rows:
                line = "\t".join(map(str, row))
                if prs_name:
                    line = prs_name + "\t" + line
                log_handle.write(line + "\n")

        # Add to PRS sum
        for arr in batch_results:
            prs += arr

        # Write batch to variant matrix if full=True
        if full and batch_results:
            batch_matrix = pd.DataFrame(
                np.column_stack(batch_results), index=samples, columns=batch_names
            )
            if not matrix_header_written:
                batch_matrix.to_csv(variant_matrix_file, sep="\t", index=True)
                matrix_header_written = True
            else:
                batch_matrix.to_csv(
                    variant_matrix_file, sep="\t", index=True, header=False
                )

        if not full:
            var_out_list.extend(batch_results)
            var_names.extend(batch_names)

    # Close files
    if log_handle:
        log_handle.close()
    if full:
        variant_matrix_file.close()

    # Build final output
    if full:
        score_out = pd.DataFrame({"sum": prs}, index=pd.Index(samples, name="IID"))
    else:
        if var_out_list:
            score_matrix = np.column_stack(var_out_list)
            score_out = pd.DataFrame(
                score_matrix, index=pd.Index(samples, name="IID"), columns=var_names
            )
        else:
            score_out = pd.DataFrame(index=pd.Index(samples, name="IID"))
        score_out.index.name = "IID"
        score_out["sum"] = prs

    logging.info(
        f"Completed with {total_genotyped} available, {total_imputed} estimated and {total_missing} missing variants."
    )
    stats = [total_genotyped, total_imputed, total_missing]
    return score_out, stats
