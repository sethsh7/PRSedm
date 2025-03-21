"""Main score_dm extension for predefined diabetes PRS."""
import argparse
import logging
import os
from dataclasses import asdict
from importlib.resources import files
import pandas as pd
from .grouped_scoring import score_grouped
from .hla_int_grs import score_int_hla
from ..core.score_bcf import score_bcf
from ..core.utilities import (
    configure_logging,
    fetch_db,
    load_meta_data,
    PRSConfig,
    normalize_bed_contigs
)


configure_logging()


def process_flag(flag, meta, db_path, config, norm=False):
    """Process a flag and add its total optionally normalized score as a new column."""
    logging.info(f"Processing flag: {flag}")
    score = fetch_db(db_path, meta[flag]['db_table'])
    method = meta[flag].get('method', 'additive')

    if method == 'grouped':
        logging.info(f"Generating grouped scores for {flag}")
        result = score_grouped(bed=score, config=config)
    elif method == 'hla_int':
        logging.info(f"Generating HLA interaction scores for {flag}")
        dq, int_df, rank = (
            fetch_db(
                db_path, meta[flag][k]) for k in (
                'db_dq', 'db_int', 'db_rank'))
        dq = normalize_bed_contigs(dq, config.bcf)
        result = score_int_hla(score, dq, int_df, rank, config)
    else:
        logging.info(f"Generating additive scores for {flag}")
        result = score_bcf(**asdict(config), bed=score)
        return pd.DataFrame(
            {f"{flag}_total": result.select_dtypes('number').sum(axis=1)})

    # Add the 'flag_total' column
    result = pd.concat([result, result.select_dtypes(
        'number').sum(axis=1).rename(f"{flag}_total")], axis=1)

    # Add normalized column if norm is True
    if norm:
        flag_total_col = f"{flag}_total"
        norm_col = f"{flag}_norm"
        min_value = meta.get(flag, {}).get('min', result[flag_total_col].min())
        max_value = meta.get(flag, {}).get('max', result[flag_total_col].max())

        result[norm_col] = (result[flag_total_col] -
                            min_value) / (max_value - min_value)
        logging.info(
            f"Added normalized column '{norm_col}' using min={min_value} and max={max_value}.")

    logging.info(f"Processed flag '{flag}' with total column '{flag}_total'.")

    return result


def gen_dm(vcf, col, prsflags, build="hg38", impute=False, refbcf=None,
           norm=False, parallel=False, ntasks=os.cpu_count(), batch_size=1):
    """Generate DM-related PRS scores using a configuration object."""
    if impute and not refbcf:
        raise ValueError(
            "Reference directory is required if imputation is enabled.")

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
        batch_size=batch_size
    )

    # Load metadata and database paths
    ext_path = files(__name__.split('.')[0]) / 'extensions'
    db_path = ext_path / 'SQL' / 'variants.db'
    meta_path = ext_path / 'JSON' / 'prs_meta.json'
    meta = load_meta_data(meta_path)

    # Process each flag and collect outputs
    outputs = []
    for flag in prsflags.split(","):
        if flag not in meta:
            logging.warning(f"Flag '{flag}' not found in metadata. Skipping.")
            continue
        outputs.append(process_flag(flag, meta, db_path, config, norm))

    # Concatenate the results or return an empty DataFrame
    result = pd.concat(outputs, axis=1) if outputs else pd.DataFrame()
    logging.info("PRS score generation completed.")
    return result


def main():
    """Main CLI entry point."""
    parser = argparse.ArgumentParser(
        description="Calculate predefined DM-related PRS.")
    parser.add_argument(
        "--vcf",
        required=True,
        help="Path to a vcf(.gz)/bcf file, or text file mapping BCF files to contigs.")
    parser.add_argument(
        "--col",
        default="GT",
        choices=[
            "GT",
            "GP"],
        help="Genotype column to score (GT=WGS or GP=Imputed array).")
    parser.add_argument(
        "--scores",
        required=True,
        help="Comma-separated list of scores (e.g., 'PRS1,PRS2').")
    parser.add_argument(
        "--build",
        default="hg38",
        choices=[
            "hg19",
            "hg38"],
        help="Genome build to use (hg19/hg38).")
    parser.add_argument(
        "--impute",
        action="store_true",
        help="Enable imputation.")
    parser.add_argument(
        "--refvcf",
        default=None,
        help="Reference directory (required if imputation is enabled).")
    parser.add_argument(
        "--norm",
        action="store_true",
        help="Perform fixed MinMax normalization (not recommended if impute is off).")
    parser.add_argument(
        "--parallel",
        action="store_true",
        help="Enable parallel processing")
    parser.add_argument(
        "--ntasks",
        type=int,
        default=os.cpu_count(),
        help="Number of tasks to use (default=max available).")
    parser.add_argument(
        "--batch-size",
        type=int,
        default=1,
        help="Number of variants per batch (default=1).")
    parser.add_argument(
        "--output",
        default="results.csv",
        help="Path to save the output file (default: 'results.csv').")

    args = parser.parse_args()
    logging.info("Parsed command line arguments.")

    if args.impute and not args.refvcf:
        parser.error("--refvcf is required if imputation is enabled.")

    try:
        # Call the gen_dm function with parsed CLI arguments
        result = gen_dm(
            vcf=args.vcf,
            col=args.col,
            prsflags=args.scores,
            build=args.build,
            impute=args.impute,
            refbcf=args.refvcf,
            norm=args.norm,
            parallel=args.parallel,
            ntasks=args.ntasks,
            batch_size=args.batch_size
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
