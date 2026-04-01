"""Main score_dm extension for predefined diabetes PRS."""

import argparse
import logging
import os
import sys
import pandas as pd

from .grouped_scoring import score_grouped
from .hla_int_grs import score_int_hla
from .SQL.get_dm_sql import get_dm_sql, get_dm_meta
from ..core.score_bcf import score_bcf
from ..core.utilities import (
    configure_logging,
    fetch_db,
    load_meta_data,
    PRSConfig,
    normalize_bed_contigs,
)

configure_logging()

_original_unraisablehook = sys.unraisablehook


def _silence_resource_tracker_errors(unraisable):
    if issubclass(unraisable.exc_type, ChildProcessError) and "ResourceTracker" in str(
        unraisable.object
    ):
        pass
    else:
        _original_unraisablehook(unraisable)


sys.unraisablehook = _silence_resource_tracker_errors


def make_variant_log_path(vcf, flag):
    """Create a sensible per-PRS variant log filename."""
    vcf_name = os.path.basename(vcf)

    for suffix in (".vcf.gz", ".bcf.gz", ".vcf", ".bcf"):
        if vcf_name.endswith(suffix):
            vcf_name = vcf_name[: -len(suffix)]
            break

    return f"{flag}_{vcf_name}_variant_log.tsv"


def load_proxy_substitutions(proxy):
    """Load proxy substitution file."""
    proxy_df = pd.read_csv(proxy, sep=r"\s+", engine="python")
    required = {
        "target_rsid",
        "target_contig_id",
        "target_position",
        "target_effect_allele",
        "sub_rsid",
        "sub_contig_id",
        "sub_position",
        "sub_effect_allele",
    }
    missing = required - set(proxy_df.columns)
    if missing:
        raise ValueError(
            f"Proxy file '{proxy}' is missing required columns: {sorted(missing)}"
        )
    return proxy_df


def apply_proxy_substitutions(score, proxy_df, build):
    """Substitute requested variants with proxy variants while keeping score weights/groups."""
    score = score.copy()
    pos_col = "position_hg38" if build == "hg38" else "position_hg19"

    n_substituted = 0
    for _, row in proxy_df.iterrows():
        mask = (
            (score["rsid"].astype(str) == str(row["target_rsid"]))
            & (score["contig_id"].astype(str) == str(row["target_contig_id"]))
            & (score[pos_col].astype(int) == int(row["target_position"]))
            & (score["effect_allele"].astype(str) == str(row["target_effect_allele"]))
        )

        n_matches = int(mask.sum())
        if n_matches:
            score.loc[mask, "rsid"] = str(row["sub_rsid"])
            score.loc[mask, "contig_id"] = str(row["sub_contig_id"])
            score.loc[mask, pos_col] = int(row["sub_position"])
            score.loc[mask, "effect_allele"] = str(row["sub_effect_allele"])
            n_substituted += n_matches

    logging.info(f"Applied {n_substituted} proxy substitutions.")
    return score


def process_flag(flag, meta, db_path, config, full=False, proxy_df=None):
    """Process a flag and add its total plus normalization if appropriate."""
    logging.info(f"Processing flag: {flag}")
    score = fetch_db(db_path, meta[flag]["db_table"])
    if proxy_df is not None:
        score = apply_proxy_substitutions(score, proxy_df, config.build)
    method = meta[flag].get("method", "additive")

    variant_log_path = make_variant_log_path(config.bcf, flag)
    if os.path.exists(variant_log_path):
        os.remove(variant_log_path)

    if method == "grouped":
        logging.info(f"Generating grouped scores for {flag}")
        result, group_total_cols, stats = score_grouped(
            bed=score,
            config=config,
            full=full,
            variant_log_path=variant_log_path,
        )
        rename_map = {col: f"{flag}__{col}" for col in result.columns}
        result.rename(columns=rename_map, inplace=True)
        group_total_cols = [f"{flag}__{col}" for col in group_total_cols]
        result[f"{flag}__total"] = result[group_total_cols].sum(axis=1)

        if not full:
            result = result[group_total_cols]

    elif method == "hla_int":
        logging.info(f"Generating HLA interaction scores for {flag}")
        dq, int_df, rank = (
            fetch_db(db_path, meta[flag][k]) for k in ("db_dq", "db_int", "db_rank")
        )
        dq = normalize_bed_contigs(dq, config.bcf)
        result, stats = score_int_hla(
            score,
            dq,
            int_df,
            rank,
            config,
            full=full,
            flag=flag,
            variant_log_path=variant_log_path,
        )

    else:
        logging.info(f"Generating additive scores for {flag}")
        result, stats = score_bcf(
            bcf=config.bcf,
            bed=score,
            col=config.col,
            build=config.build,
            estimate=config.refbcf,
            ntasks=config.ntasks,
            batch_size=config.batch_size,
            variant_log_path=variant_log_path,
            prs_name=flag,
            full=full,
        )

        rename_map = {}
        for col in result.columns:
            if col == "sum":
                rename_map[col] = f"{flag}__total"
            else:
                rename_map[col] = f"{flag}__{col}"
        result.rename(columns=rename_map, inplace=True)

        if not full:
            result = result[[f"{flag}__total"]]

    if f"{flag}__total" in result.columns:
        flag_total_col = f"{flag}__total"
        norm_col = f"{flag}__norm"

        if stats[2] > 0:
            logging.warning(
                f"Skipping normalization for '{flag}' because one or more variants were missing. "
                f"Use --estimate with a reference directory to estimate missing variants."
            )
        else:
            min_value = meta.get(flag, {}).get("min", result[flag_total_col].min())
            max_value = meta.get(flag, {}).get("max", result[flag_total_col].max())
            result[norm_col] = (result[flag_total_col] - min_value) / (
                max_value - min_value
            )
            logging.info(
                f"Added normalized column '{norm_col}' using min={min_value} and max={max_value}."
            )

    logging.info(f"Processed flag '{flag}'.")
    return result


def gen_dm(
    vcf,
    col,
    scores,
    build="hg38",
    estimate=None,
    ntasks=1,
    batch_size=5000,
    full=False,
    proxy=None,
):
    """Generate DM-related PRS scores using a configuration object."""
    logging.info("Starting PRS score generation...")

    config = PRSConfig(
        bcf=vcf,
        col=col,
        build=build,
        impute=bool(estimate),
        refbcf=estimate,
        parallel=(ntasks > 1),
        ntasks=ntasks,
        batch_size=batch_size,
    )

    db_path = get_dm_sql()
    meta_path = get_dm_meta()
    meta = load_meta_data(meta_path)
    proxy_df = load_proxy_substitutions(proxy) if proxy else None

    outputs = []
    for flag in scores.split(","):
        if flag not in meta:
            logging.warning(f"Flag '{flag}' not found in metadata. Skipping.")
            continue
        outputs.append(process_flag(flag, meta, db_path, config, full, proxy_df))

    result = pd.concat(outputs, axis=1) if outputs else pd.DataFrame()
    logging.info("PRS score generation completed.")
    return result


def main():
    """Main CLI entry point."""
    parser = argparse.ArgumentParser(description="Calculate predefined DM-related PRS.")
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
        "--scores", help="Comma-separated list of scores (e.g., 'PRS1,PRS2')."
    )
    parser.add_argument(
        "--build",
        default="hg38",
        choices=["hg19", "hg38"],
        help="Genome build to use (hg19/hg38).",
    )
    parser.add_argument(
        "--estimate",
        default=None,
        help="Reference directory used to estimate missing variants.",
    )
    parser.add_argument(
        "--proxy", default=None, help="Path to proxy substitution table."
    )
    parser.add_argument(
        "--ntasks", type=int, default=1, help="Number of tasks to use (default=1)."
    )
    parser.add_argument(
        "--batch-size",
        type=int,
        default=5000,
        help="Number of variants per batch (default=5000).",
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
        help="Download or locate the PRS SQL database (variants.db) and metadata JSON (prs_meta.json) and exit.",
    )

    args = parser.parse_args()
    logging.info("Parsed command line arguments.")

    if args.getsql:
        try:
            db_path = get_dm_sql()
            meta_path = get_dm_meta()
            print(db_path)
            print(meta_path)
            logging.info(f"SQL database available at: {db_path}")
            logging.info(f"Metadata JSON available at: {meta_path}")
        except Exception:
            logging.error(
                "Failed to locate or download SQL database and/or metadata JSON.",
                exc_info=True,
            )
            raise
        return

    if not args.vcf:
        parser.error("--vcf is required unless --getsql is used.")
    if not args.scores:
        parser.error("--scores is required unless --getsql is used.")

    try:
        result = gen_dm(
            vcf=args.vcf,
            col=args.col,
            scores=args.scores,
            build=args.build,
            estimate=args.estimate,
            ntasks=args.ntasks,
            batch_size=args.batch_size,
            full=args.full,
            proxy=args.proxy,
        )
        logging.info("Generated PRS scores successfully.")

        result.to_csv(args.output, index=True)
        logging.info(f"Results saved to {args.output}")

    except Exception as e:
        logging.error(f"Error during PRS generation: {e}", exc_info=True)
        raise


if __name__ == "__main__":
    main()
