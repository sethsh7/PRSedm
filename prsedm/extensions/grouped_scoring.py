"""Contains grouped scoring functions."""

import logging
import pandas as pd
from ..core.score_bcf import score_bcf


def score_grouped(bed, config, full=False, variant_log_path=None, prs_name=None):
    """Generate grouped scores by cluster."""
    group_total_cols = []
    total_available, total_estimated, total_missing = 0, 0, 0

    def _process_group(group, group_df, config):
        logging.info(f"Generating scores for cluster {group}...")

        if full:
            # Full mode: keep SNP-level scores + a total per group
            # Ask score_bcf to return all variant columns + 'sum'
            group_score, stats = score_bcf(
                bcf=config.bcf,
                bed=group_df,
                col=config.col,
                build=config.build,
                estimate=config.refbcf,
                ntasks=config.ntasks,
                batch_size=config.batch_size,
                variant_log_path=variant_log_path,
                prs_name=f"{prs_name}_{group}",
                full=True,
            )

            # Rename columns:
            #   - per-variant columns:  group_<variant>
            #   - 'sum' column:         group_total
            new_cols = {}
            for col in group_score.columns:
                if col == "sum":
                    new_cols[col] = f"{group}__total"
                else:
                    new_cols[col] = f"{group}__{col}"
            group_score.rename(columns=new_cols, inplace=True)

            total_col = f"{group}__total"
            group_total_cols.append(total_col)
            return group_score, stats

        else:
            # Compact mode: only keep one total column per group
            # Ask score_bcf to return just the 'sum' column
            group_score, stats = score_bcf(
                bcf=config.bcf,
                bed=group_df,
                col=config.col,
                build=config.build,
                estimate=config.refbcf,
                ntasks=config.ntasks,
                batch_size=config.batch_size,
                variant_log_path=variant_log_path,
                prs_name=group,
                full=False,
            )

            # score_bcf returns a single column named 'sum'
            # Rename it to the group name (or f"{group}__total")
            group_score.rename(columns={"sum": group}, inplace=True)
            group_total_cols.append(group)
            return group_score[[group]], stats

    # Run per-group and concatenate on columns
    grouped_scores = []
    for group, group_df in bed.groupby("group"):
        group_score, stats = _process_group(group, group_df, config)
        grouped_scores.append(group_score)
        total_available += stats[0]
        total_estimated += stats[1]
        total_missing += stats[2]

    grouped_stats = [total_available, total_estimated, total_missing]
    return pd.concat(grouped_scores, axis=1), group_total_cols, grouped_stats
