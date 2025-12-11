"""Contains grouped scoring functions."""
import logging
import pandas as pd
from ..core.score_bcf import score_bcf


def score_grouped(bed, config, full=False):
	"""Generate grouped scores by cluster."""
	group_total_cols = []

	def _process_group(group, group_df, config):
		logging.info(f"Generating scores for cluster {group}...")

		if full:
			# Full mode: keep SNP-level scores + a total per group
			# Ask score_bcf to return all variant columns + 'sum'
			group_score = score_bcf(
				bed=group_df,
				stream=False,
				**config.__dict__
			)

			# Rename columns:
			#   - per-variant columns:  group_<variant>
			#   - 'sum' column:         group_total
			new_cols = {}
			for col in group_score.columns:
				if col == "sum":
					new_cols[col] = f"{group}_total"
				else:
					new_cols[col] = f"{group}_{col}"
			group_score.rename(columns=new_cols, inplace=True)

			total_col = f"{group}_total"
			group_total_cols.append(total_col)
			return group_score

		else:
			# Compact mode: only keep one total column per group
			# Ask score_bcf to return just the 'sum' column
			group_score = score_bcf(
				bed=group_df,
				stream=True,
				**config.__dict__
			)

			# score_bcf returns a single column named 'sum'
			# Rename it to the group name (or you could use f"{group}_total")
			group_score.rename(columns={"sum": group}, inplace=True)
			group_total_cols.append(group)
			return group_score[[group]]

	# Run per-group and concatenate on columns
	grouped_scores = [
		_process_group(group, group_df, config)
		for group, group_df in bed.groupby("group")
	]

	return pd.concat(grouped_scores, axis=1), group_total_cols