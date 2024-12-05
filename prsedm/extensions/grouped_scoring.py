"""Contains grouped scoring functions."""
#region standard imports
import logging
#endregion

#region third party imports
import pandas as pd
#endregion

#region local imports
from ..core.score_bcf import score_bcf
#endregion

#region Grouped Scoring
def score_grouped(bed, config):
	"""Generate grouped scores by cluster."""
	def _process_group(group, group_df, config):
		"""Process a single group and return its score."""
		logging.info(f"Generating scores for cluster {group}...")
		group_score = score_bcf(bed=group_df, **config.__dict__)
		group_score = group_score.assign(**{group: group_score.sum(axis=1)})
		return group_score[[group]]
	grouped_scores = [
		_process_group(group, group_df, config)
		for group, group_df in bed.groupby('group')
	]
	return pd.concat(grouped_scores, axis=1)
#endregion
