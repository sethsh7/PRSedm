"""Contains HLA interaction PRS functions"""
import pandas as pd
from ..core.score_bcf import score_bcf
from .grouped_scoring import score_grouped


def generate_dosage_table(dq, config):
    """Generate and rename the dosage table."""
    # Select the correct position column based on the build from the config
    position = f"position_{config.build}"
    if position not in dq.columns:
        raise ValueError(
            f"Position column '{position}' not found in the dosage table.")
    dq['beta'] = 1
    dosages = score_bcf(bed=dq, **config.__dict__)
    dq['contig_id'] = dq['contig_id'].apply(lambda x: x.lstrip('chr'))
    dosages.columns = [col.lstrip('chr') if isinstance(
        col, str) else col for col in dosages.columns]
    reverse_map = {
        f"{c}:{p}": t
        for c, p, t in zip(dq['contig_id'], dq[position], dq['tag'])
    }
    return dosages.rename(columns=reverse_map)


def convert_to_categorical(dosages, rank):
    """Convert dosage values to categorical HLA types."""
    rank_dict = rank.set_index(rank.columns[0]).to_dict()[rank.columns[1]]
    calls = [
        (i, str(col)) for i, row in dosages.iterrows()
        for col in dosages.columns for _ in range(round(row[col])) if row[col] > 0
    ]
    cat = pd.DataFrame(index=dosages.index, columns=['a1', 'a2'], data='X')
    for i, group in pd.DataFrame(
            calls, columns=['index', 'call']).groupby('index'):
        sorted_calls = sorted(
            group['call'], key=lambda x: rank_dict.get(x, float('inf'))
        )
        if sorted_calls:  # Only if there are sorted calls
            cat.at[i, 'a1'] = sorted_calls[0]  # Assign first call
            if len(sorted_calls) > 1:
                # Assign second call if available
                cat.at[i, 'a2'] = sorted_calls[1]
    cat['HLA_type'] = cat['a1'] + '/' + cat['a2']
    return cat


def calculate_hla_scores(cat, int_df):
    """Calculate HLA scores with fallback handling for missing matches."""
    # Build dictionaries for exact and fallback beta values
    beta_dict = pd.concat([
        int_df.query('a1 != "NULL"').set_index(['a1', 'a2'])['beta'],
        int_df.query('a1 != "NULL"').set_index(['a2', 'a1'])['beta']
    ]).to_dict()
    null_beta = int_df.query('a1 == "NULL"').set_index('a2')['beta'].to_dict()

    # Compute exact scores with fallback for missing matches
    exact_scores = cat.set_index(['a1', 'a2']).index.map(beta_dict).fillna(0.0)
    fallback_scores = (
        cat['a1'].map(null_beta).fillna(0.0) +
        cat['a2'].map(null_beta).fillna(0.0)
    )

    # Create the final score series and update the DataFrame
    cat = cat.assign(
        HLA_DRDQ=exact_scores.where(
            exact_scores != 0,
            fallback_scores))
    return cat.drop(columns=['a1', 'a2'])


def score_int_hla(score, dq, int_df, rank, config):
    """Generate and combine HLA and grouped scores."""
    dosages = generate_dosage_table(dq, config)
    cat = convert_to_categorical(dosages, rank)
    hla_score = calculate_hla_scores(cat, int_df)
    grouped_score = score_grouped(score, config)
    return pd.concat([hla_score, grouped_score], axis=1)
