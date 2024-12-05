"""Contains variant processing/pysam functions."""
#region third party imports
import pandas as pd
#endregion

#region Variant Processing Functions
def fetch_variant(var_obj, r, samples):
    """Fetch variant by allele frequency."""
    contig=r['contig_id']; position=r['position']
    records = list(var_obj.fetch(contig=f"{contig}", start=position - 1, stop=position))
    best_variant, max_af = None, -1
    for rec in records:
        if rec:
            df = geno_to_df(rec, samples)
            if int(df['POS'].iloc[0]) == position:
                af = min(get_af_rec(df), 1 - get_af_rec(df))
                if af > max_af:
                    best_variant, max_af = rec, af
                    if max_af == 0.5:  # Perfect balance, no need to continue
                        break
    return best_variant

def geno_to_df(row, samples):
    """Convert a genotype row into a DataFrame."""
    columns = ['CHROM', 'POS', 'ID', 'REF', 'ALT', 
               'QUAL', 'FILTER', 'INFO', 'FORMAT'] + samples
    row_data = str(row).split()  # Split row data by whitespace
    return pd.DataFrame([row_data], columns=columns)

def get_af_rec(df):
	"""
	Extract the allele frequency (AF) from a BCF/VCF INFO column.
	If AF does not exist, calculate it manually from GT if available.
	If neither AF nor GT exists, raise an error.
	"""
	info = df['INFO'].iloc[0]
	info_fields = info.replace(';', ',').split(',')
	af_str = next((s for s in info_fields if s.startswith("AF=")), None)
	if af_str:
		return float(af_str.split("=")[1])
	if 'FORMAT' in df.columns and any(col.startswith('GT') for col in df.columns):
		allele_counts = [0, 0]
		for col in df.columns[9:]:
			gt = df[col].iloc[0].split(':')[0]
			if gt and gt != "./.":
				alleles = gt.replace('/', '|').split('|')
				for allele in alleles:
					if allele.isdigit():
						allele_counts[int(allele)] += 1
		total_alleles = sum(allele_counts)
		if total_alleles > 0:
			return allele_counts[1] / total_alleles
		else:
			raise ValueError("Error: Unable to calculate AF; no valid allele data found.")
	raise ValueError("Error: No AF in INFO column and no GT in VCF to calculate AF.")
#endregion