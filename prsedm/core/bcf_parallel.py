# bcf_parallel.py
import os
import logging
import pysam
from .scoring import (
    score_geno,
    impute_score_ref,
    resolve_site,
    match_contig_name,
    summarize_site_for_log,
)


def process_batch(batch, bcf_files, samples, col, impute, refbcf, log_handle=None):
    """Process a batch of variants and return results using single-variant fetch."""
    batch_results, batch_names, batch_log_rows, genotyped, imputed, missing = (
        [],
        [],
        [],
        0,
        0,
        0,
    )

    try:
        for contig, contig_batch in batch.groupby("contig_id", sort=False):
            contig_genotyped, contig_imputed, contig_missing = 0, 0, 0

            bcf_file = (
                bcf_files.get(contig)
                or bcf_files.get(f"chr{contig}")
                or bcf_files.get(contig.lstrip("chr"))
                or bcf_files.get("all")
            )

            if not bcf_file or not os.path.isfile(bcf_file):
                logging.warning(f"Skipping contig {contig}: BCF file not found")
                for _, variant in contig_batch.iterrows():
                    batch_log_rows.append(
                        (
                            variant["contig_id"],
                            variant["position"],
                            "NA",
                            "NA",
                            variant["effect_allele"],
                            "failed",
                            "NA",
                            "NA",
                        )
                    )
                    missing += 1
                    contig_missing += 1
                continue

            with pysam.VariantFile(bcf_file, "r") as var_obj:
                fetch_contig = match_contig_name(var_obj, contig)

                for _, variant in contig_batch.iterrows():
                    pos = int(variant["position"])
                    try:
                        recs = list(
                            var_obj.fetch(contig=fetch_contig, start=pos - 1, stop=pos)
                        )
                        site = resolve_site(recs, variant)

                        if site is not None:
                            scores = score_geno(site, variant, samples, col)
                            batch_results.append(scores)
                            batch_names.append(f"{variant['contig_id']}:{pos}")
                            genotyped += 1
                            contig_genotyped += 1

                            refs, alts_list, afs, r2s = summarize_site_for_log(site)

                            # Flatten multi-allelic sites into comma-separated strings
                            ref_str = ",".join(str(x) for x in refs)
                            alt_str = ",".join(str(x) for x in alts_list)
                            af_str = ",".join(
                                str(x) if x is not None else "NA" for x in afs
                            )
                            r2_str = ",".join(
                                str(x) if x is not None else "NA" for x in r2s
                            )

                            row = [
                                variant["contig_id"],
                                pos,
                                ref_str,
                                alt_str,
                                variant["effect_allele"],
                                "direct",
                                af_str,
                                r2_str,
                            ]
                            batch_log_rows.append(row)

                        elif impute and refbcf:
                            scores = impute_score_ref(samples, variant, refbcf)
                            batch_results.append(scores)
                            batch_names.append(f"{variant['contig_id']}:{pos}_imputed")
                            imputed += 1
                            contig_imputed += 1

                            batch_log_rows.append(
                                [
                                    variant["contig_id"],
                                    pos,
                                    "NA",
                                    "NA",
                                    variant["effect_allele"],
                                    "imputed",
                                    "NA",
                                    "NA",
                                ]
                            )

                        else:
                            batch_log_rows.append(
                                [
                                    variant["contig_id"],
                                    pos,
                                    "NA",
                                    "NA",
                                    variant["effect_allele"],
                                    "failed",
                                    "NA",
                                    "NA",
                                ]
                            )
                            missing += 1
                            contig_missing += 1

                    except Exception as e:
                        logging.error(
                            f"Exception for variant {variant['contig_id']}:{pos}: {e}",
                            exc_info=True,
                        )
                        batch_log_rows.append(
                            [
                                variant["contig_id"],
                                pos,
                                "NA",
                                "NA",
                                variant["effect_allele"],
                                "failed",
                                "NA",
                                "NA",
                            ]
                        )
                        missing += 1
                        contig_missing += 1

                logging.info(
                    f"Finished batch on {contig}: "
                    f"{len(contig_batch)} variants "
                    f"({contig_genotyped} available, {contig_imputed} estimated, {contig_missing} missing)"
                )

    except Exception as e:
        logging.error(
            f"Exception processing batch for contig {contig}: {e}", exc_info=True
        )

    # Write header if log_handle is provided and file is empty
    if log_handle and batch_log_rows:
        if log_handle.tell() == 0:
            header = [
                "contig_id",
                "position",
                "ref",
                "alt",
                "effect_allele",
                "status",
                "variant_af",
                "variant_r2",
            ]
            log_handle.write("\t".join(header) + "\n")
        for row in batch_log_rows:
            log_handle.write("\t".join(map(str, row)) + "\n")

    return batch_results, batch_names, batch_log_rows, genotyped, imputed, missing
