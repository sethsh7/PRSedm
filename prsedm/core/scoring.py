"""Contains vcf scoring functions."""

import os
import re
import numpy as np
import pysam


def score_geno(record, variant_row, samples, mode):
    effect_allele = variant_row["effect_allele"]
    beta = variant_row["beta"]

    if isinstance(record, dict):
        if record.get("type") != "synthetic_ref":
            raise ValueError(f"Unexpected synthetic site type: {record.get('type')}")

        ref = record["ref"]

        if mode == "GT":
            dosage = score_multiallelic_gt(record, samples)
        elif mode == "GP":
            dosage = score_multiallelic_gp(record, samples)
        else:
            raise ValueError(f"Invalid mode '{mode}'")

    else:
        ref = record.ref

        if mode == "GT":
            samples_data = record.samples
            gt = np.array([samples_data[s]["GT"] for s in samples], dtype=object)
            dosage = np.full(len(samples), np.nan, dtype=float)

            for i, call in enumerate(gt):
                if call is None or len(call) < 2 or call[0] is None or call[1] is None:
                    continue
                if call[0] == 0 and call[1] == 0:
                    dosage[i] = 2.0
                elif (call[0] == 0 and call[1] == 1) or (call[0] == 1 and call[1] == 0):
                    dosage[i] = 1.0
                elif call[0] == 1 and call[1] == 1:
                    dosage[i] = 0.0

        elif mode == "GP":
            samples_data = record.samples
            gp = np.array([samples_data[s]["GP"] for s in samples], dtype=object)
            dosage = np.full(len(samples), np.nan, dtype=float)

            for i, probs in enumerate(gp):
                if probs is None or len(probs) != 3 or any(p is None for p in probs):
                    continue
                dosage[i] = 2 * probs[0] + probs[1]

        else:
            raise ValueError(f"Invalid mode '{mode}'")

    scores = beta * (dosage if effect_allele == ref else 2 - dosage)
    scores[np.isnan(scores)] = 0
    return scores


def score_multiallelic_gt(site, samples):
    records = site["records"]
    dosage = np.full(len(samples), np.nan, dtype=float)

    for i, sample in enumerate(samples):
        nonref_count = 0
        missing = False

        for rec in records:
            call = rec.samples[sample].get("GT")
            if call is None or len(call) < 2 or call[0] is None or call[1] is None:
                missing = True
                break

            nonref_count += (call[0] != 0) + (call[1] != 0)

        if missing:
            continue

        if nonref_count > 2:
            nonref_count = 2

        dosage[i] = 2.0 - nonref_count

    return dosage


def score_multiallelic_gp(site, samples):
    records = site["records"]
    dosage = np.full(len(samples), np.nan, dtype=float)

    for i, sample in enumerate(samples):
        total_alt_af = 0.0
        missing = False

        for rec in records:
            probs = rec.samples[sample].get("GP")
            if probs is None or len(probs) != 3 or any(p is None for p in probs):
                missing = True
                break

            alt_dosage = probs[1] + 2.0 * probs[2]
            total_alt_af += alt_dosage

        if missing:
            continue

        if total_alt_af > 2.0:
            total_alt_af = 2.0

        dosage[i] = 2.0 - total_alt_af

    return dosage


def match_contig_name(var_obj, contig):
    contigs = var_obj.header.contigs
    if contig in contigs:
        return contig

    alt = contig.replace("chr", "") if contig.startswith("chr") else f"chr{contig}"
    if alt in contigs:
        return alt

    raise ValueError(f"Contig {contig} not found in reference VCF header.")


def impute_score_ref(samples, r, refbcf):
    contig = r["contig_id"].replace("chr", "")

    if refbcf.endswith(".txt"):
        if not os.path.exists(refbcf):
            raise FileNotFoundError(f"Mapping file not found: {refbcf}")

        vcf_paths = []
        with open(refbcf) as handle:
            for line in handle:
                parts = line.strip().split()
                if len(parts) < 2:
                    continue
                file_path, chrom = parts[0], parts[1].replace("chr", "")
                if chrom == contig:
                    vcf_paths.append(file_path)

        if not vcf_paths:
            raise ValueError(f"No file for {r['contig_id']} in mapping.")
        if len(vcf_paths) > 1:
            raise ValueError(f"Multiple files found for {r['contig_id']} in mapping.")

        vcf_path = vcf_paths[0]

    else:
        if not re.match(r".*\.vcf(\..*)?$|.*\.bcf(\..*)?$", refbcf, re.IGNORECASE):
            raise ValueError(
                "Invalid file format. Provide a .vcf/.bcf file or valid mapping."
            )

        vcf_path = refbcf

    if not os.path.exists(vcf_path):
        raise FileNotFoundError(f"VCF file not found: {vcf_path}")

    with pysam.VariantFile(vcf_path, "r") as var_obj:
        fetch_contig = match_contig_name(var_obj, r["contig_id"])

        records = list(
            var_obj.fetch(
                contig=fetch_contig, start=r["position"] - 1, stop=r["position"]
            )
        )

        site = resolve_site(records, r)

        if site is None:
            raise ValueError(f"SNP {r['contig_id']}:{r['position']} not found.")

        if isinstance(site, dict):
            ref, af = impute_multiallelic(site)

        else:
            var = site
            ref = var.ref
            if not var.alts:
                raise ValueError(
                    f"No ALT allele found for {r['contig_id']}:{r['position']}."
                )
            if len(var.alts) != 1:
                raise ValueError(
                    f"Expected biallelic site at {r['contig_id']}:{r['position']}, found {len(var.alts)} ALT alleles."
                )
            alt = var.alts[0]
            af_info = var.info.get("AF")
            if af_info is None:
                raise ValueError(
                    f"No AF field found for {r['contig_id']}:{r['position']}."
                )
            if isinstance(af_info, (tuple, list)):
                if len(af_info) != 1:
                    raise ValueError(
                        f"Expected one AF value at {r['contig_id']}:{r['position']}, got {len(af_info)}."
                    )
                af = float(af_info[0])
            else:
                af = float(af_info)
            if r["effect_allele"] not in {ref, alt}:
                raise ValueError(
                    f"Alleles don't match at {r['contig_id']}:{r['position']}: effect_allele={r['effect_allele']}, REF={ref}, ALT={alt}"
                )

    if r["effect_allele"] == ref:
        af = 1 - af

    imputed_score = r["beta"] * 2 * af
    return np.full(len(samples), imputed_score, dtype=float)


def impute_multiallelic(site):
    total_alt_af = 0.0

    for rec in site["records"]:
        af_info = rec.info.get("AF")
        if af_info is None:
            continue
        if not isinstance(af_info, (tuple, list)):
            af_info = (af_info,)
        for af_val in af_info:
            if af_val is not None:
                total_alt_af += float(af_val)

    total_alt_af = min(max(total_alt_af, 0.0), 1.0)
    return site["ref"], total_alt_af


def resolve_site(records, variant_row):
    position = variant_row["position"]
    effect_allele = variant_row["effect_allele"]

    best_alt_rec = None
    best_alt_maf = -1.0
    ref_records = []
    shared_ref = None

    for rec in records:
        if rec.pos != position:
            continue
        alts = rec.alts
        if not alts:
            continue
        ref = rec.ref
        af_info = rec.info.get("AF")
        if af_info is None:
            continue
        if not isinstance(af_info, (tuple, list)):
            af_info = (af_info,)

        for i, alt in enumerate(alts):
            if alt != effect_allele:
                continue
            if i >= len(af_info) or af_info[i] is None:
                continue
            af = float(af_info[i])
            maf = af if af <= 0.5 else 1.0 - af
            if maf > best_alt_maf:
                best_alt_maf = maf
                best_alt_rec = rec
                if maf == 0.5:
                    return rec

        if ref == effect_allele:
            ref_records.append(rec)
            shared_ref = ref

    if best_alt_rec is not None:
        return best_alt_rec
    if len(ref_records) > 1:
        ret = {"type": "synthetic_ref", "records": ref_records, "ref": shared_ref}
        return ret
    elif len(ref_records) == 1:
        return ref_records[0]
    return None


def get_info_float(info, *keys):
    for key in keys:
        val = info.get(key)
        if val is None:
            continue
        if isinstance(val, (tuple, list)):
            for x in val:
                if x is not None:
                    return float(x)
        else:
            return float(val)
    return None


def summarize_site_for_log(site):
    refs, alts_list, afs, r2s = [], [], [], []
    if isinstance(site, dict) and site.get("type") == "synthetic_ref":
        for rec in site["records"]:
            ref = rec.ref
            rec_alts = rec.alts or ["NON_REF"]
            af_info = rec.info.get("AF")
            if af_info is None:
                af_vals = [None] * len(rec_alts)
            elif not isinstance(af_info, (tuple, list)):
                af_vals = [float(af_info)]
            else:
                af_vals = [float(x) if x is not None else None for x in af_info]
            r2_info = None
            for key in ("R2", "DR2", "RSQ"):
                if key in rec.info:
                    r2_info = rec.info[key]
                    break
            for i, alt in enumerate(rec_alts):
                refs.append(ref)
                alts_list.append(alt)
                af = af_vals[i] if i < len(af_vals) else None
                if af is not None:
                    af = min(max(af, 0.0), 1.0)
                afs.append(af)
                r2s.append(r2_info)
    elif isinstance(site, pysam.libcbcf.VariantRecord):
        ref = site.ref
        site_alts = site.alts or ["NON_REF"]
        af_info = site.info.get("AF")
        if af_info is None:
            af_vals = [None] * len(site_alts)
        elif not isinstance(af_info, (tuple, list)):
            af_vals = [float(af_info)]
        else:
            af_vals = [float(x) if x is not None else None for x in af_info]
        r2_info = None
        for key in ("R2", "DR2", "RSQ"):
            if key in site.info:
                r2_info = site.info[key]
                break
        for i, alt in enumerate(site_alts):
            refs.append(ref)
            alts_list.append(alt)
            af = af_vals[i] if i < len(af_vals) else None
            if af is not None:
                af = min(max(af, 0.0), 1.0)
            afs.append(af)
            r2s.append(r2_info)
    return refs, alts_list, afs, r2s
