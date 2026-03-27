import re
import pandas as pd
from functools import lru_cache


# Tier 1 — HAMLET canonical leukemia gene panel
# Variants here are expected/known AML hits; novel ones are still relevant
# but the real discovery interest is OUTSIDE this list
HAMLET_GENES = {
    "ASXL1", "AXL", "BAALC", "BCOR", "BRAF", "CALR", "CBL", "CEBPA",
    "CSF3R", "CUX1", "DDX41", "DNMT3A", "ETV6", "EZH2", "FLT3", "GATA2",
    "IDH1", "IDH2", "JAK2", "KIT", "KMT2A", "KRAS", "MPL", "MYC", "NF1",
    "NPM1", "NRAS", "PHF6", "PTPN11", "RUNX1", "SETBP1", "SF3B1",
    "SMC1A", "SRSF2", "STAG2", "TET2", "TP53", "U2AF1", "WT1", "ZRSR2",
}

# Alias kept for backward compatibility across pages
AML_GENES = HAMLET_GENES


CSQ_FIELDS = [
    "Allele","Consequence","IMPACT","SYMBOL","Gene","Feature_type","Feature",
    "BIOTYPE","EXON","INTRON","HGVSc","HGVSp","cDNA_position","CDS_position",
    "Protein_position","Amino_acids","Codons","Existing_variation","DISTANCE",
    "STRAND","FLAGS","VARIANT_CLASS","SYMBOL_SOURCE","HGNC_ID","CANONICAL",
    "MANE","MANE_SELECT","MANE_PLUS_CLINICAL","TSL","APPRIS","CCDS","ENSP",
    "SWISSPROT","TREMBL","UNIPARC","UNIPROT_ISOFORM","GENE_PHENO","SIFT",
    "PolyPhen","DOMAINS","miRNA","HGVS_OFFSET","AF","AFR_AF","AMR_AF",
    "EAS_AF","EUR_AF","SAS_AF","gnomADe_AF","gnomADe_AFR_AF","gnomADe_AMR_AF",
    "gnomADe_ASJ_AF","gnomADe_EAS_AF","gnomADe_FIN_AF","gnomADe_MID_AF",
    "gnomADe_NFE_AF","gnomADe_REMAINING_AF","gnomADe_SAS_AF","gnomADg_AF",
    "gnomADg_AFR_AF","gnomADg_AMI_AF","gnomADg_AMR_AF","gnomADg_ASJ_AF",
    "gnomADg_EAS_AF","gnomADg_FIN_AF","gnomADg_MID_AF","gnomADg_NFE_AF",
    "gnomADg_REMAINING_AF","gnomADg_SAS_AF","MAX_AF","MAX_AF_POPS","FREQS",
    "CLIN_SIG","SOMATIC","PHENO","PUBMED","MOTIF_NAME","MOTIF_POS",
    "HIGH_INF_POS","MOTIF_SCORE_CHANGE","TRANSCRIPTION_FACTORS"
]

def _safe_float(val):
    try:
        return float(val) if val and val not in ("", ".") else None
    except:
        return None

def _parse_info(info_str):
    result = {}
    for part in info_str.split(";"):
        if "=" in part:
            k, v = part.split("=", 1)
            result[k] = v
        else:
            result[part] = True
    return result

def _best_csq(csq_str):
    """Pick the canonical transcript entry with highest impact."""
    impact_order = {"HIGH": 0, "MODERATE": 1, "LOW": 2, "MODIFIER": 3}
    entries = csq_str.split(",")
    best = None
    best_rank = 99
    for entry in entries:
        fields = entry.split("|")
        if len(fields) < 25:
            continue
        canonical = fields[24]
        impact = fields[2]
        rank = impact_order.get(impact, 99)
        if canonical == "YES" and rank < best_rank:
            best_rank = rank
            best = fields
    if best is None:
        # fallback: just take first entry
        fields = entries[0].split("|")
        best = fields
    return best

def load_vcf(vcf_path):
    rows = []
    samples = []

    with open(vcf_path) as f:
        for line in f:
            if line.startswith("##"):
                continue
            if line.startswith("#CHROM"):
                cols = line.strip().split("\t")
                samples = cols[9:]
                continue

            cols = line.strip().split("\t")
            if len(cols) < 8:
                continue

            chrom, pos, vid, ref, alt, qual, filt, info_str = cols[:8]
            fmt = cols[8] if len(cols) > 8 else ""
            sample_data = cols[9:] if len(cols) > 9 else []

            info = _parse_info(info_str)
            csq_raw = info.get("CSQ", "")

            # Parse best CSQ entry
            csq = {}
            if csq_raw:
                fields = _best_csq(csq_raw)
                for i, fname in enumerate(CSQ_FIELDS):
                    csq[fname] = fields[i] if i < len(fields) else ""

            # Variant type
            ref_len = len(ref)
            alt_alleles = alt.split(",")
            if all(len(a) == 1 and ref_len == 1 for a in alt_alleles):
                var_type = "SNP"
            else:
                var_type = "INDEL"

            # Novel = no existing_variation in canonical transcript
            existing = csq.get("Existing_variation", "")
            is_novel = (existing == "" or existing == ".")

            row = {
                "CHROM": chrom,
                "POS": int(pos),
                "ID": vid,
                "REF": ref,
                "ALT": alt,
                "QUAL": _safe_float(qual),
                "FILTER": filt,
                "VAR_TYPE": var_type,
                "IS_NOVEL": is_novel,
                # INFO fields
                "AC": info.get("AC", ""),
                "AF": _safe_float(info.get("AF", "").split(",")[0] if "," in info.get("AF","") else info.get("AF","")),
                "AN": _safe_float(info.get("AN", "")),
                "DP": _safe_float(info.get("DP", "")),
                "QD": _safe_float(info.get("QD", "")),
                "FS": _safe_float(info.get("FS", "")),
                "MQ": _safe_float(info.get("MQ", "")),
                "SOR": _safe_float(info.get("SOR", "")),
                "VQSLOD": _safe_float(info.get("VQSLOD", "")),
                "InbreedingCoeff": _safe_float(info.get("InbreedingCoeff", "")),
                "ExcessHet": _safe_float(info.get("ExcessHet", "")),
                # VEP fields
                "SYMBOL": csq.get("SYMBOL", ""),
                "CONSEQUENCE": csq.get("Consequence", ""),
                "IMPACT": csq.get("IMPACT", ""),
                "HGVSc": csq.get("HGVSc", ""),
                "HGVSp": csq.get("HGVSp", ""),
                "CANONICAL": csq.get("CANONICAL", ""),
                "BIOTYPE": csq.get("BIOTYPE", ""),
                "EXISTING_VARIATION": existing,
                "VARIANT_CLASS": csq.get("VARIANT_CLASS", ""),
                "SIFT": csq.get("SIFT", ""),
                "PolyPhen": csq.get("PolyPhen", ""),
                "gnomADe_AF": _safe_float(csq.get("gnomADe_AF", "")),
                "gnomADg_AF": _safe_float(csq.get("gnomADg_AF", "")),
                "MAX_AF": _safe_float(csq.get("MAX_AF", "")),
                "CLIN_SIG": csq.get("CLIN_SIG", ""),
                "IS_AML_GENE": csq.get("SYMBOL", "") in HAMLET_GENES,
                "GENE_TIER": (
                    "HAMLET canonical" if csq.get("SYMBOL", "") in HAMLET_GENES
                    else ("Non-canonical" if csq.get("SYMBOL", "") else "Intergenic")
                ),
                # Sample FORMAT fields
                "ROW_FORMAT": fmt,
                "ROW_SAMPLES": sample_data,
                "ROW_SAMPLE_NAMES": samples,
            }
            rows.append(row)

    df = pd.DataFrame(rows)

    # Numeric coerce
    for col in ["QUAL", "AF", "AN", "DP", "QD", "FS", "MQ", "SOR",
                "VQSLOD", "InbreedingCoeff", "ExcessHet", "gnomADe_AF",
                "gnomADg_AF", "MAX_AF"]:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce")

    return df, samples


def get_per_sample_stats(df, samples):
    """Compute per-sample GT, DP, GQ stats."""
    rows = []
    for df_row in df.itertuples():
        fmt = df_row.ROW_FORMAT
        sample_data = df_row.ROW_SAMPLES
        fmt_keys = fmt.split(":")
        gt_idx = fmt_keys.index("GT") if "GT" in fmt_keys else None
        dp_idx = fmt_keys.index("DP") if "DP" in fmt_keys else None
        gq_idx = fmt_keys.index("GQ") if "GQ" in fmt_keys else None

        for i, sname in enumerate(samples):
            if i >= len(sample_data):
                continue
            vals = sample_data[i].split(":")
            gt = vals[gt_idx] if gt_idx is not None and gt_idx < len(vals) else "./."
            dp = _safe_float(vals[dp_idx]) if dp_idx is not None and dp_idx < len(vals) else None
            gq = _safe_float(vals[gq_idx]) if gq_idx is not None and gq_idx < len(vals) else None

            # Classify genotype
            gt_clean = gt.replace("|", "/")
            alleles = gt_clean.split("/")
            if "." in alleles:
                gt_class = "missing"
            elif all(a == "0" for a in alleles):
                gt_class = "hom_ref"
            elif len(set(alleles)) == 1:
                gt_class = "hom_alt"
            else:
                gt_class = "het"

            rows.append({
                "sample": sname,
                "CHROM": df_row.CHROM,
                "POS": df_row.POS,
                "REF": df_row.REF,
                "ALT": df_row.ALT,
                "SYMBOL": df_row.SYMBOL,
                "IMPACT": df_row.IMPACT,
                "CONSEQUENCE": df_row.CONSEQUENCE,
                "IS_NOVEL": df_row.IS_NOVEL,
                "IS_AML_GENE": df_row.IS_AML_GENE,
                "FILTER": df_row.FILTER,
                "GT": gt,
                "GT_CLASS": gt_class,
                "DP": dp,
                "GQ": gq,
                "VQSLOD": df_row.VQSLOD,
                "gnomADe_AF": df_row.gnomADe_AF,
            })

    return pd.DataFrame(rows)
