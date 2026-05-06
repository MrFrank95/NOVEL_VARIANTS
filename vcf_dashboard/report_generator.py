#!/usr/bin/env python3
"""
AML Per-Sample Variant Report Generator
========================================
Reads a GATK joint-germline VCF annotated with VEP and generates one PDF
report per sample.  Each report contains:
  - Cover page
  - QC summary (cohort-level and per-sample statistics)
  - Section A: HAMLET canonical variant table (all carried variants)
  - Section B: Non-canonical novel candidates (HIGH/MODERATE, no population DB)
  - Glossary

Usage:
    python generate_reports.py <vcf_path> [--out_dir <dir>] [--samples L40 L41 ...]

Dependencies:
    pip install reportlab
"""

import argparse, os, sys, re
from datetime import datetime
from collections import defaultdict

# ── ReportLab imports ──────────────────────────────────────────────────────────
from reportlab.lib.pagesizes import A4
from reportlab.lib import colors
from reportlab.lib.units import cm, mm
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.lib.enums import TA_CENTER, TA_LEFT, TA_RIGHT, TA_JUSTIFY
from reportlab.platypus import (
    SimpleDocTemplate, Paragraph, Spacer, Table, TableStyle,
    PageBreak, HRFlowable, KeepTogether
)
from reportlab.platypus.flowables import Flowable

# ═══════════════════════════════════════════════════════════════════════════════
# CONSTANTS
# ═══════════════════════════════════════════════════════════════════════════════

HAMLET_GENES = {
    "ASXL1","AXL","BAALC","BCOR","BRAF","CALR","CBL","CEBPA",
    "CSF3R","CUX1","DDX41","DNMT3A","ETV6","EZH2","FLT3","GATA2",
    "IDH1","IDH2","JAK2","KIT","KMT2A","KRAS","MPL","MYC","NF1",
    "NPM1","NRAS","PHF6","PTPN11","RUNX1","SETBP1","SF3B1",
    "SMC1A","SRSF2","STAG2","TET2","TP53","U2AF1","WT1","ZRSR2",
}

# GATK BEST-PRACTICE HARD-FILTER THRESHOLDS (VQSR pipeline as reference)
FILTER_THRESHOLDS = {
    "VQSLOD_min":  0.0,      # variants with VQSLOD < 0 are low-confidence
    "QD_min":      2.0,      # Quality by Depth
    "FS_max_snp":  60.0,     # FisherStrand (SNP)
    "FS_max_indel":200.0,    # FisherStrand (INDEL)
    "SOR_max_snp": 3.0,      # Strand Odds Ratio (SNP)
    "MQ_min":      40.0,     # RMS Mapping Quality
    "MQRankSum_min": -12.5,  # MQ rank sum
    "ReadPosRankSum_min": -8.0,  # Read position rank sum
    "gnomAD_AF_max": 0.01,   # Max population AF to keep
    "DP_min":       10,      # Minimum depth per sample
    "GQ_min":       20,      # Minimum genotype quality
}

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

GLOSSARY = [
    ("CHROM / POS / REF / ALT",
     "Genomic coordinates on GRCh38: chromosome, base position, reference allele, and alternate allele."),
    ("VAF (%)",
     "Variant Allele Frequency for THIS sample. Calculated from the FORMAT AD field as "
     "(sum of alt allele reads / total reads at site) × 100. High VAF (~50%) suggests "
     "heterozygous germline; ~100% suggests homozygous."),
    ("Cohort AF",
     "Allele frequency across all 31 samples in this cohort (INFO/AF field). "
     "Allows you to see how many samples carry the same variant."),
    ("GT",
     "Genotype call: 0/0 = homozygous reference, 0/1 = heterozygous, "
     "1/1 = homozygous alternate, ./. = missing call."),
    ("DP",
     "Total read depth at the variant site in this sample."),
    ("GQ",
     "Genotype Quality: Phred-scaled confidence that the genotype call is correct. "
     "GQ >= 20 is generally considered reliable."),
    ("HGVSc",
     "HGVS coding DNA notation (e.g. c.123A>T). Describes the base change "
     "relative to the canonical transcript."),
    ("HGVSp",
     "HGVS protein notation (e.g. p.Lys41Asn). Describes the amino-acid consequence. "
     "Empty for synonymous or non-coding variants."),
    ("IMPACT",
     "VEP severity tier: HIGH = likely protein-disrupting (stop-gain, frameshift, "
     "splice-site); MODERATE = missense, in-frame indel; LOW = synonymous; "
     "MODIFIER = non-coding / intergenic."),
    ("Consequence",
     "The specific molecular effect predicted by VEP (e.g. missense_variant, "
     "splice_acceptor_variant, frameshift_variant)."),
    ("ClinVar",
     "Clinical significance from ClinVar: pathogenic, likely_pathogenic, "
     "uncertain_significance, benign, etc. Empty = not in ClinVar."),
    ("SIFT",
     "Predicts whether an amino-acid substitution affects protein function. "
     "deleterious (score < 0.05) / tolerated."),
    ("PolyPhen",
     "Predicts the structural/functional impact of a missense change: "
     "probably_damaging, possibly_damaging, or benign."),
    ("gnomADe AF",
     "Allele frequency in the gnomAD exome cohort (>125,000 exomes). "
     "Absent = not observed in the general population."),
    ("MAX_AF",
     "Maximum allele frequency across all gnomAD sub-populations. "
     "Useful for catching variants common in a specific ancestry."),
    ("VQSLOD",
     "Variant Quality Score Log-Odds from GATK VQSR. Higher = more likely a "
     "true variant. Negative values indicate low-confidence calls."),
    ("QD",
     "Quality by Depth: variant quality score divided by allele depth. "
     "QD < 2 may indicate a false positive."),
    ("FS",
     "FisherStrand: Phred-scaled p-value for strand bias. "
     "High values (>60 for SNPs, >200 for INDELs) suggest strand-specific artefacts."),
    ("SOR",
     "StrandOddsRatio: alternative strand bias metric. SOR > 3 for SNPs is flagged."),
    ("MQ",
     "RMS Mapping Quality of reads supporting the variant. MQ < 40 may indicate "
     "mismapped reads."),
    ("MQRankSum",
     "Mann-Whitney U test comparing mapping quality of reads with REF vs ALT alleles. "
     "Values < -12.5 are flagged."),
    ("ReadPosRankSum",
     "Mann-Whitney U test for the position of the variant within reads. "
     "Values < -8.0 suggest position bias."),
    ("Novel",
     "True if the variant has no entry in Existing_variation (dbSNP rsID / ClinVar "
     "accession) for the canonical VEP transcript. Assessed against the canonical "
     "transcript only."),
    ("HAMLET canonical gene",
     "One of the 41 genes in the HAMLET AML panel: ASXL1, AXL, BAALC, BCOR, BRAF, "
     "CALR, CBL, CEBPA, CSF3R, CUX1, DDX41, DNMT3A, ETV6, EZH2, FLT3, GATA2, IDH1, "
     "IDH2, JAK2, KIT, KMT2A, KRAS, MPL, MYC, NF1, NPM1, NRAS, PHF6, PTPN11, RUNX1, "
     "SETBP1, SF3B1, SMC1A, SRSF2, STAG2, TET2, TP53, U2AF1, WT1, ZRSR2."),
    ("Non-canonical novel candidate",
     "A novel HIGH/MODERATE impact variant in a gene OUTSIDE the HAMLET panel. "
     "These are the primary scientific discovery targets."),
]

# ═══════════════════════════════════════════════════════════════════════════════
# VCF PARSER
# ═══════════════════════════════════════════════════════════════════════════════

def _safe_float(val):
    try:
        return float(val) if val and val not in ("",".","./.") else None
    except Exception:
        return None

def _parse_info(info_str):
    d = {}
    for part in info_str.split(";"):
        if "=" in part:
            k, v = part.split("=", 1)
            d[k] = v
        else:
            d[part] = True
    return d

def _best_csq(csq_str):
    """Return the highest-impact canonical-transcript CSQ entry as a dict."""
    impact_rank = {"HIGH": 0, "MODERATE": 1, "LOW": 2, "MODIFIER": 3}
    best = None
    best_rank = 99
    for entry in csq_str.split(","):
        fields = entry.split("|")
        if len(fields) < 25:
            continue
        canonical = fields[24]
        impact    = fields[2]
        rank      = impact_rank.get(impact, 99)
        if canonical == "YES" and rank < best_rank:
            best_rank = rank
            best      = fields
    if best is None:
        best = csq_str.split(",")[0].split("|")
    csq = {}
    for i, fname in enumerate(CSQ_FIELDS):
        csq[fname] = best[i] if i < len(best) else ""
    return csq

def _calc_vaf(gt, ad_str, dp_val):
    """
    Compute VAF (0-100) for a sample from GT and AD fields.
    Returns None if missing / reference.
    """
    gt_clean = gt.replace("|","/")
    alleles  = gt_clean.split("/")
    if "." in alleles:
        return None
    alt_alleles = set([a for a in alleles if a != "0"])
    if not alt_alleles:
        return 0.0   # homozygous ref
    try:
        ad = [int(x) for x in ad_str.split(",")]
        total = sum(ad)
        if total == 0:
            return None

        alt_reads = sum(ad[int(a)] for a in alt_alleles if int(a) < len(ad))
        return round(alt_reads / total * 100, 1)
    except Exception:
        return None


def load_vcf(vcf_path):
    """
    Parse VCF and return:
      - variants: list of dicts with INFO + VEP annotation
      - samples:  ordered list of sample names
      - sample_gts: dict[sample_name] → list of per-variant dicts
                    {GT, DP, GQ, VAF, IS_CARRIER}
    """
    variants = []
    samples  = []

    with open(vcf_path, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            line = line.rstrip("\r\n")

            if line.startswith("##"):
                continue

            if line.startswith("#CHROM"):
                header_cols = line.split("\t")
                samples = [c.split("_")[0] for c in header_cols[9:]]
                continue

            cols = line.split("\t")
            if len(cols) < 8:
                continue

            chrom, pos, vid, ref, alt, qual, filt, info_str = cols[:8]
            fmt       = cols[8] if len(cols) > 8 else ""
            samp_data = cols[9:] if len(cols) > 9 else []

            info = _parse_info(info_str)
            csq  = _best_csq(info.get("CSQ","")) if "CSQ" in info else {}

            # Variant type
            ref_len   = len(ref)
            alt_alleles_all = alt.split(",")
            var_type  = "SNP" if all(len(a)==1 and ref_len==1
                                     for a in alt_alleles_all) else "INDEL"

            existing  = csq.get("Existing_variation","")
            is_novel  = (existing in ("",".",None))
            symbol    = csq.get("SYMBOL","")
            is_hamlet = symbol in HAMLET_GENES

            # Cohort AF (first value if comma-separated)
            af_raw = info.get("AF","")
            if "," in af_raw:
                af_raw = af_raw.split(",")[0]
            cohort_af = _safe_float(af_raw)

            fmt_keys = fmt.split(":")
            gt_idx = fmt_keys.index("GT") if "GT" in fmt_keys else 0
            carriers = 0
            for s in samp_data:
                vals = s.split(":")
                gt = vals[gt_idx] if gt_idx < len(vals) else "./."
                alleles = gt.replace("|", "/").split("/")
                if any(a not in ("0", ".","-") for a in alleles):
                    carriers += 1

            row = {
                "CHROM":    chrom,
                "POS":      int(pos),
                "ID":       vid,
                "REF":      ref,
                "ALT":      alt,
                "QUAL":     _safe_float(qual),
                "FILTER":   filt,
                "VAR_TYPE": var_type,
                "IS_NOVEL": is_novel,
                "IS_HAMLET":is_hamlet,
                # INFO
                "COHORT_AF": f'{carriers}/{len(samples)}',
                "AC":    info.get("AC",""),
                "AN":    _safe_float(info.get("AN","")),
                "DP":    _safe_float(info.get("DP","")),
                "QD":    _safe_float(info.get("QD","")),
                "FS":    _safe_float(info.get("FS","")),
                "MQ":    _safe_float(info.get("MQ","")),
                "SOR":   _safe_float(info.get("SOR","")),
                "VQSLOD":_safe_float(info.get("VQSLOD","")),
                "MQRankSum":        _safe_float(info.get("MQRankSum","")),
                "ReadPosRankSum":   _safe_float(info.get("ReadPosRankSum","")),
                "InbreedingCoeff":  _safe_float(info.get("InbreedingCoeff","")),
                "ExcessHet":        _safe_float(info.get("ExcessHet","")),
                # VEP
                "SYMBOL":      symbol,
                "CONSEQUENCE": csq.get("Consequence",""),
                "IMPACT":      csq.get("IMPACT",""),
                "HGVSc":       csq.get("HGVSc",""),
                "HGVSp":       csq.get("HGVSp",""),
                "BIOTYPE":     csq.get("BIOTYPE",""),
                "EXISTING_VARIATION": existing,
                "SIFT":        csq.get("SIFT",""),
                "PolyPhen":    csq.get("PolyPhen",""),
                "gnomADe_AF":  _safe_float(csq.get("gnomADe_AF","")),
                "MAX_AF":      _safe_float(csq.get("MAX_AF","")),
                "CLIN_SIG":    csq.get("CLIN_SIG",""),
                # Raw format/sample data for per-sample parsing
                "_FORMAT":     fmt,
                "_SAMPLES":    samp_data,
            }
            variants.append(row)

    # Build per-sample genotype table
    sample_gts = {s: [] for s in samples}
    fmt_keys_cache = {}

    for v in variants:
        fmt   = v["_FORMAT"]
        sdata = v["_SAMPLES"]
        if fmt not in fmt_keys_cache:
            fmt_keys_cache[fmt] = fmt.split(":")
        fmt_keys = fmt_keys_cache[fmt]
        gt_idx   = fmt_keys.index("GT") if "GT" in fmt_keys else None
        ad_idx   = fmt_keys.index("AD") if "AD" in fmt_keys else None
        dp_idx   = fmt_keys.index("DP") if "DP" in fmt_keys else None
        gq_idx   = fmt_keys.index("GQ") if "GQ" in fmt_keys else None

        for i, sname in enumerate(samples):
            if i >= len(sdata):
                sample_gts[sname].append(None)
                continue
            vals = sdata[i].split(":")
            gt  = vals[gt_idx]  if gt_idx  is not None and gt_idx  < len(vals) else "./."
            ad  = vals[ad_idx]  if ad_idx  is not None and ad_idx  < len(vals) else ""
            dp  = _safe_float(vals[dp_idx]) if dp_idx is not None and dp_idx < len(vals) else None
            gq  = _safe_float(vals[gq_idx]) if gq_idx is not None and gq_idx < len(vals) else None

            gt_clean = gt.replace("|","/").replace("./.","./.")
            alleles  = gt_clean.split("/")
            is_missing  = "." in alleles
            is_ref      = not is_missing and all(a=="0" for a in alleles)
            is_carrier  = not is_missing and not is_ref

            vaf = _calc_vaf(gt, ad, dp) if is_carrier else None

            sample_gts[sname].append({
                "GT":         gt,
                "DP":         dp,
                "GQ":         gq,
                "VAF":        vaf,
                "IS_CARRIER": is_carrier,
                "IS_MISSING": is_missing,
            })

    return variants, samples, sample_gts

# ═══════════════════════════════════════════════════════════════════════════════
# REPORT STYLES
# ═══════════════════════════════════════════════════════════════════════════════

# Palette
C_NAVY      = colors.HexColor("#1B2A4A")
C_TEAL      = colors.HexColor("#1A6B72")
C_ACCENT    = colors.HexColor("#C0392B")
C_AMBER     = colors.HexColor("#E67E22")
C_GOLD      = colors.HexColor("#F1C40F")
C_LIGHT     = colors.HexColor("#F5F7FA")
C_MID       = colors.HexColor("#DDE3ED")
C_HAMLET_BG = colors.HexColor("#FFFBEA")
C_NOVEL_BG  = colors.HexColor("#FFF0EE")
C_HIGH      = colors.HexColor("#FDDCDB")
C_MODERATE  = colors.HexColor("#FEF0DC")
C_PASS_GREEN= colors.HexColor("#E8F5E9")

def make_styles():
    ss = getSampleStyleSheet()

    cover_title = ParagraphStyle("cover_title",
        fontName="Helvetica-Bold", fontSize=28, textColor=C_NAVY,
        alignment=TA_CENTER, spaceAfter=8)

    cover_sub = ParagraphStyle("cover_sub",
        fontName="Helvetica", fontSize=14, textColor=C_TEAL,
        alignment=TA_CENTER, spaceAfter=6)

    cover_meta = ParagraphStyle("cover_meta",
        fontName="Helvetica", fontSize=10, textColor=colors.grey,
        alignment=TA_CENTER, spaceAfter=4)

    section_title = ParagraphStyle("section_title",
        fontName="Helvetica-Bold", fontSize=14, textColor=C_NAVY,
        spaceBefore=14, spaceAfter=6,
        borderPad=4, borderWidth=0)

    subsection = ParagraphStyle("subsection",
        fontName="Helvetica-Bold", fontSize=11, textColor=C_TEAL,
        spaceBefore=10, spaceAfter=4)

    body = ParagraphStyle("body",
        fontName="Helvetica", fontSize=8.5, textColor=colors.black,
        spaceBefore=3, spaceAfter=3, leading=12)

    caption = ParagraphStyle("caption",
        fontName="Helvetica-Oblique", fontSize=7.5, textColor=colors.grey,
        spaceBefore=2, spaceAfter=4, leading=10)

    cell = ParagraphStyle("cell",
        fontName="Helvetica", fontSize=7.5, textColor=colors.black,
        leading=10, wordWrap="CJK")

    cell_bold = ParagraphStyle("cell_bold",
        fontName="Helvetica-Bold", fontSize=7.5, textColor=colors.black,
        leading=10)

    th = ParagraphStyle("th",
        fontName="Helvetica-Bold", fontSize=7.5, textColor=colors.white,
        leading=10, alignment=TA_CENTER)

    small = ParagraphStyle("small",
        fontName="Helvetica", fontSize=7, textColor=colors.grey, leading=9)

    gloss_term = ParagraphStyle("gloss_term",
        fontName="Helvetica-Bold", fontSize=8, textColor=C_NAVY, leading=11)
    gloss_def = ParagraphStyle("gloss_def",
        fontName="Helvetica", fontSize=8, textColor=colors.black, leading=11,
        leftIndent=12)

    return {
        "cover_title": cover_title,
        "cover_sub":   cover_sub,
        "cover_meta":  cover_meta,
        "section":     section_title,
        "subsection":  subsection,
        "body":        body,
        "caption":     caption,
        "cell":        cell,
        "cell_bold":   cell_bold,
        "th":          th,
        "small":       small,
        "gloss_term":  gloss_term,
        "gloss_def":   gloss_def,
    }

# ═══════════════════════════════════════════════════════════════════════════════
# HEADER / FOOTER
# ═══════════════════════════════════════════════════════════════════════════════

class PageHeader:
    """Stores info used by the canvas callback."""
    def __init__(self, sample_id):
        self.sample_id = sample_id

def make_canvas_callback(sample_id, gen_date):
    def on_page(canvas, doc):
        # w, h = A4
        w = 25*cm
        h = 35*cm
        # Header bar
        canvas.setFillColor(C_NAVY)
        canvas.rect(0, h - 1.4*cm, w, 1.4*cm, fill=1, stroke=0)
        canvas.setFont("Helvetica-Bold", 9)
        canvas.setFillColor(colors.white)
        canvas.drawString(1*cm, h - 0.9*cm, "AML Variant Report")
        canvas.setFont("Helvetica", 9)
        canvas.drawRightString(w - 1*cm, h - 0.9*cm, f"Sample: {sample_id}")

        # Footer
        canvas.setFillColor(C_MID)
        canvas.rect(0, 0, w, 0.8*cm, fill=1, stroke=0)
        canvas.setFont("Helvetica", 7.5)
        canvas.setFillColor(colors.grey)
        canvas.drawString(1*cm, 0.28*cm,
            f"Generated {gen_date}  ·  GATK 4.6 · VEP 115 · GRCh38")
        canvas.drawRightString(w - 1*cm, 0.28*cm,
            f"Page {doc.page}")
        canvas.setStrokeColor(C_TEAL)
        canvas.setLineWidth(1)
        canvas.line(0, 0.8*cm, w, 0.8*cm)
    return on_page

# ═══════════════════════════════════════════════════════════════════════════════
# HELPER: coloured horizontal rule
# ═══════════════════════════════════════════════════════════════════════════════

def hrule(color=C_TEAL, thickness=1):
    return HRFlowable(width="100%", thickness=thickness,
                      color=color, spaceAfter=4, spaceBefore=4)

# ═══════════════════════════════════════════════════════════════════════════════
# HELPER: metric badge table (KPI row)
# ═══════════════════════════════════════════════════════════════════════════════

def kpi_table(metrics, styles):
    """metrics: list of (label, value, bg_color)"""
    n = len(metrics)
    col_w = (A4[0] - 2*cm) / n

    header_row = []
    value_row  = []
    for label, val, bg in metrics:
        header_row.append(Paragraph(label, ParagraphStyle("kpi_lbl",
            fontName="Helvetica", fontSize=7, textColor=colors.grey,
            alignment=TA_CENTER, leading=9)))
        value_row.append(Paragraph(str(val), ParagraphStyle("kpi_val",
            fontName="Helvetica-Bold", fontSize=14, textColor=C_NAVY,
            alignment=TA_CENTER, leading=17)))

    t = Table([header_row, value_row], colWidths=[col_w]*n, rowHeights=[14, 22])

    ts = [
        ("ALIGN",  (0,0), (-1,-1), "CENTER"),
        ("VALIGN", (0,0), (-1,-1), "MIDDLE"),
        ("TOPPADDING",    (0,0), (-1,-1), 4),
        ("BOTTOMPADDING", (0,0), (-1,-1), 4),
        ("LEFTPADDING",   (0,0), (-1,-1), 4),
        ("RIGHTPADDING",  (0,0), (-1,-1), 4),
        ("BOX",    (0,0), (-1,-1), 0.4, C_MID),
    ]
    for i, (_, _, bg) in enumerate(metrics):
        ts.append(("BACKGROUND", (i,0), (i,-1), bg))
    t.setStyle(TableStyle(ts))
    return t

# ═══════════════════════════════════════════════════════════════════════════════
# HELPER: variant table builder
# ═══════════════════════════════════════════════════════════════════════════════

def _fmt(val, decimals=2, pct=False, na="—"):
    if val is None:
        return na
    if pct:
        return f"{val:.1f}%"
    if isinstance(val, float):
        return f"{val:.{decimals}f}"
    return str(val)

def _impact_bg(impact):
    return {
        "HIGH":     C_HIGH,
        "MODERATE": C_MODERATE,
    }.get(impact, colors.white)

def variant_table(rows, sample_gt_rows, styles, include_sample_cols=True):
    """
    Build a ReportLab Table for a list of variant dicts.
    sample_gt_rows: list of per-sample GT dicts aligned with rows.
    """
    col_defs = [
        ("Gene",        2.0*cm),
        ("Location",    3.2*cm),
        ("Consequence", 2.8*cm),
        ("Impact",      1.65*cm),
        ("HGVSp",       3.4*cm),
        ("ClinVar",     1.8*cm),
        ("gnomAD-e AF", 1.5*cm),
        ("Cohort AF",   1.3*cm),
    ]
    if include_sample_cols:
        col_defs += [
            ("VAF", 1.3*cm),
            ("GT",      1.0*cm),
            ("DP",      0.9*cm),
            ("GQ",      0.9*cm),
        ]

    headers    = [h for h, _ in col_defs]
    col_widths = [w for _, w in col_defs]

    th_style = ParagraphStyle("th",
        fontName="Helvetica-Bold", fontSize=7, textColor=colors.white,
        alignment=TA_CENTER, leading=9)
    cell_style = ParagraphStyle("cell",
        fontName="Helvetica", fontSize=7, textColor=colors.black,
        leading=9, wordWrap="CJK")
    cell_bold = ParagraphStyle("cell_b",
        fontName="Helvetica-Bold", fontSize=7, textColor=C_NAVY, leading=9)

    data = [[Paragraph(h, th_style) for h in headers]]
    row_backgrounds = []

    for idx, v in enumerate(rows):
        gt_info = sample_gt_rows[idx] if sample_gt_rows and idx < len(sample_gt_rows) else None
        loc = f"{v['CHROM']}:{v['POS']}\n{v['REF']}>{v['ALT']}"
        gad = _fmt(v.get("gnomADe_AF"), 4) if v.get("gnomADe_AF") is not None else "absent"
        # caf = f"{v['COHORT_AF']:.3f}" if v.get("COHORT_AF") is not None else "—"
        caf = f"{v['COHORT_AF']}" if v.get("COHORT_AF") is not None else "—"
        clin = v.get("CLIN_SIG","") or "—"
        impact = v.get("IMPACT","")
        hgvsp  = v.get("HGVSp","") or "—"

        row_cells = [
            Paragraph(v.get("SYMBOL","—"), cell_bold),
            Paragraph(loc,                  cell_style),
            Paragraph(v.get("CONSEQUENCE","—").replace("_"," "), cell_style),
            Paragraph(impact,               cell_style),
            Paragraph(hgvsp[:40],           cell_style),
            Paragraph(clin[:22],            cell_style),
            Paragraph(gad,                  cell_style),
            Paragraph(caf,                  cell_style),
        ]

        if include_sample_cols:
            if gt_info:
                vaf_str = f"{gt_info['VAF']:.1f}%" if gt_info.get("VAF") is not None else "—"
                dp_str  = str(int(gt_info["DP"])) if gt_info.get("DP") is not None else "—"
                gq_str  = str(int(gt_info["GQ"])) if gt_info.get("GQ") is not None else "—"
                row_cells += [
                    Paragraph(vaf_str,         cell_style),
                    Paragraph(gt_info["GT"],   cell_style),
                    Paragraph(dp_str,          cell_style),
                    Paragraph(gq_str,          cell_style),
                ]
            else:
                row_cells += [Paragraph("—", cell_style)]*4

        data.append(row_cells)
        row_backgrounds.append(_impact_bg(impact))

    t = Table(data, colWidths=col_widths, repeatRows=1)

    ts = [
        ("BACKGROUND",    (0,0), (-1,0),  C_NAVY),
        ("ROWBACKGROUNDS",(0,1), (-1,-1), [C_LIGHT, colors.white]),
        ("ALIGN",         (0,0), (-1,-1), "LEFT"),
        ("VALIGN",        (0,0), (-1,-1), "TOP"),
        ("FONTSIZE",      (0,0), (-1,-1), 7),
        ("TOPPADDING",    (0,0), (-1,-1), 3),
        ("BOTTOMPADDING", (0,0), (-1,-1), 3),
        ("LEFTPADDING",   (0,0), (-1,-1), 3),
        ("RIGHTPADDING",  (0,0), (-1,-1), 3),
        ("GRID",          (0,0), (-1,-1), 0.3, C_MID),
        ("LINEBELOW",     (0,0), (-1,0),  1,   C_TEAL),
    ]
    # Apply impact colours per row
    for ridx, bg in enumerate(row_backgrounds):
        if bg != colors.white:
            ts.append(("BACKGROUND", (0, ridx+1), (-1, ridx+1), bg))

    t.setStyle(TableStyle(ts))
    return t

# ═══════════════════════════════════════════════════════════════════════════════
# QC FLAGS HELPER
# ═══════════════════════════════════════════════════════════════════════════════

def qc_flag(v):
    """Return list of QC warning strings for a variant."""
    flags = []
    vtype = v.get("VAR_TYPE","SNP")
    fs_max = FILTER_THRESHOLDS["FS_max_indel"] if vtype=="INDEL" else FILTER_THRESHOLDS["FS_max_snp"]
    if v.get("QD") is not None and v["QD"] < FILTER_THRESHOLDS["QD_min"]:
        flags.append(f"QD={v['QD']:.1f}<{FILTER_THRESHOLDS['QD_min']}")
    if v.get("FS") is not None and v["FS"] > fs_max:
        flags.append(f"FS={v['FS']:.1f}>{fs_max}")
    if vtype=="SNP" and v.get("SOR") is not None and v["SOR"] > FILTER_THRESHOLDS["SOR_max_snp"]:
        flags.append(f"SOR={v['SOR']:.2f}>{FILTER_THRESHOLDS['SOR_max_snp']}")
    if v.get("MQ") is not None and v["MQ"] < FILTER_THRESHOLDS["MQ_min"]:
        flags.append(f"MQ={v['MQ']:.1f}<{FILTER_THRESHOLDS['MQ_min']}")
    if v.get("MQRankSum") is not None and v["MQRankSum"] < FILTER_THRESHOLDS["MQRankSum_min"]:
        flags.append(f"MQRkS={v['MQRankSum']:.1f}")
    if v.get("ReadPosRankSum") is not None and v["ReadPosRankSum"] < FILTER_THRESHOLDS["ReadPosRankSum_min"]:
        flags.append(f"RPRkS={v['ReadPosRankSum']:.1f}")
    if v.get("VQSLOD") is not None and v["VQSLOD"] < FILTER_THRESHOLDS["VQSLOD_min"]:
        flags.append(f"VQSLOD={v['VQSLOD']:.2f}<0")
    return flags

# ═══════════════════════════════════════════════════════════════════════════════
# REPORT BUILDER
# ═══════════════════════════════════════════════════════════════════════════════

def build_report(sample_id, variants, sample_gts, all_samples, out_path):
    gen_date = datetime.now().strftime("%d %B %Y at %H:%M")
    styles   = make_styles()
    callback = make_canvas_callback(sample_id, gen_date)

    doc = SimpleDocTemplate(
        out_path,
        pagesize=(25*cm,35*cm),
        leftMargin=1*cm, rightMargin=1*cm,
        topMargin=1.4*cm, bottomMargin=1.2*cm,
        title=f"AML Variant Report — {sample_id}",
        author="AML Variant Dashboard",
    )

    story = []
    S     = styles

    # ── Per-sample GT lookup ─────────────────────────────────────────────────
    my_gts = sample_gts.get(sample_id, [None]*len(variants))

    # ── Collect carrier variants ─────────────────────────────────────────────
    carrier_vars   = []
    carrier_gts    = []
    for i, v in enumerate(variants):
        gt_info = my_gts[i] if i < len(my_gts) else None
        if gt_info and gt_info.get("IS_CARRIER"):
            carrier_vars.append(v)
            carrier_gts.append(gt_info)

    # PASS only
    pass_vars = [v for v in carrier_vars if v["FILTER"] == "PASS"]
    pass_gts  = [g for v, g in zip(carrier_vars, carrier_gts) if v["FILTER"]=="PASS"]

    hamlet_vars = [v for v in pass_vars if v["IS_HAMLET"]]
    hamlet_gts  = [g for v, g in zip(pass_vars, pass_gts) if v["IS_HAMLET"]]

    novel_noncan = [v for v in pass_vars
                    if v["IS_NOVEL"] and not v["IS_HAMLET"]
                    and v.get("IMPACT") in ("HIGH","MODERATE")]
    novel_nc_gts = [g for v, g in zip(pass_vars, pass_gts)
                    if v["IS_NOVEL"] and not v["IS_HAMLET"]
                    and v.get("IMPACT") in ("HIGH","MODERATE")]

    # ── Cohort stats ──────────────────────────────────────────────────────────
    total_vars   = len(variants)
    pass_total   = sum(1 for v in variants if v["FILTER"]=="PASS")
    snp_total    = sum(1 for v in variants if v["FILTER"]=="PASS" and v["VAR_TYPE"]=="SNP")
    indel_total  = sum(1 for v in variants if v["FILTER"]=="PASS" and v["VAR_TYPE"]=="INDEL")

    # Per-sample stats
    n_called   = len(carrier_vars)
    n_pass     = len(pass_vars)
    n_hamlet   = len(hamlet_vars)
    n_novel_nc = len(novel_noncan)

    het_count   = sum(1 for g in carrier_gts if "/" in g["GT"] and g["GT"] not in ("0/0","1/1","././."))
    hom_count   = sum(1 for g in carrier_gts if "1/1" in g["GT"] or "2/2" in g["GT"])
    miss_count  = sum(1 for g in my_gts if g and g.get("IS_MISSING"))
    mean_dp     = (sum(g["DP"] for g in carrier_gts if g.get("DP")) /
                   max(sum(1 for g in carrier_gts if g.get("DP")), 1))
    mean_gq     = (sum(g["GQ"] for g in carrier_gts if g.get("GQ")) /
                   max(sum(1 for g in carrier_gts if g.get("GQ")), 1))

    # ═════════════════════════════════════════════════════════════════════════
    # PAGE 1: COVER
    # ═════════════════════════════════════════════════════════════════════════

    story.append(Spacer(1, 2.5*cm))

    # Top coloured banner
    banner_data = [[Paragraph("AML VARIANT REPORT", ParagraphStyle("banner",
        fontName="Helvetica-Bold", fontSize=22, textColor=colors.white,
        alignment=TA_CENTER))]]
    banner = Table(banner_data, colWidths=[A4[0]-2*cm])
    banner.setStyle(TableStyle([
        ("BACKGROUND",    (0,0),(-1,-1), C_NAVY),
        ("TOPPADDING",    (0,0),(-1,-1), 14),
        ("BOTTOMPADDING", (0,0),(-1,-1), 14),
        ("ALIGN",         (0,0),(-1,-1), "CENTER"),
    ]))
    story.append(banner)
    story.append(Spacer(1, 0.4*cm))

    story.append(Paragraph(f"Sample: <b>{sample_id}</b>",
        ParagraphStyle("cov_samp", fontName="Helvetica-Bold", fontSize=20,
                       textColor=C_TEAL, alignment=TA_CENTER, spaceAfter=6)))

    story.append(Spacer(1, 0.5*cm))
    story.append(Paragraph(f"Generated on: {gen_date}", S["cover_meta"]))
    story.append(Paragraph(
        "Pipeline: nf-core/SAREK 3.8.1 · GATK HaplotypeCaller Joint Germline · "
        "VQSR Recalibration · VEP 115 · GRCh38",
        S["cover_meta"]))
    story.append(Paragraph(
        f"Cohort: {len(all_samples)} WES samples · HAMLET gene panel v2.3",
        S["cover_meta"]))

    story.append(Spacer(1, 0.8*cm))
    story.append(hrule(C_TEAL, 1.5))
    story.append(Spacer(1, 0.5*cm))

    # Summary KPIs on cover
    story.append(kpi_table([
        ("PASS variants\n(cohort)",    f"{pass_total:,}",  C_LIGHT),
        ("Variants carried",           f"{n_called:,}",     C_LIGHT),
        ("PASS carried",               f"{n_pass:,}",       C_PASS_GREEN),
        ("HAMLET gene hits",           f"{n_hamlet:,}",     C_HAMLET_BG),
        ("Novel non-canonical\nH/M",   f"{n_novel_nc:,}",   C_NOVEL_BG),
    ], S))

    story.append(Spacer(1, 1.5*cm))

    # Table of contents
    toc_items = [
        "1.  QC Summary",
        "2.  HAMLET Canonical Variants",
        "3.  Non-Canonical Novel Candidates",
        "4.  Glossary",
    ]
    toc_data = [[Paragraph("Table of Contents", ParagraphStyle("toc_h",
        fontName="Helvetica-Bold", fontSize=11, textColor=C_NAVY,
        alignment=TA_CENTER))]]
    for item in toc_items:
        toc_data.append([Paragraph(item, ParagraphStyle("toc_i",
            fontName="Helvetica", fontSize=9, textColor=colors.black,
            leftIndent=20))])
    toc = Table(toc_data, colWidths=[A4[0]-2*cm])
    toc.setStyle(TableStyle([
        ("BACKGROUND", (0,0),(-1,0),  C_MID),
        ("TOPPADDING", (0,0),(-1,-1), 4),
        ("BOTTOMPADDING",(0,0),(-1,-1),4),
        ("BOX", (0,0),(-1,-1), 0.5, C_NAVY),
        ("LINEBELOW",(0,0),(-1,0), 0.5, C_NAVY),
    ]))
    story.append(toc)

    story.append(PageBreak())

    # ═════════════════════════════════════════════════════════════════════════
    # PAGE 2: QC SUMMARY
    # ═════════════════════════════════════════════════════════════════════════

    story.append(Paragraph("1. Quality Control Summary", S["section"]))
    story.append(hrule())

    story.append(Paragraph("1.1  Cohort-Level Variant Statistics", S["subsection"]))
    story.append(Paragraph(
        "Statistics reflect all variants in the joint-germline VCF after VQSR recalibration, "
        "prior to per-sample filtering.",
        S["body"]))

    cohort_stats = [
        ["Metric", "Value"],
        ["Total variants (all FILTER statuses)", f"{total_vars:,}"],
        ["PASS variants", f"{pass_total:,}"],
        ["PASS SNPs", f"{snp_total:,}"],
        ["PASS INDELs", f"{indel_total:,}"],
        ["Number of samples", str(len(all_samples))],
    ]
    ct = Table(cohort_stats, colWidths=[9*cm, 4*cm])
    ct.setStyle(TableStyle([
        ("BACKGROUND",    (0,0),(-1,0),  C_NAVY),
        ("TEXTCOLOR",     (0,0),(-1,0),  colors.white),
        ("FONTNAME",      (0,0),(-1,0),  "Helvetica-Bold"),
        ("ROWBACKGROUNDS",(0,1),(-1,-1), [C_LIGHT, colors.white]),
        ("FONTSIZE",      (0,0),(-1,-1), 8),
        ("TOPPADDING",    (0,0),(-1,-1), 4),
        ("BOTTOMPADDING", (0,0),(-1,-1), 4),
        ("LEFTPADDING",   (0,0),(-1,-1), 6),
        ("GRID",          (0,0),(-1,-1), 0.3, C_MID),
    ]))
    story.append(ct)
    story.append(Spacer(1, 0.3*cm))

    story.append(Paragraph("1.2  Per-Sample Statistics", S["subsection"]))

    sample_stats = [
        ["Metric", "Value", "Threshold / Interpretation"],
        ["Variants carried (all FILTER)", f"{n_called:,}", "—"],
        ["PASS variants carried", f"{n_pass:,}", "Variants passing all GATK filters"],
        ["Heterozygous (carried)", f"{het_count:,}", "Expected: ~2× hom_alt for WES germline"],
        ["Homozygous alternate (carried)", f"{hom_count:,}", "—"],
        ["Missing genotype calls", f"{miss_count:,}", "Low values preferred"],
        ["Mean read depth (carried)", f"{mean_dp:.1f}x", "≥ 10× recommended"],
        ["Mean genotype quality (carried)", f"{mean_gq:.1f}", "≥ 20 reliable"],
        ["HAMLET canonical variants (PASS)", f"{n_hamlet:,}", "Known AML genes"],
        ["Novel non-canonical HIGH/MOD (PASS)", f"{n_novel_nc:,}", "Primary discovery targets"],
    ]
    st = Table(sample_stats, colWidths=[7*cm, 2.5*cm, 8*cm])
    st.setStyle(TableStyle([
        ("BACKGROUND",    (0,0),(-1,0),  C_NAVY),
        ("TEXTCOLOR",     (0,0),(-1,0),  colors.white),
        ("FONTNAME",      (0,0),(-1,0),  "Helvetica-Bold"),
        ("ROWBACKGROUNDS",(0,1),(-1,-1), [C_LIGHT, colors.white]),
        ("FONTSIZE",      (0,0),(-1,-1), 8),
        ("TOPPADDING",    (0,0),(-1,-1), 4),
        ("BOTTOMPADDING", (0,0),(-1,-1), 4),
        ("LEFTPADDING",   (0,0),(-1,-1), 6),
        ("GRID",          (0,0),(-1,-1), 0.3, C_MID),
    ]))
    story.append(st)
    story.append(Spacer(1, 0.3*cm))

    story.append(Paragraph("1.3  Quality Filter Thresholds Applied", S["subsection"]))
    story.append(Paragraph(
        "The following thresholds conform to GATK Best Practices for germline short-variant discovery. "
        "Variants in Sections 2 and 3 pass the GATK FILTER=PASS status (VQSR recalibration). "
        "Individual per-variant QC warnings are shown in the tables below where applicable.",
        S["body"]))

    thresh_data = [
        ["Metric", "Threshold", "Rationale"],
        ["FILTER", "PASS", "Variant passed GATK VQSR recalibration"],
        ["VQSLOD", f"≥ {FILTER_THRESHOLDS['VQSLOD_min']:.1f}", "Bayesian variant quality score; negative = low confidence"],
        ["QD (Quality by Depth)", f"≥ {FILTER_THRESHOLDS['QD_min']:.1f}", "Normalises quality score by depth; QD<2 → likely FP"],
        ["FS (FisherStrand SNP)", f"≤ {FILTER_THRESHOLDS['FS_max_snp']:.0f}", "Strand bias test; high FS → artefactual strand enrichment"],
        ["FS (FisherStrand INDEL)", f"≤ {FILTER_THRESHOLDS['FS_max_indel']:.0f}", "Same metric, more permissive threshold for indels"],
        ["SOR (StrandOddsRatio)", f"≤ {FILTER_THRESHOLDS['SOR_max_snp']:.1f}", "Symmetric strand bias; SOR>3 for SNPs is flagged"],
        ["MQ (Mapping Quality)", f"≥ {FILTER_THRESHOLDS['MQ_min']:.0f}", "RMS mapping quality; MQ<40 suggests mismapping"],
        ["MQRankSum", f"≥ {FILTER_THRESHOLDS['MQRankSum_min']:.1f}", "MQ comparison ref vs alt reads; <-12.5 is suspicious"],
        ["ReadPosRankSum", f"≥ {FILTER_THRESHOLDS['ReadPosRankSum_min']:.1f}", "Read position bias; <-8 suggests end-of-read artefact"],
        ["gnomAD exome AF", f"≤ {FILTER_THRESHOLDS['gnomAD_AF_max']:.2f}", "Exclude common population variants (>1%)"],
        ["Per-sample DP", f"≥ {FILTER_THRESHOLDS['DP_min']}", "Minimum read depth for reliable genotype"],
        ["Per-sample GQ", f"≥ {FILTER_THRESHOLDS['GQ_min']}", "Genotype quality ≥ 20 is considered reliable"],
    ]
    tt = Table(thresh_data, colWidths=[5*cm, 3*cm, 9.5*cm])
    tt.setStyle(TableStyle([
        ("BACKGROUND",    (0,0),(-1,0),  C_NAVY),
        ("TEXTCOLOR",     (0,0),(-1,0),  colors.white),
        ("FONTNAME",      (0,0),(-1,0),  "Helvetica-Bold"),
        ("ROWBACKGROUNDS",(0,1),(-1,-1), [C_LIGHT, colors.white]),
        ("FONTSIZE",      (0,0),(-1,-1), 7.5),
        ("TOPPADDING",    (0,0),(-1,-1), 3),
        ("BOTTOMPADDING", (0,0),(-1,-1), 3),
        ("LEFTPADDING",   (0,0),(-1,-1), 5),
        ("GRID",          (0,0),(-1,-1), 0.3, C_MID),
    ]))
    story.append(tt)
    story.append(Paragraph(
        "Note: Novel status is assessed against the Existing_variation field of the canonical VEP "
        "transcript (dbSNP rsIDs and ClinVar accessions). Variants absent from this field are "
        "classified as novel.",
        S["caption"]))

    story.append(PageBreak())

    # ═════════════════════════════════════════════════════════════════════════
    # PAGE 3+: HAMLET CANONICAL VARIANTS
    # ═════════════════════════════════════════════════════════════════════════

    story.append(Paragraph("2. HAMLET Canonical Variants", S["section"]))
    story.append(hrule(C_AMBER))
    story.append(Paragraph(
        "All PASS variants carried by this sample in the 41 HAMLET AML gene panel. "
        "These genes are the established targets for AML genomic profiling. "
        "Novel variants in this section may represent previously undescribed alleles of "
        "canonical AML drivers.",
        S["body"]))
    story.append(Paragraph(
        f"Total HAMLET variants carried: <b>{n_hamlet}</b>",
        S["body"]))
    story.append(Spacer(1, 0.15*cm))

    if not hamlet_vars:
        story.append(Paragraph(
            "No PASS variants in HAMLET canonical genes carried by this sample.",
            S["body"]))
    else:
        # Sort: novel first, then by impact rank, then VQSLOD desc
        impact_rank = {"HIGH":0,"MODERATE":1,"LOW":2,"MODIFIER":3}
        sorted_pairs = sorted(zip(hamlet_vars, hamlet_gts),
            key=lambda x: (not x[0]["IS_NOVEL"],
                           impact_rank.get(x[0]["IMPACT"],9),
                           -(x[0]["VQSLOD"] or -999)))
        hamlet_vars_s = [p[0] for p in sorted_pairs]
        hamlet_gts_s  = [p[1] for p in sorted_pairs]

        story.append(variant_table(hamlet_vars_s, hamlet_gts_s, S))
        story.append(Paragraph(
            "Columns: Impact colour — RED = HIGH (protein-disrupting), ORANGE = MODERATE (missense/in-frame). "
            "VAF = variant allele frequency for this sample. Cohort AF = fraction of cohort alleles. "
            "Novel variants (absent from dbSNP/ClinVar) are sorted first.",
            S["caption"]))

        # QC warning table for HAMLET variants
        qc_rows = [(v, g) for v, g in zip(hamlet_vars_s, hamlet_gts_s) if qc_flag(v)]
        if qc_rows:
            story.append(Spacer(1, 0.2*cm))
            story.append(Paragraph("QC Warnings (HAMLET variants):", S["subsection"]))
            qc_data = [["Location", "Gene", "Flags"]]
            for v, g in qc_rows:
                qc_data.append([
                    f"{v['CHROM']}:{v['POS']}",
                    v["SYMBOL"],
                    ", ".join(qc_flag(v)),
                ])
            qct = Table(qc_data, colWidths=[4*cm, 3*cm, 10.5*cm])
            qct.setStyle(TableStyle([
                ("BACKGROUND",   (0,0),(-1,0), C_AMBER),
                ("TEXTCOLOR",    (0,0),(-1,0), colors.white),
                ("FONTNAME",     (0,0),(-1,0), "Helvetica-Bold"),
                ("FONTSIZE",     (0,0),(-1,-1), 7.5),
                ("TOPPADDING",   (0,0),(-1,-1), 3),
                ("BOTTOMPADDING",(0,0),(-1,-1), 3),
                ("LEFTPADDING",  (0,0),(-1,-1), 5),
                ("GRID",         (0,0),(-1,-1), 0.3, C_MID),
                ("ROWBACKGROUNDS",(0,1),(-1,-1), [colors.HexColor("#FFF8E1"), colors.white]),
            ]))
            story.append(qct)

    story.append(PageBreak())

    # ═════════════════════════════════════════════════════════════════════════
    # PAGE 4+: NON-CANONICAL NOVEL CANDIDATES
    # ═════════════════════════════════════════════════════════════════════════

    story.append(Paragraph("3. Non-Canonical Novel Candidates", S["section"]))
    story.append(hrule(C_ACCENT))
    story.append(Paragraph(
        "Novel HIGH/MODERATE impact variants in genes <b>outside</b> the HAMLET panel. "
        "A variant is classified as novel if it has no entry in dbSNP or ClinVar for "
        "the canonical VEP transcript. These are the primary scientific discovery targets "
        "— potential new AML-associated genes or rare pathogenic alleles.",
        S["body"]))
    story.append(Paragraph(
        f"Total non-canonical novel HIGH/MOD carried: <b>{n_novel_nc}</b>",
        S["body"]))
    story.append(Spacer(1, 0.15*cm))

    if not novel_noncan:
        story.append(Paragraph(
            "No novel non-canonical HIGH/MODERATE impact variants carried by this sample.",
            S["body"]))
    else:
        impact_rank = {"HIGH":0,"MODERATE":1,"LOW":2,"MODIFIER":3}
        sorted_pairs = sorted(zip(novel_noncan, novel_nc_gts),
            key=lambda x: (impact_rank.get(x[0]["IMPACT"],9),
                           -(x[0]["VQSLOD"] or -999)))
        nc_vars_s = [p[0] for p in sorted_pairs]
        nc_gts_s  = [p[1] for p in sorted_pairs]

        story.append(variant_table(nc_vars_s, nc_gts_s, S))
        story.append(Paragraph(
            "All variants shown are: FILTER=PASS · Novel (absent from dbSNP/ClinVar) · "
            "Outside HAMLET gene panel · Impact = HIGH or MODERATE.",
            S["caption"]))

        # QC warnings for non-canonical
        qc_rows = [(v, g) for v, g in zip(nc_vars_s, nc_gts_s) if qc_flag(v)]
        if qc_rows:
            story.append(Spacer(1, 0.2*cm))
            story.append(Paragraph("QC Warnings (Non-canonical variants):", S["subsection"]))
            qc_data = [["Location", "Gene", "Consequence", "Flags"]]
            for v, g in qc_rows:
                qc_data.append([
                    f"{v['CHROM']}:{v['POS']}",
                    v["SYMBOL"],
                    v["CONSEQUENCE"].replace("_"," "),
                    ", ".join(qc_flag(v)),
                ])
            qct = Table(qc_data, colWidths=[4*cm, 2.5*cm, 4*cm, 7*cm])
            qct.setStyle(TableStyle([
                ("BACKGROUND",   (0,0),(-1,0), C_ACCENT),
                ("TEXTCOLOR",    (0,0),(-1,0), colors.white),
                ("FONTNAME",     (0,0),(-1,0), "Helvetica-Bold"),
                ("FONTSIZE",     (0,0),(-1,-1), 7.5),
                ("TOPPADDING",   (0,0),(-1,-1), 3),
                ("BOTTOMPADDING",(0,0),(-1,-1), 3),
                ("LEFTPADDING",  (0,0),(-1,-1), 5),
                ("GRID",         (0,0),(-1,-1), 0.3, C_MID),
                ("ROWBACKGROUNDS",(0,1),(-1,-1), [C_NOVEL_BG, colors.white]),
            ]))
            story.append(qct)

    story.append(PageBreak())

    # ═════════════════════════════════════════════════════════════════════════
    # LAST PAGE: GLOSSARY
    # ═════════════════════════════════════════════════════════════════════════

    story.append(Paragraph("4. Glossary", S["section"]))
    story.append(hrule())
    story.append(Paragraph(
        "Definitions of all fields and metrics used in this report.",
        S["body"]))
    story.append(Spacer(1, 0.2*cm))

    for term, defn in GLOSSARY:
        story.append(KeepTogether([
            Paragraph(term, S["gloss_term"]),
            Paragraph(defn, S["gloss_def"]),
            Spacer(1, 0.15*cm),
        ]))

    # ── Build ─────────────────────────────────────────────────────────────────
    doc.build(story, onFirstPage=callback, onLaterPages=callback)
    return out_path

# ═══════════════════════════════════════════════════════════════════════════════
# ENTRY POINT
# ═══════════════════════════════════════════════════════════════════════════════

def main():
    parser = argparse.ArgumentParser(
        description="Generate one PDF report per sample from a GATK VEP-annotated VCF.")
    parser.add_argument("vcf", help="Path to VCF file")
    parser.add_argument("--out_dir", default="reports",
                        help="Output directory for PDF reports (default: ./reports)")
    parser.add_argument("--samples", nargs="*",
                        help="Subset of samples to process (default: all)")
    args = parser.parse_args()

    if not os.path.isfile(args.vcf):
        print(f"ERROR: VCF not found: {args.vcf}", file=sys.stderr)
        sys.exit(1)

    os.makedirs(args.out_dir, exist_ok=True)

    print(f"Parsing VCF: {args.vcf}")
    variants, all_samples, sample_gts = load_vcf(args.vcf)
    print(f"  Loaded {len(variants):,} variants · {len(all_samples)} samples")

    target = args.samples if args.samples else all_samples
    invalid = [s for s in target if s not in all_samples]
    if invalid:
        print(f"WARNING: unknown samples: {invalid}", file=sys.stderr)
    target = [s for s in target if s in all_samples]

    for i, sname in enumerate(target, 1):
        out_path = os.path.join(args.out_dir, f"AML_report_{sname}.pdf")
        print(f"  [{i}/{len(target)}] Generating report for {sname} → {out_path}")
        build_report(sname, variants, sample_gts, all_samples, out_path)

    print(f"\nDone. {len(target)} reports written to: {args.out_dir}/")

if __name__ == "__main__":
    main()