# 🧬 AML Novel Variant Dashboard

Streamlit dashboard for exploring joint germline VCF output from nf-core/Sarek,
focused on **novel variant discovery for Acute Myeloid Leukaemia (AML)**.

## Setup

```bash
# 1. Install dependencies
pip install -r requirements.txt

# 2. Place your VCF file
cp /path/to/joint_germline_recalibrated_VEP_ann.vcf data/

# 3. Run the app
streamlit run app.py
```

## Pages

| Page | Description |
|------|-------------|
| 🏠 Overview | Cohort-level summary cards, variant counts, chromosome distribution |
| 🔴 Novel AML Candidates | Prioritised novel HIGH/MODERATE impact variants — the discovery page |
| 📊 Quality Metrics | VQSLOD, QD, FS, MQ, SOR distributions (novel vs known) |
| 👥 Per-Sample Burden | Per-sample variant counts, het/hom ratio, AML gene heatmap |
| 🔍 Variant Browser | Searchable/filterable full variant table with genotype drill-down |

## Novel Variant Definition

A variant is flagged as **novel** when the `Existing_variation` field in the canonical
VEP transcript annotation is empty — i.e. the variant is absent from dbSNP and ClinVar.

## AML Gene Panel

The dashboard includes a curated panel of ~54 AML driver genes including:
FLT3, NPM1, DNMT3A, IDH1/2, TET2, RUNX1, TP53, ASXL1/2, EZH2, RAD21, SETD2, and more.

## Data

- **VCF**: `joint_germline_recalibrated_VEP_ann.vcf`
- **Pipeline**: nf-core/Sarek 3.8.1, GATK 4.6.1, VEP 115, GRCh38
- **Samples**: 31 samples, WES
- **Variants**: ~11,597 total, ~11,461 PASS
