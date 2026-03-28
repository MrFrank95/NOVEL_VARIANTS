import streamlit as st
import pandas as pd

df = st.session_state["df_view"]
samples = st.session_state["samples"]

st.title("🔍 Variant Browser")
st.markdown("Full variant table with interactive filters. Expand any row for per-sample genotype breakdown.")

# ── Glossary ──────────────────────────────────────────────────────────────────
with st.expander("📖 Glossary — what do these fields mean?"):
    st.markdown("""
    **Variant identifiers**
    - **CHROM / POS / REF / ALT** — Genomic coordinates: chromosome, position, reference allele, and alternate allele.
    - **ID** — Variant identifier (e.g. rsID from dbSNP, or `.` if none).
    - **HGVSc** — HGVS coding DNA notation (e.g. `c.123A>T`): describes the base change relative to the transcript.
    - **HGVSp** — HGVS protein notation (e.g. `p.Lys41Asn`): describes the amino-acid consequence.
    - **EXISTING_VARIATION** — Known identifiers for this site (dbSNP rsIDs, ClinVar accessions). Empty = novel.
    - **VAR_TYPE** — Variant type: SNV, insertion, deletion, MNV, etc.

    **Gene & annotation flags**
    - **SYMBOL** — HGNC gene symbol of the overlapping gene.
    - **IS_NOVEL** — `True` if the variant has no entry in `Existing_variation` for the canonical transcript (not in dbSNP / ClinVar).
    - **IS_AML_GENE** — `True` if the gene belongs to the curated AML gene panel.
    - **FILTER** — GATK filter status. `PASS` = variant passed all quality filters; other values (e.g. `LowQual`) indicate a flagged call.

    **Impact & consequence**
    - **IMPACT** — VEP severity tier: `HIGH` (likely protein-disrupting: stop-gain, frameshift, splice-site), `MODERATE` (missense, in-frame indel), `LOW` (synonymous), `MODIFIER` (non-coding / intergenic).
    - **CONSEQUENCE** — The specific molecular effect predicted by VEP (e.g. `missense_variant`, `splice_acceptor_variant`).
    - **SIFT** — Predicts whether an amino-acid substitution affects protein function. Score < 0.05 = **deleterious**; higher = tolerated.
    - **PolyPhen** — Predicts the structural/functional impact of a missense change. Values: `probably_damaging`, `possibly_damaging`, `benign`.

    **Quality scores**
    - **VQSLOD** — Variant Quality Score Log-Odds from GATK VQSR. Higher = more likely a true variant; negative values indicate low-confidence calls.
    - **QD** — Quality by Depth: variant quality divided by read depth. QD < 2 can indicate a false positive.
    - **FS** — FisherStrand: strand bias test (Fisher's exact). High values (> 60 for SNPs) suggest strand-specific artefacts.
    - **MQ** — RMS Mapping Quality. Low MQ (< 40) may reflect mismapped reads.

    **Allele frequencies**
    - **Cohort AF** — Allele frequency within this study cohort.
    - **gnomAD-e AF** — Allele frequency in the gnomAD exome cohort (> 125 000 exomes). Absent = not observed, which supports novelty.
    - **MAX_AF** — Maximum allele frequency across all gnomAD sub-populations.

    **Clinical**
    - **CLIN_SIG** — ClinVar clinical significance: `pathogenic`, `likely_pathogenic`, `uncertain_significance`, `benign`, etc. Empty = not in ClinVar.

    **Per-sample fields (drill-down table)**
    - **GT** — Genotype call (e.g. `0/1` = heterozygous, `1/1` = homozygous alt, `0/0` = homozygous ref).
    - **GT_CLASS** — Simplified label: `het`, `hom_alt`, `hom_ref`, `missing`.
    - **DP** — Total read depth at the variant site in that sample.
    - **GQ** — Genotype Quality: Phred-scaled confidence in the genotype call. GQ ≥ 20 is generally considered reliable.
    """)

# ── Filters ──────────────────────────────────────────────────────────────────
with st.expander("🔧 Filters", expanded=True):
    col1,col2,col3,col4 = st.columns(4)
    with col1:
        chroms = ["All"] + sorted(df["CHROM"].unique().tolist(),
                  key=lambda x: int(x.replace("chr","").replace("X","23").replace("Y","24")))
        sel_chrom = st.selectbox("Chromosome", chroms)
    with col2:
        impacts = ["All"] + sorted(df["IMPACT"].dropna().unique().tolist())
        sel_impact = st.selectbox("VEP Impact", impacts,
                                  help="HIGH = protein-disrupting · MODERATE = missense/in-frame indel · LOW = synonymous · MODIFIER = non-coding")
    with col3:
        novel_sel = st.selectbox("Novelty", ["All","Novel only","Known only"],
                                 help="Novel = absent from dbSNP / ClinVar in the canonical transcript")
    with col4:
        aml_sel = st.checkbox("AML genes only", value=False,
                              help="Restrict to variants overlapping the curated AML gene panel")

    col5,col6 = st.columns(2)
    with col5:
        gene_search = st.text_input("Gene symbol search", placeholder="e.g. EZH2, TP53")
    with col6:
        consq_search = st.text_input("Consequence search", placeholder="e.g. frameshift, missense")

filt = df.copy()
if sel_chrom != "All":
    filt = filt[filt["CHROM"] == sel_chrom]
if sel_impact != "All":
    filt = filt[filt["IMPACT"] == sel_impact]
if novel_sel == "Novel only":
    filt = filt[filt["IS_NOVEL"]]
elif novel_sel == "Known only":
    filt = filt[~filt["IS_NOVEL"]]
if aml_sel:
    filt = filt[filt["IS_AML_GENE"]]
if gene_search:
    filt = filt[filt["SYMBOL"].str.contains(gene_search, case=False, na=False)]
if consq_search:
    filt = filt[filt["CONSEQUENCE"].str.contains(consq_search, case=False, na=False)]

st.caption(f"Showing {len(filt):,} of {len(df):,} variants")

# ── Table ─────────────────────────────────────────────────────────────────────
display_cols = ["CHROM","POS","ID","REF","ALT","SYMBOL","CONSEQUENCE","IMPACT",
                "FILTER","VAR_TYPE","IS_NOVEL","IS_AML_GENE",
                "HGVSc","HGVSp","AF","VQSLOD","QD","FS","MQ",
                "gnomADe_AF","MAX_AF","CLIN_SIG","SIFT","PolyPhen",
                "EXISTING_VARIATION"]
show = filt[[c for c in display_cols if c in filt.columns]].reset_index(drop=True)

st.caption(
    "**Column guide:** IMPACT — VEP severity tier · HGVSc — coding DNA change · "
    "HGVSp — protein change · VQSLOD — variant confidence (higher = better) · "
    "Cohort AF — allele frequency in this cohort · gnomAD-e AF — population frequency · "
    "SIFT/PolyPhen — missense pathogenicity predictors · FILTER: PASS = high-quality call"
)

st.dataframe(show, width='stretch', height=450,
             column_config={
                 "IS_NOVEL":   st.column_config.CheckboxColumn("Novel",
                                help="True = absent from dbSNP / ClinVar (canonical transcript)"),
                 "IS_AML_GENE": st.column_config.CheckboxColumn("AML Gene",
                                help="True = gene is part of the curated AML panel"),
                 "VQSLOD":     st.column_config.NumberColumn("VQSLOD", format="%.2f",
                                help="GATK variant confidence score — higher is better; negative = low confidence"),
                 "gnomADe_AF": st.column_config.NumberColumn("gnomAD-e AF", format="%.4f",
                                help="Allele frequency in gnomAD exomes (>125 000 samples); absent = not observed"),
                 "AF":         st.column_config.NumberColumn("Cohort AF", format="%.3f",
                                help="Allele frequency within this study cohort"),
             })

csv = show.to_csv(index=False).encode()
st.download_button("⬇️ Download filtered table as CSV", csv,
                   "variants_filtered.csv", mime="text/csv")

# ── Per-variant sample genotype drill-down ────────────────────────────────────
st.markdown("---")
st.subheader("🔎 Per-Variant Genotype Drill-down")
st.caption(
    "Select a variant below to see genotype calls across all samples. "
    "**GT**: genotype (0/1 = heterozygous, 1/1 = homozygous alt, 0/0 = hom ref) · "
    "**DP**: read depth · **GQ**: genotype quality (≥20 reliable)"
)

sdf = st.session_state["sample_df_view"]

variant_options = show.apply(
    lambda r: f"{r['CHROM']}:{r['POS']} {r['REF']}>{r['ALT']}  [{r.get('SYMBOL','')}]",
    axis=1
).tolist()

if variant_options:
    sel_var = st.selectbox("Select variant", variant_options[:200])
    idx = variant_options.index(sel_var)
    row = show.iloc[idx]

    var_sample_data = sdf[(sdf["CHROM"]==row["CHROM"]) & (sdf["POS"]==row["POS"])]
    if len(var_sample_data):
        gt_class_colors = {"hom_ref":"#d5e8d4","het":"#ffe6cc","hom_alt":"#f8cecc","missing":"#f0f0f0"}
        var_sample_data = var_sample_data[["sample","GT","GT_CLASS","DP","GQ"]].copy()
        st.dataframe(
            var_sample_data.reset_index(drop=True).style.applymap(
                lambda v: f"background-color:{gt_class_colors.get(v,'')}", subset=["GT_CLASS"]
            ),
            width='stretch'
        )
    else:
        st.info("No per-sample data available for this variant.")
else:
    st.info("No variants to display with current filters.")