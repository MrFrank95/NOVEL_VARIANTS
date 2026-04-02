import streamlit as st
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go

df        = st.session_state["df_view"]
sample_df = st.session_state["sample_df_view"]
samples   = st.session_state["samples"]

from vcf_loader import HAMLET_GENES

st.title("🔴 Novel AML Variant Candidates")

st.markdown("""
Variants are classified into **two discovery tiers**:

| Tier | Definition | Interest level |
|------|-----------|---------------|
| 🟡 **HAMLET canonical** | Novel variant in one of the 41 HAMLET leukemia genes | Expected territory — confirms known biology |
| 🔴 **Non-canonical** | Novel variant in a gene **outside** the HAMLET panel | **Primary discovery target** — potential new AML genes |

*Novel = absent from `Existing_variation` in the canonical VEP transcript (not in dbSNP / ClinVar).*
""")

# ── Glossary ──────────────────────────────────────────────────────────────────
with st.expander("📖 Glossary — what do these fields mean?"):
    st.markdown("""
    **Variant identifiers**
    - **HGVSc** — Human Genome Variation Society notation at the *coding DNA* level (e.g. `c.123A>T`). Describes the base change relative to the transcript.
    - **HGVSp** — HGVS notation at the *protein* level (e.g. `p.Lys41Asn`). Shows the amino-acid consequence.
    - **VARIANT_CLASS** — Structural type of the variant: SNV (single nucleotide variant), insertion, deletion, MNV, etc.
    - **BIOTYPE** — Ensembl biotype of the overlapping transcript (e.g. `protein_coding`, `lncRNA`).

    **Impact & consequence**
    - **IMPACT** — VEP severity tier: `HIGH` (likely protein-disrupting: stop-gain, frameshift, splice-site), `MODERATE` (missense, in-frame indel), `LOW` (synonymous), `MODIFIER` (non-coding / intergenic).
    - **CONSEQUENCE** — The specific molecular effect predicted by VEP (e.g. `missense_variant`, `splice_acceptor_variant`).
    - **SIFT** — Predicts whether an amino-acid substitution affects protein function. Score < 0.05 = **deleterious** (tolerated otherwise).
    - **PolyPhen** — Predicts the effect of a missense change on protein structure/function. Values: `probably_damaging`, `possibly_damaging`, `benign`.

    **Quality scores**
    - **VQSLOD** — Variant Quality Score Log-Odds from GATK's VQSR (Variant Quality Score Recalibration). Higher = more likely a true variant; negative values indicate low-confidence calls. A common threshold is ≥ 0.
    - **QD** — Quality by Depth: variant quality score divided by read depth. Low QD (< 2) can indicate false positives. 
    - **FS** — FisherStrand: measures strand bias using Fisher's exact test. High values (> 60 for SNPs, > 200 for indels) suggest strand-specific artefacts.
    - **MQ** — RMS Mapping Quality of reads supporting the variant. Low MQ (< 40) may reflect mismapped reads.

    **Population frequency**
    - **gnomADe AF** — Allele frequency in the gnomAD *exome* cohort (> 125 000 exomes). Absent = not observed in gnomAD, which strengthens novelty.
    - **MAX_AF** — Maximum allele frequency across all gnomAD sub-populations. Useful for catching variants common in a specific ancestry even if rare overall.

    **Clinical / database**
    - **ClinSig (CLIN_SIG)** — ClinVar clinical significance: `pathogenic`, `likely_pathogenic`, `uncertain_significance`, `benign`, etc. Empty = not in ClinVar.
    - **IS_NOVEL** — `True` if the variant has no entry in `Existing_variation` (dbSNP rsID / ClinVar accession) for the canonical transcript.

    **Sample-level fields**
    - **GT** — Genotype call (e.g. `0/1` = heterozygous, `1/1` = homozygous alt).
    - **GT_CLASS** — Simplified genotype label: `het`, `hom_alt`, `hom_ref`, `unknown`.
    - **DP** — Total read depth at the variant site in that sample.
    - **GQ** — Genotype Quality: Phred-scaled confidence in the genotype call. GQ ≥ 20 is generally considered reliable.
    """)

# ── Sidebar filters ───────────────────────────────────────────────────────────
with st.sidebar:
    st.markdown("### 🔬 Candidate Filters")
    novel_only  = st.checkbox("Novel variants only", value=True,
                              help="Keep only variants absent from dbSNP / ClinVar (IS_NOVEL = True)")
    impacts     = st.multiselect("VEP Impact", ["HIGH","MODERATE","LOW","MODIFIER"],
                                 default=["HIGH","MODERATE"],
                                 help="HIGH = protein-disrupting; MODERATE = missense / in-frame indel")
    max_gnomad  = st.slider("Max gnomAD exome AF", 0.0, 0.05, 0.01,
                            step=0.001, format="%.3f",
                            help="Exclude variants seen at higher frequency in the gnomAD exome cohort")
    min_vqslod  = st.slider("Min VQSLOD", -5.0, 20.0, 2.0, step=0.5,
                            help="VQSLOD: GATK variant confidence score — higher is better. Negative values indicate low-confidence calls.")
    tier_filter = st.multiselect(
        "Gene Tier", ["HAMLET canonical","Non-canonical","Intergenic"],
        default=["HAMLET canonical","Non-canonical"])
    st.markdown("---")

# ── Apply filters ─────────────────────────────────────────────────────────────
cand = df.copy()
if novel_only:
    cand = cand[cand["IS_NOVEL"]]
if impacts:
    cand = cand[cand["IMPACT"].isin(impacts)]
if tier_filter:
    cand = cand[cand["GENE_TIER"].isin(tier_filter)]
cand = cand[cand["VQSLOD"] >= min_vqslod]
cand = cand[(cand["gnomADe_AF"].isna()) | (cand["gnomADe_AF"] <= max_gnomad)]

hamlet_cand = cand[cand["GENE_TIER"] == "HAMLET canonical"]
noncan_cand = cand[cand["GENE_TIER"] == "Non-canonical"]

# ── KPI cards ─────────────────────────────────────────────────────────────────
c1, c2, c3, c4, c5 = st.columns(5)
c1.metric("Total candidates",      len(cand),
          help="All variants passing the current sidebar filters")
c2.metric("🟡 HAMLET canonical",    len(hamlet_cand),
          help="Novel variants inside the 41 HAMLET leukemia genes")
c3.metric("🔴 Non-canonical",       len(noncan_cand),
          help="Novel variants OUTSIDE HAMLET — primary discovery targets")
c4.metric("HIGH impact",           (cand["IMPACT"] == "HIGH").sum(),
          help="Likely protein-disrupting: stop-gain, frameshift, splice-site")
c5.metric("MODERATE impact",       (cand["IMPACT"] == "MODERATE").sum(),
          help="Missense substitutions and in-frame indels")

st.markdown("---")

# ═══════════════════════════════════════════════════════════════════════════════
# SECTION 1 — Non-canonical (primary discovery)
# ═══════════════════════════════════════════════════════════════════════════════
st.subheader("🔴 Non-Canonical Novel Candidates — Primary Discovery Targets")
st.caption("Genes outside the HAMLET panel that carry novel HIGH/MODERATE impact variants. "
           "These are the most scientifically interesting findings.")

BIOLOGICAL_CONTEXT = {
    "ASXL2":   "ASXL family (like ASXL1) — Polycomb repressor complex",
    "ASXL3":   "ASXL family — Polycomb repressor complex",
    "TET3":    "TET family (like TET2) — DNA demethylation",
    "SETD2":   "H3K36 methyltransferase — epigenetic regulator, mutated in AML",
    "RAD21":   "Cohesin complex (like STAG2) — chromosome cohesion",
    "NUP98":   "Nucleoporin — well-known AML translocation partner (NUP98 fusions)",
    "DOT1L":   "H3K79 methyltransferase — therapeutic target in MLL-rearranged AML",
    "BRD4":    "BET bromodomain — epigenetic reader, active therapeutic target in AML",
    "KDM2B":   "H3K36 demethylase — Polycomb recruiter, implicated in leukaemia",
    "IKZF3":   "Ikaros family transcription factor — lymphoid/myeloid development",
    "SMARCA2":  "SWI/SNF chromatin remodelling complex (like SMARCB1)",
    "ARID1B":  "SWI/SNF chromatin remodelling complex",
    "DHX15":   "RNA helicase — pre-mRNA splicing factor",
    "FOXP1":   "Transcription factor — B-cell and myeloid development",
    "NOTCH2":  "NOTCH signalling — mutated in B-cell malignancies",
    "KLF2":    "Transcription factor — haematopoietic stem cell regulation",
    "AURKB":   "Aurora kinase B — mitosis regulator, therapeutic target",
    "AXL":     "Receptor tyrosine kinase — drug resistance in AML",
    "PTPRO":   "Protein tyrosine phosphatase — tumour suppressor",
    "MGA":     "MAX dimerisation protein — MYC network regulator",
    "CHD8":    "Chromodomain helicase — chromatin remodeller",
}

if len(noncan_cand) == 0:
    st.info("No non-canonical candidates with current filters.")
else:
    color_map = {"HIGH":"#C0392B","MODERATE":"#E67E22","LOW":"#F1C40F","MODIFIER":"#95A5A6"}
    gene_nc = noncan_cand.groupby(["SYMBOL","IMPACT"]).size().reset_index(name="Count")
    fig = px.bar(gene_nc, x="SYMBOL", y="Count", color="IMPACT",
                 color_discrete_map=color_map,
                 title="Non-Canonical Genes with Novel Variants (by VEP Impact)",
                 labels={"SYMBOL":"Gene","Count":"# Variants","IMPACT":"VEP Impact"})
    fig.update_layout(xaxis_tickangle=-45, margin=dict(t=40, b=10))
    st.plotly_chart(fig, width='stretch')

    # Expandable cards, HIGH first
    for _, row in noncan_cand.sort_values(["IMPACT","VQSLOD"],
                                           ascending=[True, False]).iterrows():
        impact_icon = "🔴" if row["IMPACT"] == "HIGH" else "🟠"
        ctx = BIOLOGICAL_CONTEXT.get(row["SYMBOL"], "")
        ctx_badge = f"  ·  _{ctx}_" if ctx else ""
        with st.expander(
            f"{impact_icon} **{row['SYMBOL']}** · {row['CONSEQUENCE']} · "
            f"`{row['CHROM']}:{row['POS']} {row['REF']}>{row['ALT']}`{ctx_badge}"
        ):
            col1, col2, col3 = st.columns(3)
            with col1:
                st.markdown("**Variant**")
                st.write(f"HGVSc (coding change): `{row['HGVSc']}`")
                st.write(f"HGVSp (protein change): `{row['HGVSp']}`")
                st.write(f"Class: `{row['VARIANT_CLASS']}`")
                st.write(f"Biotype: `{row['BIOTYPE']}`")
            with col2:
                st.markdown("**Quality** *(higher VQSLOD / QD = more confident call)*")
                st.write(f"VQSLOD (variant confidence): `{row['VQSLOD']}`")
                st.write(f"QD (quality/depth): `{row['QD']}`")
                st.write(f"FS (strand bias, lower = better): `{row['FS']}`")
                st.write(f"MQ (mapping quality, ≥40 preferred): `{row['MQ']}`")
            with col3:
                st.markdown("**Population / Clinical**")
                gnomad_str = (f"`{row['gnomADe_AF']:.6f}`"
                              if pd.notna(row["gnomADe_AF"])
                              else "**absent from gnomAD** ⚠️")
                st.markdown(f"gnomADe AF (population freq.): {gnomad_str}")
                st.write(f"MAX_AF (max across populations): `{row['MAX_AF'] if pd.notna(row['MAX_AF']) else 'N/A'}`")
                st.write(f"ClinSig (ClinVar significance): `{row['CLIN_SIG'] or 'not in ClinVar'}`")
                st.write(f"SIFT (<0.05 = deleterious): `{row['SIFT']}`")
                st.write(f"PolyPhen (missense pathogenicity): `{row['PolyPhen']}`")

            # Which samples carry this variant?
            carriers = row['CARRIERS']

            if carriers is not None:
                st.markdown(
                    f"**Carriers: {carriers}/{len(samples)} samples** "
                    #f"— GT: genotype (0/1 = het, 1/1 = hom alt) · "
                    #f"DP: read depth · GQ: genotype quality (≥20 reliable)"
                )
            else:
                st.info("No carriers found in sample data.")

            if ctx:
                st.info(f"💡 **Biological context:** {ctx}")

    # Download
    nc_cols = ["CHROM","POS","REF","ALT","SYMBOL","CONSEQUENCE","IMPACT",
               "HGVSc","HGVSp","VQSLOD","AF","gnomADe_AF","MAX_AF",
               "CLIN_SIG","SIFT","PolyPhen","EXISTING_VARIATION"]
    nc_dl = noncan_cand[[c for c in nc_cols if c in noncan_cand.columns]]
    st.download_button("⬇️ Download non-canonical candidates CSV",
                       nc_dl.to_csv(index=False).encode(),
                       "noncanonical_novel_candidates.csv", mime="text/csv")

st.markdown("---")

# ═══════════════════════════════════════════════════════════════════════════════
# SECTION 2 — HAMLET canonical (confirmatory)
# ═══════════════════════════════════════════════════════════════════════════════
st.subheader("🟡 HAMLET Canonical Novel Variants — Confirmatory")
st.caption("Novel variants inside the 41 HAMLET leukemia genes. These confirm known AML gene "
           "biology and may represent undescribed alleles of canonical drivers.")

if len(hamlet_cand) == 0:
    st.info("No HAMLET canonical candidates with current filters.")
else:
    color_map = {"HIGH":"#C0392B","MODERATE":"#E67E22","LOW":"#F1C40F","MODIFIER":"#95A5A6"}
    gene_h = hamlet_cand.groupby(["SYMBOL","IMPACT"]).size().reset_index(name="Count")
    fig2 = px.bar(gene_h, x="SYMBOL", y="Count", color="IMPACT",
                  color_discrete_map=color_map,
                  title="HAMLET Canonical Genes with Novel Variants (by VEP Impact)",
                  labels={"SYMBOL":"Gene","Count":"# Variants","IMPACT":"VEP Impact"})
    fig2.update_layout(xaxis_tickangle=-45, margin=dict(t=40, b=10))
    st.plotly_chart(fig2, width='stretch')

    disp_cols = ["CHROM","POS","SYMBOL","CONSEQUENCE","IMPACT","HGVSc","HGVSp",
                 "VQSLOD","AF","gnomADe_AF","CLIN_SIG","SIFT","PolyPhen"]
    show_h = hamlet_cand[[c for c in disp_cols if c in hamlet_cand.columns]].copy()
    st.caption(
        "**Column guide:** IMPACT — VEP severity tier · HGVSc — coding DNA change · "
        "HGVSp — protein change · VQSLOD — variant confidence (higher = better) · "
        "gnomADe_AF — population allele frequency · SIFT/PolyPhen — pathogenicity predictors"
    )
    show_h.reset_index(drop=True, inplace=True)
    show_h.index += 1
    st.dataframe(
        show_h.style.map(
            lambda v: {"HIGH": "background-color:#ffd6d6",
                       "MODERATE": "background-color:#ffe5c0"}.get(v, ""),
            subset=["IMPACT"]
        ),
        width='stretch', height=350
    )
    st.download_button("⬇️ Download HAMLET canonical candidates CSV",
                       show_h.to_csv(index=False).encode(),
                       "hamlet_canonical_candidates.csv", mime="text/csv")

st.markdown("---")

# ═══════════════════════════════════════════════════════════════════════════════
# SECTION 3 — Combined download
# ═══════════════════════════════════════════════════════════════════════════════
st.subheader("⬇️ Download All Candidates (both tiers)")
all_cols = ["CHROM","POS","REF","ALT","SYMBOL","GENE_TIER","CONSEQUENCE","IMPACT",
            "HGVSc","HGVSp","VQSLOD","AF","gnomADe_AF","MAX_AF",
            "CLIN_SIG","SIFT","PolyPhen","IS_NOVEL","EXISTING_VARIATION"]
all_dl = cand[[c for c in all_cols if c in cand.columns]]
st.download_button("⬇️ Download all candidates CSV",
                   all_dl.to_csv(index=False).encode(),
                   "all_novel_candidates.csv", mime="text/csv")