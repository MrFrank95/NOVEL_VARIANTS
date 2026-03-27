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

# ── Sidebar filters ───────────────────────────────────────────────────────────
with st.sidebar:
    st.markdown("### 🔬 Candidate Filters")
    novel_only  = st.checkbox("Novel variants only", value=True)
    impacts     = st.multiselect("VEP Impact", ["HIGH","MODERATE","LOW","MODIFIER"],
                                 default=["HIGH","MODERATE"])
    max_gnomad  = st.slider("Max gnomAD exome AF", 0.0, 0.05, 0.01,
                            step=0.001, format="%.3f")
    min_vqslod  = st.slider("Min VQSLOD", -5.0, 20.0, 2.0, step=0.5)
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
c1.metric("Total candidates",      len(cand))
c2.metric("🟡 HAMLET canonical",    len(hamlet_cand),
          help="Novel variants inside the 41 HAMLET genes")
c3.metric("🔴 Non-canonical",       len(noncan_cand),
          help="Novel variants OUTSIDE HAMLET — primary discovery targets")
c4.metric("HIGH impact",           (cand["IMPACT"] == "HIGH").sum())
c5.metric("MODERATE impact",       (cand["IMPACT"] == "MODERATE").sum())

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
                 title="Non-Canonical Genes with Novel Variants (by Impact)",
                 labels={"SYMBOL":"Gene","Count":"# Variants"})
    fig.update_layout(xaxis_tickangle=-45, margin=dict(t=40, b=10))
    st.plotly_chart(fig, use_container_width=True)

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
                st.write(f"HGVSc: `{row['HGVSc']}`")
                st.write(f"HGVSp: `{row['HGVSp']}`")
                st.write(f"Class: `{row['VARIANT_CLASS']}`")
                st.write(f"Biotype: `{row['BIOTYPE']}`")
            with col2:
                st.markdown("**Quality**")
                st.write(f"VQSLOD: `{row['VQSLOD']}`")
                st.write(f"QD: `{row['QD']}`")
                st.write(f"FS: `{row['FS']}`")
                st.write(f"MQ: `{row['MQ']}`")
            with col3:
                st.markdown("**Population / Clinical**")
                gnomad_str = (f"`{row['gnomADe_AF']:.6f}`"
                              if pd.notna(row["gnomADe_AF"])
                              else "**absent from gnomAD** ⚠️")
                st.markdown(f"gnomADe AF: {gnomad_str}")
                st.write(f"MAX_AF: `{row['MAX_AF'] if pd.notna(row['MAX_AF']) else 'N/A'}`")
                st.write(f"ClinSig: `{row['CLIN_SIG'] or 'not in ClinVar'}`")
                st.write(f"SIFT: `{row['SIFT']}`")
                st.write(f"PolyPhen: `{row['PolyPhen']}`")

            # Which samples carry this variant?
            carriers = sample_df[
                (sample_df["CHROM"] == row["CHROM"]) &
                (sample_df["POS"]   == row["POS"]) &
                (sample_df["GT_CLASS"].isin(["het","hom_alt"]))
            ][["sample","GT","GT_CLASS","DP","GQ"]].copy()

            if len(carriers):
                st.markdown(f"**Carriers: {len(carriers)}/{len(samples)} samples**")
                st.dataframe(carriers.reset_index(drop=True), use_container_width=True)
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
                  title="HAMLET Canonical Genes with Novel Variants",
                  labels={"SYMBOL":"Gene","Count":"# Variants"})
    fig2.update_layout(xaxis_tickangle=-45, margin=dict(t=40, b=10))
    st.plotly_chart(fig2, use_container_width=True)

    disp_cols = ["CHROM","POS","SYMBOL","CONSEQUENCE","IMPACT","HGVSc","HGVSp",
                 "VQSLOD","AF","gnomADe_AF","CLIN_SIG","SIFT","PolyPhen"]
    show_h = hamlet_cand[[c for c in disp_cols if c in hamlet_cand.columns]].copy()
    st.dataframe(
        show_h.reset_index(drop=True).style.applymap(
            lambda v: {"HIGH": "background-color:#ffd6d6",
                       "MODERATE": "background-color:#ffe5c0"}.get(v, ""),
            subset=["IMPACT"]
        ),
        use_container_width=True, height=350
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
