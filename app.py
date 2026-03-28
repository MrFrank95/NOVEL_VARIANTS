import streamlit as st
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
import sys, os

sys.path.insert(0, os.path.dirname(__file__))
from vcf_loader import load_vcf, get_per_sample_stats, HAMLET_GENES

st.set_page_config(
    page_title="AML Variant Dashboard",
    page_icon="🧬",
    layout="wide",
)

# ── Data loading ──────────────────────────────────────────────────────────────
@st.cache_data(show_spinner="Parsing VCF and computing sample statistics… (~30s first run)")
def load_data():
    vcf_path = os.path.join(os.path.dirname(__file__), "data",
                            "joint_germline_recalibrated_VEP_ann.vcf")
    df, samples = load_vcf(vcf_path)
    sdf = get_per_sample_stats(df, samples)
    return df, sdf, samples

df_all, sdf_all, samples = load_data()

# ── Global sidebar filter ─────────────────────────────────────────────────────
with st.sidebar:
    st.title("🧬 AML Variant Dashboard")
    st.caption("nf-core/Sarek 3.8.1 · GATK 4.6 · VEP 115 · GRCh38\n31 samples · WES · Joint germline")
    st.markdown("---")
    st.markdown("### Global Filter")
    pass_only = st.checkbox("PASS variants only", value=True)
    st.markdown("---")
    st.markdown(f"**HAMLET panel:** {len(HAMLET_GENES)} genes")
    st.markdown(", ".join(sorted(HAMLET_GENES)))

df  = df_all[df_all["FILTER"] == "PASS"].copy() if pass_only else df_all.copy()
sdf = sdf_all[sdf_all["FILTER"] == "PASS"].copy() if pass_only else sdf_all.copy()

# ── Tabs ──────────────────────────────────────────────────────────────────────
tab1, tab2, tab3 = st.tabs([
    "🏠  Overview",
    "🔴  Novel AML Candidates",
    "🔍  Variant Browser",
])

# ══════════════════════════════════════════════════════════════════════════════
# TAB 1 — OVERVIEW
# ══════════════════════════════════════════════════════════════════════════════
with tab1:
    st.header("Cohort Overview")

    total   = len(df)
    novel   = df["IS_NOVEL"].sum()
    snp     = (df["VAR_TYPE"] == "SNP").sum()
    indel   = (df["VAR_TYPE"] == "INDEL").sum()
    hamlet  = (df["GENE_TIER"] == "HAMLET canonical").sum()
    noncan  = (df["GENE_TIER"] == "Non-canonical").sum()

    c1,c2,c3,c4,c5,c6 = st.columns(6)
    c1.metric("Total variants",         f"{total:,}")
    c2.metric("Novel (no pop. DB hit)", f"{novel:,}", f"{novel/total*100:.1f}%")
    c3.metric("SNPs",                   f"{snp:,}")
    c4.metric("INDELs",                 f"{indel:,}")
    c5.metric("In HAMLET genes",        f"{hamlet:,}")
    c6.metric("Non-canonical genes",    f"{noncan:,}")

    st.markdown("---")
    col1, col2 = st.columns(2)

    with col1:
        # FILTER breakdown (always from full data)
        fc = df_all["FILTER"].value_counts().reset_index()
        fc.columns = ["FILTER","Count"]
        fig = px.pie(fc, names="FILTER", values="Count", title="FILTER Status",
                     hole=0.4, color_discrete_sequence=px.colors.qualitative.Set2)
        fig.update_traces(textposition="outside", textinfo="percent+label")
        fig.update_layout(showlegend=False, margin=dict(t=40,b=0,l=0,r=0))
        st.plotly_chart(fig, use_container_width=True)

    with col2:
        cats = pd.DataFrame({
            "Category": ["SNP Known","SNP Novel","INDEL Known","INDEL Novel"],
            "Count": [
                ((df["VAR_TYPE"]=="SNP")   & ~df["IS_NOVEL"]).sum(),
                ((df["VAR_TYPE"]=="SNP")   &  df["IS_NOVEL"]).sum(),
                ((df["VAR_TYPE"]=="INDEL") & ~df["IS_NOVEL"]).sum(),
                ((df["VAR_TYPE"]=="INDEL") &  df["IS_NOVEL"]).sum(),
            ],
            "Status": ["Known","Novel","Known","Novel"],
        })
        fig2 = px.bar(cats, x="Category", y="Count", color="Status",
                      color_discrete_map={"Novel":"#E74C3C","Known":"#85C1E9"},
                      title="SNP / INDEL × Novel / Known", barmode="stack")
        fig2.update_layout(margin=dict(t=40,b=0))
        st.plotly_chart(fig2, use_container_width=True)

    # Variants per chromosome
    st.subheader("Variant Density per Chromosome")
    AML_CHROMS = {"chr5","chr7","chr11","chr17","chr8","chr21"}
    chrom_order = [f"chr{i}" for i in list(range(1,23))+["X"]]
    cc = df.groupby(["CHROM","VAR_TYPE"]).size().reset_index(name="Count")
    cc = cc[cc["CHROM"].isin(chrom_order)]
    cc["CHROM"] = pd.Categorical(cc["CHROM"], categories=chrom_order, ordered=True)
    cc = cc.sort_values("CHROM")
    fig3 = px.bar(cc, x="CHROM", y="Count", color="VAR_TYPE",
                  color_discrete_map={"SNP":"#5DADE2","INDEL":"#F39C12"},
                  title="Variants per Chromosome",
                  labels={"CHROM":"Chromosome"})
    for ch in AML_CHROMS:
        if ch in chrom_order:
            fig3.add_vrect(x0=chrom_order.index(ch)-.5, x1=chrom_order.index(ch)+.5,
                           fillcolor="rgba(231,76,60,0.12)", line_width=0,
                           annotation_text="AML", annotation_position="top",
                           annotation_font_size=9)
    fig3.update_layout(xaxis_tickangle=-45, margin=dict(t=50,b=0))
    st.plotly_chart(fig3, use_container_width=True)

    # VEP Impact
    col3, col4 = st.columns(2)
    color_map_imp = {"HIGH":"#C0392B","MODERATE":"#E67E22","LOW":"#F1C40F","MODIFIER":"#95A5A6"}
    with col3:
        id1 = df["IMPACT"].value_counts().reset_index()
        id1.columns = ["IMPACT","Count"]
        f4 = px.pie(id1, names="IMPACT", values="Count", color="IMPACT",
                    color_discrete_map=color_map_imp, title="All Variants by Impact", hole=0.4)
        f4.update_traces(textposition="outside", textinfo="percent+label")
        f4.update_layout(showlegend=False, margin=dict(t=40,b=0,l=0,r=0))
        st.plotly_chart(f4, use_container_width=True)
    with col4:
        id2 = df[df["IS_NOVEL"]]["IMPACT"].value_counts().reset_index()
        id2.columns = ["IMPACT","Count"]
        f5 = px.pie(id2, names="IMPACT", values="Count", color="IMPACT",
                    color_discrete_map=color_map_imp, title="Novel Variants by Impact", hole=0.4)
        f5.update_traces(textposition="outside", textinfo="percent+label")
        f5.update_layout(showlegend=False, margin=dict(t=40,b=0,l=0,r=0))
        st.plotly_chart(f5, use_container_width=True)


# ══════════════════════════════════════════════════════════════════════════════
# TAB 2 — NOVEL AML CANDIDATES
# ══════════════════════════════════════════════════════════════════════════════
with tab2:
    st.header("Novel AML Variant Candidates")
    st.markdown("""
Variants classified into two tiers:
- 🔴 **Non-canonical** — novel variant in a gene **outside** the HAMLET panel → *primary discovery target*
- 🟡 **HAMLET canonical** — novel variant in one of the 40 HAMLET leukemia genes → *confirmatory*

*Novel = absent from `Existing_variation` in the canonical VEP transcript.*
""")

    # ── Candidate filters ─────────────────────────────────────────────────────
    with st.expander("🔧 Candidate Filters", expanded=True):
        fc1, fc2, fc3, fc4 = st.columns(4)
        with fc1:
            novel_only = st.checkbox("Novel variants only", value=True)
        with fc2:
            sel_impacts = st.multiselect("VEP Impact",
                ["HIGH","MODERATE","LOW","MODIFIER"], default=["HIGH","MODERATE"])
        with fc3:
            max_gnomad = st.slider("Max gnomAD-e AF", 0.0, 0.05, 0.01,
                                   step=0.001, format="%.3f")
        with fc4:
            min_vqslod = st.slider("Min VQSLOD", -5.0, 20.0, 2.0, step=0.5)

    cand = df.copy()
    if novel_only:
        cand = cand[cand["IS_NOVEL"]]
    if sel_impacts:
        cand = cand[cand["IMPACT"].isin(sel_impacts)]
    cand = cand[cand["VQSLOD"] >= min_vqslod]
    cand = cand[(cand["gnomADe_AF"].isna()) | (cand["gnomADe_AF"] <= max_gnomad)]
    cand = cand[cand["GENE_TIER"].isin(["HAMLET canonical","Non-canonical"])]

    hamlet_cand = cand[cand["GENE_TIER"] == "HAMLET canonical"]
    noncan_cand = cand[cand["GENE_TIER"] == "Non-canonical"]

    # ── Summary numbers ───────────────────────────────────────────────────────
    k1,k2,k3,k4,k5,k6 = st.columns(6)
    k1.metric("Total candidates",          len(cand))
    k2.metric("🔴 Non-canonical",          len(noncan_cand))
    k3.metric("  └ Unique genes",          noncan_cand["SYMBOL"].nunique())
    k4.metric("🟡 HAMLET canonical",       len(hamlet_cand))
    k5.metric("  └ Unique HAMLET genes",   hamlet_cand["SYMBOL"].nunique())
    k6.metric("HIGH impact total",         (cand["IMPACT"]=="HIGH").sum())

    st.markdown("---")

    # ══════════════════════════════════════════════════════════════════════════
    # NON-CANONICAL SECTION
    # ══════════════════════════════════════════════════════════════════════════
    st.subheader("🔴 Non-Canonical Novel Candidates — Primary Discovery Targets")

    if len(noncan_cand) == 0:
        st.info("No non-canonical candidates with current filters.")
    else:
        # Gene-level summary table
        st.markdown("#### Gene-Level Summary")
        gene_summary = (
            noncan_cand.groupby("SYMBOL")
            .agg(
                Variants      = ("POS", "count"),
                HIGH          = ("IMPACT", lambda x: (x=="HIGH").sum()),
                MODERATE      = ("IMPACT", lambda x: (x=="MODERATE").sum()),
                LOW           = ("IMPACT", lambda x: (x=="LOW").sum()),
                MODIFIER      = ("IMPACT", lambda x: (x=="MODIFIER").sum()),
                In_gnomAD     = ("gnomADe_AF", lambda x: x.notna().sum()),
                Absent_gnomAD = ("gnomADe_AF", lambda x: x.isna().sum()),
                Best_VQSLOD   = ("VQSLOD", "max"),
            )
            .reset_index()
            .sort_values(["HIGH","MODERATE","Variants"], ascending=False)
        )
        gene_summary.columns = [
            "Gene","Variants","HIGH","MODERATE","LOW","MODIFIER",
            "In gnomAD","Not in gnomAD","Best VQSLOD"
        ]

        # Biological context annotations
        BIO_CTX = {
            "ASXL2":   "ASXL family (like ASXL1) — Polycomb complex",
            "ASXL3":   "ASXL family — Polycomb complex",
            "TET3":    "TET family (like TET2) — DNA demethylation",
            "SETD2":   "H3K36 methyltransferase — epigenetic regulator",
            "RAD21":   "Cohesin complex (like STAG2)",
            "NUP98":   "Nucleoporin — known AML translocation partner",
            "DOT1L":   "H3K79 methyltransferase — MLL-AML therapeutic target",
            "BRD4":    "BET bromodomain — epigenetic reader / AML target",
            "KDM2B":   "H3K36 demethylase — Polycomb recruiter",
            "IKZF3":   "Ikaros TF — haematopoietic development",
            "SMARCA2":  "SWI/SNF chromatin remodelling",
            "ARID1B":  "SWI/SNF chromatin remodelling",
            "DHX15":   "RNA helicase — splicing factor",
            "FOXP1":   "TF — myeloid/B-cell development",
            "NOTCH2":  "NOTCH signalling",
            "KLF2":    "TF — HSC regulation",
            "AURKB":   "Aurora kinase B — mitosis",
            "AXL":     "Receptor TK — drug resistance in AML",
            "MGA":     "MYC network — MAX binding partner",
            "CHD8":    "Chromodomain helicase remodeller",
            "PTPRO":   "Protein tyrosine phosphatase — tumour suppressor",
        }
        gene_summary["Biological Context"] = gene_summary["Gene"].map(BIO_CTX).fillna("")
        st.dataframe(gene_summary, use_container_width=True, height=300,
                     column_config={
                         "HIGH":     st.column_config.NumberColumn("HIGH",     format="%d"),
                         "MODERATE": st.column_config.NumberColumn("MODERATE", format="%d"),
                         "Best VQSLOD": st.column_config.NumberColumn("Best VQSLOD", format="%.2f"),
                     })
        st.caption(f"**{len(noncan_cand):,} novel variants** across **{noncan_cand['SYMBOL'].nunique()} unique non-canonical genes** "
                   f"({noncan_cand[noncan_cand['IMPACT'].isin(['HIGH','MODERATE'])]['SYMBOL'].nunique()} with HIGH/MOD impact)")

        # Chart
        col_ch1, col_ch2 = st.columns([2,1])
        with col_ch1:
            hi_mod_genes = noncan_cand[noncan_cand["IMPACT"].isin(["HIGH","MODERATE"])]
            if len(hi_mod_genes):
                gc = hi_mod_genes.groupby(["SYMBOL","IMPACT"]).size().reset_index(name="Count")
                color_map = {"HIGH":"#C0392B","MODERATE":"#E67E22"}
                fig = px.bar(gc, x="SYMBOL", y="Count", color="IMPACT",
                             color_discrete_map=color_map,
                             title="Non-Canonical Genes — HIGH / MODERATE Impact Variants",
                             labels={"SYMBOL":"Gene"})
                fig.update_layout(xaxis_tickangle=-45, margin=dict(t=40,b=0))
                st.plotly_chart(fig, use_container_width=True)

        with col_ch2:
            absent = noncan_cand["gnomADe_AF"].isna().sum()
            rare   = ((noncan_cand["gnomADe_AF"].notna()) & (noncan_cand["gnomADe_AF"] < 0.001)).sum()
            common = (noncan_cand["gnomADe_AF"] >= 0.001).sum()
            pop_df = pd.DataFrame({
                "Category": ["Absent from gnomAD","Rare (AF<0.001)","AF≥0.001"],
                "Count":    [absent, rare, common],
            })
            fig_p = px.pie(pop_df, names="Category", values="Count",
                           color="Category",
                           color_discrete_map={
                               "Absent from gnomAD":"#C0392B",
                               "Rare (AF<0.001)":   "#E67E22",
                               "AF≥0.001":          "#85C1E9",
                           },
                           title="gnomAD-e Frequency Category", hole=0.4)
            fig_p.update_traces(textposition="outside", textinfo="percent+label")
            fig_p.update_layout(showlegend=False, margin=dict(t=40,b=0,l=0,r=0))
            st.plotly_chart(fig_p, use_container_width=True)

        # Full variant table
        st.markdown("#### Complete Non-Canonical Variant Table")
        nc_cols = ["CHROM","POS","REF","ALT","SYMBOL","CONSEQUENCE","IMPACT",
                   "HGVSc","HGVSp","VQSLOD","AF","gnomADe_AF","MAX_AF",
                   "CLIN_SIG","SIFT","PolyPhen","EXISTING_VARIATION"]
        nc_show = noncan_cand[[c for c in nc_cols if c in noncan_cand.columns]].copy()
        nc_show = nc_show.sort_values(["IMPACT","VQSLOD"], ascending=[True,False]).reset_index(drop=True)
        st.dataframe(
            nc_show.style.applymap(
                lambda v: {"HIGH":"background-color:#ffd6d6",
                           "MODERATE":"background-color:#ffe5c0"}.get(v,""),
                subset=["IMPACT"]
            ),
            use_container_width=True, height=400,
            column_config={
                "gnomADe_AF": st.column_config.NumberColumn("gnomAD-e AF", format="%.5f"),
                "VQSLOD":     st.column_config.NumberColumn("VQSLOD",      format="%.2f"),
                "AF":         st.column_config.NumberColumn("Cohort AF",   format="%.3f"),
            }
        )
        st.download_button("⬇️ Download non-canonical candidates CSV",
                           nc_show.to_csv(index=False).encode(),
                           "noncanonical_novel_candidates.csv", mime="text/csv")

        # Expandable per-variant cards (HIGH only to keep it clean)
        high_nc = noncan_cand[noncan_cand["IMPACT"] == "HIGH"].sort_values("VQSLOD", ascending=False)
        if len(high_nc):
            st.markdown("#### 🔬 HIGH Impact — Detail Cards")
            for _, row in high_nc.iterrows():
                ctx = BIO_CTX.get(row["SYMBOL"], "")
                with st.expander(
                    f"**{row['SYMBOL']}** · {row['CONSEQUENCE']} · "
                    f"`{row['CHROM']}:{row['POS']} {row['REF']}>{row['ALT']}`"
                    + (f"  — _{ctx}_" if ctx else "")
                ):
                    col1, col2, col3 = st.columns(3)
                    with col1:
                        st.markdown("**Variant**")
                        st.write(f"HGVSc: `{row['HGVSc']}`")
                        st.write(f"HGVSp: `{row['HGVSp']}`")
                        st.write(f"Class: `{row['VARIANT_CLASS']}`")
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
                                      else "**⚠️ absent from gnomAD**")
                        st.markdown(f"gnomADe AF: {gnomad_str}")
                        st.write(f"ClinSig: `{row['CLIN_SIG'] or 'not in ClinVar'}`")
                        st.write(f"SIFT: `{row['SIFT']}`")
                        st.write(f"PolyPhen: `{row['PolyPhen']}`")

                    # Carriers
                    carriers = sdf[
                        (sdf["CHROM"] == row["CHROM"]) &
                        (sdf["POS"]   == row["POS"]) &
                        (sdf["GT_CLASS"].isin(["het","hom_alt"]))
                    ][["sample","GT","GT_CLASS","DP","GQ"]].copy()
                    if len(carriers):
                        st.markdown(f"**Carriers: {len(carriers)}/{len(samples)} samples**")
                        st.dataframe(carriers.reset_index(drop=True), use_container_width=True)
                    else:
                        st.info("No carriers in sample data.")
                    if ctx:
                        st.info(f"💡 {ctx}")

    st.markdown("---")

    # ══════════════════════════════════════════════════════════════════════════
    # HAMLET CANONICAL SECTION
    # ══════════════════════════════════════════════════════════════════════════
    st.subheader("🟡 HAMLET Canonical Novel Variants — Confirmatory")
    st.caption("Novel variants inside the 40 HAMLET genes. Confirms known biology; "
               "may represent new alleles of established AML drivers.")

    if len(hamlet_cand) == 0:
        st.info("No HAMLET canonical candidates with current filters.")
    else:
        # Gene-level summary
        st.markdown("#### Gene-Level Summary")
        hamlet_gene_sum = (
            hamlet_cand.groupby("SYMBOL")
            .agg(
                Variants  = ("POS","count"),
                HIGH      = ("IMPACT", lambda x: (x=="HIGH").sum()),
                MODERATE  = ("IMPACT", lambda x: (x=="MODERATE").sum()),
                LOW       = ("IMPACT", lambda x: (x=="LOW").sum()),
                MODIFIER  = ("IMPACT", lambda x: (x=="MODIFIER").sum()),
                Best_VQSLOD = ("VQSLOD","max"),
            )
            .reset_index()
            .sort_values(["HIGH","MODERATE","Variants"], ascending=False)
        )
        hamlet_gene_sum.columns = [
            "Gene","Variants","HIGH","MODERATE","LOW","MODIFIER","Best VQSLOD"
        ]
        st.dataframe(hamlet_gene_sum, use_container_width=True, height=280,
                     column_config={
                         "Best VQSLOD": st.column_config.NumberColumn("Best VQSLOD", format="%.2f"),
                     })
        st.caption(f"**{len(hamlet_cand):,} novel variants** across "
                   f"**{hamlet_cand['SYMBOL'].nunique()} of 40 HAMLET genes** "
                   f"({hamlet_cand[hamlet_cand['IMPACT'].isin(['HIGH','MODERATE'])]['SYMBOL'].nunique()} with HIGH/MOD impact)")

        # Bar chart
        hg = hamlet_cand.groupby(["SYMBOL","IMPACT"]).size().reset_index(name="Count")
        color_map = {"HIGH":"#C0392B","MODERATE":"#E67E22","LOW":"#F1C40F","MODIFIER":"#95A5A6"}
        fig2 = px.bar(hg, x="SYMBOL", y="Count", color="IMPACT",
                      color_discrete_map=color_map,
                      title="HAMLET Canonical Genes — Novel Variants by Impact",
                      labels={"SYMBOL":"Gene"})
        fig2.update_layout(xaxis_tickangle=-45, margin=dict(t=40,b=0))
        st.plotly_chart(fig2, use_container_width=True)

        # Full table
        st.markdown("#### Complete HAMLET Canonical Variant Table")
        h_cols = ["CHROM","POS","REF","ALT","SYMBOL","CONSEQUENCE","IMPACT",
                  "HGVSc","HGVSp","VQSLOD","AF","gnomADe_AF","CLIN_SIG","SIFT","PolyPhen"]
        h_show = hamlet_cand[[c for c in h_cols if c in hamlet_cand.columns]].copy()
        h_show = h_show.sort_values(["SYMBOL","IMPACT"], ascending=[True,True]).reset_index(drop=True)
        st.dataframe(
            h_show.style.applymap(
                lambda v: {"HIGH":"background-color:#ffd6d6",
                           "MODERATE":"background-color:#ffe5c0"}.get(v,""),
                subset=["IMPACT"]
            ),
            use_container_width=True, height=380,
            column_config={
                "gnomADe_AF": st.column_config.NumberColumn("gnomAD-e AF", format="%.5f"),
                "VQSLOD":     st.column_config.NumberColumn("VQSLOD",      format="%.2f"),
            }
        )
        st.download_button("⬇️ Download HAMLET canonical candidates CSV",
                           h_show.to_csv(index=False).encode(),
                           "hamlet_canonical_candidates.csv", mime="text/csv")

    st.markdown("---")

    # ── Combined download ─────────────────────────────────────────────────────
    all_cols = ["CHROM","POS","REF","ALT","SYMBOL","GENE_TIER","CONSEQUENCE","IMPACT",
                "HGVSc","HGVSp","VQSLOD","AF","gnomADe_AF","MAX_AF",
                "CLIN_SIG","SIFT","PolyPhen","IS_NOVEL","EXISTING_VARIATION"]
    all_dl = cand[[c for c in all_cols if c in cand.columns]]
    st.download_button("⬇️ Download ALL candidates (both tiers) CSV",
                       all_dl.to_csv(index=False).encode(),
                       "all_novel_candidates.csv", mime="text/csv")


# ══════════════════════════════════════════════════════════════════════════════
# TAB 3 — VARIANT BROWSER
# ══════════════════════════════════════════════════════════════════════════════
with tab3:
    st.header("Variant Browser")
    st.markdown("Full searchable variant table. Expand any variant below for per-sample genotype detail.")

    # Filters
    with st.expander("🔧 Filters", expanded=True):
        bc1,bc2,bc3,bc4 = st.columns(4)
        with bc1:
            chrom_opts = ["All"] + sorted(
                df["CHROM"].unique(),
                key=lambda x: int(x.replace("chr","").replace("X","23").replace("Y","24"))
            )
            sel_chrom = st.selectbox("Chromosome", chrom_opts)
        with bc2:
            sel_impact = st.selectbox("Impact", ["All"] + sorted(df["IMPACT"].dropna().unique()))
        with bc3:
            sel_novel = st.selectbox("Novelty", ["All","Novel only","Known only"])
        with bc4:
            sel_tier = st.selectbox("Gene Tier",
                ["All","HAMLET canonical","Non-canonical","Intergenic"])
        bc5, bc6 = st.columns(2)
        with bc5:
            gene_q = st.text_input("Gene symbol", placeholder="e.g. EZH2")
        with bc6:
            consq_q = st.text_input("Consequence", placeholder="e.g. frameshift")

    bdf = df.copy()
    if sel_chrom != "All":
        bdf = bdf[bdf["CHROM"] == sel_chrom]
    if sel_impact != "All":
        bdf = bdf[bdf["IMPACT"] == sel_impact]
    if sel_novel == "Novel only":
        bdf = bdf[bdf["IS_NOVEL"]]
    elif sel_novel == "Known only":
        bdf = bdf[~bdf["IS_NOVEL"]]
    if sel_tier != "All":
        bdf = bdf[bdf["GENE_TIER"] == sel_tier]
    if gene_q:
        bdf = bdf[bdf["SYMBOL"].str.contains(gene_q, case=False, na=False)]
    if consq_q:
        bdf = bdf[bdf["CONSEQUENCE"].str.contains(consq_q, case=False, na=False)]

    st.caption(f"Showing **{len(bdf):,}** of {len(df):,} variants")

    disp = ["CHROM","POS","ID","REF","ALT","SYMBOL","GENE_TIER","CONSEQUENCE","IMPACT",
            "FILTER","IS_NOVEL","HGVSc","HGVSp","AF","VQSLOD","QD","FS","MQ",
            "gnomADe_AF","MAX_AF","CLIN_SIG","SIFT","PolyPhen","EXISTING_VARIATION"]
    show_b = bdf[[c for c in disp if c in bdf.columns]].reset_index(drop=True)

    st.dataframe(show_b, use_container_width=True, height=450,
                 column_config={
                     "IS_NOVEL":   st.column_config.CheckboxColumn("Novel"),
                     "GENE_TIER":  st.column_config.TextColumn("Gene Tier"),
                     "VQSLOD":     st.column_config.NumberColumn("VQSLOD",    format="%.2f"),
                     "gnomADe_AF": st.column_config.NumberColumn("gnomAD-e AF", format="%.4f"),
                     "AF":         st.column_config.NumberColumn("Cohort AF", format="%.3f"),
                 })

    st.download_button("⬇️ Download filtered table CSV",
                       show_b.to_csv(index=False).encode(),
                       "variants_filtered.csv", mime="text/csv")

    # Per-variant genotype drill-down
    st.markdown("---")
    st.subheader("🔎 Per-Variant Genotype Drill-down")
    variant_labels = show_b.apply(
        lambda r: f"{r['CHROM']}:{r['POS']}  {r['REF']}>{r['ALT']}  [{r.get('SYMBOL','')}]  {r.get('CONSEQUENCE','')}",
        axis=1
    ).tolist()

    if variant_labels:
        sel_var = st.selectbox("Select variant", variant_labels[:300])
        idx = variant_labels.index(sel_var)
        row = show_b.iloc[idx]
        var_sdf = sdf[(sdf["CHROM"] == row["CHROM"]) & (sdf["POS"] == row["POS"])]
        if len(var_sdf):
            gt_colors = {"hom_ref":"#d5e8d4","het":"#ffe6cc","hom_alt":"#f8cecc","missing":"#f0f0f0"}
            vsd = var_sdf[["sample","GT","GT_CLASS","DP","GQ"]].reset_index(drop=True)
            st.dataframe(
                vsd.style.applymap(lambda v: f"background-color:{gt_colors.get(v,'')}", subset=["GT_CLASS"]),
                use_container_width=True
            )
        else:
            st.info("No per-sample data for this variant.")
    else:
        st.info("No variants to display with current filters.")
