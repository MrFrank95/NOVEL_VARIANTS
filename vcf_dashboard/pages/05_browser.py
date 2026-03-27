import streamlit as st
import pandas as pd

df = st.session_state["df_view"]
samples = st.session_state["samples"]

st.title("🔍 Variant Browser")
st.markdown("Full variant table with interactive filters. Expand any row for per-sample genotype breakdown.")

# ── Filters ──────────────────────────────────────────────────────────────────
with st.expander("🔧 Filters", expanded=True):
    col1,col2,col3,col4 = st.columns(4)
    with col1:
        chroms = ["All"] + sorted(df["CHROM"].unique().tolist(),
                  key=lambda x: int(x.replace("chr","").replace("X","23").replace("Y","24")))
        sel_chrom = st.selectbox("Chromosome", chroms)
    with col2:
        impacts = ["All"] + sorted(df["IMPACT"].dropna().unique().tolist())
        sel_impact = st.selectbox("VEP Impact", impacts)
    with col3:
        novel_sel = st.selectbox("Novelty", ["All","Novel only","Known only"])
    with col4:
        aml_sel = st.checkbox("AML genes only", value=False)

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

st.dataframe(show, use_container_width=True, height=450,
             column_config={
                 "IS_NOVEL":  st.column_config.CheckboxColumn("Novel"),
                 "IS_AML_GENE": st.column_config.CheckboxColumn("AML Gene"),
                 "VQSLOD":    st.column_config.NumberColumn("VQSLOD", format="%.2f"),
                 "gnomADe_AF":st.column_config.NumberColumn("gnomAD-e AF", format="%.4f"),
                 "AF":        st.column_config.NumberColumn("Cohort AF", format="%.3f"),
             })

csv = show.to_csv(index=False).encode()
st.download_button("⬇️ Download filtered table as CSV", csv,
                   "variants_filtered.csv", mime="text/csv")

# ── Per-variant sample genotype drill-down ────────────────────────────────────
st.markdown("---")
st.subheader("🔎 Per-Variant Genotype Drill-down")

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
            use_container_width=True
        )
    else:
        st.info("No per-sample data available for this variant.")
else:
    st.info("No variants to display with current filters.")
