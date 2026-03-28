import streamlit as st
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go

df = st.session_state["df_view"]

st.title("🏠 Cohort Overview")
st.caption(f"31 samples · WES · GATK HaplotypeCaller joint germline · VQSR recalibrated · VEP 115")

# ── KPI Cards ────────────────────────────────────────────────────────────────
c1, c2, c3, c4, c5, c6 = st.columns(6)
total = len(df)
pass_n = (df["FILTER"] == "PASS").sum()
novel = df["IS_NOVEL"].sum()
snp = (df["VAR_TYPE"] == "SNP").sum()
indel = (df["VAR_TYPE"] == "INDEL").sum()
aml_hits = df["IS_AML_GENE"].sum()

c1.metric("Total Variants", f"{total:,}")
c2.metric("PASS", f"{pass_n:,}", f"{pass_n/total*100:.1f}%")
c3.metric("Novel", f"{novel:,}", f"{novel/total*100:.1f}%")
c4.metric("SNPs", f"{snp:,}")
c5.metric("INDELs", f"{indel:,}")
c6.metric("In AML Genes", f"{aml_hits:,}")

st.markdown("---")

col1, col2 = st.columns(2)

# ── FILTER breakdown ─────────────────────────────────────────────────────────
with col1:
    all_df = st.session_state["df"]
    filter_counts = all_df["FILTER"].value_counts().reset_index()
    filter_counts.columns = ["FILTER", "Count"]
    fig = px.pie(filter_counts, names="FILTER", values="Count",
                 title="FILTER Status",
                 color_discrete_sequence=px.colors.qualitative.Set2,
                 hole=0.4)
    fig.update_traces(textposition="outside", textinfo="percent+label")
    fig.update_layout(showlegend=False, margin=dict(t=40,b=0,l=0,r=0))
    st.plotly_chart(fig, use_container_width=True)

# ── SNP/INDEL + Novel/Known ──────────────────────────────────────────────────
with col2:
    cats = pd.DataFrame({
        "Category": ["SNP Known", "SNP Novel", "INDEL Known", "INDEL Novel"],
        "Count": [
            ((df["VAR_TYPE"]=="SNP") & ~df["IS_NOVEL"]).sum(),
            ((df["VAR_TYPE"]=="SNP") &  df["IS_NOVEL"]).sum(),
            ((df["VAR_TYPE"]=="INDEL") & ~df["IS_NOVEL"]).sum(),
            ((df["VAR_TYPE"]=="INDEL") &  df["IS_NOVEL"]).sum(),
        ],
        "Type": ["SNP","SNP","INDEL","INDEL"],
        "Status": ["Known","Novel","Known","Novel"],
    })
    fig2 = px.bar(cats, x="Type", y="Count", color="Status",
                  title="SNP / INDEL × Novel / Known",
                  color_discrete_map={"Novel":"#E74C3C","Known":"#85C1E9"},
                  barmode="stack")
    fig2.update_layout(margin=dict(t=40,b=0,l=0,r=0))
    st.plotly_chart(fig2, use_container_width=True)

# ── Variants per chromosome ──────────────────────────────────────────────────
st.subheader("Variant Density per Chromosome")
AML_CHROMS = {"chr5","chr7","chr11","chr17","chr8","chr21"}  # classic AML cytogenetics

chrom_order = [f"chr{i}" for i in list(range(1,23))+["X","Y"]]
chr_counts = df.groupby(["CHROM","VAR_TYPE"]).size().reset_index(name="Count")
chr_counts = chr_counts[chr_counts["CHROM"].isin(chrom_order)]
chr_counts["CHROM"] = pd.Categorical(chr_counts["CHROM"], categories=chrom_order, ordered=True)
chr_counts = chr_counts.sort_values("CHROM")

fig3 = px.bar(chr_counts, x="CHROM", y="Count", color="VAR_TYPE",
              color_discrete_map={"SNP":"#5DADE2","INDEL":"#F39C12"},
              title="Variants per Chromosome (SNP / INDEL)",
              labels={"CHROM":"Chromosome","Count":"# Variants"})

# Highlight AML-relevant chromosomes
for chrom in AML_CHROMS:
    if chrom in chrom_order:
        fig3.add_vrect(
            x0=chrom_order.index(chrom)-0.5,
            x1=chrom_order.index(chrom)+0.5,
            fillcolor="rgba(231,76,60,0.12)", line_width=0,
            annotation_text="AML", annotation_position="top",
            annotation_font_size=9,
        )
fig3.update_layout(margin=dict(t=50,b=0), xaxis_tickangle=-45)
st.plotly_chart(fig3, use_container_width=True)

# ── IMPACT breakdown ─────────────────────────────────────────────────────────
st.subheader("Variant Impact Distribution")
col3, col4 = st.columns(2)

with col3:
    impact_df = df["IMPACT"].value_counts().reset_index()
    impact_df.columns = ["IMPACT","Count"]
    color_map = {"HIGH":"#C0392B","MODERATE":"#E67E22","LOW":"#F1C40F","MODIFIER":"#95A5A6"}
    fig4 = px.pie(impact_df, names="IMPACT", values="Count",
                  color="IMPACT", color_discrete_map=color_map,
                  title="All Variants by VEP Impact", hole=0.4)
    fig4.update_traces(textposition="outside", textinfo="percent+label")
    fig4.update_layout(showlegend=False, margin=dict(t=40,b=0,l=0,r=0))
    st.plotly_chart(fig4, use_container_width=True)

with col4:
    novel_impact = df[df["IS_NOVEL"]]["IMPACT"].value_counts().reset_index()
    novel_impact.columns = ["IMPACT","Count"]
    fig5 = px.pie(novel_impact, names="IMPACT", values="Count",
                  color="IMPACT", color_discrete_map=color_map,
                  title="Novel Variants by VEP Impact", hole=0.4)
    fig5.update_traces(textposition="outside", textinfo="percent+label")
    fig5.update_layout(showlegend=False, margin=dict(t=40,b=0,l=0,r=0))
    st.plotly_chart(fig5, use_container_width=True)
