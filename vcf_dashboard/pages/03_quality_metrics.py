import streamlit as st
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go

df = st.session_state["df_view"]

st.title("📊 Quality Metrics (VQSR)")
st.markdown("Distribution of GATK/VQSR quality annotations. Novel variants should have comparable quality to known variants — if they skew worse, they may be artefacts.")

# ── Controls ─────────────────────────────────────────────────────────────────
col_ctrl1, col_ctrl2 = st.columns(2)
with col_ctrl1:
    split_by = st.selectbox("Colour/split by", ["IS_NOVEL","VAR_TYPE","IMPACT"])
with col_ctrl2:
    var_type_filter = st.multiselect("Variant Type", ["SNP","INDEL"], default=["SNP","INDEL"])

plot_df = df[df["VAR_TYPE"].isin(var_type_filter)].copy()
plot_df["IS_NOVEL"] = plot_df["IS_NOVEL"].map({True:"Novel",False:"Known"})

color_maps = {
    "IS_NOVEL":  {"Novel":"#E74C3C","Known":"#85C1E9"},
    "VAR_TYPE":  {"SNP":"#5DADE2","INDEL":"#F39C12"},
    "IMPACT":    {"HIGH":"#C0392B","MODERATE":"#E67E22","LOW":"#F1C40F","MODIFIER":"#95A5A6"},
}

cmap = color_maps.get(split_by, {})

def hist(col, title, xrange=None, log_x=False):
    sub = plot_df[plot_df[col].notna()]
    fig = px.histogram(sub, x=col, color=split_by,
                       nbins=60, barmode="overlay", opacity=0.7,
                       color_discrete_map=cmap,
                       title=title,
                       log_x=log_x)
    if xrange:
        fig.update_xaxes(range=xrange)
    fig.update_layout(margin=dict(t=40,b=0), legend_title_text=split_by)
    return fig

st.markdown("---")

# Row 1
c1, c2 = st.columns(2)
with c1:
    st.plotly_chart(hist("VQSLOD", "VQSLOD — VQSR log-odds score (higher = better)"), use_container_width=True)
with c2:
    st.plotly_chart(hist("QD", "QD — Quality by Depth (expect > 2)", xrange=[0,40]), use_container_width=True)

# Row 2
c3, c4 = st.columns(2)
with c3:
    st.plotly_chart(hist("FS", "FS — FisherStrand bias (expect < 60 SNP / < 200 INDEL)", xrange=[0,100]), use_container_width=True)
with c4:
    st.plotly_chart(hist("MQ", "MQ — Mapping Quality (expect ~60)", xrange=[40,70]), use_container_width=True)

# Row 3
c5, c6 = st.columns(2)
with c5:
    st.plotly_chart(hist("SOR", "SOR — Strand Odds Ratio (expect < 3)", xrange=[0,10]), use_container_width=True)
with c6:
    st.plotly_chart(hist("InbreedingCoeff", "InbreedingCoeff (expect near 0)"), use_container_width=True)

# ── VQSLOD scatter Novel vs Known ─────────────────────────────────────────────
st.markdown("---")
st.subheader("VQSLOD: Novel vs Known (violin)")
sub2 = plot_df[plot_df["VQSLOD"].notna()].copy()
fig_vio = px.violin(sub2, x="IS_NOVEL", y="VQSLOD", color="IS_NOVEL",
                    color_discrete_map={"Novel":"#E74C3C","Known":"#85C1E9"},
                    box=True, points="outliers",
                    title="VQSLOD Distribution: Novel vs Known Variants")
fig_vio.update_layout(showlegend=False, margin=dict(t=40,b=0))
st.plotly_chart(fig_vio, use_container_width=True)

# ── Culprit annotation ────────────────────────────────────────────────────────
st.markdown("---")
st.subheader("VQSR Culprit Annotation")
all_df = st.session_state["df"]
if "culprit" in all_df.columns:
    culp = all_df["culprit"].value_counts().reset_index()
    culp.columns = ["Culprit","Count"]
    fig_c = px.bar(culp, x="Culprit", y="Count",
                   title="Most common VQSR culprit annotation",
                   color="Count", color_continuous_scale="Reds")
    fig_c.update_layout(margin=dict(t=40,b=0))
    st.plotly_chart(fig_c, use_container_width=True)
else:
    # derive from INFO
    import re
    culprit_vals = all_df.apply(lambda _: None, axis=1)  # placeholder
    st.info("Culprit field not parsed — add 'culprit' to vcf_loader.py INFO fields if needed.")
