import streamlit as st
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go

df      = st.session_state["df_view"]
sdf     = st.session_state["sample_df_view"]
samples = st.session_state["samples"]

from vcf_loader import HAMLET_GENES

st.title("👥 Per-Sample Burden Analysis")
st.markdown("How many (novel) variants does each sample carry, split by gene tier?")

# ── Per-sample summary ────────────────────────────────────────────────────────
summary = []
for sname in samples:
    s      = sdf[sdf["sample"] == sname]
    called = s[s["GT_CLASS"].isin(["het","hom_alt"])]
    het    = s[s["GT_CLASS"] == "het"]
    hom    = s[s["GT_CLASS"] == "hom_alt"]

    novel       = called[called["IS_NOVEL"]]
    hamlet_var  = called[called["IS_AML_GENE"]]
    noncan_nov  = called[(called["IS_NOVEL"]) &
                         (~called["SYMBOL"].isin(HAMLET_GENES)) &
                         (called["SYMBOL"] != "")]
    hi_mod      = called[called["IMPACT"].isin(["HIGH","MODERATE"])]

    summary.append({
        "Sample":                     sname,
        "Total called":               len(called),
        "Het":                        len(het),
        "Hom alt":                    len(hom),
        "Het/Hom ratio":              round(len(het) / max(len(hom), 1), 2),
        "Novel total":                len(novel),
        "🔴 Novel non-canonical H/M": len(noncan_nov[noncan_nov["IMPACT"].isin(["HIGH","MODERATE"])]),
        "🟡 In HAMLET genes":          len(hamlet_var),
        "HIGH/MOD impact":            len(hi_mod),
        "Mean GQ":                    round(s["GQ"].mean(), 1) if s["GQ"].notna().any() else None,
        "Mean DP":                    round(s["DP"].mean(), 1) if s["DP"].notna().any() else None,
    })

sum_df = pd.DataFrame(summary).set_index("Sample")

st.subheader("Per-Sample Summary Table")
st.dataframe(
    sum_df.style.background_gradient(
        subset=["🔴 Novel non-canonical H/M","🟡 In HAMLET genes"], cmap="YlOrRd"
    ),
    use_container_width=True,
)
st.download_button("⬇️ Download table",
                   sum_df.reset_index().to_csv(index=False).encode(),
                   "per_sample_summary.csv", mime="text/csv")

st.markdown("---")

col1, col2 = st.columns(2)
with col1:
    fig = px.bar(sum_df.reset_index(), x="Sample",
                 y="🔴 Novel non-canonical H/M",
                 color="🔴 Novel non-canonical H/M",
                 color_continuous_scale="Reds",
                 title="🔴 Novel Non-Canonical HIGH/MOD Variants per Sample")
    fig.update_layout(xaxis_tickangle=-60, margin=dict(t=40,b=0))
    st.plotly_chart(fig, use_container_width=True)

with col2:
    fig2 = px.bar(sum_df.reset_index(), x="Sample",
                  y="🟡 In HAMLET genes",
                  color="🟡 In HAMLET genes",
                  color_continuous_scale="YlOrBr",
                  title="🟡 Variants in HAMLET Canonical Genes per Sample")
    fig2.update_layout(xaxis_tickangle=-60, margin=dict(t=40,b=0))
    st.plotly_chart(fig2, use_container_width=True)

# ── Het/Hom ratio ─────────────────────────────────────────────────────────────
st.subheader("Het / Hom Alt Ratio (QC — expect ~2.0–2.5 for WES germline)")
fig3 = px.bar(sum_df.reset_index(), x="Sample", y="Het/Hom ratio",
              color="Het/Hom ratio", color_continuous_scale="RdYlGn")
fig3.add_hline(y=2.0, line_dash="dash", line_color="grey",
               annotation_text="expected lower bound")
fig3.add_hline(y=2.5, line_dash="dash", line_color="grey",
               annotation_text="expected upper bound")
fig3.update_layout(xaxis_tickangle=-60, margin=dict(t=40,b=0))
st.plotly_chart(fig3, use_container_width=True)

st.markdown("---")

# ── Heatmap: sample × gene ────────────────────────────────────────────────────
st.subheader("🧬 Sample × Gene Heatmap")

tab1, tab2 = st.tabs(["🔴 Non-Canonical Novel", "🟡 HAMLET Canonical"])

def make_heatmap(data, title):
    if len(data) == 0:
        st.info("No data for current filters.")
        return
    impact_rank = {"HIGH":3,"MODERATE":2,"LOW":1,"MODIFIER":0}
    data = data.copy()
    data["impact_rank"] = data["IMPACT"].map(impact_rank).fillna(0)
    pivot = data.pivot_table(index="sample", columns="SYMBOL",
                             values="impact_rank", aggfunc="max").fillna(-1)
    gene_order = pivot.max().sort_values(ascending=False).index.tolist()
    pivot = pivot[gene_order]
    fig_hm = go.Figure(go.Heatmap(
        z=pivot.values, x=pivot.columns.tolist(), y=pivot.index.tolist(),
        colorscale=[
            [0.0,  "#f0f0f0"],[0.25, "#f0f0f0"],
            [0.26, "#F1C40F"],[0.5,  "#F1C40F"],
            [0.51, "#E67E22"],[0.75, "#E67E22"],
            [0.76, "#C0392B"],[1.0,  "#C0392B"],
        ],
        zmin=-1, zmax=3,
        colorbar=dict(tickvals=[-1,0,1,2,3],
                      ticktext=["none","MODIFIER","LOW","MODERATE","HIGH"]),
        hoverongaps=False,
    ))
    fig_hm.update_layout(title=title, xaxis_tickangle=-45,
                         margin=dict(t=50,b=0),
                         height=max(400, len(pivot)*18))
    st.plotly_chart(fig_hm, use_container_width=True)

with tab1:
    hm_nc = sdf[
        (sdf["GT_CLASS"].isin(["het","hom_alt"])) &
        (sdf["IS_NOVEL"]) &
        (sdf["IMPACT"].isin(["HIGH","MODERATE"])) &
        (~sdf["SYMBOL"].isin(HAMLET_GENES)) &
        (sdf["SYMBOL"] != "")
    ]
    make_heatmap(hm_nc, "Sample × Non-Canonical Gene (Novel HIGH/MOD)")

with tab2:
    hm_h = sdf[
        (sdf["GT_CLASS"].isin(["het","hom_alt"])) &
        (sdf["IS_AML_GENE"])
    ]
    with st.expander("Heatmap filters"):
        hm_nov = st.checkbox("Novel only", value=False, key="hm_h_novel")
        hm_imp = st.multiselect("Impacts", ["HIGH","MODERATE","LOW","MODIFIER"],
                                default=["HIGH","MODERATE"], key="hm_h_imp")
    if hm_nov:
        hm_h = hm_h[hm_h["IS_NOVEL"]]
    if hm_imp:
        hm_h = hm_h[hm_h["IMPACT"].isin(hm_imp)]
    make_heatmap(hm_h, "Sample × HAMLET Gene")
