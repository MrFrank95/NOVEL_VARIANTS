"""
06_lollyplot.py  —  Interactive Variant Lollyplot
==================================================
Lollypop chart where:
  • Y-axis (stem height) = Cohort AF  (how common is the variant across the 31 samples)
  • Circle diameter       = Mean VAF   (how much of each carrier's reads carry the alt)
  • Click a lolly         → filters the sample table below

Performance notes
-----------------
* VAF is pre-computed in vcf_loader.get_per_sample_stats() — no re-parsing here.
* All heavy data lives in st.session_state, loaded once by app.py.
* Interactivity uses st.plotly_chart(on_select="rerun") — no external package needed.

Glossary
--------
Cohort AF   : INFO/AF — fraction of all alleles in the 31-sample cohort that
              carry the alternate allele.  High value → common across cohort.
VAF         : Variant Allele Frequency for a single sample, computed as
              alt_reads / (ref_reads + alt_reads) from the AD FORMAT field.
              Shown as a percentage (0–100 %).
Mean VAF    : Average VAF across all samples that carry the variant
              (het or hom-alt).  Used to size the lolly circle.
HGVSp       : Protein-level HGVS notation (e.g. p.Asp835Tyr). Used as the
              x-axis label; falls back to genomic position when absent.
"""

import streamlit as st
import pandas as pd
import plotly.graph_objects as go

# ─────────────────────────────────────────────────────────────────────────────
# DATA — read from session_state populated by app.py
# ─────────────────────────────────────────────────────────────────────────────

if "df" not in st.session_state:
    st.error("Data not loaded. Please start from the main app page.")
    st.stop()

df_full    = st.session_state["df"]         # variant-level DataFrame (one row per variant)
sample_df  = st.session_state["sample_df"]  # per-sample DataFrame  (variants × samples)
samples    = st.session_state["samples"]

# ─────────────────────────────────────────────────────────────────────────────
# SIDEBAR FILTERS
# ─────────────────────────────────────────────────────────────────────────────

st.sidebar.markdown("### Lollyplot Filters")

pass_only = st.sidebar.checkbox("PASS only", value=True, key='chk_pass_lolly')
novel_only = st.sidebar.checkbox(
    "Novel variants only",
    value=False,
    help="Restrict to variants absent from dbSNP / ClinVar (canonical transcript)",
)

min_cohort_af = st.sidebar.slider(
    "Min Cohort AF",
    min_value=0.0, max_value=1.0, value=0.0, step=0.01,
    help="Only show variants with Cohort AF ≥ this value",
)

impact_opts = ["HIGH", "MODERATE", "LOW", "MODIFIER"]
sel_impacts = st.sidebar.multiselect(
    "VEP Impact",
    options=impact_opts,
    default=["HIGH", "MODERATE","LOW","MODIFIER"],
    help="HIGH = protein-disrupting · MODERATE = missense / in-frame indel",
)

tier_opts = sorted(df_full["GENE_TIER"].dropna().unique().tolist())
sel_tiers = st.sidebar.multiselect(
    "Gene category",
    options=tier_opts,
    default=tier_opts,
    help="HAMLET canonical = 41-gene AML panel · Non-canonical = everything else",
)

# ─────────────────────────────────────────────────────────────────────────────
# APPLY FILTERS TO THE VARIANT-LEVEL DF
# ─────────────────────────────────────────────────────────────────────────────

df_v = df_full.copy()
if pass_only:
    df_v = df_v[df_v["FILTER"] == "PASS"]
if sel_impacts:
    df_v = df_v[df_v["IMPACT"].isin(sel_impacts)]
if sel_tiers:
    df_v = df_v[df_v["GENE_TIER"].isin(sel_tiers)]
if novel_only:
    df_v = df_v[df_v["IS_NOVEL"]]
df_v = df_v[df_v["SYMBOL"].notna() & (df_v["SYMBOL"] != "")]
df_v = df_v[df_v["AF"].notna() & (df_v["AF"] >= min_cohort_af)]

# Narrow sample_df to matching variants (CHROM+POS join is most reliable)
variant_keys = set(zip(df_v["CHROM"], df_v["POS"]))
sdf = sample_df[
    sample_df.apply(lambda r: (r["CHROM"], r["POS"]) in variant_keys, axis=1)
].copy()

# Keep only carriers
sdf_carriers = sdf[sdf["GT_CLASS"].isin(["Heterozygous", "Homozygous Alternate"])].copy()

# ─────────────────────────────────────────────────────────────────────────────
# PAGE HEADER
# ─────────────────────────────────────────────────────────────────────────────

st.title("🍭 Variant Lollyplot Explorer")
st.markdown(
    "**Stem height** = Cohort AF &nbsp;·&nbsp; "
    "**Circle diameter** = Mean VAF across carriers &nbsp;·&nbsp; "
    "**Click a lolly** to filter the sample table below."
)

with st.expander("📖 Glossary", expanded=False):
    st.markdown("""
| Term | Definition |
|------|-----------|
| **Cohort AF** | INFO/AF — fraction of alleles in this 31-sample cohort carrying the alt allele. |
| **VAF (%)** | Variant Allele Frequency for a single sample: `alt_reads / total_reads × 100`, from the AD FORMAT field. |
| **Mean VAF** | Average VAF across all carrier samples (het or hom-alt). Sizes the lolly circle. |
| **HGVSp** | Protein-level HGVS notation (e.g. `p.Asp835Tyr`). X-axis label; falls back to genomic position. |
| **HAMLET canonical** | One of the 41 AML driver genes in the HAMLET panel. |
| **Non-canonical** | Any gene outside the HAMLET panel — primary discovery targets. |
| **Novel** | Absent from dbSNP / ClinVar for the canonical VEP transcript. |
""")

if df_v.empty:
    st.warning("No variants match the current filters.")
    st.stop()

# ─────────────────────────────────────────────────────────────────────────────
# GENE SELECTOR
# ─────────────────────────────────────────────────────────────────────────────

# Build gene summary for the macro view and selector ordering
gene_summary = (
    df_v.groupby(["SYMBOL", "GENE_TIER"])
    .agg(Max_Cohort_AF=("AF", "max"), N_Variants=("POS", "nunique"))
    .reset_index()
    .sort_values("Max_Cohort_AF", ascending=False)
)

col_gene, col_spacer = st.columns([2, 3])
with col_gene:
    gene_list = gene_summary["SYMBOL"].tolist()
    sel_gene = st.selectbox(
        "Select gene",
        options=gene_list,
        format_func=lambda g: (
            f"{g}  (max cohort AF: "
            f"{gene_summary.loc[gene_summary['SYMBOL']==g, 'Max_Cohort_AF'].values[0]:.3f})"
        ),
    )

# ─────────────────────────────────────────────────────────────────────────────
# BUILD LOLLYPLOT DATA FOR SELECTED GENE
# ─────────────────────────────────────────────────────────────────────────────

gene_variants = df_v[df_v["SYMBOL"] == sel_gene].copy()
gene_carriers = sdf_carriers[sdf_carriers["SYMBOL"] == sel_gene].copy()

if gene_variants.empty:
    st.info(f"No variants found for **{sel_gene}** with current filters.")
    st.stop()

# X-axis label: HGVSp if available, else chr:pos ref>alt
def make_label(row):
    hp = str(row.get("HGVSp", "")).strip()
    # Strip transcript prefix if present (e.g. ENSP00000241453.7:p.Asp835Tyr → p.Asp835Tyr)
    if ":" in hp:
        hp = hp.split(":")[-1]
    if hp and hp != "." and hp != "nan":
        return hp
    return f"{row['CHROM']}:{row['POS']} {row['REF']}>{row['ALT']}"

gene_variants["_label"] = gene_variants.apply(make_label, axis=1)
gene_variants = gene_variants.sort_values("POS")  # biological order on x-axis

# Per-variant mean VAF from carriers
mean_vaf = (
    gene_carriers
    .groupby(["CHROM", "POS"])["VAF"]
    .mean()
    .reset_index()
    .rename(columns={"VAF": "Mean_VAF"})
)

lolly_df = gene_variants.merge(mean_vaf, on=["CHROM", "POS"], how="left")
lolly_df["Mean_VAF"] = lolly_df["Mean_VAF"].fillna(0)
lolly_df["N_Carriers"] = (
    gene_carriers.groupby(["CHROM", "POS"])["sample"]
    .nunique()
    .reindex(pd.MultiIndex.from_arrays([lolly_df["CHROM"], lolly_df["POS"]]))
    .values
)
lolly_df["N_Carriers"] = lolly_df["N_Carriers"].fillna(0).astype(int)

# Build a unique stable key for each lolly (used for click matching)
lolly_df["_key"] = lolly_df["_label"]

# ─────────────────────────────────────────────────────────────────────────────
# LOLLYPLOT FIGURE
# ─────────────────────────────────────────────────────────────────────────────

IMPACT_COLORS = {
    "HIGH":     "#C0392B",
    "MODERATE": "#E67E22",
    "LOW":      "#F1C40F",
    "MODIFIER": "#95A5A6",
}

# Size scaling: map Mean_VAF (0–100) to circle area
# sizeref formula: 2 * max_val / (desired_max_px ** 2)
max_vaf   = max(lolly_df["Mean_VAF"].max(), 1.0)
size_ref  = 2.0 * max_vaf / (40.0 ** 2)

fig = go.Figure()

# ── Stems (one vertical line per variant) ────────────────────────────────────
for _, row in lolly_df.iterrows():
    fig.add_shape(
        type="line",
        x0=row["_key"], y0=0,
        x1=row["_key"], y1=row["AF"],
        line=dict(color="#B0BEC5", width=1.5),
        layer="below",
    )

# ── Lolly circles (one scatter trace per impact so legend works) ─────────────
for impact, grp in lolly_df.groupby("IMPACT"):
    color  = IMPACT_COLORS.get(impact, "#95A5A6")
    novel_markers = ["star" if n else "circle" for n in grp["IS_NOVEL"]]

    # Build hover lines: add dbSNP / ClinVar links if available
    hover_texts = []
    for _, r in grp.iterrows():
        ex  = str(r.get("EXISTING_VARIATION", "") or "")
        clin = str(r.get("CLIN_SIG", "") or "")

        # Parse rsIDs from Existing_variation (may be comma/& separated)
        rs_ids = [x.strip() for x in ex.replace("&", ",").split(",")
                  if x.strip().startswith("rs")]
        db_links = " | ".join(
            f'<a href="https://www.ncbi.nlm.nih.gov/snp/{rs}">{rs}</a>'
            for rs in rs_ids
        ) if rs_ids else "—"

        hover_texts.append(
            f"<b>{r['_key']}</b><br>"
            f"Cohort AF: {r['AF']:.4f}<br>"
            f"Mean VAF: {r['Mean_VAF']:.1f}%<br>"
            f"Carriers: {r['N_Carriers']}<br>"
            f"Consequence: {r['CONSEQUENCE'].replace('_',' ')}<br>"
            f"ClinVar: {clin or '—'}<br>"
            f"dbSNP: {db_links}"
        )

    fig.add_trace(go.Scatter(
        x=grp["_key"],
        y=grp["AF"],
        mode="markers",
        name=impact,
        marker=dict(
            size=grp["Mean_VAF"],
            sizemode="area",
            sizeref=size_ref,
            sizemin=10,
            color=color,
            opacity=0.85,
            line=dict(width=1.5, color="white"),
            symbol=novel_markers,
        ),
        customdata=grp[["N_Carriers", "Mean_VAF", "CLIN_SIG",
                         "CONSEQUENCE", "EXISTING_VARIATION"]].values,
        text=hover_texts,
        hovertemplate="%{text}<extra></extra>",
        hoverlabel=dict(bgcolor="white", font_size=12),
    ))

fig.update_layout(
    title=dict(
        text=f"<b>{sel_gene}</b> — Variant Distribution",
        font=dict(size=16, color="#1B2A4A"),
    ),
    xaxis=dict(
        title="Protein consequence (HGVSp) / Genomic position",
        tickangle=-40,
        showgrid=False,
        zeroline=False,
    ),
    yaxis=dict(
        title="Cohort AF",
        gridcolor="rgba(200,200,200,0.4)",
        zeroline=True,
        zerolinecolor="#B0BEC5",
    ),
    plot_bgcolor="rgba(245,247,250,0.6)",
    paper_bgcolor="white",
    legend=dict(title="VEP Impact", orientation="v"),
    height=520,
    margin=dict(l=60, r=20, t=60, b=120),
    clickmode="event+select",
)

# ─────────────────────────────────────────────────────────────────────────────
# RENDER + CLICK HANDLING
# Uses st.plotly_chart(on_select="rerun") — built into Streamlit ≥ 1.33
# No external package required.
# ─────────────────────────────────────────────────────────────────────────────

st.caption(
    "🍭 **Circle size** = mean VAF across carriers · "
    "**Stem height** = cohort AF · "
    "**★** = novel variant · "
    "**Click a circle** to filter the table below"
)

event = st.plotly_chart(
    fig,
    width='stretch',
    on_select="rerun",         # triggers a rerun on point click
    selection_mode="points",
    key=f"lolly_{sel_gene}",
)

# Extract the clicked variant label
clicked_label = None
if event and event.selection and event.selection.points:
    clicked_label = event.selection.points[0].get("x")

# ─────────────────────────────────────────────────────────────────────────────
# SAMPLE TABLE (filtered by click or showing all carriers)
# ─────────────────────────────────────────────────────────────────────────────

st.markdown("---")
st.subheader("📋 Sample-Level Detail")

# Merge carrier data with the lolly labels
detail_df = gene_carriers.merge(
    gene_variants[["CHROM", "POS", "REF", "ALT", "_label",
                   "HGVSc", "CLIN_SIG", "EXISTING_VARIATION", "AF"]],
    on=["CHROM", "POS"],
    how="inner",
    suffixes=("", "_v"),
)

# If something was clicked, filter to that variant
if clicked_label:
    filtered = detail_df[detail_df["_label"] == clicked_label]
    st.success(f"📌 Filtered to: **{clicked_label}** — {len(filtered)} carrier samples")
    if st.button("✖ Clear selection"):
        st.rerun()
else:
    filtered = detail_df
    st.caption(
        f"Showing all {len(filtered)} carrier observations across {sel_gene}. "
        "Click a lolly above to filter."
    )

# Build display table
display_cols_map = {
    "sample":      "Sample",
    "_label":      "Variant (HGVSp)",
    "GT_CLASS":    "GT Class",
    "GT":          "GT",
    "VAF":         "VAF (%)",
    "AF":          "Cohort AF",
    "DP":          "DP",
    "GQ":          "GQ",
    "CLIN_SIG":    "ClinVar",
    "EXISTING_VARIATION": "dbSNP / ID",
}

show = (
    filtered[[c for c in display_cols_map if c in filtered.columns]]
    .rename(columns=display_cols_map)
    .sort_values("VAF (%)", ascending=False)
    .reset_index(drop=True)
)
show.index += 1

# GT_CLASS colour map
GT_COLORS = {
    "Heterozygous":       "#FFF2CC",
    "Homozygous Alternate":"#FFCDD2",
    "Homozygous Reference":"#C8E6C9",
    "Missing":            "#F5F5F5",
}

def _color_gt(val):
    bg = GT_COLORS.get(val, "")
    return f"background-color: {bg}; color: black" if bg else ""

fmt_dict = {}
if "VAF (%)" in show.columns:
    fmt_dict["VAF (%)"] = "{:.1f}"
if "Cohort AF" in show.columns:
    fmt_dict["Cohort AF"] = "{:.4f}"
if "DP" in show.columns:
    fmt_dict["DP"] = "{:.0f}"
if "GQ" in show.columns:
    fmt_dict["GQ"] = "{:.0f}"

styled = show.style
if "GT Class" in show.columns:
    styled = styled.map(_color_gt, subset=["GT Class"])
if fmt_dict:
    styled = styled.format(fmt_dict, na_rep="—")

st.dataframe(styled, width='stretch')

# Download button
csv_bytes = show.to_csv(index=False).encode()
st.download_button(
    "⬇️ Download table as CSV",
    csv_bytes,
    file_name=f"lollyplot_{sel_gene}_carriers.csv",
    mime="text/csv",
)