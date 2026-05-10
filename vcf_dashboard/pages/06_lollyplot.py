"""
06_lollyplot.py  —  Lollyplot con posición proteica numérica + dominios UniProt
================================================================================
Ejes:
  • Eje X  = posición aminoácido (numérico, de Protein_position en VEP/CSQ)
  • Eje Y  = N° de muestras portadoras de esa variante
  • Tamaño = VAF promedio de los portadores (%)
  • Color  = tipo de variante (missense, frameshift, nonsense, deleción, etc.)
  • Borde  = impacto VEP (HIGH, MODERATE, LOW, MODIFIER)

Panel inferior:
  • Dominios proteicos obtenidos de la API REST de UniProt (entradas revisadas).
"""

import streamlit as st
import pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import requests

# ─────────────────────────────────────────────────────────────────────────────
# DATOS  — desde session_state (cargados por app.py)
# ─────────────────────────────────────────────────────────────────────────────

if "df" not in st.session_state:
    st.error("Datos no cargados. Inicia desde la página principal (app.py).")
    st.stop()

df_full   = st.session_state["df"]
sample_df = st.session_state["sample_df"]
samples   = st.session_state["samples"]

# ─────────────────────────────────────────────────────────────────────────────
# PALETAS DE COLORES
# ─────────────────────────────────────────────────────────────────────────────

# Mapeo consecuencia VEP → tipo de variante legible
CONSEQUENCE_TYPE_MAP = {
    "missense_variant":       "missense",
    "stop_gained":            "nonsense",
    "stop_lost":              "nonsense",
    "frameshift_variant":     "frameshift",
    "inframe_insertion":      "inframe_indel",
    "inframe_deletion":       "inframe_indel",
    "splice_donor_variant":   "splicing",
    "splice_acceptor_variant":"splicing",
    "splice_region_variant":  "splicing",
    "start_lost":             "nonsense",
    "synonymous_variant":     "sinónimo",
    "coding_sequence_variant":"otra_codificante",
    "protein_altering_variant":"otra_codificante",
}

# Colores de relleno de los círculos (tipo de variante)
VARTYPE_COLORS = {
    "missense":         "#E8613C",
    "nonsense":         "#C0392B",
    "frameshift":       "#F39C12",
    "inframe_indel":    "#8E44AD",
    "splicing":         "#2980B9",
    "sinónimo":         "#27AE60",
    "otra_codificante": "#95A5A6",
    "other":            "#BDC3C7",
}

# Colores de borde de los círculos (impacto VEP) — no repiten con VARTYPE_COLORS
IMPACT_BORDER_COLOR = {
    "HIGH":     "#840032",   # rojo oscuro
    "MODERATE": "#002642",   # naranja quemado
    "LOW":      "#e59500",   # verde oscuro
    "MODIFIER": "#6f5e53",   # gris pizarra
}

# Ancho de borde según impacto
IMPACT_BORDER_WIDTH = {
    "HIGH":     3.0,
    "MODERATE": 2.0,
    "LOW":      1.5,
    "MODIFIER": 1.2,
}

# Colores de fondo para la tabla (impacto VEP)
IMPACT_TEXT = {}

# Colores de background para la tabla (impacto VEP)
IMPACT_BG = {
    "HIGH":     "#d00000",
    "MODERATE": "#81a4cd",
    "LOW":      "#ffba08",
    "MODIFIER": "#7b6b43",
}

# Colores de fondo para la tabla (genotipo)
GT_COLORS = {
    "Heterocigoto":      "#FFF2CC",
    "Homocigoto Alt.":   "#FFCDD2",
    "Homocigoto Ref.":   "#C8E6C9",
    "Sin llamada":       "#F5F5F5",
}

# ─────────────────────────────────────────────────────────────────────────────
# HELPERS
# ─────────────────────────────────────────────────────────────────────────────

def classify_consequence(csq_str: str) -> str:
    """Devuelve el tipo canónico del primer término de consecuencia."""
    if not csq_str:
        return "other"
    for term in csq_str.split("&"):
        if term in CONSEQUENCE_TYPE_MAP:
            return CONSEQUENCE_TYPE_MAP[term]
    return "other"


@st.cache_data(show_spinner="Consultando UniProt…", ttl=86_400)
def fetch_uniprot_domains(gene_symbol: str) -> tuple[list[dict], int]:
    """
    Retorna (lista_de_dominios, longitud_proteína).
    Cada dominio: {"name": str, "start": int, "end": int, "type": str}
    Usa la API REST v2 de UniProt buscando la entrada Swiss-Prot humana canónica.
    """
    base   = "https://rest.uniprot.org/uniprotkb/search"
    params = {
        "query":  f"gene_exact:{gene_symbol} AND reviewed:true AND organism_id:9606",
        "fields": "accession,ft_domain,ft_region,ft_motif,length",
        "format": "json",
        "size":   1,
    }
    try:
        r = requests.get(base, params=params, timeout=10)
        r.raise_for_status()
        data    = r.json()
        results = data.get("results", [])
        if not results:
            return [], 0

        entry   = results[0]
        seq_len = entry.get("sequence", {}).get("length", 0)
        features = entry.get("features", [])

        domains = []
        for feat in features:
            ftype = feat.get("type", "")
            if ftype not in ("Domain", "Region", "Motif"):
                continue
            loc   = feat.get("location", {})
            start = loc.get("start", {}).get("value")
            end   = loc.get("end",   {}).get("value")
            name  = feat.get("description", ftype)
            if start is None or end is None:
                continue
            domains.append({
                "name":  name,
                "start": int(start),
                "end":   int(end),
                "type":  ftype,
            })
        return domains, seq_len

    except Exception as e:
        st.warning(f"No se pudieron obtener dominios UniProt para {gene_symbol}: {e}")
        return [], 0


def make_hgvsp_short(hp: str) -> str:
    """p.Asp835Tyr → p.D835Y  (notación de 1 letra)."""
    AA3 = {
        "Ala":"A","Arg":"R","Asn":"N","Asp":"D","Cys":"C","Gln":"Q",
        "Glu":"E","Gly":"G","His":"H","Ile":"I","Leu":"L","Lys":"K",
        "Met":"M","Phe":"F","Pro":"P","Ser":"S","Thr":"T","Trp":"W",
        "Tyr":"Y","Val":"V","Ter":"*",
    }
    if not hp or hp in (".", "nan"):
        return hp
    if ":" in hp:
        hp = hp.split(":")[-1]
    for aa3, aa1 in AA3.items():
        hp = hp.replace(aa3, aa1)
    return hp

# ─────────────────────────────────────────────────────────────────────────────
# FILTROS SIDEBAR
# ─────────────────────────────────────────────────────────────────────────────

st.sidebar.markdown("### Filtros Lollyplot")

pass_only  = st.sidebar.checkbox("Solo PASS", value=True, key="chk_pass_lolly")
novel_only = st.sidebar.checkbox("Solo variantes nuevas", value=False)
min_carriers = st.sidebar.slider(
    "Mín. muestras portadoras",
    min_value=1, max_value=len(samples), value=1, step=1,
    help="Solo mostrar variantes presentes en ≥ N muestras",
)

impact_opts = ["HIGH", "MODERATE", "LOW", "MODIFIER"]
sel_impacts = st.sidebar.multiselect(
    "Impacto VEP",
    options=impact_opts,
    default=["HIGH", "MODERATE"],
)

tier_opts = sorted(df_full["GENE_TIER"].dropna().unique().tolist())
sel_tiers = st.sidebar.multiselect(
    "Categoría de gen",
    options=tier_opts,
    default=tier_opts,
)

# ─────────────────────────────────────────────────────────────────────────────
# APLICAR FILTROS
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

variant_keys = set(zip(df_v["CHROM"], df_v["POS"]))
sdf = sample_df[
    sample_df.apply(lambda r: (r["CHROM"], r["POS"]) in variant_keys, axis=1)
].copy()
sdf_carriers = sdf[sdf["GT_CLASS"].isin(["Heterozygous", "Homozygous Alternate"])].copy()

# ─────────────────────────────────────────────────────────────────────────────
# ENCABEZADO
# ─────────────────────────────────────────────────────────────────────────────

st.title("Lollyplot de Variantes")

with st.expander("📖 Glosario — ¿qué significan estos campos?"):
    st.markdown("""
    **Identificadores de variante**
    - **HGVSc** — Notación HGVS a nivel de *ADN codificante* (ej. `c.123A>T`). Describe el cambio de base respecto al transcrito.
    - **HGVSp** — Notación HGVS a nivel de *proteína* (ej. `p.Lys41Asn`). Muestra la consecuencia aminoácida.
    - **VARIANT_CLASS** — Tipo estructural de la variante: SNV (sustitución de un solo nucleótido), inserción, deleción, MNV, etc.
    - **BIOTYPE** — Biotype de Ensembl del transcrito solapado (ej. `protein_coding`, `lncRNA`).

    **Impacto y consecuencia**
    - **IMPACT** — Nivel de severidad de VEP: `HIGH` (probable disrupción proteica: stop-gain, frameshift, sitio de splicing), `MODERATE` (missense, indel en marco), `LOW` (sinónimo), `MODIFIER` (no codificante / intergénico).
    - **CONSECUENCIA** — Efecto molecular específico predicho por VEP (ej. `missense_variant`, `splice_acceptor_variant`).
    - **SIFT** — Predice si una sustitución aminoácida afecta la función proteica. Puntuación < 0.05 = **deletéreo** (tolerado en otro caso).
    - **PolyPhen** — Predice el efecto de un cambio missense en la estructura/función proteica. Valores: `probably_damaging`, `possibly_damaging`, `benign`.

    **Puntuaciones de calidad**
    - **VQSLOD** — Log-odds de calidad de variante de VQSR de GATK. Más alto = más probable variante verdadera; valores negativos indican baja confianza.
    - **QD** — Calidad por profundidad: puntuación de calidad de variante dividida por profundidad de lectura. QD bajo (< 2) puede indicar falsos positivos.
    - **FS** — FisherStrand: mide sesgo de hebra usando el test exacto de Fisher. Valores altos (> 60 en SNPs, > 200 en indels) sugieren artefactos de hebra.
    - **MQ** — Calidad de mapeo RMS de lecturas que soportan la variante. MQ bajo (< 40) puede reflejar lecturas mal mapeadas.

    **Frecuencia poblacional**
    - **gnomADe AF** — Frecuencia alélica en la cohorte de *exomas* de gnomAD (> 125 000 exomas). Ausente = no observado en gnomAD, lo que refuerza la novedad.
    - **MAX_AF** — Frecuencia alélica máxima en todas las sub-poblaciones de gnomAD. Útil para detectar variantes comunes en una ascendencia específica aunque sean raras en general.

    **Clínica / bases de datos**
    - **ClinSig (CLIN_SIG)** — Significado clínico de ClinVar: `pathogenic`, `likely_pathogenic`, `uncertain_significance`, `benign`, etc. Vacío = no está en ClinVar.
    - **IS_NOVEL** — `True` si la variante no tiene entrada en `Existing_variation` (rsID de dbSNP / accesión de ClinVar) para el transcrito canónico.

    **Campos a nivel de muestra**
    - **GT** — Genotipo (ej. `0/1` = heterocigoto, `1/1` = homocigoto alt.).
    - **GT_CLASS** — Etiqueta de genotipo simplificada: `het`, `hom_alt`, `hom_ref`, `desconocido`.
    - **DP** — Profundidad total de lectura en el sitio de la variante en esa muestra.
    - **GQ** — Calidad de genotipo: confianza en la llamada de genotipo escalada en Phred. GQ ≥ 20 se considera generalmente fiable.
    """)

st.markdown(
    "**Eje X** = Posición del aminoácido &nbsp;·&nbsp; "
    "**Eje Y** = N° de muestras portadoras &nbsp;·&nbsp; "
    "**Tamaño** = VAF promedio (%) &nbsp;·&nbsp; "
    "**Color de relleno** = tipo de variante &nbsp;·&nbsp; "
    "**Color de borde** = impacto VEP"
)

if df_v.empty:
    st.warning("Ninguna variante coincide con los filtros actuales.")
    st.stop()

# ─────────────────────────────────────────────────────────────────────────────
# SELECTOR DE GEN
# ─────────────────────────────────────────────────────────────────────────────

# ─────────────────────────────────────────────────────────────────────────────
# SELECTOR DE GEN
# ─────────────────────────────────────────────────────────────────────────────
gene_summary = (
    df_v.groupby(["SYMBOL", "GENE_TIER"])
    .agg(N_Variantes=("POS", "nunique"), Max_AF=("AF", "max"))
    .reset_index()
    .sort_values("N_Variantes", ascending=False)
)

gene_list = gene_summary["SYMBOL"].tolist()

col_gene, _ = st.columns([2, 3])
with col_gene:
    busqueda = st.text_input(
        "Buscar gen",
        placeholder="Ej: TP53, BRCA1…",
        help="Escribe parte del nombre del gen para filtrar la lista",
    )

    # Filtrar lista según texto ingresado (insensible a mayúsculas)
    if busqueda.strip():
        gene_list_filtrada = [
            g for g in gene_list
            if busqueda.strip().upper() in g.upper()
        ]
    else:
        gene_list_filtrada = gene_list

    if not gene_list_filtrada:
        st.warning(f"Ningún gen coincide con **{busqueda}**.")
        st.stop()

    sel_gene = st.selectbox(
        "Seleccionar gen",
        options=gene_list_filtrada,
        format_func=lambda g: (
            f"{g}  "
            f"({gene_summary.loc[gene_summary['SYMBOL']==g, 'N_Variantes'].values[0]} variantes)"
        ),
        help=f"{len(gene_list_filtrada)} gen(es) encontrados",
    )
# ─────────────────────────────────────────────────────────────────────────────
# CONSTRUIR DATOS DEL LOLLYPLOT PARA EL GEN SELECCIONADO
# ─────────────────────────────────────────────────────────────────────────────

gene_variants = df_v[df_v["SYMBOL"] == sel_gene].copy()
gene_carriers = sdf_carriers[sdf_carriers["SYMBOL"] == sel_gene].copy()

if gene_variants.empty:
    st.info(f"No se encontraron variantes para **{sel_gene}** con los filtros actuales.")
    st.stop()

# Clasificar tipo de variante
gene_variants["VAR_CLASS"] = gene_variants["CONSEQUENCE"].apply(classify_consequence)

# Etiqueta HGVSp corta (para hover)
def make_label(row):
    hp = str(row.get("HGVSp", "")).strip()
    if ":" in hp:
        hp = hp.split(":")[-1]
    if hp and hp not in (".", "nan"):
        return make_hgvsp_short(hp)
    return f"{row['CHROM']}:{row['POS']} {row['REF']}>{row['ALT']}"

gene_variants["_label"] = gene_variants.apply(make_label, axis=1)
gene_variants = gene_variants.sort_values("Protein_position")

# Conteo de portadores por variante (CHROM + POS)
carrier_counts = (
    gene_carriers.groupby(["CHROM", "POS"])["sample"]
    .nunique()
    .reset_index()
    .rename(columns={"sample": "N_Carriers"})
)

# VAF promedio por variante
mean_vaf = (
    gene_carriers.groupby(["CHROM", "POS"])["VAF"]
    .mean()
    .reset_index()
    .rename(columns={"VAF": "Mean_VAF"})
)

lolly_df = (
    gene_variants
    .merge(carrier_counts, on=["CHROM", "POS"], how="left")
    .merge(mean_vaf,       on=["CHROM", "POS"], how="left")
)
lolly_df["N_Carriers"] = lolly_df["N_Carriers"].fillna(0).astype(int)
lolly_df["Mean_VAF"]   = lolly_df["Mean_VAF"].fillna(0)

# Usar posición proteica; si no existe, caer a posición genómica
lolly_df["_xpos"] = lolly_df["Protein_position"]
no_prot = lolly_df["_xpos"].isna()
if no_prot.any():
    lolly_df.loc[no_prot, "_xpos"] = lolly_df.loc[no_prot, "POS"]

# Filtrar por mínimo de portadores
lolly_df = lolly_df[lolly_df["N_Carriers"] >= min_carriers]

if lolly_df.empty:
    st.info(f"Ninguna variante de **{sel_gene}** tiene ≥ {min_carriers} portadores.")
    st.stop()

# ─────────────────────────────────────────────────────────────────────────────
# DOMINIOS UNIPROT
# ─────────────────────────────────────────────────────────────────────────────

domains, prot_len = fetch_uniprot_domains(sel_gene)

if prot_len == 0:
    prot_len = int(lolly_df["_xpos"].max() * 1.05) if not lolly_df.empty else 500

x_min = 0
x_max = max(prot_len, int(lolly_df["_xpos"].max()) + 10) if not lolly_df.empty else prot_len

# ─────────────────────────────────────────────────────────────────────────────
# FIGURA: dos filas — lollyplot + dominios
# ─────────────────────────────────────────────────────────────────────────────

has_domains  = len(domains) > 0
row_heights  = [0.85, 0.15] if has_domains else [1.0]
n_rows       = 2 if has_domains else 1

fig = make_subplots(
    rows=n_rows, cols=1,
    row_heights=row_heights,
    shared_xaxes=True,
    vertical_spacing=0.04,
    subplot_titles=[
        f"Lollyplot de variantes en {sel_gene}",
        "Dominios proteicos (UniProt)" if has_domains else "",
    ],
)

# ── Tallos (uno por variante) ──────────────────────────────────────────────────
for _, row in lolly_df.iterrows():
    fig.add_shape(
        type="line",
        x0=row["_xpos"], y0=0,
        x1=row["_xpos"], y1=row["N_Carriers"],
        line=dict(color="#90A4AE", width=1.2),
        row=1, col=1,
    )

# ── Círculos (una traza por tipo de variante) ──────────────────────────────────
max_vaf  = max(lolly_df["Mean_VAF"].max(), 1.0)
size_ref = 2.0 * max_vaf / (38.0 ** 2)

for vclass, grp in lolly_df.groupby("VAR_CLASS"):
    fill_color   = VARTYPE_COLORS.get(vclass, VARTYPE_COLORS["other"])
    novel_symbol = ["star" if n else "circle" for n in grp["IS_NOVEL"]]

    # Borde del círculo según impacto VEP
    border_colors = grp["IMPACT"].map(
        lambda v: IMPACT_BORDER_COLOR.get(v, "#888888")
    ).tolist()
    border_widths = grp["IMPACT"].map(
        lambda v: IMPACT_BORDER_WIDTH.get(v, 1.2)
    ).tolist()

    # Textos de hover
    # Textos de hover
    hover_texts = []
    for _, r in grp.iterrows():
        ex   = str(r.get("EXISTING_VARIATION", "") or "")
        clin = str(r.get("CLIN_SIG", "") or "")
        rs_ids = [
            x.strip() for x in ex.replace("&", ",").split(",")
            if x.strip().startswith("rs")
        ]
        db_links = " | ".join(
            f'<a href="https://www.ncbi.nlm.nih.gov/snp/{rs}">{rs}</a>'
            for rs in rs_ids
        ) if rs_ids else "—"
 
        impact_val = str(r.get("IMPACT", "") or "—")
 
        hover_texts.append(
            f"<b>{r['_label']}</b><br>"
            f"Posición aa: {int(r['_xpos']) if pd.notna(r['_xpos']) else '?'}<br>"
            f"Portadoras: {r['N_Carriers']}<br>"
            f"VAF promedio: {r['Mean_VAF']:.1f}%<br>"
            f"AF Cohorte: {r['AF']:.4f}%<br>"
            f"Consecuencia: {r['CONSEQUENCE'].replace('_',' ')}<br>"
            f"Impacto VEP: <b>{impact_val}</b><br>"
            f"ClinVar: {clin or '—'}<br>"
            f"dbSNP: {db_links}"
        )

    fig.add_trace(
        go.Scatter(
            x=grp["_xpos"],
            y=grp["N_Carriers"],
            mode="markers",
            name=vclass,
            legendgroup="tipo",
            legendgrouptitle_text="Tipo de variante",
            marker=dict(
                size=grp["Mean_VAF"],
                sizemode="area",
                sizeref=size_ref,
                sizemin=8,
                color=fill_color,
                opacity=0.88,
                line=dict(
                    width=border_widths,
                    color=border_colors,
                ),
                symbol=novel_symbol,
            ),
            text=hover_texts,
            hovertemplate="%{text}<extra></extra>",
            hoverlabel=dict(bgcolor="white", font_size=12),
        ),
        row=1, col=1,
    )

# ── Leyenda de impacto (puntos fantasma para la leyenda de bordes) ─────────────
for imp_label, imp_color in IMPACT_BORDER_COLOR.items():
    fig.add_trace(
        go.Scatter(
            x=[None], y=[None],
            mode="markers",
            name=imp_label,
            legendgroup="impacto",
            legendgrouptitle_text="Impacto VEP",
            marker=dict(
                size=14,
                color="rgba(200,200,200,0.35)",
                line=dict(
                    width=IMPACT_BORDER_WIDTH.get(imp_label, 1.5),
                    color=imp_color,
                ),
                symbol="circle",
            ),
            showlegend=True,
            hoverinfo="skip",
        ),
        row=1, col=1,
    )

# ── Track de dominios ──────────────────────────────────────────────────────────
DOMAIN_COLOR_MAP = {
    "Domain": "rgba(58, 117, 175, 0.80)",
    "Region": "rgba(39, 120, 77,  0.70)",
    "Motif":  "rgba(195, 90, 22,  0.70)",
}

if has_domains:
    # Backbone proteína completa
    fig.add_trace(
        go.Scatter(
            x=[0, prot_len, prot_len, 0, 0],
            y=[0.42, 0.42, 0.58, 0.58, 0.42],
            fill="toself",
            mode="lines",
            line=dict(width=0),
            fillcolor="rgba(189,195,199,0.35)",
            showlegend=False,
            hoverinfo="skip",
        ),
        row=2, col=1,
    )

    # Un rectángulo por dominio
    for dom in domains:
        color = DOMAIN_COLOR_MAP.get(dom["type"], "rgba(58,117,175,0.75)")
        fig.add_trace(
            go.Scatter(
                x=[dom["start"], dom["end"], dom["end"], dom["start"], dom["start"]],
                y=[0.15, 0.15, 0.85, 0.85, 0.15],
                fill="toself",
                mode="lines",
                line=dict(width=0.5, color="white"),
                fillcolor=color,
                name=dom["name"],
                showlegend=False,
                hovertemplate=(
                    f"<b>{dom['name']}</b><br>"
                    f"Tipo: {dom['type']}<br>"
                    f"Posición: {dom['start']}–{dom['end']} aa<extra></extra>"
                ),
                hoveron="fills",
            ),
            row=2, col=1,
        )

# ── Layout global ──────────────────────────────────────────────────────────────
fig.update_layout(
    height=600 if has_domains else 480,
    plot_bgcolor="rgba(245,247,250,0.7)",
    paper_bgcolor="white",
    legend=dict(
        title=None,
        orientation="v",
        x=1.02, y=1,
        tracegroupgap=16,
    ),
    margin=dict(l=60, r=190, t=60, b=40),
)

# Eje X compartido
if has_domains:
    fig.update_xaxes(
        range=[x_min - prot_len * 0.01, x_max * 1.01],
        showgrid=False,
        row=2, col=1,
    )

# Eje Y lollyplot
fig.update_yaxes(
    title_text="N° de muestras portadoras",
    gridcolor="rgba(200,200,200,0.35)",
    zeroline=True,
    zerolinecolor="#B0BEC5",
    rangemode="tozero",
    row=1, col=1,
)

# Eje Y dominios (oculto, rango fijo [0,1])
if has_domains:
    fig.update_yaxes(
        visible=False,
        range=[0, 1],
        row=2, col=1,
    )

# ── Sub-leyenda VAF dentro del plot ───────────────────────────────────────────
y_max      = int(lolly_df["N_Carriers"].max())
vaf_legend = [10, 25, 50, 75, 100]
x_leg      = x_max * 1.2
y_spacing  = max(y_max * 0.12, 1)
y_base     = y_max * 0.97

fig.add_trace(
    go.Scatter(
        x=[x_leg] * len(vaf_legend),
        y=[y_base - i * y_spacing for i in range(len(vaf_legend))],
        mode="markers+text",
        text=[f" {v}%" for v in vaf_legend],
        textposition="middle right",
        textfont=dict(size=10, color="#555"),
        marker=dict(
            size=vaf_legend,
            sizemode="area",
            sizeref=size_ref,
            sizemin=8,
            color="rgba(120,120,120,0.35)",
            line=dict(width=1, color="#888"),
        ),
        hoverinfo="skip",
        showlegend=False,
    ),
    row=1, col=1,
)
fig.add_annotation(
    x=x_leg, y=y_base + y_spacing * 0.8,
    text="<b>VAF media (%)</b>",
    showarrow=False,
    font=dict(size=10, color="#333"),
    xanchor="center",
    row=1, col=1,
)

# ─────────────────────────────────────────────────────────────────────────────
# RENDER + SELECCIÓN POR CLIC
# ─────────────────────────────────────────────────────────────────────────────

st.caption(
    "**Tamaño del círculo** = VAF promedio de portadores &nbsp;·&nbsp; "
    "**Color de relleno** = tipo de variante &nbsp;·&nbsp; "
    "**Color de borde** = impacto VEP &nbsp;·&nbsp; "
    "**★** = Variante nueva &nbsp;·&nbsp; "
    "**Haz clic** en un círculo para filtrar la tabla inferior"
)

event = st.plotly_chart(
    fig,
    width="stretch",
    on_select="rerun",
    selection_mode="points",
    key=f"lolly_{sel_gene}",
)

clicked_x = None
if event and event.selection and event.selection.points:
    clicked_x = event.selection.points[0].get("x")

# ─────────────────────────────────────────────────────────────────────────────
# TABLA DE MUESTRAS
# ─────────────────────────────────────────────────────────────────────────────

st.markdown("---")
st.subheader("📋 Detalle por muestra")
# Enriquecer carriers con etiqueta y datos de variante
detail_df = gene_carriers.merge(
    gene_variants[[
        "CHROM", "POS", "REF", "ALT", "_label",
        "Protein_position", "HGVSc", "CLIN_SIG",
        "EXISTING_VARIATION", "AF", "VAR_CLASS", "IMPACT",  # ← IMPACT incluido
    ]],
    on=["CHROM", "POS"],
    how="inner",
    suffixes=("", "_v"),
)

if clicked_x is not None:
    tol    = 0
    mask   = (detail_df["Protein_position"] - clicked_x).abs() == tol
    filt   = detail_df[mask]
    label_ = filt["_label"].iloc[0] if not filt.empty else str(clicked_x)
    st.success(f"📌 Filtrado a: **{label_}** — {len(filt)} muestras portadoras")
    #if st.button("✖ Limpiar selección"):
    #    st.rerun()
else:
    filt = detail_df
    st.caption(
        f"Mostrando las {len(filt)} observaciones portadoras de **{sel_gene}**. "
        "Haz clic en un lolly para filtrar."
    )

display_cols_map = {
    "sample":           "Muestra",
    "_label":           "Variante (HGVSp)",
    "Protein_position": "Posición aa",
    "IMPACT":           "Impacto VEP",      # ← nueva columna
    "GT_CLASS":         "Genotipo",
    "GT":               "GT",
    "VAF":              "VAF (%)",
    "AF":               "AF Cohorte",
    "DP":               "DP",
    "GQ":               "GQ",
    "CLIN_SIG":         "ClinVar",
}

show = (
    filt[[c for c in display_cols_map if c in filt.columns]]
    .rename(columns=display_cols_map)
    .sort_values(["Posición aa", "VAF (%)"], ascending=[True, False])
    .reset_index(drop=True)
)
show.index += 1

# ── Traducir valores de genotipo ───────────────────────────────────────────
GT_TRANSLATE = {
    "Heterozygous":         "Heterocigoto",
    "Homozygous Alternate": "Homocigoto Alt.",
    "Homozygous Reference": "Homocigoto Ref.",
    "Missing":              "Sin llamada",
}
if "Genotipo" in show.columns:
    show["Genotipo"] = show["Genotipo"].map(lambda v: GT_TRANSLATE.get(v, v))

# ── Construir URL de dbSNP (primer rsID disponible) ───────────────────────
def _first_rs_url(ex_str):
    if not ex_str or pd.isna(ex_str):
        return None
    rs_ids = [
        x.strip() for x in str(ex_str).replace("&", ",").split(",")
        if x.strip().startswith("rs")
    ]
    return f"https://www.ncbi.nlm.nih.gov/snp/{rs_ids[0]}" if rs_ids else None

show["dbSNP"] = filt["EXISTING_VARIATION"].values[: len(show)]
show["dbSNP"] = show["dbSNP"].apply(_first_rs_url)

# ── Funciones de estilo ────────────────────────────────────────────────────
def _color_gt(val):
    bg = GT_COLORS.get(val, "")
    return f"background-color: {bg}; color: black" if bg else ""

def _color_impact(val):
    bg   = IMPACT_BG.get(val, "")
    text = IMPACT_TEXT.get(val, "black")
    return f"background-color: {bg}; color: {text}" if bg else ""

# ── Formato numérico ───────────────────────────────────────────────────────
fmt_dict = {}
for col, fmt in [
    ("VAF (%)", "{:.1f}"),
    ("AF Cohorte", "{:.4f}"),
    ("DP", "{:.0f}"),
    ("GQ", "{:.0f}"),
    ("Posición aa", "{:.0f}"),
]:
    if col in show.columns:
        fmt_dict[col] = fmt

# ── Aplicar estilos ────────────────────────────────────────────────────────
styled = show.style
if "Genotipo" in show.columns:
    styled = styled.map(_color_gt, subset=["Genotipo"])
if "Impacto VEP" in show.columns:
    styled = styled.map(_color_impact, subset=["Impacto VEP"])
if fmt_dict:
    styled = styled.format(fmt_dict, na_rep="—")

# ── Renderizar tabla ───────────────────────────────────────────────────────
st.dataframe(
    styled,
    width="stretch",
    column_config={
        "dbSNP": st.column_config.LinkColumn(
            "dbSNP",
            display_text=r"rs\d+",
            help="Enlace a la entrada en dbSNP (NCBI)",
        ),
        "VAF (%)":      st.column_config.NumberColumn(format="%.1f"),
        "AF Cohorte":   st.column_config.NumberColumn(format="%.4f"),
        "DP":           st.column_config.NumberColumn(format="%d"),
        "GQ":           st.column_config.NumberColumn(format="%d"),
        "Posición aa":  st.column_config.NumberColumn(format="%d"),
    },
    hide_index=False,
)

csv_bytes = show.to_csv(index=False).encode()
st.download_button(
    "⬇️ Descargar tabla CSV",
    csv_bytes,
    file_name=f"lollyplot_{sel_gene}_portadores.csv",
    mime="text/csv",
)
