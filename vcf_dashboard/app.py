import streamlit as st
import pandas as pd
import os, sys

st.set_page_config(
    page_title="AML Variant Dashboard",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded",
)

# ── Cached data loader ──────────────────────────────────────────────────────
@st.cache_data(show_spinner="Parsing VCF… this takes ~30s the first time")
def load_data():
    sys.path.insert(0, os.path.dirname(__file__))
    from vcf_loader import load_vcf, get_per_sample_stats
    vcf_path = os.path.join(os.path.dirname(__file__), "data", "joint_germline_recalibrated_VEP_ann.vcf")
    df, samples = load_vcf(vcf_path)
    sample_df = get_per_sample_stats(df, samples)
    return df, sample_df, samples

df, sample_df, samples = load_data()
st.session_state["df"] = df
st.session_state["sample_df"] = sample_df
st.session_state["samples"] = samples

# ── Sidebar ─────────────────────────────────────────────────────────────────
st.sidebar.image("https://upload.wikimedia.org/wikipedia/commons/thumb/1/16/DNA_orbit_animated.gif/200px-DNA_orbit_animated.gif", width=60)
st.sidebar.title("🧬 AML Variant Dashboard")
st.sidebar.caption("Joint Germline · GATK 4.6 · VEP 115 · GRCh38")
st.sidebar.markdown("---")

pages = {
    "🏠 Overview": "pages/01_overview.py",
    "🔴 Novel AML Candidates": "pages/02_novel_candidates.py",
    #"📊 Quality Metrics": "pages/03_quality_metrics.py",
    #"👥 Per-Sample Burden": "pages/04_per_sample.py",
    "🔍 Variant Browser": "pages/05_browser.py",
}

st.sidebar.markdown("### Navigation")
page = st.sidebar.radio("", list(pages.keys()), label_visibility="collapsed")

# ── Global filter ────────────────────────────────────────────────────────────
st.sidebar.markdown("---")
st.sidebar.markdown("### Global Filters")
filter_pass = st.sidebar.checkbox("PASS only", value=True)
if filter_pass:
    df_view = df[df["FILTER"] == "PASS"].copy()
    sample_df_view = sample_df[sample_df["FILTER"] == "PASS"].copy()
else:
    df_view = df.copy()
    sample_df_view = sample_df.copy()

st.session_state["df_view"] = df_view
st.session_state["sample_df_view"] = sample_df_view

n_total = len(df_view)
n_novel = df_view["IS_NOVEL"].sum()
n_aml = df_view["IS_AML_GENE"].sum()
st.sidebar.metric("Variants shown", f"{n_total:,}")
st.sidebar.metric("Novel (no population DB)", f"{n_novel:,}")
st.sidebar.metric("In AML gene panel", f"{n_aml:,}")

# ── Routing ──────────────────────────────────────────────────────────────────
import importlib.util, pathlib

page_file = pathlib.Path(__file__).parent / pages[page]
spec = importlib.util.spec_from_file_location("page_module", page_file)
mod = importlib.util.module_from_spec(spec)
spec.loader.exec_module(mod)
