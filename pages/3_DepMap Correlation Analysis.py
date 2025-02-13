import streamlit as st
import pandas as pd
import plotly.express as px
import scipy.stats as stats
import numpy as np

# st.set_page_config(layout="wide")
# --- Data Loading ---
@st.cache_data
def load_data(csv_path="data/Depmap_AML_subset.csv"): # You might need to update the CSV path if your CNV data is in a different file
    return pd.read_csv(csv_path)

# --- Prefix Stripping Utility ---
def strip_prefix(col_name: str):
    prefixes = [
        "Hotspot Mutations: ", "Damaging Mutations: ", "CN Gene Public 24Q4: ", # ADDED CNV prefix
        "Hotspot Mutations", "Damaging Mutations", "Omics Absolute CN Gene Public 24Q4 ", # ADDED CNV prefix (without colon and space, for completeness)
        "RNAi: ", "CRISPR: ",
        "RNAi (Achilles+DRIVE+Marcotte, DEMETER2) ",
        "CRISPR (DepMap Public 24Q4+Score, Chronos)"
    ]
    for prefix in prefixes:
        if col_name.startswith(prefix):
            return col_name.replace(prefix, "").strip()
    return col_name

# --- Main Application ---
def main():
    st.title("AML DepMap Gene Ranking")

    st.markdown("""
    **Gene Ranking based on Mutation and Dependency Scores from DepMap:**
    This app ranks genes based on their correlation with a selected mutation type and gene,
    using gene dependency scores (RNAi and CRISPR) from the DepMap project, focusing on AML cell lines.
    """)

    # --- Data Loading ---
    df = load_data()

    # --- Sidebar for User Selections ---
    st.sidebar.header("Selections")

    mutation_type = st.sidebar.selectbox(
        "Select Mutation Type",
        ["Hotspot Mutations", "Damaging Mutations", "CNV Data", "All"], # UPDATED dropdown options
        index = 3, # UPDATED index to default to "All"
        help="""
        **Mutation Type:** Choose the type of mutations or genomic alteration to analyze.
        - **Hotspot Mutations:** Mutations frequently found in cancer.
        - **Damaging Mutations:** Mutations predicted to disrupt gene function.
        - **CNV Data:** Copy Number Variation data from Omics Absolute CN Gene Public 24Q4.
        - **All:** Analyze all mutation and CNV types.
        """
    )

    # --- Column Filtering and Dictionary Creation for Mutations---
    if mutation_type != "All": # UPDATED from "Both" to "All"
        raw_mutation_cols = []
        if mutation_type == "Hotspot Mutations":
            raw_mutation_cols = [
                col for col in df.columns if col.startswith(mutation_type)
            ]
        elif mutation_type == "Damaging Mutations":
            raw_mutation_cols = [
                col for col in df.columns if col.startswith(mutation_type)
            ]
        elif mutation_type == "CNV Data": # ADDED CNV Data condition
            raw_mutation_cols = [
                col for col in df.columns if col.startswith("Omics Absolute CN Gene Public 24Q4 ")
            ]
        mutation_col_dict = {strip_prefix(c): c for c in raw_mutation_cols}
    else: # mutation_type == "All" # UPDATED from "Both" to "All"
        raw_mutation_cols = [
            col for col in df.columns
            if col.startswith("Hotspot Mutations") or
               col.startswith("Damaging Mutations") or
               col.startswith("Omics Absolute CN Gene Public 24Q4 ") # ADDED CNV Data
        ]
        mutation_col_dict = {}
        for col in raw_mutation_cols:
            gene_name = strip_prefix(col).replace("Hotspot ", "").replace("Damaging ", "").replace("Omics Absolute CN Gene Public 24Q4  ", "") # UPDATED strip for CNV
            mutation_col_dict.setdefault(gene_name, {})
            if "Hotspot" in col:
                mutation_col_dict[gene_name]['Hotspot'] = col
            elif "Damaging" in col:
                mutation_col_dict[gene_name]['Damaging'] = col
            elif "Omics Absolute CN Gene Public 24Q4 " in col: # ADDED CNV Data
                mutation_col_dict[gene_name]['CNV'] = col

    mutation_options = list(mutation_col_dict.keys())
    default_mutation_index = mutation_options.index("TP53") if "TP53" in mutation_options else 0
    selected_mut_display = st.sidebar.selectbox(
        "Select Mutation Gene", mutation_options, index=default_mutation_index,
        help="Choose the gene for mutation analysis." # Updated help text to include CNV implicitly through "mutation analysis"
    )

    binary_option = st.sidebar.checkbox("Convert mutation column to binary", help="Convert mutation data into binary (Mutation Present/Absent).")

    # --- NA Cutoff Controls ---
    na_cutoff_rnai = st.sidebar.number_input(
        "Max NA values for RNAi", min_value=0, value=30, # reasonable default, adjust as needed
        help="Maximum number of NA values allowed for RNAi score. Genes exceeding this will be skipped in RNAi analysis."
    )
    na_cutoff_crispr = st.sidebar.number_input(
        "Max NA values for CRISPR", min_value=0, value=30, # reasonable default, adjust as needed
        help="Maximum number of NA values allowed for CRISPR score. Genes exceeding this will be skipped in CRISPR analysis."
    )

    submit_button = st.sidebar.button("Submit")

    raw_score_cols_rnai = [col for col in df.columns if col.startswith("RNAi")]
    raw_score_cols_crispr = [col for col in df.columns if col.startswith("CRISPR")]
    score_col_dict_rnai = {strip_prefix(c): c for c in raw_score_cols_rnai}
    score_col_dict_crispr = {strip_prefix(c): c for c in raw_score_cols_crispr}
    combined_score_options = {}
    for gene, rnai_col in score_col_dict_rnai.items():
        combined_score_options.setdefault(gene, {})['RNAi'] = rnai_col
    for gene, crispr_col in score_col_dict_crispr.items():
        combined_score_options.setdefault(gene, {})['CRISPR'] = crispr_col

    if submit_button:
        st.subheader("Gene Ranking Results")

        # --- Determine Mutation Column based on User Input ---
        if mutation_type == "All": # UPDATED from "Both" to "All"
            hotspot_col = mutation_col_dict[selected_mut_display].get('Hotspot')
            damaging_col = mutation_col_dict[selected_mut_display].get('Damaging')
            cnv_col = mutation_col_dict[selected_mut_display].get('CNV') # ADDED CNV column

            if hotspot_col and damaging_col and cnv_col: # UPDATED condition to include CNV
                df["combined_mutation_col"] = df[hotspot_col].fillna(0) + df[damaging_col].fillna(0) + df[cnv_col].fillna(0) # UPDATED combination to include CNV
                x_col = "combined_mutation_col"
            elif hotspot_col and damaging_col:
                df["combined_mutation_col"] = df[hotspot_col].fillna(0) + df[damaging_col].fillna(0)
                x_col = "combined_mutation_col"
            elif hotspot_col and cnv_col:
                df["combined_mutation_col"] = df[hotspot_col].fillna(0) + df[cnv_col].fillna(0)
                x_col = "combined_mutation_col"
            elif damaging_col and cnv_col:
                df["combined_mutation_col"] = df[damaging_col].fillna(0) + df[cnv_col].fillna(0)
                x_col = "combined_mutation_col"
            elif hotspot_col:
                x_col = hotspot_col
            elif damaging_col:
                x_col = damaging_col
            elif cnv_col: # ADDED CNV condition
                x_col = cnv_col
            else:
                st.error("No Hotspot, Damaging, or CNV columns found for selected gene.") # UPDATED error message to include CNV
                return
        elif mutation_type == "CNV Data": # ADDED CNV Data condition
            x_col = mutation_col_dict[selected_mut_display]
        else: # Hotspot or Damaging Mutations (no change needed here)
            x_col = mutation_col_dict[selected_mut_display]

        if binary_option:
            df["binary_mutation_col"] = df[x_col].astype(bool).astype(int)
            x_col = "binary_mutation_col"

        ranking_data = []
        score_gene_options = list(combined_score_options.keys())

        with st.spinner("Calculating gene rankings..."):
            for score_gene in score_gene_options:
                selected_score_cols = combined_score_options[score_gene]
                selected_score_col_rnai = selected_score_cols.get('RNAi')
                selected_score_col_crispr = selected_score_cols.get('CRISPR')

                rnai_corr, rnai_pval = np.nan, np.nan
                crispr_corr, crispr_pval = np.nan, np.nan

                # --- NA Filtering ---
                if selected_score_col_rnai:
                    na_count_rnai = df[selected_score_col_rnai].isna().sum()
                    if na_count_rnai > na_cutoff_rnai:
                        selected_score_col_rnai = None # Skip RNAi analysis for this gene

                if selected_score_col_crispr:
                    na_count_crispr = df[selected_score_col_crispr].isna().sum()
                    if na_count_crispr > na_cutoff_crispr:
                        selected_score_col_crispr = None # Skip CRISPR analysis for this gene
                # --- End NA Filtering ---


                if selected_score_col_rnai:
                    df_rnai_clean = df[[x_col, selected_score_col_rnai]].dropna()
                    if not df_rnai_clean.empty and len(df_rnai_clean[x_col].unique()) > 1: # avoid error when all mutation values are the same
                        rnai_corr, rnai_pval = stats.pearsonr(df_rnai_clean[x_col], df_rnai_clean[selected_score_col_rnai])

                if selected_score_col_crispr:
                    df_crispr_clean = df[[x_col, selected_score_col_crispr]].dropna()
                    if not df_crispr_clean.empty and len(df_crispr_clean[x_col].unique()) > 1: # avoid error when all mutation values are the same
                        crispr_corr, crispr_pval = stats.pearsonr(df_crispr_clean[x_col], df_crispr_clean[selected_score_col_crispr])

                ranking_data.append({
                    "Gene": score_gene,
                    "RNAi Correlation": rnai_corr,
                    "RNAi P-value": rnai_pval,
                    "CRISPR Correlation": crispr_corr,
                    "CRISPR P-value": crispr_pval
                })

        ranking_df = pd.DataFrame(ranking_data)

        # Rank by RNAi Correlation, then CRISPR if RNAi is NaN
        ranking_df['RNAi Rank'] = ranking_df['RNAi Correlation'].rank(ascending=False, na_option='bottom').astype(int)
        ranking_df['CRISPR Rank'] = ranking_df['CRISPR Correlation'].rank(ascending=False, na_option='bottom').astype(int)
        ranking_df['Combined Rank'] = np.where(ranking_df['RNAi Rank'].isna(), ranking_df['CRISPR Rank'], ranking_df['RNAi Rank']) # prioritize RNAi rank


        ranked_df = ranking_df.sort_values(by='Combined Rank', ascending=True) # sort by combined rank

        st.dataframe(ranked_df, column_config={ # improve display and formatting
            "RNAi Correlation": st.column_config.NumberColumn(format="%.2f"),
            "RNAi P-value": st.column_config.NumberColumn(format="%.3f"),
            "CRISPR Correlation": st.column_config.NumberColumn(format="%.2f"),
            "CRISPR P-value": st.column_config.NumberColumn(format="%.3f"),
        }, use_container_width=True)


if __name__ == "__main__":
    main()
