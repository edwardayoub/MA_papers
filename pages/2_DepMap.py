import streamlit as st
import pandas as pd
import plotly.express as px
import scipy.stats as stats
import numpy as np

# --- Data Loading ---
@st.cache_data
def load_data(csv_path="Depmap_AML_subset.csv"):
    return pd.read_csv(csv_path)

# --- Prefix Stripping Utility ---
def strip_prefix(col_name: str):
    prefixes = [
        "Hotspot Mutations: ", "Damaging Mutations: ",
        "Hotspot Mutations", "Damaging Mutations",
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
    st.title("AML DepMap Explorer")

    st.markdown("""
    **Gene Scores from DepMap:** This app explores gene scores from the DepMap project, focusing on AML (Acute Myeloid Leukemia) cell lines. 
    Gene scores, in this context, represent the **dependency** of cancer cells on specific genes. 
    Lower scores (more negative) generally indicate that a gene is more essential for cell survival.
    """)

    st.markdown("""
    **Data Source:** The gene scores are derived from large-scale **RNA interference (RNAi)** and **CRISPR-Cas9** gene knockout screens conducted as part of the DepMap project.
    Specifically, the RNAi scores are from the **Achilles, DRIVE, and Marcotte** datasets (DEMETER2 algorithm), and CRISPR scores are from **DepMap Public 24Q4+Score** (Chronos algorithm).
    These screens systematically assess how knocking down (RNAi) or knocking out (CRISPR) each gene affects the viability of various cancer cell lines.
    """)

    # --- Data Loading ---
    df = load_data()

    # --- Sidebar for User Selections ---
    st.sidebar.header("Selections") # Adding a header to the sidebar for clarity

    mutation_type = st.sidebar.selectbox(
        "Select Mutation Type",
        ["Hotspot Mutations", "Damaging Mutations", "Both"],
        help="""
        **Mutation Type:** Choose the type of mutations to analyze.
        - **Hotspot Mutations:** Mutations frequently found in cancer.
        - **Damaging Mutations:** Mutations predicted to disrupt gene function.
        - **Both:** Analyze both Hotspot and Damaging mutations.
        """
    )


    # --- Column Filtering and Dictionary Creation ---
    if mutation_type != "Both":
        raw_mutation_cols = [
            col for col in df.columns if col.startswith(mutation_type)
        ]
        mutation_col_dict = {strip_prefix(c): c for c in raw_mutation_cols}
    else: # mutation_type == "Both"
        raw_mutation_cols = [
            col for col in df.columns
            if col.startswith("Hotspot Mutations") or col.startswith("Damaging Mutations")
        ]
        mutation_col_dict = {}
        for col in raw_mutation_cols:
            gene_name = strip_prefix(col).replace("Hotspot ", "").replace("Damaging ", "")
            mutation_col_dict.setdefault(gene_name, {})
            if "Hotspot" in col:
                mutation_col_dict[gene_name]['Hotspot'] = col
            elif "Damaging" in col:
                mutation_col_dict[gene_name]['Damaging'] = col

    raw_score_cols_rnai = [col for col in df.columns if col.startswith("RNAi")]
    raw_score_cols_crispr = [col for col in df.columns if col.startswith("CRISPR")]

    if not raw_mutation_cols:
        st.warning("No columns found matching the selected Mutation Type.")
        return
    if not raw_score_cols_rnai and not raw_score_cols_crispr:
        st.warning("No RNAi or CRISPR score columns found.")
        return

    score_col_dict_rnai = {strip_prefix(c): c for c in raw_score_cols_rnai}
    score_col_dict_crispr = {strip_prefix(c): c for c in raw_score_cols_crispr}

    combined_score_options = {}
    for gene, rnai_col in score_col_dict_rnai.items():
        combined_score_options.setdefault(gene, {})['RNAi'] = rnai_col
    for gene, crispr_col in score_col_dict_crispr.items():
        combined_score_options.setdefault(gene, {})['CRISPR'] = crispr_col

    score_gene_options = list(combined_score_options.keys())
    selected_score_gene = st.sidebar.selectbox(
        "Select Score Gene", score_gene_options,
        help="""
        **Score Gene:** Choose the gene for which you want to examine the gene dependency scores (RNAi and CRISPR scores).
        The app will then display how mutations in another gene (Mutation Gene selection below) relate to the dependency of this selected Score Gene.
        """
    )

    mutation_options = list(mutation_col_dict.keys())
    default_mutation_index = mutation_options.index("TP53") if "TP53" in mutation_options else 0
    selected_mut_display = st.sidebar.selectbox(
        "Select Mutation Gene", mutation_options, index=default_mutation_index,
        help="""
        **Mutation Gene:** Choose the gene for which you want to analyze mutations.
        The app will show how mutations in this gene correlate with the gene dependency scores of the Score Gene selected above.
        """
    )


    # --- Column Selection based on User Input ---
    selected_score_cols = combined_score_options[selected_score_gene]
    selected_score_col_rnai = selected_score_cols.get('RNAi')
    selected_score_col_crispr = selected_score_cols.get('CRISPR')
    selected_mut_col = mutation_col_dict[selected_mut_display]

    binary_option = st.sidebar.checkbox("Convert mutation column to binary", help="Convert mutation data into binary (Mutation Present/Absent) for simpler analysis.")

    if mutation_type == "Both":
        hotspot_col = selected_mut_col.get('Hotspot')
        damaging_col = selected_mut_col.get('Damaging')

        if hotspot_col and damaging_col:
            df["combined_mutation_col"] = df[hotspot_col].fillna(0) + df[damaging_col].fillna(0)
            x_col = "combined_mutation_col"
            x_title_base = f"{selected_mut_display} (Hotspot + Damaging)"
        elif hotspot_col:
            x_col = hotspot_col
            x_title_base = f"{selected_mut_display} (Hotspot)"
        elif damaging_col:
            x_col = damaging_col
            x_title_base = f"{selected_mut_display} (Damaging)"
        else:
            st.error("No Hotspot or Damaging columns found for selected gene.")
            return
    else:
        x_col = selected_mut_col
        x_title_base = selected_mut_display

    if binary_option:
        df["binary_mutation_col"] = df[x_col].astype(bool) # More concise binary conversion
        x_col = "binary_mutation_col"
        x_title = f"{x_title_base} (Binary)"
    else:
        x_title = x_title_base

    # --- Plotting ---
    num_cols = sum(1 for col in [selected_score_col_rnai, selected_score_col_crispr] if col) # Count valid score cols
    plot_cols = st.columns(num_cols) if num_cols > 1 else [st.columns(1)[0]] # Handle single or double column layout
    col_index = 0 # Index to track column for plotting

    if selected_score_col_rnai:
        with plot_cols[col_index]:
            col_index += 1 # Increment for next plot
            st.subheader(f"RNAi Scores for {selected_score_gene}")
            fig_rnai = px.box(
                df, x=x_col, y=selected_score_col_rnai,
                color=x_col if binary_option else None,
                color_discrete_sequence=px.colors.qualitative.Set1 if binary_option else None,
                labels={x_col: x_title, selected_score_col_rnai: "RNAi Score"}
            )
            fig_rnai.update_layout(showlegend=False)
            st.plotly_chart(fig_rnai, use_container_width=True)

    if selected_score_col_crispr:
        with plot_cols[col_index]: # Use the incremented index
            st.subheader(f"CRISPR Scores for {selected_score_gene}")
            fig_crispr = px.box(
                df, x=x_col, y=selected_score_col_crispr,
                color=x_col if binary_option else None,
                color_discrete_sequence=px.colors.qualitative.Set1 if binary_option else None,
                labels={x_col: x_title, selected_score_col_crispr: "CRISPR Score"}
            )
            fig_crispr.update_layout(showlegend=False)
            st.plotly_chart(fig_crispr, use_container_width=True)

    # --- Filtered Data Display ---
    with st.expander("Show Filtered Data"):
        cols_to_show = ["cell_line_display_name"]
        mutation_cols = []
        score_cols = []

        # Add mutation columns based on mutation_type and binary_option
        if mutation_type == "Both":
            if 'Hotspot' in selected_mut_col:
                mutation_cols.append(selected_mut_col['Hotspot'])
            if 'Damaging' in selected_mut_col:
                mutation_cols.append(selected_mut_col['Damaging'])
            if 'combined_mutation_col' in df.columns:
                mutation_cols.append("combined_mutation_col")
        else:
            mutation_cols.append(x_col) # x_col is already set to the correct mutation column

        if binary_option and 'binary_mutation_col' in df.columns: # Only add if binary option is checked and column exists
            mutation_cols.append('binary_mutation_col')

        # Always add score columns if they are selected
        if selected_score_col_rnai:
            score_cols.append(selected_score_col_rnai)
        if selected_score_col_crispr:
            score_cols.append(selected_score_col_crispr)

        # Construct ordered list of columns
        cols_to_show_ordered = ["cell_line_display_name"] + mutation_cols + score_cols
        # Remove potential duplicates, maintaining order as much as possible
        cols_to_show_unique_ordered = []
        seen_cols = set()
        for col in cols_to_show_ordered:
            if col not in seen_cols:
                cols_to_show_unique_ordered.append(col)
                seen_cols.add(col)


        df_filtered = df[cols_to_show_unique_ordered].copy() # Copy to avoid modifying original dataframe
        st.dataframe(df_filtered)

if __name__ == "__main__":
    main()
