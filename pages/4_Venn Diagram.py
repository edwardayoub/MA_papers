import streamlit as st
import pandas as pd
import matplotlib.pyplot as plt
from venn import venn, pseudovenn  # Import venn library

def create_gene_presence_dataframe(gene_sets):
    """
    Creates a DataFrame showing gene presence in each gene set.

    Args:
        gene_sets (dict): Dictionary of gene sets, where keys are set names and values are sets of genes.

    Returns:
        pandas.DataFrame: DataFrame with genes as rows and gene sets as columns,
                          showing boolean presence (True/False).
    """
    all_genes = set()
    for gene_set in gene_sets.values():
        all_genes.update(gene_set)

    gene_presence_data = []
    sorted_genes = sorted(list(all_genes))  # Sort genes alphabetically
    sorted_gene_set_names = sorted(list(gene_sets.keys())) # Sort gene set names

    for gene in sorted_genes:
        presence = [gene in gene_sets[gene_set_name] for gene_set_name in sorted_gene_set_names]
        gene_presence_data.append(presence)

    return pd.DataFrame(gene_presence_data, index=sorted_genes, columns=sorted_gene_set_names)


def create_geneset_overlap_dataframe(gene_sets):
    """
    Creates a DataFrame showing the number of common genes between each pair of gene sets.

    Args:
        gene_sets (dict): Dictionary of gene sets.

    Returns:
        pandas.DataFrame: DataFrame showing the overlap counts between gene sets.
    """
    gene_set_names = sorted(list(gene_sets.keys())) # Sort gene set names for consistent order
    overlap_data = []

    for i, set1_name in enumerate(gene_set_names):
        overlap_row = []
        for j, set2_name in enumerate(gene_set_names):
            common_genes = gene_sets[set1_name].intersection(gene_sets[set2_name])
            overlap_row.append(len(common_genes))  # Count of common genes
        overlap_data.append(overlap_row)

    return pd.DataFrame(overlap_data, index=gene_set_names, columns=gene_set_names)


def main():
    st.title("Interactive Venn Diagram for Gene Sets")

    st.markdown("""
    Upload a CSV file where columns represent gene sets.
    Select up to 6 columns from the sidebar to create a Venn Diagram.
    The app will generate an interactive Venn diagram using the 'venn' library.
    Detailed gene lists for intersections are not provided in this version for more than 3 sets.
    """)

    with st.sidebar: # Move file uploader to sidebar
        st.header("Data Upload and Selection")
        uploaded_file = st.file_uploader("Upload CSV file", type=["csv"])

    if uploaded_file is not None:
        try:
            df = pd.read_csv(uploaded_file)
            available_columns = df.columns.tolist()

            with st.sidebar: # Keep column selection in sidebar
                selected_columns = st.multiselect(
                    "Select up to 6 columns for Venn diagram:",
                    available_columns,
                    max_selections=6
                )
                submit_button = st.button("Submit") # Keep submit button in sidebar

            if submit_button:
                if not selected_columns:
                    st.warning("Please select at least 2 columns to generate a Venn diagram.")
                elif len(selected_columns) > 6:
                    st.warning("Please select a maximum of 6 columns for the Venn diagram.")
                else:
                    gene_set_columns = selected_columns # Use user selected columns
                    gene_sets = {}
                    for col in gene_set_columns:
                        gene_list = df[col].dropna().astype(str).tolist()
                        gene_sets[col] = set(gene_list)

                    if gene_sets:
                        set_names = list(gene_sets.keys())
                        dataset_dict = {name: gene_sets[name] for name in set_names} # create dataset_dict for venn

                        num_sets = len(gene_sets)

                        plt.figure(figsize=(8, 8)) # Create a matplotlib figure
                        if num_sets >= 2 and num_sets <= 5:
                            venn(dataset_dict) # Use venn for 2-5 sets
                        elif num_sets == 6:
                            pseudovenn(dataset_dict) # Use pseudovenn for 6 sets
                        else:
                            st.error("Venn diagram is only supported for 2 to 6 sets.")

                        plt.title("Gene Set Venn Diagram") # Add a title to the matplotlib plot
                        st.pyplot(plt) # Display the matplotlib figure in Streamlit
                        plt.close() # Close the matplotlib figure to prevent display issues

                        # Gene Presence Dataframe in expander
                        with st.expander("Gene Presence in Gene Sets"):
                            gene_presence_df = create_gene_presence_dataframe(gene_sets)
                            st.dataframe(gene_presence_df)

                        # Gene Set Overlap Dataframe in expander
                        with st.expander("Gene Set Overlap Count"):
                            overlap_df = create_geneset_overlap_dataframe(gene_sets)
                            st.dataframe(overlap_df)


        except pd.errors.ParserError:
            st.error("Error: Could not parse CSV file. Please ensure it is a valid CSV format.")
        except Exception as e:
            st.error(f"An error occurred: {e}")

if __name__ == "__main__":
    main()
