import streamlit as st
import pandas as pd
import matplotlib.pyplot as plt
from venn import venn, pseudovenn  # Import venn library
import re

def create_gene_presence_dataframe(gene_sets):
    """
    Creates a DataFrame showing gene presence in each gene set, including dataset count.

    Args:
        gene_sets (dict): Dictionary of gene sets, where keys are set names and values are sets of genes.

    Returns:
        pandas.DataFrame: DataFrame with genes as rows (now a 'Gene' column) and gene sets as columns,
                          showing boolean presence (True/False), and a 'Dataset Count' column.
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

    gene_presence_df = pd.DataFrame(gene_presence_data, index=sorted_genes, columns=sorted_gene_set_names)

    # Add 'Dataset Count' column
    dataset_counts = []
    for index, row in gene_presence_df.iterrows():
        dataset_counts.append(sum(row)) # Sum boolean values to get count of True values
    gene_presence_df['Dataset Count'] = dataset_counts

    # Reset index to make 'Gene' a regular column
    gene_presence_df = gene_presence_df.reset_index()
    gene_presence_df = gene_presence_df.rename(columns={'index': 'Gene'}) # Rename 'index' column to 'Gene'
    gene_presence_df['Gene'] = gene_presence_df['Gene'].apply(lambda x: re.sub(r'\'', '', x))
    gene_presence_df = gene_presence_df.sort_values(by='Dataset Count', ascending=False)

    return gene_presence_df


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


def parse_gene_string(gene_string):
    """
    Parses the gene string from the CSV to extract a list of genes.

    Args:
        gene_string (str): String containing genes, e.g., "[FCN1, NFKBIA, CD14, ...]".

    Returns:
        list: List of gene names.
    """
    if isinstance(gene_string, str):
        gene_string = gene_string.strip('[]') # Remove brackets
        genes = [gene.strip() for gene in gene_string.split(',') if gene.strip()] # Split by comma and remove empty strings
        return genes
    return [] # Return empty list if not a string


def main():
    st.title("Interactive Venn Diagram for Gene Sets (File & Direction)")

    st.markdown("""
    Upload a CSV file. Select up to 6 file/direction combinations from the sidebar to create a Venn Diagram.
    The app will generate an interactive Venn diagram using the 'venn' library based on genes listed
    in the 'genes' column for each selected file/direction combination.
    Detailed gene lists for intersections are not provided in this version for more than 3 sets.
    The 'Gene Presence in Gene Sets' table now includes a 'Dataset Count' column showing the number of selected datasets each gene is present in.
    The genes in the 'Gene Presence in Gene Sets' table are displayed without quotes.
    """)

    with st.sidebar: # Move file uploader to sidebar
        st.header("Data Upload and Selection")
        uploaded_file = st.file_uploader("Upload CSV file", type=["csv"])

    if uploaded_file is not None:
        try:
            df = pd.read_csv(uploaded_file)

            # Create 'file_direction' column for selection
            df['file_direction'] = df['file'] + " - " + df['direction']
            available_file_directions = df['file_direction'].unique().tolist()

            with st.sidebar: # Keep selection in sidebar
                selected_file_directions = st.multiselect(
                    "Select up to 6 file/direction combinations for Venn diagram:",
                    available_file_directions,
                    max_selections=6
                )
                submit_button = st.button("Submit") # Keep submit button in sidebar

            if submit_button:
                if not selected_file_directions:
                    st.warning("Please select at least 2 file/direction combinations to generate a Venn diagram.")
                elif len(selected_file_directions) > 6:
                    st.warning("Please select a maximum of 6 file/direction combinations for the Venn diagram.")
                else:
                    gene_sets = {}
                    for file_direction_name in selected_file_directions:
                        file_direction_df = df[df['file_direction'] == file_direction_name]
                        if not file_direction_df.empty: # Check if combination exists in df
                            gene_string = file_direction_df['genes'].iloc[0] # Take the first row
                            gene_list = parse_gene_string(gene_string)
                            gene_sets[file_direction_name] = set(gene_list)
                        else:
                            st.warning(f"File/Direction combination '{file_direction_name}' not found in the file.")
                            continue # Skip to the next combination if not found

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
        except KeyError as e:
            st.error(f"Error: Column '{e}' not found in CSV file. Please ensure your CSV file has the required columns ('file', 'direction', 'genes').")
        except Exception as e:
            st.error(f"An error occurred: {e}")

if __name__ == "__main__":
    main()
