# MA_papers

This repository contains a set of Streamlit applications and related data used to explore AML (acute myeloid leukemia) datasets and perform literature searches.

## Contents
- `Home.py` – main Streamlit application used to chat with OpenAI assistants.
- `pages/` – additional Streamlit pages:
  - `2_DepMap` – interactive visualization of gene dependency scores from DepMap.
  - `3_DepMap Correlation Analysis` – correlation-based ranking of genes.
  - `4_Venn Diagram` – venn diagram comparisons of gene sets.
  - `PubMed.py` – PubMed query helper using Biopython.
- `Depmap_AML_subset.csv` – subset of DepMap data focused on AML cell lines.
- `requirements.txt` – Python dependencies for running the apps.

## Setup
1. Create a Python environment (Python 3.11 or later recommended).
2. Install dependencies:
   ```bash
   pip install -r requirements.txt
   ```
3. Obtain an OpenAI API key and store it in Streamlit secrets (`.streamlit/secrets.toml`) as `OPENAI_API_KEY`.

## Running
Run the main Streamlit app:
```bash
streamlit run Home.py
```
This will open a web browser with the available pages where you can explore the DepMap data or search PubMed.

## Notes
The repository uses a large CSV dataset (`Depmap_AML_subset.csv`). Ensure this file remains in the root directory when running the apps.
