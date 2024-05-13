import streamlit as st
import os
from Bio import Entrez
import pandas as pd

# Set your email here to let NCBI know who you are
Entrez.email = "your_email@example.com"

def fetch_articles(pubmed_ids):
    handle = Entrez.efetch(db="pubmed", id=pubmed_ids, retmode="xml")
    records = Entrez.read(handle)['PubmedArticle']
    articles = []
    for record in records:
        article = record['MedlineCitation']['Article']
        try:
            title = article['ArticleTitle']
        except KeyError:
            title = "No title available"
        try:
            abstract = article['Abstract']['AbstractText'][0]
        except KeyError:
            abstract = "No abstract available"
        try:
            authors_list = article['AuthorList']
            authors = ", ".join([f"{author['LastName']}, {author['ForeName']}" for author in authors_list if 'LastName' in author and 'ForeName' in author])
        except KeyError:
            authors = "No authors available"
        
        pmid = record['MedlineCitation']['PMID']
        link = f"https://pubmed.ncbi.nlm.nih.gov/{pmid}"
        articles.append({
            "PMID": pmid,
            "Title": title,
            "Authors": authors,
            "Abstract": abstract,
            "Link": link
        })
    return articles

def search_pubmed(query):
    handle = Entrez.esearch(db="pubmed", term=query, retmax=500)
    record = Entrez.read(handle)
    pubmed_ids = record['IdList']
    return pubmed_ids

# Streamlit application
st.title("PubMed Article Search")

query = st.text_input("Enter search keywords:")

if st.button("Search"):
    pubmed_ids = search_pubmed(query)
    if pubmed_ids:
        articles = fetch_articles(pubmed_ids)
        df = pd.DataFrame(articles)
        st.dataframe(df)
        st.success("Articles have been fetched and saved.")
    else:
        st.warning("No articles found for the given keywords.")
