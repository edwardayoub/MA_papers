import streamlit as st
import os
from Bio import Entrez
import csv
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

def save_to_csv(articles, filename):
    with open(filename, mode='w', newline='', encoding='utf-8') as file:
        fieldnames = ["PMID", "Title", "Authors", "Abstract", "Link"]
        writer = csv.DictWriter(file, fieldnames=fieldnames)
        writer.writeheader()
        for article in articles:
            writer.writerow(article)

def save_articles_as_txt(articles):
    folder_path = 'articles'
    os.makedirs(folder_path, exist_ok=True)  # Create the folder if it does not exist
    
    for article in articles:
        file_path = os.path.join(folder_path, f"{article['PMID']}.txt")
        if not os.path.exists(file_path):  # Check if file already exists
            with open(file_path, 'w', encoding='utf-8') as file:
                file.write(f"Title: {article['Title']}\n")
                file.write(f"Authors: {article['Authors']}\n")
                file.write(f"Abstract: {article['Abstract']}\n")
                file.write(f"Link: {article['Link']}\n")

# Streamlit application
st.title("PubMed Article Search")

query = st.text_input("Enter search keywords:")

if st.button("Search"):
    pubmed_ids = search_pubmed(query)
    if pubmed_ids:
        articles = fetch_articles(pubmed_ids)
        df = pd.DataFrame(articles)
        st.dataframe(df)

        save_to_csv(articles, 'MA_articles.csv')
        save_articles_as_txt(articles)
        st.success("Articles have been fetched and saved.")
    else:
        st.warning("No articles found for the given keywords.")
