from flask import Flask
from flask import request
from Bio import Entrez

email = ''

app = Flask(__name__)

@app.route("/")
def entrance():
    doi = request.args.getlist('doi')
    abstract_dict = {}
    without_abstract = []
    Entrez.email = email
    handle = Entrez.efetch(db="pubmed", id=''.join(doi),rettype="xml", retmode="text") 
    records = Entrez.read(handle)
    for pubmed_article in records['PubmedArticle']:
        pmid = int(str(pubmed_article['MedlineCitation']['PMID']))
        article = pubmed_article['MedlineCitation']['Article']
        if 'Abstract' in article:
            abstract = article['Abstract']['AbstractText'][0]
            abstract_dict[pmid] = abstract
        else:
            without_abstract.append(pmid)
    handle.close()
    return f'{abstract_dict}\n{without_abstract}'